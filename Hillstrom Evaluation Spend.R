# This script evaluates several different models on the Hillstrom data set which is contained in this repository.
# The models evaluated are Contextual Treatment Selection (CTS), Causal Forest, Separate Model Approach with 
# Random Forest and a custom Tree/Random Forest with two different gain functions ("Simple" and "Frac). More 
# information about the custom Tree and Random Forest can be found under ModelImplementations/DecisionTreeImplementation.R
# The user can specify the parameter n_predictions. If it is set to one, each model is trained once on the 
# original data set and then evaluated. If n_predictions is greater than 1, there will be n_predictions iterations.
# For each iteration a bootstrap sample is taken as the new data and then the models are build. After the models
# have been built and the predictions have been made, the predictions are evaluated and the results ploted.

library(ggplot2)
library(caret)
library(plyr)
library(dplyr)
library(reshape2)

source('ModelImplementations/DecisionTreeImplementation.R')
source('ModelImplementations/RzepakowskiTree.R')
source('Evaluation Methods.R')
source('ModelImplementations/CausalTree.R')
source('ModelImplementations/Separate Model Approach.R')
source('ModelImplementations/ContextualTreatmentSelection.R')
source('ModelImplementations/VisualizationHelper.R')
source("ModelImplementations/PredictionFunctions.R")


set.seed(1234)
n_predictions <- 25
#Data import and preprocessing
email <- read.csv('Data/Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$conversion <- email$segment <- NULL

response <- 'spend'
control <- "control"
treatment_list <- c('men_treatment','women_treatment')

original_email <- email

bootstrap_df <- read.csv("bootstrap.csv")
test_df <- read.csv("test.csv")
# The training and prediction part
for(f in 1:n_predictions){
  
  # If n_predictions is > 1 as bootstrap sample is created
  if(n_predictions > 1){
    email <- original_email[bootstrap_df[,f],]
  }
  idx <- test_df[,f]
  train <- email[-idx, ]
  
  test <- email[idx, ]
  
  test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                     "newbie","channel")],TRUE, max_cases = 10)
  
  start_time <- Sys.time()
  for(c in c("simple","frac","max")){
    print(c)
    #Random Forest
    forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 100,n_features = 3,
                                           criterion = c,remain_cores = 2)
    pred <- predict_forest_df(forest,test)
    write.csv(pred, paste("Predictions/Hillstrom/random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  }

  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control,ntree = 1000,
                                                s_rule = "TOT", s_true = T)
  write.csv(causal_forest_pred, paste("Predictions/Hillstrom/causal_forest",as.character(f),".csv",sep = ""),
            row.names = FALSE)

  # Separate Model Approach
  pred_sma_rf <- dt_models(train, response, "anova",treatment_list,control,test,"rf", mtry = 3, ntree = 300)
  write.csv(pred_sma_rf, paste("Predictions/Hillstrom/sma rf",as.character(f),".csv",sep = ""),
            row.names = FALSE)

  # CTS
  cts_forest <- build_cts(response, control, treatment_list, train, 100, nrow(train), 3, 0.15, 100, parallel = TRUE,
                          remain_cores = 1)
  pred <- predict_forest_df(cts_forest, test)
  write.csv(pred, paste("Predictions/Hillstrom/cts",as.character(f),".csv",sep = ""), row.names = FALSE)
  end_time <- Sys.time()
  print(difftime(end_time,start_time,units = "mins"))
}


# Here the predictions are evaluated. Additionally we look at the treatment distribution, to see which treatments
# are assigned how often by the models.
start_time <- Sys.time()
folder <- "Predictions/Hillstrom/"
outcomes <- c()
result_qini <- c()
decile_treated <- c()
for(model in c("random_forest","cts","sma rf","causal_forest")){
  if(sum(model == c("tree","random_forest")) > 0){
    for(c in c("frac")){
      for(f in 1:n_predictions){
        pred <- read.csv(paste(folder,model,"_",c,as.character(f),".csv",sep = ""))
        if(length(outcomes) == 0){
          outcomes <- c(new_expected_quantile_response(response,control,treatment_list,pred),
                        paste(model,"_",c,sep = ""))
          decile_treated <- cbind(decile_perc_treated(pred,treatment_list),
                                  rep(paste(model,"_",c,sep = ""),11*length(treatment_list)))
          result_qini <- cbind(qini_curve(pred,control,treatment_list),
                               paste(model,"_",c,sep = ""))
        } else{
          outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),
                                       paste(model,"_",c,sep = "")))
          decile_treated <- rbind(decile_treated,
                                  cbind(decile_perc_treated(pred,treatment_list),
                                        rep(paste(model,"_",c,sep = ""),11*length(treatment_list))))
          result_qini <- rbind(result_qini,cbind(qini_curve(pred,control,treatment_list),
                                                 paste(model,"_",c,sep = "")))
        }
      }
    }
  } else{
    colnames(result_qini) <- c("Percentile","Values","Treatment","model")
    for(f in 1:n_predictions){
      pred <- read.csv(paste(folder,model,as.character(f),".csv",sep = ""))
      outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),model))
      decile_treated <- rbind(decile_treated,
                              cbind(decile_perc_treated(pred,treatment_list),
                                    rep(model,11*length(treatment_list))))
      result_qini <- rbind(result_qini,cbind(qini_curve(pred,control,treatment_list),
                                             model))    
    }
  }
}
outcome_df <- data.frame(outcomes)
decile_treated_df <- data.frame(decile_treated)
colnames(outcome_df) <- c(0,10,20,30,40,50,60,70,80,90,100,"Model")
colnames(decile_treated_df) <- c("PercTreated","Treatment","Decile","Model")
rownames(outcome_df) <- 1:nrow(outcome_df)
rownames(decile_treated_df) <- 1:nrow(decile_treated_df)
for(c in 1:11){
  outcome_df[,c] <- as.numeric(as.character(outcome_df[,c]))
}
outcome_df[,12] <- as.character(outcome_df[,12])
decile_treated_df[,1] <- as.numeric(as.character(decile_treated_df[,1]))
decile_treated_df[,3] <- as.numeric(as.character(decile_treated_df[,3]))
colnames(result_qini) <- c("percentile","values","treatment","model")
start <- mean(result_qini[result_qini$percentile == 0.0,"values"])
finish <- mean(result_qini[result_qini$percentile == 1.0,"values"])
qini_random <- seq(start,finish,by = (finish-start)/10)
random_df <- cbind(seq(0,1,by=0.1),qini_random,"random","random")
colnames(random_df) <- c("percentile","values","treatment","model")
result_qini <- rbind(result_qini,random_df)
result_qini[,2] <- as.numeric(result_qini[,2])
result_qini[,1] <- as.numeric(result_qini[,1])
# result_qini[,4] <- as.character(result_qini[,4])
print(difftime(Sys.time(),start_time,units = "mins"))


#Visualize the results
visualize_qini_uplift(result_qini,type = "qini")
visualize_qini_uplift(result_qini,type = "qini",errorbars = F,multiplot = F)
visualize(outcome_df,n_treated = decile_treated_df,multiplot = T)
visualize(outcome_df,multiplot = F,errorbars = F)
