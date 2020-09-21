# This script evaluates several different models on the Hillstrom data set which is contained in this repository.
# The models evaluated are Contextual Treatment Selection (CTS), Causal Forest, Separate Model Approach with 
# Random Forest and a custom Tree/Random Forest with two different gain functions ("Simple" and "Frac). More 
# information about the custom Tree and Random Forest can be found under Algorithm Implementations/DecisionTreeImplementation.R
# The user can specify the parameter n_predictions. If it is set to one, each model is trained once on the 
# original data set and then evaluated. If n_predictions is greater than 1, there will be n_predictions iterations.
# For each iteration a bootstrap sample is taken as the new data and then the models are build. After the models
# have been built and the predictions have been made, the predictions are evaluated and the results ploted.

library(ggplot2)
library(caret)
library(plyr)
library(dplyr)
library(reshape2)

source('src/Algorithm Implementations/DOM.R')
source('src/Algorithm Implementations/RzepakowskiTree.R')
source('src/Algorithm Implementations/CausalTree.R')
source('src/Algorithm Implementations/Separate Model Approach.R')
source('src/Algorithm Implementations/ContextualTreatmentSelection.R')
source('src/Algorithm Implementations/PredictionFunctions.R')

source('src/Helper Functions/VisualizationHelper.R')
source('src/Helper Functions/Evaluation Methods.R')

set.seed(1234)
n_predictions <- 25
remain_cores <- 1
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

#Create and save the bootrap samples and train test splits. This is done so, if we want to change something
#on one model we can retrain and test it on the sample bootstrap samples in order for fair comparison with
#the other models
if(!file.exists("bootstrap.csv")){
  bootstrap_idx <- c()
  for(f in 1:n_predictions){
    bootstrap_idx <- cbind(bootstrap_idx,sample(nrow(original_email),nrow(original_email),replace = TRUE))
  }
  bootstrap_df <- data.frame(bootstrap_idx)
  write.csv(bootstrap_idx,"bootstrap.csv")
} else{
  bootstrap_df <- read.csv("bootstrap.csv")
}
if(!file.exists("test.csv")){
  test_idx <- c()
  for(f in 1:n_predictions){
    email <- original_email[bootstrap_df[,f],]
    test_idx <- cbind(test_idx,createDataPartition(y = email[ , response], p=0.2, list = FALSE))
  }
  test_df <- data.frame(test_idx)
  write.csv(test_idx,"test.csv")
} else{
  test_df <- read.csv("test.csv")
}
folder <- "Predictions/Hillstrom/"
# The training and prediction part
for(f in order(1:25,decreasing = T)){
  
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
  for(c in c("frac","max")){
    print(c)
    #Random Forest
    forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 500,n_features = 3,
                                           criterion = c,remain_cores = remain_cores)
    pred <- predict_forest_df(forest,test, treatment_list, control,remain_cores = remain_cores)
    write.csv(pred, paste(folder,"random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  }

  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control,ntree = 1000,
                                                s_rule = "TOT", s_true = T)
  write.csv(causal_forest_pred, paste(folder,"causal_forest",as.character(f),".csv",sep = ""),
            row.names = FALSE)

  # Separate Model Approach
  pred_sma_rf <- dt_models(train, response, "anova",treatment_list,control,test,"rf", mtry = 3, ntree = 500)
  write.csv(pred_sma_rf, paste(folder,"sma rf",as.character(f),".csv",sep = ""),
            row.names = FALSE)

  # CTS
  cts_forest <- build_cts(response, control, treatment_list, train, ntree = 500, nrow(train), m_try = 4,
                          n_reg = 4, min_split = 10, parallel = TRUE, remain_cores = remain_cores)
  pred <- predict_forest_df(forest = cts_forest,test_data = test, treatment_list =  treatment_list,
                            control =  control, remain_cores =  remain_cores,parallel_pred = F)
  write.csv(pred, paste(folder,"CTS",as.character(f),".csv",sep = ""), row.names = FALSE)
  end_time <- Sys.time()
  print(difftime(end_time,start_time,units = "mins"))
}

# Here the predictions are evaluated. Additionally we look at the treatment distribution, to see which treatments
# are assigned how often by the models.
start_time <- Sys.time()
outcomes <- c()
result_qini <- c()
decile_treated <- c()
for(model in c("random_forest","cts","sma rf","causal_forest","new_cts")){
  if(sum(model == c("random_forest")) > 0){
    for(c in c("frac",'max')){
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
    colnames(result_qini) <- c("Percentile","Values","model")
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
colnames(result_qini) <- c("percentile","values","model")
result_qini[,2] <- as.numeric(result_qini[,2])
result_qini[,1] <- as.numeric(result_qini[,1])
result_qini$model <- as.character(result_qini$model)
decile_treated_df$Model <- as.character(decile_treated_df$Model)
decile_treated_df$PercTreated <- decile_treated_df$PercTreated * 100
decile_treated_df$Decile <- decile_treated_df$Decile * 10
result_qini$percentile <- result_qini$percentile*100
outcome_df[outcome_df$Model == "cts","Model"] <- "CTS"
outcome_df[outcome_df$Model == "causal_forest","Model"] <- "Causal Forest"
outcome_df[outcome_df$Model == "random_forest_frac","Model"] <- "DOM"
outcome_df[outcome_df$Model == "sma rf","Model"] <- "SMA"
result_qini[result_qini$model == "cts","model"] <- "CTS"
result_qini[result_qini$model == "causal_forest","model"] <- "Causal Forest"
result_qini[result_qini$model == "random_forest_frac","model"] <- "DOM"
result_qini[result_qini$model == "sma rf","model"] <- "SMA"
decile_treated_df[decile_treated_df$Model == "cts","Model"] <- "CTS"
decile_treated_df[decile_treated_df$Model == "causal_forest","Model"] <- "Causal Forest"
decile_treated_df[decile_treated_df$Model == "random_forest_frac","Model"] <- "DOM"
decile_treated_df[decile_treated_df$Model == "sma rf","Model"] <- "SMA"
print(difftime(Sys.time(),start_time,units = "mins"))


new_qini <- result_qini[!(result_qini$model %in% c("random_forest_max")),]
new_qini <- new_qini[order(new_qini$model),]
colnames(new_qini) <- c("percentile","values","Model")
new_outcome <- outcome_df[!(outcome_df$Model %in% c("random_forest_max")),]
new_outcome <- new_outcome[order(new_outcome$Model),]

#Visualize the results
visualize_qini_uplift(new_qini,type = "qini",errorbars = F,multiplot = F,ylabel = "Cumulative Gained Spend")
visualize(new_outcome,ylabel = "Expected Amount Spend per Person",n_treated = decile_treated_df[!(decile_treated_df$Model %in% c("random_forest_max")),],multiplot = T)
visualize(new_outcome,ylabel = "Expected Amount Spend per Person",multiplot = F,errorbars = F)
outcome_boxplot(new_outcome[,2:12],"Expected Amount Spend per Customer")