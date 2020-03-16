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
source("ModelImplementations/RzepakowskiTree.R")


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

email$spend <- email$visit <- email$segment <- NULL

response <- 'conversion'
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
folder <- "Predictions/Conversion/"
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
  for(c in c("frac","max")){
    print(c)
    #Random Forest
    forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 500,n_features = 3,
                                           criterion = c, remain_cores = remain_cores)
    pred <- predict_forest_df(forest,test, treatment_list, control,remain_cores = remain_cores)
    write.csv(pred, paste(folder,"random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  }
  
  #Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control,ntree = 1000,
                                                s_rule = "TOT", s_true = T)
  write.csv(causal_forest_pred, paste(folder,"causal_forest",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  # Separate Model Approach
  pred_sma_rf <- dt_models(train, response, "class",treatment_list,control,test,"rf")
  write.csv(pred_sma_rf, paste(folder,"sma rf",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  # CTS
  cts_forest <- build_cts(response, control, treatment_list, train, 500, nrow(train), 5, 2, 100, parallel = TRUE,
                          remain_cores = remain_cores)
  pred <- predict_forest_df(cts_forest, test, treatment_list, control,remain_cores = remain_cores)
  write.csv(pred, paste(folder,"cts",as.character(f),".csv",sep = ""), row.names = FALSE)
  
  #Rzp
  for(div in c("binary_KL_divergence")){
    rzp_forest <- build_random_rzp_forest(train_data = train, treatment_list = treatment_list,
                                                   response = response,control = control, n_trees = 500, n_features = 3,
                                                   normalize = T, max_depth = 100,test_list = test_list, remain_cores = 12,
                                                   divergence = div)
    pred <- predict_forest_df(rzp_forest, test, treatment_list, control,remain_cores = 12)
    write.csv(pred, paste(folder,"rzp",div,as.character(f),".csv",sep = ""), row.names = FALSE)
  }
  
  end_time <- Sys.time()
  print(difftime(end_time,start_time,units = "mins"))
}


# Here the predictions are evaluated. Additionally we look at the treatment distribution, to see which treatments
# are assigned how often by the models.
start_time <- Sys.time()
outcomes <- c()
decile_treated <- c()
result_qini <- c()
for(model in c("random_forest","cts","sma rf","causal_forest","rzp")){
  if(sum(model == c("random_forest")) > 0){
    for(c in c("frac","max")){
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
                                        rep(paste(model,"_",c,sep = ""),11*length(treatment_list),11)))
          result_qini <- rbind(result_qini,cbind(qini_curve(pred,control,treatment_list),
                                                 paste(model,"_",c,sep = "")))
        }
      }
    }
  } else if(model == "rzp"){
    for(div in c("EucDistance","binary_KL_divergence")){
      for(f in 1:n_predictions){
        pred <- read.csv(paste(folder,model,div,as.character(f),".csv",sep = ""))
        outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),
                                     paste(model,div,sep = "_")))
        decile_treated <- rbind(decile_treated,
                                cbind(decile_perc_treated(pred,treatment_list),
                                      rep(paste(model,div,sep = "_"),11*length(treatment_list))))
        result_qini <- rbind(result_qini,setNames(cbind(qini_curve(pred,control,treatment_list),
                                                        paste(model,div,sep = "_")),names(result_qini)))    
      }
    }
  } 
  else{
    for(f in 1:n_predictions){
      pred <- read.csv(paste(folder,model,as.character(f),".csv",sep = ""))
      outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),model))
      decile_treated <- rbind(decile_treated,
                              cbind(decile_perc_treated(pred,treatment_list),
                                    rep(model,11*length(treatment_list))))
      result_qini <- rbind(result_qini,setNames(cbind(qini_curve(pred,control,treatment_list),
                                                      model),names(result_qini)))    
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
result_qini[,2] <- as.numeric(result_qini[,2])
result_qini[,1] <- as.numeric(result_qini[,1])
result_qini$model <- as.character(result_qini$model)
decile_treated_df$Model <- as.character(decile_treated_df$Model)
result_qini$percentile <- result_qini$percentile*100
outcome_df[outcome_df$Model == "cts","Model"] <- "CTS"
outcome_df[outcome_df$Model == "causal_forest","Model"] <- "Causal Forest"
outcome_df[outcome_df$Model == "random_forest_frac","Model"] <- "DOM"
outcome_df[outcome_df$Model == "sma rf","Model"] <- "SMA"
outcome_df[outcome_df$Model == "rzp_EucDistance","Model"] <- "Rzp-ED"
result_qini[result_qini$model == "cts","model"] <- "CTS"
result_qini[result_qini$model == "causal_forest","model"] <- "Causal Forest"
result_qini[result_qini$model == "random_forest_frac","model"] <- "DOM"
result_qini[result_qini$model == "sma rf","model"] <- "SMA"
result_qini[result_qini$model == "rzp_EucDistance","model"] <- "Rzp-ED"
decile_treated_df[decile_treated_df$Model == "cts","Model"] <- "CTS"
decile_treated_df[decile_treated_df$Model == "causal_forest","Model"] <- "Causal Forest"
decile_treated_df[decile_treated_df$Model == "random_forest_frac","Model"] <- "DOM"
decile_treated_df[decile_treated_df$Model == "sma rf","Model"] <- "SMA"
decile_treated_df[decile_treated_df$Model == "rzp_EucDistance","Model"] <- "Rzp-ED"
print(difftime(Sys.time(),start_time,units = "mins"))


new_qini <- result_qini[!(result_qini$model %in% c("random","random_forest_max","rzp_binary_KL_divergence")),]
new_qini <- new_qini[order(new_qini$model),]
colnames(new_qini) <- c("percentile","values","Model")
new_outcome <- outcome_df[!(outcome_df$Model %in% c("random","random_forest_max","rzp_binary_KL_divergence")),]
new_outcome <- new_outcome[order(new_outcome$Model),]

visualize_qini_uplift(new_qini,type = "qini")
visualize_qini_uplift(new_qini,type = "qini",errorbars = F,multiplot = F,ylabel = "Cummulated Gained Conversion")
visualize(new_outcome,ylabel = "Expected Converion Probability per Person",n_treated = decile_treated_df[!(decile_treated_df$Model %in% c("random","random_forest_max","rzp_binary_KL_divergence")),],multiplot = T)
visualize(new_outcome,ylabel = "Expected Converion Probability per Person",multiplot = F,errorbars = F)

outcome_boxplot(new_outcome,"Expected Conversion Probability per Customer")

temp_data <- outcome_df
values <- c()
percentile <- c()
model <- c()
for(f in 1:nrow(temp_data)){
  if(length(values) == 0){
    values <- temp_data[f,1:11]
    percentile <- colnames(temp_data)[1:11]
    model <- rep(temp_data[f,12],11)
  } else{
    values <- c(values,temp_data[f,1:11])
    percentile <- c(percentile, colnames(temp_data)[1:11])
    model <- c(model,rep(temp_data[f,12],11))
  }
}
temp_df <- data.frame(cbind(values,percentile,model))
rownames(temp_df) <- 1:nrow(temp_df)
colnames(temp_df) <- c("values","percentile","model")
for(c in 1:2){
  temp_df[,c] <- as.numeric(as.character(temp_df[,c]))
}
temp_df[,3] <- as.character(temp_df[,3]) 
tgc <- summarySE(data=temp_df, measurevar="values", groupvars=c("percentile","model"))

result_outcome <- tgc[,c("percentile","model","mean","sd")] 
result_outcome[,c("mean","sd")] <- round(result_outcome[,c("mean","sd")]*100,2)

mean_sd <- c()
for(x in 1:nrow(result_outcome)){
  temp_string <- paste(as.character(result_outcome[x,]$mean), " (",as.character(result_outcome[x,]$sd),")",sep = "")
  mean_sd <- c(mean_sd,temp_string)
}
result_outcome[,c("mean","sd")] <- NULL
result_outcome[,"mean_sd"] <- mean_sd
test_df <- result_outcome[order(result_outcome$percentile,decreasing = T),][order(result_outcome$model),]
write.csv(test_df,"ResultGraphs/Conversion/results_conversion.csv")



temp_df <- result_qini[result_qini$model != "random",]
tgc <- summarySE(data=temp_df, measurevar="values", groupvars=c("percentile","model"))

result_outcome <- tgc[,c("percentile","model","mean","sd")] 
result_outcome[,c("mean","sd")] <- round(result_outcome[,c("mean","sd")],2)

mean_sd <- c()
for(x in 1:nrow(result_outcome)){
  temp_string <- paste(as.character(result_outcome[x,]$mean), " (",as.character(result_outcome[x,]$sd),")",sep = "")
  mean_sd <- c(mean_sd,temp_string)
}
result_outcome[,c("mean","sd")] <- NULL
result_outcome[,"mean_sd"] <- mean_sd
test_df <- result_outcome[order(result_outcome$percentile,decreasing = T),][order(result_outcome$model),]
write.csv(test_df,"ResultGraphs/Conversion/qini_conversion.csv")
