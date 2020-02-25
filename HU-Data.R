# HU-Data Evaluation
# This script runs with data which is not public. Please see Hillstrom Evaluation Spend.R for an example.

library(ggplot2)
library(caret)
library(plyr)
library(dplyr)
library(reshape2)
library(gbm)

source('ModelImplementations/DecisionTreeImplementation.R')
source('ModelImplementations/RzepakowskiTree.R')
source('Evaluation Methods.R')
source('ModelImplementations/CausalTree.R')
source('ModelImplementations/Separate Model Approach.R')
source('ModelImplementations/ContextualTreatmentSelection.R')
source('ModelImplementations/VisualizationHelper.R')
source("ModelImplementations/PredictionFunctions.R")

set.seed(1234)
#Preprocessing---- 
if(!file.exists("Data/hu-data.csv")){
  hu_data <- read.csv("Data/explore_mt.csv",sep = ";")
  tempfunction <- function(x){
    templist <- strsplit(x,",")
    new_string <- ""
    r <- 1
    for(s in templist[[1]]){
      if (r == 2) {
        new_string <- paste(new_string,s,sep = ".")
      } else{
        new_string <- paste(new_string,s,sep="")
      }
      r <- r+1
    }
    return(new_string)
  }
  
  hu_data$checkoutAmount <- as.numeric(lapply(as.character(hu_data$checkoutAmount), tempfunction))
  hu_data$Number.of.seconds.between.last.and.previous.views <- 
    as.numeric(lapply(as.character(hu_data$Number.of.seconds.between.last.and.previous.views), tempfunction))
  for (x in grep("^log.of",colnames(hu_data))) {
    hu_data[,x] <- as.numeric(lapply(as.character(hu_data[,x]), tempfunction))
  }
  for (x in colnames(hu_data[,16:155])) {
    if(is.na(mean(hu_data[[x]]))){
      print(x)
    }
  }
  
  remove_names <- c("X","epochSecond","converted","confirmed","aborted","dropOff")
  for(name in setdiff(colnames(hu_data)[-(1:15)],"DeviceCategory")){
    if(var(hu_data[,name])==0){
      remove_names <- c(remove_names,name)
    }
  }
  
  hu_data <- hu_data[,setdiff(colnames(hu_data),remove_names)]
  
  tmp_data <- hu_data[,-(1:8)]
  tmp_data$DeviceCategory <- NULL
  tmp <- cor(tmp_data)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  tmp_data <- tmp_data[,!apply(tmp,2,function(x) any(x > 0.9))]
  new_hu_data <- cbind(hu_data[,7:8],tmp_data)
  new_hu_data$DeviceCategory <- hu_data$DeviceCategory
  
  control <- trainControl(method="repeatedcv", number=5, repeats=1)
  # train the model
  model <- train(checkoutAmount~., data=new_hu_data[,-c(1)], method="gbm", trControl=control)
  # estimate variable importance
  importance <- varImp(model, scale=FALSE)
  
  importance_df <- importance$importance
  importance_df$Temp <- 1
  importance_df <- importance_df[importance_df$Overall > 0,]
  importance_df$Temp <- NULL
  importance$importance <- importance_df
  plot(importance)
  new_hu_data <- hu_data[,c(colnames(new_hu_data)[1:2],rownames(importance$importance)[1:26],"DeviceCategory")]
  
  for(x in levels(new_hu_data$multi_treat)){
    new_hu_data[x] <- ifelse(new_hu_data$multi_treat == x ,1,0)
  }
  
  write.csv(new_hu_data, "Data/hu-data.csv", row.names = FALSE)
} else{
  new_hu_data <- read.csv("Data/hu-data.csv")
}

# Model Building----
response <- 'checkoutAmount'
control <- 'X0'

for(t in c(treatment_list,control)){
  print(t)
  print(sum(new_hu_data[,t]==1))
  print(mean(new_hu_data[new_hu_data[,t]==1,response]))
}

n_predictions <- 15
treatment_list <- levels(new_hu_data$multi_treat)[2:7]
n_treatments <- length(treatment_list)
new_hu_data$multi_treat <- NULL
feature_list <- setdiff(colnames(new_hu_data),c(treatment_list,control,response))


for(f in 1:n_predictions){
  
  hu_data <- new_hu_data[sample(nrow(new_hu_data),nrow(new_hu_data),replace = TRUE),]
  idx <- createDataPartition(y = hu_data[ , response], p=0.2, list = FALSE)
  train <- hu_data[-idx, ]
  
  test <- hu_data[idx, ]
  
  test_list <- set_up_tests(train[,colnames(train[,feature_list])],TRUE,max_cases = 10)
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  
  val <- train[p_idx,]
  train_val <- train[-p_idx,]
  
  start_time <- Sys.time()
  for(c in c("simple","max","frac")){
    print(c)
    # Single Tree
    raw_tree <- build_tree(train_val,0,100,treatment_list,response,control,test_list,criterion = c)
    pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control,criterion = c)
    pred <- predict.dt.as.df(pruned_tree, test)
    write.csv(pred, paste("Predictions/HU-Data/tree_",c,as.character(f),".csv",sep = ""), row.names = FALSE)

    #Random Forest
    forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 300,n_features = 5,
                                           criterion = c)
    pred <- predict_forest_df(forest,test)
    write.csv(pred, paste("Predictions/HU-Data/random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  }

  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control)
  write.csv(causal_forest_pred, paste("Predictions/HU-Data/causal_forest",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  #Separate Model Approach
  pred <- dt_models(train, response, "anova",treatment_list,control,test,"rf")
  write.csv(pred, paste("Predictions/HU-Data/sma rf",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  # CTS
  cts_forest <- build_cts(response, control, treatment_list, train, 300, nrow(train), 5, 1, 100, parallel = TRUE,
                          remain_cores = 1)
  pred <- predict_forest_df(cts_forest, test)
  write.csv(pred, paste("Predictions/HU-Data/cts",as.character(f),".csv",sep = ""), row.names = FALSE)
  end_time <- Sys.time()
  print(difftime(end_time,start_time,units = "mins"))
}



start_time <- Sys.time()
folder <- "Predictions/HU-Data/"
outcomes <- c()
decile_treated <- c()
result_qini <- c()
result_uplift <- c()

for(model in c("random_forest","cts","sma rf","causal_forest")){
  if(sum(model == c("tree","random_forest")) > 0){
    for(c in c("max")){
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
    colnames(result_uplift) <- c("Percentile","combined","model")
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
colnames(result_uplift) <- c("percentile","values","model")
start <- mean(result_qini[result_qini$percentile == 0.0,"values"])
finish <- mean(result_qini[result_qini$percentile == 1.0,"values"])
qini_random <- seq(start,finish,by = (finish-start)/10)
random_df <- cbind(seq(0,1,by=0.1),qini_random,"random","random")
colnames(random_df) <- c("percentile","values","treatment","model")
result_qini <- rbind(result_qini,random_df)
result_qini[,2] <- as.numeric(result_qini[,2])
result_qini[,1] <- as.numeric(result_qini[,1])
print(difftime(Sys.time(),start_time,units = "mins"))


visualize_qini_uplift(result_qini,type = "qini")
visualize_qini_uplift(result_qini,type = "qini",errorbars = F,multiplot = F)
visualize(outcome_df,n_treated = decile_treated_df,multiplot = T)
visualize(outcome_df,multiplot = F,errorbars = F)
