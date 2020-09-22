# HU-Data Evaluation
# This script runs with data which is not public. Please see Hillstrom Evaluation Spend.R for an example.

library(caret)


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
remain_cores <- 2

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
  new_hu_data <- read.csv("Data/hu-data.csv")
} else{
  new_hu_data <- read.csv("Data/hu-data.csv")
}

# Model Building----
response <- 'checkoutAmount'
control <- 'X0'


treatment_list <- levels(new_hu_data$multi_treat)[2:7]
n_treatments <- length(treatment_list)
new_hu_data$multi_treat <- NULL
feature_list <- setdiff(colnames(new_hu_data),c(treatment_list,control,response))

#Create and save the bootrap samples and train test splits. This is done so, if we want to change something
#on one model we can retrain and test it on the sample bootstrap samples in order for fair comparison with
#the other models
if(!file.exists("bootstrap_hu.csv")){
  bootstrap_idx <- c()
  for(f in 1:n_predictions){
    bootstrap_idx <- cbind(bootstrap_idx,sample(nrow(new_hu_data),nrow(new_hu_data),replace = TRUE))
  }
  bootstrap_df <- data.frame(bootstrap_idx)
  write.csv(bootstrap_idx,"bootstrap_hu.csv")
} else{
  bootstrap_df <- read.csv("bootstrap_hu.csv")
}
if(!file.exists("test_hu.csv")){
  test_idx <- c()
  for(f in 1:n_predictions){
    hu_data <- new_hu_data[bootstrap_df[,f],]
    test_idx <- cbind(test_idx,createDataPartition(y = hu_data[ , response], p=0.2, list = FALSE))
  }
  test_df <- data.frame(test_idx)
  write.csv(test_idx,"test_hu.csv")
} else{
  test_df <- read.csv("test_hu.csv")
}

folder <- "Predictions/HU-Data/"

for(f in 1:n_predictions){
  hu_data <- new_hu_data[bootstrap_df[,f],]
  train <- hu_data[-test_df[,f], ]
  
  test <- hu_data[test_df[,f], ]
  
  test_list <- set_up_tests(train[,colnames(train[,feature_list])],TRUE,max_cases = 10)
  
  start_time <- Sys.time()
  for(c in c("max","frac")){
    print(c)
    #Random Forest
    forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 500,n_features = 5,
                                           criterion = c, min_split = 100, remain_cores = remain_cores)
    pred <- predict_forest_df(forest,test, treatment_list, control, remain_cores = remain_cores)
    write.csv(pred, paste(folder,"random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  }

  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control,ntree = 1000,
                                                s_rule = "TOT", s_true = T)
  write.csv(causal_forest_pred, paste(folder,"causal_forest",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  #Separate Model Approach
  pred <- dt_models(train, response, "anova",treatment_list,control,test,"rf")
  write.csv(pred, paste(folder,"sma rf",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  # CTS
  cts_forest <- build_cts(response, control, treatment_list, train, 500, nrow(train), 5, 2, 10, parallel = TRUE,
                          remain_cores = remain_cores)
  pred <- predict_forest_df(cts_forest, test, treatment_list, control, remain_cores = remain_cores)
  write.csv(pred, paste(folder,"cts",as.character(f),".csv",sep = ""), row.names = FALSE)
  end_time <- Sys.time()
  print(difftime(end_time,start_time,units = "mins"))
}



start_time <- Sys.time()
outcomes <- c()
decile_treated <- c()
result_qini <- c()
for(model in c("random_forest","cts","sma rf","causal_forest")){
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
colnames(random_df) <- c("percentile","values","model")
result_qini[,2] <- as.numeric(result_qini[,2])
result_qini[,1] <- as.numeric(result_qini[,1])
result_qini$model <- as.character(result_qini$model)
decile_treated_df$Model <- as.character(decile_treated_df$Model)
decile_treated_df$PercTreated <- decile_treated_df$PercTreated*100
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
decile_treated_df$Decile <- decile_treated_df$Decile*10
print(difftime(Sys.time(),start_time,units = "mins"))


new_qini <- result_qini[!(result_qini$model %in% c("random_forest_max")),]
new_qini <- new_qini[order(new_qini$model),]
colnames(new_qini) <- c("percentile","values","Model")
new_outcome <- outcome_df[!(outcome_df$Model %in% c("random_forest_max")),]
new_outcome <- new_outcome[order(new_outcome$Model),]

#Visualize the results
visualize_qini_uplift(new_qini,type = "qini",errorbars = F,multiplot = F,ylabel = "Cumulative Gained  Checkout Amount")
visualize(new_outcome,ylabel = "Expected Checkout Amount per Person",n_treated = decile_treated_df[!(decile_treated_df$Model %in% c("random_forest_max")),],multiplot = T)
visualize(new_outcome,ylabel = "Expected Checkout Amount per Person",multiplot = F,errorbars = F)
outcome_boxplot(new_outcome[,2:12],"Expected Checkout Amount per Customer")