#Test Script

library(ggplot2)
library(caret)
library(plyr)
library(dplyr)
library(reshape2)
####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('DecisionTreeImplementation.R')
source('RzepakowskiTree.R')
source('Evaluation Methods.R')
source('CausalTree.R')
source('Causal Forest.R')
source('Separate Model Approach.R')
source('ContextualTreatmentSelection.R')
source('VisualizationHelper.R')
source("PredictionFunctions.R")


set.seed(1234)

#Data import
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


#email$spend <- email$spend/max(email$spend)

folds <- createFolds(email$spend, k = 5, list = TRUE, returnTrain = FALSE)

original_email <- email

for(f in 1:25){
  
  email <- original_email[sample(nrow(original_email),nrow(original_email),replace = TRUE),]
  idx <- createDataPartition(y = email[ , response], p=0.2, list = FALSE)
  train <- email[-idx, ]
  
  test <- email[idx, ]
  
  test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                     "newbie","channel")],TRUE, max_cases = 10)
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  
  val <- train[p_idx,]
  train_val <- train[-p_idx,]
  
  
  start_time <- Sys.time()
  # for(c in c("simple","max","frac")){
  #   print(c)
  #   # Single Tree
  #   raw_tree <- build_tree(train_val,0,100,treatment_list,response,control,test_list,criterion = c)
  #   pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control,criterion = c)
  #   pred <- predict.dt.as.df(pruned_tree, test)
  #   write.csv(pred, paste("Predictions/Hillstrom/tree_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  #   
  #   
  #   #Forest
  #   forest <- parallel_build_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3,
  #                                   pruning = F, criterion = c)
  #   pred <- predict_forest_df(forest,test)
  #   write.csv(pred, paste("Predictions/Hillstrom/forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  # 
  #   #Random Forest
  #   forest <- parallel_build_random_forest(train,treatment_list,response,control,n_trees = 100,n_features = 3, 
  #                                          criterion = c)
  #   pred <- predict_forest_df(forest,test)
  #   write.csv(pred, paste("Predictions/Hillstrom/random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  # }
  
  start_time <- Sys.time()
  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)
  
  # Causal Tree
  causal_pred <- causalTreePredicitons(train, test, treatment_list, response)
  
  print(difftime(Sys.time(),start_time,units = "mins"))
  
  # Separate Model Approach
  pred_sma_rf <- dt_models(train, response, "anova",treatment_list,control,test,"rf")
  
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


# for(f in 1:5){
#   train <- email[-folds[[f]], ]
#   
#   test <- email[folds[[f]], ]
#   
#   test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
#                                      "newbie","channel")],TRUE, max_cases = 10)
#   
#   causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)
#   
#   causal_forest_pred[ , "uplift_men_treatment"] <- causal_forest_pred[ , 1] - causal_forest_pred[ , 3]
#   causal_forest_pred[ , "uplift_women_treatment"] <- causal_forest_pred[ , 2] - causal_forest_pred[ , 3]
#   causal_forest_pred[ , "Treatment"] <- colnames(causal_forest_pred)[apply(causal_forest_pred[, 1:3], 1, which.max)]
#   
#   causal_forest_pred[ , "Outcome"] <- test[, response]
#   # get the actual assignment from test data
#   causal_forest_pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
#   
#   write.csv(causal_forest_pred, paste("Predictions/causal forest spend pred",as.character(f),".csv",sep = ""),
#             row.names = FALSE)
#   causal_forest_pred <- read.csv(paste("Predictions/causal forest spend pred",as.character(f),".csv",sep = ""))
#   if(sum(causal_forest_pred$men_treatment == 0 & causal_forest_pred$women_treatment == 0)>0){
#     causal_forest_pred[causal_forest_pred$uplift_men_treatment==0 & causal_forest_pred$uplift_women_treatment ==0,]$Treatment <- control
#   }
#   assign(paste("new_exp_inc_outcome_c_forest",as.character(f),sep = ""),
#          new_expected_quantile_response(response,control,treatment_list,causal_forest_pred))
# }


for(f in 1:5){
  train <- email[-folds[[f]], ]
  
  test <- email[folds[[f]], ]
  
  test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                     "newbie","channel")],TRUE, max_cases = 10)
  
  #Separate Model Approach
  pred_sma_rf <- dt_models(train, response, "anova",treatment_list,control,test,"rf")
  
  write.csv(pred_sma_rf, paste("Predictions/Hillstrom/sma rf",as.character(f),".csv",sep = ""),
            row.names = FALSE)
  
  # CTS
  cts_forest <- build_cts(response, control, treatment_list, train, 100, nrow(train), 3, 0.15, 100, parallel = TRUE,
                          remain_cores = 1)
  
  pred <- predict_forest_df(cts_forest, test)
  
  pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
  pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
  pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)
  
  pred[ , "Outcome"] <- test[, response]
  # get the actual assignment from test data
  pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)
  
  write.csv(pred, paste("Predictions/Hillstrom/cts_",as.character(f),".csv",sep = ""), row.names = FALSE)
}

start_time <- Sys.time()
folder <- "Predictions/Hillstrom/"
outcomes <- c()
for(model in c("tree","forest","random_forest","cts","sma rf")){
  if(sum(model == c("tree","forest","random_forest")) > 0){
    for(c in c("simple","max","frac")){
      for(f in 1:25){
        pred <- read.csv(paste(folder,model,"_",c,as.character(f),".csv",sep = ""))
        if(length(outcomes) == 0){
          outcomes <- c(new_expected_quantile_response(response,control,treatment_list,pred),
                        paste(model,"_",c,sep = ""))
        } else{
          outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),
                                       paste(model,"_",c,sep = "")))
        }
      }
    }
  } else{
    for(f in 1:25){
      pred <- read.csv(paste(folder,model,as.character(f),".csv",sep = ""))
      outcomes <- rbind(outcomes,c(new_expected_quantile_response(response,control,treatment_list,pred),model))
    }
  }
}
outcome_df <- data.frame(outcomes)
colnames(outcome_df) <- c(0,10,20,30,40,50,60,70,80,90,100,"Model")
rownames(outcome_df) <- 1:nrow(outcome_df)
for(c in 1:11){
  outcome_df[,c] <- as.numeric(as.character(outcome_df[,c]))
}
outcome_df[,12] <- as.character(outcome_df[,12])
print(difftime(Sys.time(),start_time,units = "mins"))

for(model in unique(outcome_df$Model)){
  temp_data <- outcome_df[outcome_df$Model == model,]
  visualize(temp_data)
}



for(model in unique(outcome_df$Model)){
  temp_matrix <- t(data.matrix(outcome_df[outcome_df$Model == model,1:11]))
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:25, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome for different Folds: ",model, sep = "")) #plot
}


for(f in 1:25){ 
  
  
  #Load original predictions
  pred <- read.csv(paste("Predictions/Hillstrom/tree_simple",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_tree_simple",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/Hillstrom/forest_simple",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_forest_simple",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/Hillstrom/random_forest_simple",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_random_forest_simple",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  #Load max predictions
  pred <- read.csv(paste("Predictions/Hillstrom/tree_max",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_tree_max",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/Hillstrom/forest_max",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_forest_max",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/Hillstrom/random_forest_max",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_random_forest_max",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  #Load frac pred
  pred <- read.csv(paste("Predictions/Hillstrom/tree_frac",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_tree_frac",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/Hillstrom/forest_frac",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_forest_frac",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/Hillstrom/random_forest_frac",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_random_forest_frac",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))

  pred_sma_rf <- read.csv(paste("Predictions/Hillstrom/sma rf",as.character(f),".csv",sep = ""))
  assign(paste("new_exp_inc_sma_rf",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred_sma_rf))

  # causal_forest_pred <- read.csv(paste("Predictions/causal forest spend pred",as.character(f),".csv",sep = ""))
  # if(sum(causal_forest_pred$men_treatment == 0 & causal_forest_pred$women_treatment == 0)>0){
  #   causal_forest_pred[causal_forest_pred$uplift_men_treatment==0 & causal_forest_pred$uplift_women_treatment ==0,]$Treatment <- control
  # }
  # assign(paste("new_exp_inc_outcome_c_forest",as.character(f),sep = ""),
  #        new_expected_quantile_response(response,control,treatment_list,causal_forest_pred))
  # 
  pred <- read.csv(paste("Predictions/Hillstrom/cts_",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_cts",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
}  

  
for(c in c("simple","max","frac")){
  temp_matrix <- matrix(cbind(eval(as.name(paste("new_exp_inc_tree_",c,"1",sep = ""))),
                              eval(as.name(paste("new_exp_inc_tree_",c,"2",sep = ""))),
                              eval(as.name(paste("new_exp_inc_tree_",c,"3",sep = ""))),
                              eval(as.name(paste("new_exp_inc_tree_",c,"4",sep = ""))),
                              eval(as.name(paste("new_exp_inc_tree_",c,"5",sep = "")))),nrow = 11, ncol = 5)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome for different Folds: Tree ",c, sep = "")) #plot
  legend("topleft", legend = 1:5, col=1:5, pch=1)
  
  temp_matrix <- matrix(cbind(eval(as.name(paste("new_exp_inc_forest_",c,"1",sep = ""))),
                              eval(as.name(paste("new_exp_inc_forest_",c,"2",sep = ""))),
                              eval(as.name(paste("new_exp_inc_forest_",c,"3",sep = ""))),
                              eval(as.name(paste("new_exp_inc_forest_",c,"4",sep = ""))),
                              eval(as.name(paste("new_exp_inc_forest_",c,"5",sep = "")))),nrow = 11, ncol = 5)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome for different Folds: Forest ",c, sep = "")) #plot
  legend("topleft", legend = 1:5, col=1:5, pch=1)
  
  temp_matrix <- matrix(cbind(eval(as.name(paste("new_exp_inc_random_forest_",c,"1",sep = ""))),
                              eval(as.name(paste("new_exp_inc_random_forest_",c,"2",sep = ""))),
                              eval(as.name(paste("new_exp_inc_random_forest_",c,"3",sep = ""))),
                              eval(as.name(paste("new_exp_inc_random_forest_",c,"4",sep = ""))),
                              eval(as.name(paste("new_exp_inc_random_forest_",c,"5",sep = "")))),nrow = 11, ncol = 5)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome for different Folds: Random Forest ",c, sep = "")) #plot
  legend("topleft", legend = 1:5, col=1:5, pch=1)
}

for(c in c("simple","max","frac")){
  mean_tree <- eval(as.name(paste("new_exp_inc_tree_",c,"1",sep = ""))) +
    eval(as.name(paste("new_exp_inc_tree_",c,"2",sep = ""))) +
    eval(as.name(paste("new_exp_inc_tree_",c,"3",sep = ""))) +
    eval(as.name(paste("new_exp_inc_tree_",c,"4",sep = ""))) +
    eval(as.name(paste("new_exp_inc_tree_",c,"5",sep = ""))) 
  mean_tree <- mean_tree/5
  
  mean_forest <- eval(as.name(paste("new_exp_inc_forest_",c,"1",sep = ""))) +
    eval(as.name(paste("new_exp_inc_forest_",c,"2",sep = ""))) +
    eval(as.name(paste("new_exp_inc_forest_",c,"3",sep = ""))) +
    eval(as.name(paste("new_exp_inc_forest_",c,"4",sep = ""))) +
    eval(as.name(paste("new_exp_inc_forest_",c,"5",sep = ""))) 
  mean_forest <- mean_forest/5
  
  mean_random_forest <- eval(as.name(paste("new_exp_inc_random_forest_",c,"1",sep = ""))) +
    eval(as.name(paste("new_exp_inc_random_forest_",c,"2",sep = ""))) +
    eval(as.name(paste("new_exp_inc_random_forest_",c,"3",sep = ""))) +
    eval(as.name(paste("new_exp_inc_random_forest_",c,"4",sep = ""))) +
    eval(as.name(paste("new_exp_inc_random_forest_",c,"5",sep = ""))) 
  mean_random_forest <- mean_random_forest/5
  
  # mean_causal_forest <- eval(as.name(paste("new_exp_inc_outcome_c_forest","1",sep = ""))) +
  #   eval(as.name(paste("new_exp_inc_outcome_c_forest","2",sep = ""))) +
  #   eval(as.name(paste("new_exp_inc_outcome_c_forest","3",sep = ""))) +
  #   eval(as.name(paste("new_exp_inc_outcome_c_forest","4",sep = ""))) +
  #   eval(as.name(paste("new_exp_inc_outcome_c_forest","5",sep = ""))) 
  # mean_causal_forest <- mean_causal_forest/5
  # 
  mean_sma_rf <- eval(as.name(paste("new_exp_inc_sma_rf","1",sep = ""))) +
    eval(as.name(paste("new_exp_inc_sma_rf","2",sep = ""))) +
    eval(as.name(paste("new_exp_inc_sma_rf","3",sep = ""))) +
    eval(as.name(paste("new_exp_inc_sma_rf","4",sep = ""))) +
    eval(as.name(paste("new_exp_inc_sma_rf","5",sep = "")))
  mean_sma_rf <- mean_sma_rf/5
  
  mean_cts <- eval(as.name(paste("exp_inc_cts","1",sep = ""))) +
    eval(as.name(paste("exp_inc_cts","2",sep = ""))) +
    eval(as.name(paste("exp_inc_cts","3",sep = ""))) +
    eval(as.name(paste("exp_inc_cts","4",sep = ""))) +
    eval(as.name(paste("exp_inc_cts","5",sep = "")))
  mean_cts <- mean_cts/5
  
  # temp_matrix <- matrix(cbind(mean_tree,mean_forest,mean_random_forest),
  #                       nrow=11,ncol=3)
  # matplot(temp_matrix, type = c("b"),pch=1,col = 1:3, ylab = "Expected Outcome per Person",
  #         xlab = "Percent Treated", main = paste("Expected Outcome Comparison ",c, sep = "")) #plot
  # legend("bottomright", legend = c("Tree","Forest","Random Forest"), col=1:3, pch=1,
  #        cex = 0.75)
  
  temp_matrix <- matrix(cbind(mean_tree,mean_forest,mean_random_forest,mean_cts,mean_sma_rf),
                        nrow=11,ncol=5)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome Comparison ",c, sep = "")) #plot
  legend("bottomright", legend = c("Tree","Forest","Random Forest","CTS","SMA-RF"), col=1:5, pch=1,
         cex = 0.75)
}

for(m in c("tree","forest","random_forest")){
  for(c in c("simple","max","frac")){
    mean_old <- eval(as.name(paste("exp_inc_",m,"_",c,"1",sep = ""))) +
      eval(as.name(paste("exp_inc_",m,"_",c,"2",sep = ""))) +
      eval(as.name(paste("exp_inc_",m,"_",c,"3",sep = ""))) +
      eval(as.name(paste("exp_inc_",m,"_",c,"4",sep = ""))) +
      eval(as.name(paste("exp_inc_",m,"_",c,"5",sep = ""))) 
    mean_old <- mean_old/5
    
    mean_new <- eval(as.name(paste("new_exp_inc_",m,"_",c,"1",sep = ""))) +
      eval(as.name(paste("new_exp_inc_",m,"_",c,"2",sep = ""))) +
      eval(as.name(paste("new_exp_inc_",m,"_",c,"3",sep = ""))) +
      eval(as.name(paste("new_exp_inc_",m,"_",c,"4",sep = ""))) +
      eval(as.name(paste("new_exp_inc_",m,"_",c,"5",sep = ""))) 
    mean_new <- mean_new/5
    
    temp_matrix <- matrix(cbind(mean_old,mean_new), nrow=11,ncol=2)
    matplot(temp_matrix, type = c("b"),pch=1,col = 1:2, ylab = "Expected Outcome per Person",
            xlab = "Percent Treated", main = paste("Expected Outcome Comparison",m,c, sep = " ")) #plot
    legend("bottomright", legend = c("Old", "New"), col=1:2, pch=1,
           cex = 0.75)
  }  
}

  
  
  
  
  
  
temp_matrix <- matrix(cbind(exp_inc_tree_frac1,exp_inc_tree_frac2,exp_inc_tree_frac3,exp_inc_tree_frac4,
                            exp_inc_tree_frac5),nrow = 11, ncol = 5)
matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
        xlab = "Percent Treated", main = "Expected Outcome for different Folds") #plot
legend("topleft", legend = 1:5, col=1:5, pch=1)

temp_matrix <- matrix(cbind(exp_inc_tree_max1,exp_inc_tree_max2,exp_inc_tree_max3,exp_inc_tree_max4,
                            exp_inc_tree_max5),nrow = 11, ncol = 5)
matplot(temp_matrix, type = c("b"),pch=1,col = 1:4) #plot
legend("topleft", legend = 1:5, col=1:5, pch=1)

temp_matrix <- matrix(cbind(exp_inc_tree_simple1,exp_inc_tree_simple2,exp_inc_tree_simple3,exp_inc_tree_simple4,
                            exp_inc_tree_simple5),nrow = 11, ncol = 5)
matplot(temp_matrix, type = c("b"),pch=1,col = 1:4) #plot
legend("topleft", legend = 1:5, col=1:5, pch=1)





temp_vec <- rep(seq(0,1,0.1),8)
name_vec <- c(rep("Simple Forest",11),rep("Simple Random Forest",11),
              rep("Max Tree",11), rep("Max Forest",11),rep("Max Random Forest",11),
              rep("Frac Tree",11),rep("Frac Forest",11),rep("Frac Random Forest",11))
temp_results <- cbind(c(exp_inc_forest_simple,exp_inc_random_forest_simple,
                        exp_inc_tree_max,exp_inc_forest_max,exp_inc_random_forest_max,
                        exp_inc_tree_frac,exp_inc_forest_frac,exp_inc_random_forest_frac))
temp_df <- data.frame(cbind(temp_vec,temp_results,name_vec))
colnames(temp_df) <- c("Percentile","Expected_Outcome","Model")
temp_df$Expected_Outcome <- as.numeric(as.character(temp_df$Expected_Outcome))

p <-  melt(temp_df, id.vars = c("Percentile","Model")) %>% ggplot(aes(x = Percentile)) +
  geom_line(aes(y = value, group= Model, color= Model), size=0.5 ) +
  labs(
    color="Base Learner",
    title = "Model Comparison",
    y = "Avg. Spending per Person",
    x ="Amount of Treated"
  ) +
  scale_colour_brewer(palette = "Dark2") +
  theme_light()
p

assign("testassign",1)

















#Benchmark old prediction vs new prediction

start_time <- Sys.time()
forest_pred <- parallel_predict_forest_df(forest, test)
forest_time1 <- difftime(Sys.time(),start_time,units = "secs")
start_time <- Sys.time()
forest_pred1 <- parallel_predict_forest_average_apply(forest,test)
forest_time2 <- difftime(Sys.time(),start_time,units = "secs")

if(sum(forest_pred==forest_pred1)==nrow(test)*3){
  print("New and old predictions are identical.")
  print(paste(c("The old time was:",forest_time1,"and the new time was:",forest_time2),collapse = " "))
  print(paste(c("The difference was:", as.double(forest_time1)-as.double(forest_time2),"seconds"),collapse = " "))
  paste(c("The new method takes ",round(as.double(forest_time2)/as.double(forest_time1),4)*100,
          "% of the time."),collapse = "")
}


start_time <- Sys.time() 
tree_pred1 <- predict.dt.as.df(pruned_tree,test)
tree_time1 <- difftime(Sys.time(),start_time,units = "secs")
start_time <- Sys.time()
tree_pred2 <- predict.dt.as.df_apply(pruned_tree,test)
tree_time2 <- difftime(Sys.time(),start_time,units = "secs")

if(sum(tree_pred1==tree_pred2)==nrow(test)*3){
  print("New and old predictions are identical.")
  print(paste(c("The old time was:",tree_time1,"and the new time was:",tree_time2),collapse = " "))
  print(paste(c("The difference was:", as.double(tree_time1)-as.double(tree_time2),"seconds"),collapse = " "))
  paste(c("The new method takes ",round(as.double(tree_time2)/as.double(tree_time1),4)*100,
          "% of the time."),collapse = "")
}











repeat_fct <- function(){
  if(depth == max_depth){
    return(final_node(data,treatment_list,target,control))
  }
  #Create current node
  if(random){
    retain_cols <- c(treatment_list,control,target)
    sample_cols <- setdiff(colnames(data),retain_cols)
    temp_cols <- sample(sample_cols,n_features,replace = F)
    chosen_cols <- c(temp_cols,retain_cols)
    test_list<- set_up_tests(data[,temp_cols],TRUE)
  }
  
  node <- list()
  #Select split with maximum gain
  temp_split <- select_split(test_list = test_list, treatment = treatment_list, control, target,data,criterion)
  #Return a leaf, if there is no split with gain > 0
  if(temp_split == -1){
    return(final_node(data,treatment_list,target,control))
  }
  #Construct current node
  node[['type']] = 'node'
  if(depth == 0){
    node[['type']] <- 'root'
  }
  #Number of training samples in current node
  node[['n_samples']] <- nrow(data)
  #The estimated effects for an observation in the current node, used for pruning
  treatment_names <- c(treatment_list,control)
  effects <- c()
  for(t in treatment_names){
    effects <- c(effects,mean(data[data[t]==1,target]))
  }
  names(effects) <- treatment_names
  node[['results']] <- effects
  #The current split
  node[['split']] <- temp_split
  return(temp_split)
}