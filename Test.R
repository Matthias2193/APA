#Test Script

library(ggplot2)
library(caret)
library(dplyr)
library(reshape2)
####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('DecisionTreeImplementation.R')
source('RzepakowskiTree.R')
source('Evaluation Methods.R')


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

idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)
folds <- createFolds(email$spend, k = 5, list = TRUE, returnTrain = FALSE)
for(f in 4:5){
  train <- email[-folds[[f]], ]
  
  test <- email[folds[[f]], ]
  
  # Partition training data for pruning
  p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)
  
  val <- train[p_idx,]
  train_val <- train[-p_idx,]
  
  
  test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                     "newbie","channel")],TRUE, max_cases = 10)
  start_time <- Sys.time()
  for(c in c("simple","frac","max")){
    print(c)
    # Single Tree
    raw_tree <- build_tree(train_val,0,100,treatment_list,response,control,test_list,criterion = c)
    pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control,criterion = c)
    
    # add to the result df the outcome, assignment and calculate uplift for each T
    pred <- predict.dt.as.df(pruned_tree, test)
    
    
    ### Results Preparation to bring into equal format
    # Calculate Uplift for each T
    pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
    pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
    pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)
    
    pred[ , "Outcome"] <- test[, response]
    # get the actual assignment from test data
    pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)
    
    write.csv(pred, paste("Predictions/tree_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
    
    
    #Forest
    forest <- parallel_build_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3,
                                    pruning = F, criterion = c)

    # add to the result df the outcome, assignment and calculate uplift for each T
    pred <- predict_forest_df(forest,test)

    ### Results Preparation to bring into equal format
    # Calculate Uplift for each T
    pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
    pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
    pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)

    pred[ , "Outcome"] <- test[, response]
    # get the actual assignment from test data
    pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)

    write.csv(pred, paste("Predictions/forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)

    #Random Forest
    forest <- parallel_build_random_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3,
                                           pruning = F, criterion = c)

    # add to the result df the outcome, assignment and calculate uplift for each T
    pred <- predict_forest_df(forest,test)

    ### Results Preparation to bring into equal format
    # Calculate Uplift for each T
    pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
    pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
    pred[ , "Treatment"] <- predictions_to_treatment(pred, treatment_list, control)

    pred[ , "Outcome"] <- test[, response]
    # get the actual assignment from test data
    pred[ , "Assignment"] <- predictions_to_treatment(test, treatment_list, control)

    write.csv(pred, paste("Predictions/random_forest_",c,as.character(f),".csv",sep = ""), row.names = FALSE)
  }
}


for(f in 1:5){  
  train <- email[-folds[[f]], ]
  
  test <- email[folds[[f]], ]
  #Load original predictions
  pred <- read.csv(paste("Predictions/tree_simple",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_tree_simple",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/forest_simple",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_forest_simple",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/random_forest_simple",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_random_forest_simple",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  #Load max predictions
  pred <- read.csv(paste("Predictions/tree_max",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_tree_max",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/forest_max",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_forest_max",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/random_forest_max",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_random_forest_max",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  #Load frac pred
  pred <- read.csv(paste("Predictions/tree_frac",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_tree_frac",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/forest_frac",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_forest_frac",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  pred <- read.csv(paste("Predictions/random_forest_frac",as.character(f),".csv",sep = ""))
  assign(paste("exp_inc_random_forest_frac",as.character(f),sep = ""),
         new_expected_quantile_response(response,control,treatment_list,pred))
  end_time <- Sys.time()
}  

  
for(c in c("simple","max","frac")){
  temp_matrix <- matrix(cbind(eval(as.name(paste("exp_inc_tree_",c,"1",sep = ""))),
                              eval(as.name(paste("exp_inc_tree_",c,"2",sep = ""))),
                              eval(as.name(paste("exp_inc_tree_",c,"3",sep = ""))),
                              eval(as.name(paste("exp_inc_tree_",c,"4",sep = ""))),
                              eval(as.name(paste("exp_inc_tree_",c,"5",sep = "")))),nrow = 11, ncol = 5)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome for different Folds: Tree ",c, sep = "")) #plot
  legend("topleft", legend = 1:5, col=1:5, pch=1)
  
  temp_matrix <- matrix(cbind(eval(as.name(paste("exp_inc_forest_",c,"1",sep = ""))),
                              eval(as.name(paste("exp_inc_forest_",c,"2",sep = ""))),
                              eval(as.name(paste("exp_inc_forest_",c,"3",sep = ""))),
                              eval(as.name(paste("exp_inc_forest_",c,"4",sep = ""))),
                              eval(as.name(paste("exp_inc_forest_",c,"5",sep = "")))),nrow = 11, ncol = 5)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome for different Folds: Forest ",c, sep = "")) #plot
  legend("topleft", legend = 1:5, col=1:5, pch=1)
  
  temp_matrix <- matrix(cbind(eval(as.name(paste("exp_inc_random_forest_",c,"1",sep = ""))),
                              eval(as.name(paste("exp_inc_random_forest_",c,"2",sep = ""))),
                              eval(as.name(paste("exp_inc_random_forest_",c,"3",sep = ""))),
                              eval(as.name(paste("exp_inc_random_forest_",c,"4",sep = ""))),
                              eval(as.name(paste("exp_inc_random_forest_",c,"5",sep = "")))),nrow = 11, ncol = 5)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:5, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome for different Folds: Random Forest ",c, sep = "")) #plot
  legend("topleft", legend = 1:5, col=1:5, pch=1)
}

for(c in c("simple","max","frac")){
  mean_tree <- eval(as.name(paste("exp_inc_tree_",c,"1",sep = ""))) +
    eval(as.name(paste("exp_inc_tree_",c,"2",sep = ""))) +
    eval(as.name(paste("exp_inc_tree_",c,"3",sep = ""))) +
    eval(as.name(paste("exp_inc_tree_",c,"4",sep = ""))) +
    eval(as.name(paste("exp_inc_tree_",c,"5",sep = ""))) 
  mean_tree <- mean_tree/5
  
  mean_forest <- eval(as.name(paste("exp_inc_forest_",c,"1",sep = ""))) +
    eval(as.name(paste("exp_inc_forest_",c,"2",sep = ""))) +
    eval(as.name(paste("exp_inc_forest_",c,"3",sep = ""))) +
    eval(as.name(paste("exp_inc_forest_",c,"4",sep = ""))) +
    eval(as.name(paste("exp_inc_forest_",c,"5",sep = ""))) 
  mean_forest <- mean_forest/5
  
  mean_random_forest <- eval(as.name(paste("exp_inc_random_forest_",c,"1",sep = ""))) +
    eval(as.name(paste("exp_inc_random_forest_",c,"2",sep = ""))) +
    eval(as.name(paste("exp_inc_random_forest_",c,"3",sep = ""))) +
    eval(as.name(paste("exp_inc_random_forest_",c,"4",sep = ""))) +
    eval(as.name(paste("exp_inc_random_forest_",c,"5",sep = ""))) 
  mean_random_forest <- mean_random_forest/5
  
  temp_matrix <- matrix(cbind(mean_tree,mean_forest,mean_random_forest),nrow=11,ncol=3)
  matplot(temp_matrix, type = c("b"),pch=1,col = 1:3, ylab = "Expected Outcome per Person",
          xlab = "Percent Treated", main = paste("Expected Outcome Comparison ",c, sep = "")) #plot
  legend("topleft", legend = c("Tree","Forest","Random Forest"), col=1:3, pch=1)
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