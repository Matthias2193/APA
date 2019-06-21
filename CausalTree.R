# install.packages("devtools")
# library(devtools) 
# install_github("susanathey/causalTree")
 library(causalTree)
# 
# email <- read.csv('Email.csv')
# email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
# email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
# email$control <- ifelse(email$segment=='No E-Mail',1,0)
# email$segment <- NULL
# email$mens <- as.factor(email$mens)
# email$womens <- as.factor(email$womens)
# email$newbie <- as.factor(email$newbie)
# 
# treatment_list <- c('men_treatment','women_treatment')
# 
# smp_size <- floor(0.75 * nrow(email))
# set.seed(123)
# train_ind <- sample(seq_len(nrow(email)), size = smp_size)
# train <- email[train_ind, ]
# test <- email[-train_ind, ]

causalTreePredicitons <- function(train, test,treatment_list, response){
  for(t in treatment_list){
    train_data <- train
    train_data <- train_data[train_data[,setdiff(treatment_list,t)] == 0,]
    train_data_new <- train_data[,1:8]
    train_data_new[, response] <- train_data[, response]
    train_data <- train_data_new
    
    test_data <- test[,1:8]
    test_data[,t] <- test[,t]
    
    if(file.exists(paste(paste('models/tree',t,sep = '_'),'rda',sep='.'))){
      load(paste(paste('models/tree',t,sep = '_'),'rda',sep='.'))
    }
    else{
      tree <- causalTree(as.formula(paste(response, "~.")), data = train_data, treatment = train[,t], split.Rule = "CT", cv.option = "CT", split.Honest = T,
                         cv.Honest = T, split.Bucket = F, xval = 5, cp = 0, minsize = 20, propensity = 0.5)
      save(tree, file = paste(paste('models/tree',t,sep = '_'),'rda',sep='.'))
    }
    
    
    #opcp <- tree$cptable[,1][which.min(tree$cptable[,3])]
    
    #pruned_tree <- prune(tree, opcp)
    
    assign(paste('predictions',t,sep = '_'),predict(tree,newdata = test_data))
  }
  
  pred <- data.frame(cbind(predictions_men_treatment,predictions_women_treatment))
  colnames(pred) <- treatment_list
  pred$control <- 0
  pred[ , "Uplift - Mens E-Mail"] <- pred[ , 1] - pred[ , 3]
  pred[ , "Uplift - Womens E-Mail"] <- pred[ , 2] - pred[ , 3]
  
  pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]
  
  pred[ , "Outcome"] <- test[, response]
  
  pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
  
  return(pred)
}

