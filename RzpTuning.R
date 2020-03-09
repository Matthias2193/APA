# Basic script for tuning Rzp-Tree

library(caret)

source("ModelImplementations/PredictionFunctions.R")
source("ModelImplementations/RzepakowskiTree.R")
source("Evaluation Methods.R")
source('ModelImplementations/VisualizationHelper.R')


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

folds <- createFolds(email[,response], k = 5, list = TRUE, returnTrain = FALSE)

alpha_candidates <- c(0,0.5,1)
lambda_candidates <- list(c(0.5,0.5),c(0.1,0.9))
max_depth_candidates <- c(5,10,50)


for(f in 1:length(folds)){
  train <- email[-folds[[1]],]
  test <- email[folds[[1]],]
  idx <- createDataPartition(train[,response],p=0.2)[[1]]
  val <- train[idx,]
  train <- train[-idx,]
  test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                     "newbie","channel")],TRUE, max_cases = 10)
  for(divergence in c('EucDistance','binary_KL_divergence')){
    for(alpha in alpha_candidates){
      for(lambda in lambda_candidates){
        for(m_depth in max_depth_candidates){
          temp_tree <- build_tree_rzp(data = train, depth = 0, max_depth = m_depth, treatment_list = treatment_list,
                                      target = response,control = control,test_list = test_list, alpha = alpha,
                                      l = lambda, divergence = divergence,normalize = T)
          pruned_tree <- prune_tree(temp_tree,val, treatment_list, test_list, response,control)
          temp_pred <- predict.dt.as.df(pruned_tree,test,treatment_list,control,additional_info = T)
          exp_out <- expected_outcome(temp_pred,response,control,treatment_list)
          print(exp_out)
        }
      }
    }
  }
}


#Visualize Tree
temp_tree <- build_tree_rzp(data = train, depth = 0, max_depth = 4, treatment_list = treatment_list,
                            target = response,control = control,test_list = test_list,normalize = F)
pruned_tree <- prune_tree(temp_tree,val, treatment_list, test_list, response,control)
visualize_tree(pruned_tree,result_digits = 2,result_multiplicator = 100)
visualize_tree(pruned_tree,T)
