source('Preprocessing.R')
source('DecisionTreeImplementation.R')
library('caret')

remain_cols <- c("female","age","voted","hh_size","treatment_CivicDuty","treatment_Self",
                 "treatment_Control","treatment_Hawthorne","treatment_Neighbors")
individuals <- individuals[,remain_cols]
individuals$female <- as.factor(individuals$female)
individuals$hh_size <- as.factor(individuals$hh_size)
individuals$voted <- as.numeric(individuals$voted) - 1

response <- 'voted'
control <- 'treatment_Control'


treatments <- c("treatment_CivicDuty","treatment_Self","treatment_Hawthorne","treatment_Neighbors")
treatment_list <- c("treatment_CivicDuty","treatment_Self")
test_list <- set_up_tests(individuals[,c("female","age","hh_size")],TRUE)
n_treatments <- length(treatment_list)
# idx <- createDataPartition(y = individuals[ , response], p=0.3, list = FALSE)
# 
# train <- individuals[-idx, ]
# test <- individuals[idx, ]
# 
# p_idx <- createDataPartition(y = train[ , response], p=0.2, list = FALSE)
# val <- train[p_idx,]
# train <- train[-p_idx,]
# 
# 
# 
# start_time <- Sys.time()
# test_tree <- create_node(train,0,100,treatment_list,response,control,test_list,
#                          normalize  = TRUE, l = rep(1/n_treatments,n_treatments),
#                          g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
# pruned_tree <- prune_tree(test_tree,val,treatment_list,test_list,response,control)
# pred <- predict.dt.as.df(pruned_tree, test)
# print(Sys.time()-start_time)
# 
# start_time <- Sys.time()
# forest <- build_forest(train,val,treatment_list,response,control,n_trees = 20,n_features = 3,criterion = 2,
#                        pruning = F)
# pred_forest <- predict_forest_df(forest, test)
# print(Sys.time()-start_time)


#Benchmark
size_vector <- c(50000,100000,200000)
time_vector_rzp <- c() 
time_vector_simple <- c() 
for(i in 2:4){
  treatment_list <- treatments[1:i]
  n_treatments <- length(treatment_list)
  temp_time_vec_rzp <- c()
  temp_time_vec_simple <- c()
  for(s in size_vector){
    temp_data <- individuals[1:s ,]
    idx <- createDataPartition(y = temp_data[ , response], p=0.2, list = FALSE)
    val <- temp_data[idx,]
    train <- temp_data[-idx,]
    start_time <- Sys.time()
    test_tree <- create_node(train,0,100,treatment_list,response,control,test_list,
                             normalize  = TRUE, l = rep(1/n_treatments,n_treatments),
                             g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
    temp_time_vec_rzp <- c(temp_time_vec_rzp,Sys.time()-start_time)
    start_time <- Sys.time()
    test_tree <- create_node(train,0,100,treatment_list,response,control,test_list,
                             normalize  = TRUE, l = rep(1/n_treatments,n_treatments),
                             g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments),criterion = 2)
    temp_time_vec_simple <- c(temp_time_vec_simple,Sys.time()-start_time)
  }
  time_vector_rzp <- c(time_vector_rzp,temp_time_vec_rzp)
  time_vector_simple <- c(time_vector_simple,temp_time_vec_simple)
}
rzp_df <- data.frame(cbind(time_vector_rzp[1],cbind(time_vector_rzp[2],time_vector_rzp[3])))
simple_df <- data.frame(cbind(time_vector_simple[1],cbind(time_vector_simple[2],time_vector_simple[3])))

