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


treatment_list <- c("treatment_CivicDuty","treatment_Self","treatment_Hawthorne","treatment_Neighbors")
test_list <- set_up_tests(individuals[,c("female","age","hh_size")],TRUE)
#individuals_reduced <- individuals[1:15000,]
idx <- createDataPartition(y = individuals[ , response], p=0.3, list = FALSE)

train <- individuals[-idx, ]

test <- individuals[idx, ]

n_treatments <- length(treatment_list)

start_time <- Sys.time()
test_tree <- create_node(train,0,100,treatment_list,response,control,test_list,
                         normalize  = TRUE, l = rep(1/n_treatments,n_treatments),
                         g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
print(Sys.time()-start_time)