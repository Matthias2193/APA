#Test Script


library(ggplot2)
library(caret)
library(dplyr)
library(reshape2)
####################################################
# Uplift DT Rzepakowski et. al 2012
####################################################

source('DecisionTreeImplementation.R')
source('Evaluation Methods.R')
source('Separate Model Approach.R')
source('CausalTree.R')
source('Causal Forest.R')

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

idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)

train <- email[-idx, ]

test <- email[idx, ]

# Partition training data for pruning
p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)

val <- train[p_idx,]
train <- train[-p_idx,]

treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE, max_cases = 10)



raw_tree <- build_tree(train,0,100,treatment_list,response,control,test_list)

pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control)

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)


### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]


write.csv(pred, 'Predictions/simple tree spend pred bench.csv', row.names = FALSE)




# Expected Response per targeted customers
exp_outcome_simple1 <- new_expected_outcome(pred,response,control,treatment_list)
exp_inc_outcome_simple1 <- new_expected_quantile_response(response,control,treatment_list,pred)
overall_matched_simple <- sum(pred$Treatment==pred$Assignment)/nrow(pred)
perc_matched_simple <- perc_matched(pred)

temp_vec <- c(exp_outcome_simple1,exp_outcome_simple)
names(temp_vec) <- c("new","benchmark")
barplot(temp_vec)
plot(exp_inc_outcome_simple)
lines(exp_inc_outcome_simple1)














