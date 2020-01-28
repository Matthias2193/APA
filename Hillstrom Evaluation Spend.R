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
source('RzepakowskiTree.R')
source("ContextualTreatmentSelection.R")
set.seed(1234)

#Data import
email <- read.csv('Data/Email.csv')

email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$conversion <- email$visit <- email$segment <- NULL

response <- 'spend'
control <- 'control'

# Split into test and train data
idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)

train <- email[-idx, ]

#train$spend <- train$spend/max(train$spend)

test <- email[idx, ]

# Partition training data for pruning
p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)

val <- train[p_idx,]
train_val <- train[-p_idx,]

treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train[,c("recency","history_segment","history","mens","womens","zip_code",
                                   "newbie","channel")],TRUE,max_cases = 20)


### Causal Forest

causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response)

causal_forest_pred[ , "uplift_men_treatment"] <- causal_forest_pred[ , 1] - causal_forest_pred[ , 3]
causal_forest_pred[ , "uplift_women_treatment"] <- causal_forest_pred[ , 2] - causal_forest_pred[ , 3]
causal_forest_pred[ , "Treatment"] <- colnames(causal_forest_pred)[apply(causal_forest_pred[, 1:3], 1, which.max)]

causal_forest_pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
causal_forest_pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(causal_forest_pred, "Predictions/scaled causal forest spend pred.csv", row.names = FALSE)


### Simple Criterion
raw_tree <- build_tree(train_val,0,100,treatment_list,response,control,test_list,criterion = "simple")

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

write.csv(pred, 'Predictions/scaled simple tree spend pred old.csv', row.names = FALSE)

### Forest
forest <- parallel_build_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3, pruning = F)

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict_forest_df(forest, test)

### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(pred, 'Predictions/scaled simple forest spend old.csv', row.names = FALSE)


forest <- parallel_build_random_forest(train,val,treatment_list,response,control,n_trees = 100,n_features = 3, pruning = F)

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict_forest_df(forest, test)

### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(pred, 'Predictions/scaled random forest spend old.csv', row.names = FALSE)


# CTS
cts_forest <- build_cts(response, control, treatment_list, train, 10, nrow(train), 3, 0.15, 100, parallel = TRUE,
                        remain_cores = 1)

pred <- predict_forest_df(cts_forest, test)

### Results Preparation to bring into equal format
# Calculate Uplift for each T
pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

pred[ , "Outcome"] <- test[, response]
# get the actual assignment from test data
pred[ , "Assignment"] <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]

write.csv(pred, 'Predictions/cts spend.csv', row.names = FALSE)



### Separate Model Approach
#################
## Decision Tree
pred_sma_dt <- dt_models(train, response, "anova",treatment_list,control,test,"dt")
#################
## Random Forest
pred_sma_rf <- dt_models(train, response, "anova",treatment_list,control,test,"rf")


# Evaluate pre saved Causal Tree results
causal_forest_pred <- read.csv("Predictions/scaled causal forest spend pred.csv")

exp_inc_outcome_c_forest <- new_expected_quantile_response(response,control,treatment_list,causal_forest_pred)


pred <- read.csv('Predictions/scaled simple tree spend pred old.csv')

# Expected Response per targeted customers
exp_inc_outcome_simple <- new_expected_quantile_response(response,control,treatment_list,pred)


pred <- read.csv('Predictions/scaled simple forest spend old.csv')

# Expected Response per targeted customers
exp_inc_outcome_simple_forest <- new_expected_quantile_response(response,control,treatment_list,pred)


pred <- read.csv('Predictions/scaled random forest spend old.csv')

# Expected Response per targeted customers
exp_inc_outcome_random_forest <- new_expected_quantile_response(response,control,treatment_list,pred)



exp_inc_outcome_sma_dt <- new_expected_quantile_response(response,control,treatment_list,pred_sma_dt)



exp_inc_outcome_sma_rf <- new_expected_quantile_response(response,control,treatment_list,pred_sma_rf)


naive_predictions <- data.frame(cbind(rep("men_treatment",nrow(test)),sample(nrow(test),replace = T)))
colnames(naive_predictions) <- c("Treatment","random_uplift")
naive_predictions$Assignment <- colnames(test)[apply(test[, 10:12], 1, which.max) + 9]
naive_predictions$Outcome <- test[,response]
exp_inc_outcome_naive_men <- new_expected_quantile_response(response,control,treatment_list,naive_predictions)



temp_vec <- rep(seq(0,1,0.1),7)
name_vec <- c(rep("C_Forest",11),rep("Simple_Tree",11),rep("Simple_Forest",11),rep("Naive Men",11),
              rep("Random Forest",11),rep("SMA-DT",11),rep("SMA-RF",11))
temp_results <- cbind(c(exp_inc_outcome_c_forest,exp_inc_outcome_simple,exp_inc_outcome_simple_forest,exp_inc_outcome_naive_men,exp_inc_outcome_random_forest,exp_inc_outcome_sma_dt,exp_inc_outcome_sma_rf))
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


