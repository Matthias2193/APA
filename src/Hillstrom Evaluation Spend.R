# This script evaluates several different models on the Hillstrom data set which is contained in this repository.
# The models evaluated are Contextual Treatment Selection (CTS), Causal Forest, Separate Model Approach with
# Random Forest and a custom Tree/Random Forest with two different gain functions ("Simple" and "Frac). More
# information about the custom Tree and Random Forest can be found under Algorithm Implementations/DecisionTreeImplementation.R
# The user can specify the parameter n_predictions. If it is set to one, each model is trained once on the
# original data set and then evaluated. If n_predictions is greater than 1, there will be n_predictions iterations.
# For each iteration a bootstrap sample is taken as the new data and then the models are build. After the models
# have been built and the predictions have been made, the predictions are evaluated and the results ploted.

library(caret)

source("src/Algorithm Implementations/DOM.R")
source("src/Algorithm Implementations/RzepakowskiTree.R")
source("src/Algorithm Implementations/CausalTree.R")
source("src/Algorithm Implementations/Separate Model Approach.R")
source("src/Algorithm Implementations/ContextualTreatmentSelection.R")
source("src/Algorithm Implementations/PredictionFunctions.R")

source("src/Helper Functions/VisualizationHelper.R")
source("src/Helper Functions/Evaluation Methods.R")

set.seed(1234)
n_predictions <- 2
remain_cores <- 2
# Data import and preprocessing
email <- read.csv("Data/Email.csv")

email$men_treatment <- ifelse(email$segment == "Mens E-Mail", 1, 0)
email$women_treatment <- ifelse(email$segment == "Womens E-Mail", 1, 0)
email$control <- ifelse(email$segment == "No E-Mail", 1, 0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$conversion <- email$segment <- NULL

response <- "spend"
control <- "control"
treatment_list <- c("men_treatment", "women_treatment")

original_email <- email

# Create and save the bootrap samples and train test splits. This is done so, if we want to change something
# on one model we can retrain and test it on the sample bootstrap samples in order for fair comparison with
# the other models
if (!file.exists("bootstrap.csv")) {
  bootstrap_idx <- c()
  for (f in 1:n_predictions) {
    bootstrap_idx <- cbind(bootstrap_idx, sample(nrow(original_email), nrow(original_email), replace = TRUE))
  }
  bootstrap_df <- data.frame(bootstrap_idx)
  write.csv(bootstrap_idx, "bootstrap.csv")
} else {
  bootstrap_df <- read.csv("bootstrap.csv")
}
if (!file.exists("test.csv")) {
  test_idx <- c()
  for (f in 1:n_predictions) {
    email <- original_email[bootstrap_df[, f], ]
    test_idx <- cbind(test_idx, createDataPartition(y = email[, response], p = 0.2, list = FALSE))
  }
  test_df <- data.frame(test_idx)
  write.csv(test_idx, "test.csv")
} else {
  test_df <- read.csv("test.csv")
}
folder <- "Predictions/Hillstrom/"
# The training and prediction part
for (f in 1:n_predictions) {

  # If n_predictions is > 1 as bootstrap sample is created
  if (n_predictions > 1) {
    email <- original_email[bootstrap_df[, f], ]
  }
  idx <- test_df[, f]
  train <- email[-idx, ]

  test <- email[idx, ]

  start_time <- Sys.time()
  for (c in c("frac", "max")) {
    print(c)
    # Random Forest
    forest <- dom_train(train, treatment_list, response, control,
      n_trees = 500, n_features = 3,
      criterion = c, remain_cores = remain_cores
    )
    pred <- predict_forest_df(forest, test, treatment_list, control, remain_cores = remain_cores)
    write.csv(pred, paste(folder, "DOM", c, as.character(f), ".csv", sep = ""), row.names = FALSE)
  }

  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control,
    ntree = 1000,
    s_rule = "TOT", s_true = T
  )
  write.csv(causal_forest_pred, paste(folder, "Causal Forest", as.character(f), ".csv", sep = ""),
    row.names = FALSE
  )

  # Separate Model Approach
  pred_sma_rf <- dt_models(train, response, "anova", treatment_list, control, test, "rf", mtry = 3, ntree = 500)
  write.csv(pred_sma_rf, paste(folder, "SMA", as.character(f), ".csv", sep = ""),
    row.names = FALSE
  )

  # CTS
  cts_forest <- cts_train(response, control, treatment_list, train,
    ntree = 500, nrow(train), m_try = 4,
    n_reg = 4, min_split = 10, parallel = TRUE, remain_cores = remain_cores
  )
  pred <- predict_forest_df(
    forest = cts_forest, test_data = test, treatment_list = treatment_list,
    control = control, remain_cores = remain_cores, parallel_pred = F
  )
  write.csv(pred, paste(folder, "CTS", as.character(f), ".csv", sep = ""), row.names = FALSE)
  end_time <- Sys.time()
  print(difftime(end_time, start_time, units = "mins"))
}

# Here the predictions are evaluated. Additionally we look at the treatment distribution, to see which treatments
# are assigned how often by the models.
start_time <- Sys.time()
outcomes <- c()
decile_treated <- c()
result_qini <- c()
for (model in c("DOM", "SMA", "Causal Forest", "CTS")) {
  if (model == "DOM") {
    for (c in c("frac", "max")) {
      for (f in 1:n_predictions) {
        pred <- read.csv(paste(folder, model, c, as.character(f), ".csv", sep = ""))
        if (length(outcomes) == 0) {
          outcomes <- expected_outcome_curve(response, control, treatment_list, pred, 
                                             model_name = paste(model, c, sep = " "))
          decile_treated <- decile_perc_treated(pred, treatment_list, paste(model, c, sep = " "))
          result_qini <- qini_curve(pred, control, treatment_list, model_name = paste(model, c, sep = " "))
        } else {
          outcomes <- rbind(outcomes, 
                            expected_outcome_curve(response, control, treatment_list, pred, 
                                                   model_name = paste(model, c, sep = " "))
          )
          decile_treated <- rbind(
            decile_treated, decile_perc_treated(pred, treatment_list, paste(model, c, sep = " "))
          )
          result_qini <- rbind(result_qini, 
                               qini_curve(pred, control, treatment_list, 
                                          model_name = paste(model, c, sep = " ")))
        }
      }
    }
  } 
  else {
    for (f in 1:n_predictions) {
      pred <- read.csv(paste(folder, model, as.character(f), ".csv", sep = ""))
      outcomes <- rbind(outcomes, 
                        expected_outcome_curve(response, control, treatment_list, pred, model_name = model))
      decile_treated <- rbind(
        decile_treated, decile_perc_treated(pred, treatment_list, model)
      )
      result_qini <- rbind(result_qini, qini_curve(pred, control, treatment_list, model_name = model))
    }
  }
}
print(difftime(Sys.time(), start_time, units = "mins"))

new_qini <- result_qini[!(result_qini$Model %in% c("DOM max")), ]
new_qini <- new_qini[order(new_qini$Model), ]

new_outcome <- outcomes[!(outcomes$Model %in% c("DOM max")), ]
new_outcome <- new_outcome[order(new_outcome$Model), ]

# Visualize the results
visualize_qini(new_qini, type = "qini", ylabel = "Cumulative Gained Spend")
visualize_outcome(new_outcome, ylabel = "Expected Amount Spend per Person", n_treated = decile_treated[!(decile_treated$Model %in% c("DOM max")), ], multiplot = T)
visualize_outcome(new_outcome, ylabel = "Expected Amount Spend per Person")
outcome_boxplot(new_outcome[, 2:12], "Expected Amount Spend per Customer")
