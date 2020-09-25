# HU-Data Evaluation
# This script runs with data which is not public. Please see Hillstrom Evaluation Spend.R for an example.

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

# Preprocessing----
if (!file.exists("Data/hu-data.csv")) {
  hu_data <- read.csv("Data/explore_mt.csv", sep = ";")
  tempfunction <- function(x) {
    templist <- strsplit(x, ",")
    new_string <- ""
    r <- 1
    for (s in templist[[1]]) {
      if (r == 2) {
        new_string <- paste(new_string, s, sep = ".")
      } else {
        new_string <- paste(new_string, s, sep = "")
      }
      r <- r + 1
    }
    return(new_string)
  }

  hu_data$checkoutAmount <- as.numeric(lapply(as.character(hu_data$checkoutAmount), tempfunction))
  hu_data$Number.of.seconds.between.last.and.previous.views <-
    as.numeric(lapply(as.character(hu_data$Number.of.seconds.between.last.and.previous.views), tempfunction))
  for (x in grep("^log.of", colnames(hu_data))) {
    hu_data[, x] <- as.numeric(lapply(as.character(hu_data[, x]), tempfunction))
  }
  for (x in colnames(hu_data[, 16:155])) {
    if (is.na(mean(hu_data[[x]]))) {
      print(x)
    }
  }

  remove_names <- c("X", "epochSecond", "converted", "confirmed", "aborted", "dropOff")
  for (name in setdiff(colnames(hu_data)[-(1:15)], "DeviceCategory")) {
    if (var(hu_data[, name]) == 0) {
      remove_names <- c(remove_names, name)
    }
  }

  hu_data <- hu_data[, setdiff(colnames(hu_data), remove_names)]

  tmp_data <- hu_data[, -(1:8)]
  tmp_data$DeviceCategory <- NULL
  tmp <- cor(tmp_data)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  tmp_data <- tmp_data[, !apply(tmp, 2, function(x) any(x > 0.9))]
  new_hu_data <- cbind(hu_data[, 7:8], tmp_data)
  new_hu_data$DeviceCategory <- hu_data$DeviceCategory

  control <- trainControl(method = "repeatedcv", number = 5, repeats = 1)
  # train the model
  model <- train(checkoutAmount ~ ., data = new_hu_data[, -c(1)], method = "gbm", trControl = control)
  # estimate variable importance
  importance <- varImp(model, scale = FALSE)

  importance_df <- importance$importance
  importance_df$Temp <- 1
  importance_df <- importance_df[importance_df$Overall > 0, ]
  importance_df$Temp <- NULL
  importance$importance <- importance_df
  plot(importance)
  new_hu_data <- hu_data[, c(colnames(new_hu_data)[1:2], rownames(importance$importance)[1:26], "DeviceCategory")]

  for (x in levels(new_hu_data$multi_treat)) {
    new_hu_data[x] <- ifelse(new_hu_data$multi_treat == x, 1, 0)
  }

  write.csv(new_hu_data, "Data/hu-data.csv", row.names = FALSE)
  new_hu_data <- read.csv("Data/hu-data.csv")
} else {
  new_hu_data <- read.csv("Data/hu-data.csv")
}

# Model Building----
response <- "checkoutAmount"
control <- "X0"


treatment_list <- levels(new_hu_data$multi_treat)[2:7]
n_treatments <- length(treatment_list)
new_hu_data$multi_treat <- NULL
feature_list <- setdiff(colnames(new_hu_data), c(treatment_list, control, response))

# Create and save the bootrap samples and train test splits. This is done so, if we want to change something
# on one model we can retrain and test it on the sample bootstrap samples in order for fair comparison with
# the other models
if (!file.exists("bootstrap_hu.csv")) {
  bootstrap_idx <- c()
  for (f in 1:n_predictions) {
    bootstrap_idx <- cbind(bootstrap_idx, sample(nrow(new_hu_data), nrow(new_hu_data), replace = TRUE))
  }
  bootstrap_df <- data.frame(bootstrap_idx)
  write.csv(bootstrap_idx, "bootstrap_hu.csv")
} else {
  bootstrap_df <- read.csv("bootstrap_hu.csv")
}
if (!file.exists("test_hu.csv")) {
  test_idx <- c()
  for (f in 1:n_predictions) {
    hu_data <- new_hu_data[bootstrap_df[, f], ]
    test_idx <- cbind(test_idx, createDataPartition(y = hu_data[, response], p = 0.2, list = FALSE))
  }
  test_df <- data.frame(test_idx)
  write.csv(test_idx, "test_hu.csv")
} else {
  test_df <- read.csv("test_hu.csv")
}

folder <- "Predictions/HU-Data/"

for (f in 1:n_predictions) {
  hu_data <- new_hu_data[bootstrap_df[, f], ]
  train <- hu_data[-test_df[, f], ]

  test <- hu_data[test_df[, f], ]

  start_time <- Sys.time()
  for (c in c("max", "frac")) {
    print(c)
    # Random Forest
    forest <- dom_train(train, treatment_list, response, control,
      n_trees = 10, n_features = 5,
      criterion = c, min_split = 100, remain_cores = remain_cores
    )
    pred <- predict_forest_df(forest, test, treatment_list, control, remain_cores = remain_cores)
    write.csv(pred, paste(folder, "DOM", c, as.character(f), ".csv", sep = ""), row.names = FALSE)
  }

  # Causal Forest
  causal_forest_pred <- causalForestPredicitons(train, test, treatment_list, response, control,
    ntree = 10,
    s_rule = "TOT", s_true = T
  )
  write.csv(causal_forest_pred, paste(folder, "Causal Forest", as.character(f), ".csv", sep = ""),
    row.names = FALSE
  )

  # Separate Model Approach
  pred <- dt_models(train, response, "anova", treatment_list, control, test, "rf")
  write.csv(pred, paste(folder, "SMA", as.character(f), ".csv", sep = ""),
    row.names = FALSE
  )

  # CTS
  cts_forest <- cts_train(response, control, treatment_list, train, 10, nrow(train), 5, 2, 10,
    parallel = TRUE,
    remain_cores = remain_cores
  )
  pred <- predict_forest_df(cts_forest, test, treatment_list, control, remain_cores = remain_cores)
  write.csv(pred, paste(folder, "CTS", as.character(f), ".csv", sep = ""), row.names = FALSE)
  end_time <- Sys.time()
  print(difftime(end_time, start_time, units = "mins"))
}



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
visualize_qini_uplift(new_qini, type = "qini", errorbars = F, multiplot = F, ylabel = "Cumulative Gained  Checkout Amount")
visualize(new_outcome, ylabel = "Expected Checkout Amount per Person", n_treated = decile_treated_df[!(decile_treated_df$Model %in% c("random_forest_max")), ], multiplot = T)
visualize(new_outcome, ylabel = "Expected Checkout Amount per Person", multiplot = F, errorbars = F)
outcome_boxplot(new_outcome[, 2:12], "Expected Checkout Amount per Customer")
