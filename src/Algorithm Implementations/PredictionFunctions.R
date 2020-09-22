# The predictions functions

# Prediction----
# Tree
predict.dt.as.df <- function(tree, new_data, treatment_list, control, additional_info = TRUE) {
  type_list <- sapply(new_data, class)
  names(type_list) <- colnames(new_data)
  temp_function <- function(x, node) {
    type <- node[["type"]]
    while (type != "leaf") {
      split <- node[["split"]]
      if (type_list[[names(split)]] == "factor") {
        if (x[names(split)] == split[[1]]) {
          node <- node[["left"]]
          type <- node[["type"]]
        } else {
          node <- node[["right"]]
          type <- node[["type"]]
        }
      } else {
        if (as.numeric(x[names(split)]) < split[[1]]) {
          node <- node[["left"]]
          type <- node[["type"]]
        } else {
          node <- node[["right"]]
          type <- node[["type"]]
        }
      }
    }
    return(node[["results"]])
  }
  results <- data.frame(t(apply(new_data, 1, temp_function, node = tree)))
  if (additional_info) {
    colnames(results) <- c(treatment_list, control)
    results[, "Treatment"] <- predictions_to_treatment(results, treatment_list, control)
    results[, "Assignment"] <- predictions_to_treatment(new_data, treatment_list, control)
    results[, "Outcome"] <- new_data[, response]
    for (t in treatment_list) {
      results[, paste("uplift", t, sep = "_")] <- results[t] - results[control]
    }
    return(results)
  } else {
    return(results)
  }
}


# Takes predictions as input and returns just the name of the best treatment for each prediction

predictions_to_treatment <- function(pred, treatment_list, control) {
  temp_list <- c()
  results <- colnames(pred[, c(treatment_list, control)])[apply(pred[, c(treatment_list, control)], 1, which.max)]
  # If all uplifts are 0, we assign control
  for (c in treatment_list) {
    if (length(temp_list) == 0) {
      temp_list <- pred[c]
    } else {
      temp_list <- temp_list + pred[, c]
    }
  }
  if (sum(temp_list == 0) > 0) {
    results[temp_list == 0] <- control
  }
  return(results)
}


# Forest

# Majority Vote
predict_forest_majority <- function(forest, test_data, treatment_list, control) {
  predictions <- predictions_to_treatment(
    predict.dt(forest[[1]], test_data, treatment_list, control),
    treatment_list, control
  )
  for (x in 2:length(forest)) {
    predictions <- cbind(predictions, predictions_to_treatment(predict.dt(forest[[x]], test_data, treatment_list, control)))
  }
  return(apply(data.frame(predictions), MARGIN = 1, FUN = forest_predictions_helper))
}

forest_predictions_helper <- function(preds) {
  temp_names <- unique(unlist(preds))
  if (length(temp_names) == 1) {
    return(temp_names)
  } else {
    return(names(which.max(table(unlist(preds)))))
  }
}


# Average
# Sequentially
predict_forest_average <- function(forest, test_data, treatment_list, control, additional_info = TRUE) {
  predictions <- list()
  for (x in 1:length(forest)) {
    predictions[[x]] <- predict.dt.as.df(forest[[x]], test_data, additional_info = FALSE)
  }
  final_predictions <- predictions[[1]]
  for (x in 2:length(predictions)) {
    final_predictions <- final_predictions + predictions[[x]]
  }
  final_predictions <- final_predictions / length(predictions)
  if (additional_info) {
    colnames(final_predictions) <- c(treatment_list, control)
    final_predictions[, "Treatment"] <- predictions_to_treatment(final_predictions, treatment_list, control)
    final_predictions[, "Assignment"] <- predictions_to_treatment(test_data, treatment_list, control)
    final_predictions[, "Outcome"] <- test_data[, response]
    for (t in treatment_list) {
      final_predictions[, paste("uplift", t, sep = "_")] <- final_predictions[t] - final_predictions[control]
    }
    return(final_predictions)
  } else {
    return(final_predictions)
  }
}

# Parallel
parallel_predict_forest_average <- function(forest, test_data, treatment_list, control, remain_cores = 1, additional_info = TRUE) {
  predictions <- list()
  numCores <- detectCores()
  cl <- makePSOCKcluster(numCores - remain_cores)
  registerDoParallel(cl)
  predictions <- foreach(x = 1:length(forest)) %dopar% {
    source("src/Algorithm Implementations/DOM.R")
    tree <- forest[[x]]
    new_data <- test_data
    type_list <- sapply(new_data, class)
    names(type_list) <- colnames(new_data)
    temp_function <- function(x, node) {
      type <- node[["type"]]
      while (type != "leaf") {
        split <- node[["split"]]
        if (type_list[[names(split)]] == "factor") {
          if (x[names(split)] == split[[1]]) {
            node <- node[["left"]]
            type <- node[["type"]]
          } else {
            node <- node[["right"]]
            type <- node[["type"]]
          }
        } else {
          if (as.numeric(x[names(split)]) < split[[1]]) {
            node <- node[["left"]]
            type <- node[["type"]]
          } else {
            node <- node[["right"]]
            type <- node[["type"]]
          }
        }
      }
      return(node[["results"]])
    }
    results <- data.frame(t(apply(new_data, 1, temp_function, node = tree)))
    return(results)
  }
  stopCluster(cl)
  final_predictions <- predictions[[1]]
  for (x in 2:length(predictions)) {
    final_predictions <- final_predictions + predictions[[x]]
  }
  final_predictions <- final_predictions / length(predictions)
  if (additional_info) {
    colnames(final_predictions) <- c(treatment_list, control)
    final_predictions[, "Treatment"] <- predictions_to_treatment(final_predictions, treatment_list, control)
    final_predictions[, "Assignment"] <- predictions_to_treatment(test_data, treatment_list, control)
    final_predictions[, "Outcome"] <- test_data[, response]
    for (t in treatment_list) {
      final_predictions[, paste("uplift", t, sep = "_")] <- final_predictions[t] - final_predictions[control]
    }
    return(final_predictions)
  } else {
    return(final_predictions)
  }
}





predict_forest_df <- function(forest, test_data, treatment_list, control, parallel_pred = TRUE, remain_cores = 1, additiona_info = TRUE) {
  if (parallel_pred) {
    return(parallel_predict_forest_average(forest, test_data, treatment_list, control, remain_cores, additional_info = additiona_info))
  } else {
    return(predict_forest_average(forest, test_data, treatment_list, control, additional_info = additiona_info))
  }
}
