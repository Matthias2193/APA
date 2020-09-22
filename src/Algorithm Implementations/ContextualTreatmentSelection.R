library(dplyr)
library(foreach)
library(doParallel)
source("src/Algorithm Implementations/DOM.R")

# The main method to build a cts forest.
# There are 5 parameters which can be tuned to impact the mdoel performance:
# ntree: The number of trees in the forest
# B: The number of samples in each bootstrap sample. B should be smaller or equal to the number of samples in the original data set.
# m_try: The number of covariates selected at each split.
# min_split: The minimum number of observations for each treatment.
# If parallel = TRUE the process will be parallalized using the number of cores available - remain_cores.
# For further information about the algorithm please read the paper.
build_cts <- function(response, control, treatments, data, ntree, B, m_try, n_reg, min_split, parallel = TRUE,
                      remain_cores = 1) {
  if (parallel) {
    # Set up the parallelization
    numCores <- detectCores()
    cl <- makePSOCKcluster(numCores - remain_cores)
    registerDoParallel(cl)
    trees <- foreach(x = 1:ntree) %dopar% {
      source("src/Algorithm Implementations/ContextualTreatmentSelection.R")
      set.seed(x)
      # In CTS we sample according to the treatment distribution in the train set instead of completely random.
      for (t in c(treatments, control)) {
        if (t == treatments[1]) {
          temp_train_data <- sample_n(data[data[, t] == 1, ], round(B * (sum(data[, t] == 1) / nrow(data)), 0), replace = TRUE)
        } else {
          temp_train_data <- rbind(temp_train_data, sample_n(data[data[, t] == 1, ],
            round(B * (sum(data[, t] == 1) / nrow(data)), 0),
            replace = TRUE
          ))
        }
      }
      temp_train_data <- na.omit(temp_train_data)
      # Once the new training data is sampled the tree is build.
      return(build_cts_tree(response, control, treatments, temp_train_data, m_try, n_reg,
        min_split,
        parent_predictions = NA
      ))
    }
    stopCluster(cl)
    return(trees)
  }
  else {
    trees <- list()
    for (x in 1:ntree) {
      for (t in c(treatments, control)) {
        if (t == treatments[1]) {
          temp_train_data <- sample_n(data[data[, t] == 1, ], B * (sum(data[, t] == 1) / nrow(data)), replace = TRUE)
        } else {
          temp_train_data <- rbind(temp_train_data, sample_n(data[data[, t] == 1, ], B * (sum(data[, t] == 1) / nrow(data)),
            replace = TRUE
          ))
        }
      }
      temp_train_data <- na.omit(temp_train_data)
      trees[[x]] <- build_cts_tree(response, control, treatments, temp_train_data, m_try, n_reg,
        min_split,
        parent_predictions = NA
      )
    }
    return(trees)
  }
}


# The main method to build a single tree for CTS
build_cts_tree <- function(response, control, treatments, data, m_try, n_reg, min_split,
                           min_gain = 0, parent_predictions = NA, depth = 0) {
  # We sample mtry covariates to use for the next split
  retain_cols <- c(treatments, control, response)
  sample_cols <- setdiff(colnames(data), retain_cols)
  temp_cols <- sample(sample_cols, m_try, replace = F)
  chosen_cols <- c(temp_cols, retain_cols)
  test_list <- set_up_tests(data[, temp_cols], TRUE)


  # We create the current node and add the current predictions and the number of samples for each treatment
  node <- list()
  results <- c()
  n_treatment_samples <- c()
  if (is.na(parent_predictions)) {
    # If the are no predictions from the parent node we are in the root. In the root the predictions are simply
    # the sample average.
    node[["type"]] <- "root"
    for (t in c(treatments, control)) {
      results <- c(results, mean(data[data[, t] == 1, response]))
      n_treatment_samples <- c(n_treatment_samples, nrow(data[data[, t] == 1, ]))
    }
    names(results) <- c(treatments, control)
    node[["results"]] <- results
  } else {
    # If there are parent predictions we calculate the predictions for the current node according to the paper.
    node[["type"]] <- "node"
    for (t in c(treatments, control)) {
      n_samples <- nrow(data[data[, t] == 1, ])
      n_treatment_samples <- c(n_treatment_samples, n_samples)
      if (n_samples < min_split) {
        # If the number of samples with treatment t is smaller than min_split we just take the parent prediction
        # for treatment t.
        results <- c(results, parent_predictions[[t]])
      } else {
        # Else with calculate the current prediction for treatment t according to the formula from the paper.
        results <- c(
          results,
          (sum(data[data[, t] == 1, response]) + parent_predictions[[t]] * n_reg) / (sum(data[, t] == 1) + n_reg)
        )
      }
    }
  }
  # Once we calculated the the predictions we add the new information to the current node.
  names(results) <- c(treatments, control)
  names(n_treatment_samples) <- c(treatments, control)
  node[["results"]] <- results
  node[["n_samples"]] <- nrow(data)
  node[["n_samplse_treatments"]] <- n_treatment_samples

  # Next we check if we terminate the algorithm.
  terminate <- TRUE
  # If there are fewer than min_split samples in the data set for each treatment we terminate.
  for (t in c(treatments, control)) {
    if (sum(data[, t] == 1) >= min_split) {
      terminate <- FALSE
    }
  }
  # If the response is the same for all samples we terminate.
  if (sum(data[, response] == data[1, response]) == nrow(data)) {
    terminate <- TRUE
  }
  # If either of the previous conditions is tree we simply change the type of the current node to "leaf" and
  # return it.
  if (terminate) {
    node[["type"]] <- "leaf"
    return(node)
  }

  # If we have not terminated we select the next split.
  temp_split <- select_cts_split(
    node[["results"]], data, treatments, response, control, n_reg, min_split,
    test_list, min_gain
  )

  # Return a leaf, if there is no split with gain > 0
  if (temp_split == -1) {
    node[["type"]] <- "leaf"
    return(node)
  }

  # If there is a split with gain > 0 the information is added to the node.
  node[["split"]] <- temp_split
  # The data is split according to the selected split and two new trees are grown.
  if (names(temp_split) %in% names(test_list$categorical)) {
    node[["left"]] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)] == temp_split[[1]], ],
      m_try, n_reg, min_split, min_gain, node[["results"]],
      depth = depth + 1
    )
    node[["right"]] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)] != temp_split[[1]], ],
      m_try, n_reg, min_split, min_gain, node[["results"]],
      depth = depth + 1
    )
  } else {
    node[["left"]] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)] < temp_split[[1]], ],
      m_try, n_reg, min_split, min_gain, node[["results"]],
      depth = depth + 1
    )
    node[["right"]] <- build_cts_tree(response, control, treatments, data[data[names(temp_split)] >= temp_split[[1]], ],
      m_try, n_reg, min_split, min_gain, node[["results"]],
      depth = depth + 1
    )
  }
  return(node)
}


# This function goes through all the possible splits and calculates the gain. Then it returns the split with
# the highest gain. If no split has a gain > 0 the the function returns -1.
select_cts_split <- function(parent_predictions, data, treatments, response, control, n_reg, min_split,
                             test_list, min_gain) {
  gain_list <- c()
  name_list <- c()
  if (length(test_list$categorical) > 0) {
    for (x in 1:length(test_list$categorical)) {
      col_name <- names(test_list$categorical[x])
      for (y in 1:length(test_list$categorical[[x]])) {
        test_case <- test_list$categorical[[x]][y]
        new_name <- paste(col_name, as.character(test_case), sep = "@@")
        gain_list <- c(gain_list, cts_gain(
          test_case, treatments, control, response, data, "categorical", col_name,
          parent_predictions, n_reg, min_split
        ))
        name_list <- c(name_list, new_name)
      }
    }
  }
  if (length(test_list$numerical) > 0) {
    for (x in 1:length(test_list$numerical)) {
      col_name <- names(test_list$numerical[x])
      for (y in 1:length(test_list$numerical[[x]])) {
        test_case <- test_list$numerical[[x]][y]
        new_name <- paste(col_name, as.character(test_case), sep = "@@")
        gain_list <- c(gain_list, cts_gain(
          test_case, treatments, control, response, data, "numerical", col_name,
          parent_predictions, n_reg, min_split
        ))
        name_list <- c(name_list, new_name)
      }
    }
  }
  if (is.na(max(gain_list))) {
    return(-1)
  }
  if (max(gain_list) > min_gain) {
    temp_string <- name_list[match(max(gain_list), gain_list)]
    temp_result <- strsplit(temp_string, split = "@@", fixed = TRUE)
    options(warn = -1)
    if (is.na(as.numeric(temp_result[[1]][2]))) {
      result <- c(temp_result[[1]][2])
      names(result) <- temp_result[[1]][1]
      options(warn = 0)
      return(result)
    } else {
      options(warn = 0)
      result <- c(as.numeric(temp_result[[1]][2]))
      names(result) <- temp_result[[1]][1]
      return(result)
    }
  } else {
    return(-1)
  }
}


# This function calculates the gain for a given split.
cts_gain <- function(test_case, treatments, control, response, data, test_type, col_name, parent_predictions,
                     n_reg, min_split) {
  # Initialize the variables for the maximum outcomes.
  max_left <- 0
  max_right <- 0
  max_root <- max(parent_predictions)
  # Split the data according to the current split.
  if (test_type == "categorical") {
    data_left <- data[data[, col_name] == test_case, ]
    data_right <- data[data[, col_name] != test_case, ]
  } else {
    data_left <- data[data[, col_name] < test_case, ]
    data_right <- data[data[, col_name] >= test_case, ]
  }
  # Calculate the percentage of samples in the left and right data set.
  p_left <- nrow(data_left) / nrow(data)
  p_right <- nrow(data_right) / nrow(data)
  if (p_left == 0 || p_right == 0) {
    return(-1)
  }
  # For each treatment and control calculate the expected outcome. If the expected outcome is greater than the
  # current max, it becomes the new max.
  for (t in c(treatments, control)) {
    n_left <- nrow(data_left[data_left[, t] == 1, ])
    n_right <- nrow(data_right[data_right[, t] == 1, ])
    if (n_left < min_split) {
      exp_left <- parent_predictions[[t]]
    } else {
      exp_left <- (sum(data_left[data_left[, t] == 1, response]) + parent_predictions[[t]] * n_reg) / (n_left + n_reg)
    }
    if (n_right < min_split) {
      exp_right <- parent_predictions[[t]]
    } else {
      exp_right <- (sum(data_right[data_right[, t] == 1, response]) + parent_predictions[[t]] * n_reg) / (n_right + n_reg)
    }
    max_left <- max(max_left, exp_left)
    max_right <- max(max_right, exp_right)
  }
  # Return the calculated gain.
  return(p_left * max_left + p_right * max_right - max_root)
}
