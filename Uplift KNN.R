# Tetsing of the uKNN implementation in the uplift package for Multiple Treatments

list.of.packages <- c("dplyr", "caret", "uplift")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) install.packages(new.packages)

library(dplyr)
library(caret)
library(uplift)

# Due in inconsitencies with loaded packages
select <- dplyr::select

# Multiple Treatment
data <- read.csv('Email.csv')

# First only look at conversion
data <- data %>% select(-visit, -spend)

# Only select first 10,000 entries as the uKNN model cannot handle the full dataset
data <- data[1:10000 , ]

# Split into Test and Training set
idx <- createDataPartition(y = data$segment, p=0.3, list = FALSE)
train <- data[-idx, ] 
test <- data[idx, ]


# Data Exploration
# explore only works with binary treatment / control variable
data$treated <- 0
data <- data %>% mutate(treated=replace(treated, segment != 'No E-Mail', 1)) %>% as.data.frame()
data <- data %>% select(-segment)


data %>% colnames()

explore(conversion~channel+newbie+history+history_segment+mens+zip_code+trt(treated), data)


# Create the model
knn <- upliftKNN(train[, 1:8], test[, 1:8], train$conversion, train$segment, k = 1, dist.method = "euclidean",
                 p = 2, ties.meth = "min", agg.method = "mean")

# Output prediction of treatment for test group
knn
# treatments only rarely assigned ...
# but No-Email also as treatment should be assigned..
