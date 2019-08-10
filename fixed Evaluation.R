
source('Evaluation Methods.R')

pred <- read.csv('rzp tree pred.csv')

head(pred)

## TAKEN FROM 'Models Training.R'

# ### Results Preparation to bring into equal format
# # Calculate Uplift for each T
# pred[ , "uplift_men_treatment"] <- pred[ , 1] - pred[ , 3]
# pred[ , "uplift_women_treatment"] <- pred[ , 2] - pred[ , 3]
# pred[ , "Treatment"] <- colnames(pred)[apply(pred[, 1:3], 1, which.max)]

colnames(pred)[4] <- 'uplift_men_treatment'
colnames(pred)[5] <- 'uplift_women_treatment'


expected_percentile_response(pred)

matching_evaluation(pred,control_level = 'control')


