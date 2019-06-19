## Model tuning

# for pre set data split
# use train data to tune parameters of DT and RF Ridge

# best param for each T needed ...


# load train data for each Treatment and use for CV fitting
train_m <- read.csv('Mens Train.csv')
train_w <- read.csv('Womens Train.csv')
train_c <- read.csv('Control Train.csv')

response <- 'conversion'

## Fitting Logit not useful...


#################
# RF
train <- train_m

train[, response] <- as.factor(train[, response])

## TODO how to set / adjust parameters ??
rf <- randomForest(as.formula(paste(response, "~.")), data = train, mtry=3, ntree = 350)

# only mtry tuned...
# how to find best ntree ???
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
set.seed(1234)
mtry <- 7
rf_random <- train(as.formula(paste(response, "~.")), data=train, method="rf", metric='Accuracy', tuneLength=15, trControl=control)
print(rf_random)
plot(rf_random)


# best mtry e (1,2,3)
