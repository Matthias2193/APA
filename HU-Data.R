source("DecisionTreeImplementation.R")
library(caret)
hu_data <- read.csv("Data/explore_mt.csv",sep = ";")
for(x in levels(hu_data$multi_treat)){
  hu_data[x] <- ifelse(hu_data$multi_treat == x ,1,0)
}
tempfunction <- function(x){
  templist <- strsplit(x,",")
  new_string <- ""
  r <- 1
  for(s in templist[[1]]){
    if (r == 2) {
      new_string <- paste(new_string,s,sep = ".")
    } else{
      new_string <- paste(new_string,s,sep="")
    }
    r <- r+1
  }
  return(new_string)
}

hu_data$checkoutAmount <- as.numeric(lapply(as.character(hu_data$checkoutAmount), tempfunction))
hu_data$Number.of.seconds.between.last.and.previous.views <- 
  as.numeric(lapply(as.character(hu_data$Number.of.seconds.between.last.and.previous.views), tempfunction))
for (x in grep("^log.of",colnames(hu_data))) {
  hu_data[,x] <- as.numeric(lapply(as.character(hu_data[,x]), tempfunction))
}

for (x in colnames(train[,16:155])) {
  if(is.na(mean(hu_data[[x]]))){
    print(x)
  }
}

response <- 'checkoutAmount'
control <- '0'




set.seed(1234)
# Split into test and train data
idx <- createDataPartition(y = hu_data[ , response], p=0.3, list = FALSE)

train <- hu_data[-idx, ]

test <- hu_data[idx, ]

# Partition training data for pruning
p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)

val <- train[p_idx,]
train <- train[-p_idx,]

treatment_list <- levels(hu_data$multi_treat)[2:7]
test_list <- set_up_tests(train[,colnames(train[,16:155])],TRUE)


raw_tree <- build_tree(train,0,100,treatment_list,response,control,test_list,criterion = 2)
pruned_tree <- simple_prune_tree(raw_tree,val,treatment_list,test_list,response,control)

forest <- parallel_build_forest(train,val,treatment_list,response,'0',n_trees = 50,n_features = 100,criterion = 2, pruning = F)

mean(hu_data[[treatment_list[1]]])
typeof(hu_data[treatment_list[1]][0])
mean(hu_data["log.of.SecondsSinceOn.about."])
