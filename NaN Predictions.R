source('DecisionTreeImplementation.R')

set.seed(212121)

#also contains NA values
## 2121
#212121

#####################################################################################
### Conversion Prediction
#####################################################################################

#Data import
email <- read.csv('Email.csv')

response <- 'conversion'
email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)
email$mens <- as.factor(email$mens)
email$womens <- as.factor(email$womens)
email$newbie <- as.factor(email$newbie)

email$visit <- email$spend <- email$segment <- NULL

email$spend_bins <- cut(email$spend, 20, include.lowest = TRUE)
levels(email$spend_bins) <- 1:20
email$spend_bins <- as.numeric(email$spend_bins)

data1 <- email[email$mens == 1,]
data2 <- email[email$womens==1,]

dist1 <- c()
dist2 <- c()

for(x in 1:20){
  dist1 <- c(dist1,nrow(data1[data1$spend_bins == x,])/nrow(data1))
  dist2 <- c(dist2,nrow(data2[data2$spend_bins == x,])/nrow(data2))
}


dist1 <- fitdist(data1$spend, distr = "gamma", method = "mle")
dist2 <- fitdistr(data2$spend,'Poisson')

temp1 <- rpois(100,dist1$estimate)
temp2 <- rpois(100,dist2$estimate)

temp_data <- cbind(data1$visit,data2$visit[1:10])

temp <-KLdiv(temp_data)
temp[2,1]# Split into test and train data
idx <- createDataPartition(y = email[ , response], p=0.3, list = FALSE)

train <- email[-idx, ]

test <- email[idx, ]

# Partition training data for pruning
p_idx <- createDataPartition(y = train[ , response], p=0.3, list = FALSE)

val <- train[p_idx,]
train_tree <- train[-p_idx,]

treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(train_tree[,c("recency","history_segment","history","mens","womens","zip_code",
                                        "newbie","channel")],TRUE)

######
# Tree simple Criterion
raw_tree <- create_node(train_tree,0,100,treatment_list,response,'control',test_list,criterion = 2)

pruned_tree <- raw_tree

# add to the result df the outcome, assignment and calculate uplift for each T
pred <- predict.dt.as.df(pruned_tree, test)

sum(is.na(pred))

