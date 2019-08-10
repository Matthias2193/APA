source('Preprocessing.R')
source('DecisionTreeImplementation.R')
library('caret')
source('CausalTree.R')
source('Causal Forest.R')

remain_cols <- c("female","age","voted","hh_size","treatment_CivicDuty","treatment_Self",
                 "treatment_Control","treatment_Hawthorne","treatment_Neighbors")
individuals <- individuals[,remain_cols]
individuals$female <- as.factor(individuals$female)
individuals$hh_size <- as.factor(individuals$hh_size)
individuals$voted <- as.numeric(individuals$voted) - 1

response <- 'voted'
control <- 'treatment_Control'


treatments <- c("treatment_CivicDuty","treatment_Self","treatment_Hawthorne","treatment_Neighbors")
treatment_list <- c("treatment_CivicDuty","treatment_Self")
test_list <- set_up_tests(individuals[,c("female","age","hh_size")],TRUE)
n_treatments <- length(treatment_list)
individuals$temp_feature1 <- rnorm(nrow(individuals),0,1)
individuals$temp_feature2 <- rbinom(n=nrow(individuals), size=1, prob=0.5)
individuals$temp_feature2 <- as.factor(individuals$temp_feature2)

#Benchmark
size_vector <- c(5000,15000,30000)
time_vector_rzp <- list()
time_vector_simple <- list()
time_vector_forest <- list()
time_vector_causal_tree <- list()
remain_cols_all <- c("female","age","voted","hh_size","treatment_Control","temp_feature1","temp_feature2")
for(r in 5:7){
  remain_cols <- remain_cols_all[1:r]
  for(i in 2:4){
    treatment_list <- treatments[1:i]
    n_treatments <- length(treatment_list)
    temp_time_vec_rzp <- c()
    temp_time_vec_simple <- c()
    temp_time_vec_forest <- c()
    temp_time_vec_causal_tree <- c()
    for(s in size_vector){
      idx <- createDataPartition(y = individuals[ , response], p=round(s/nrow(individuals),4), list = FALSE)
      temp_data <- individuals[idx,c(remain_cols,treatment_list)]
      idx <- createDataPartition(y = temp_data[ , response], p=0.2, list = FALSE)
      val <- temp_data[idx,]
      train <- temp_data[-idx,]
      start_time <- Sys.time()
      test_tree <- create_node(train,0,100,treatment_list,response,control,test_list,
                               normalize  = TRUE, l = rep(1/n_treatments,n_treatments),
                               g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
      temp_time_vec_rzp <- c(temp_time_vec_rzp,difftime(Sys.time(), start_time, units='mins'))
      start_time <- Sys.time()
      test_tree <- create_node(train,0,100,treatment_list,response,control,test_list,
                               normalize  = TRUE, l = rep(1/n_treatments,n_treatments),
                               g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments),criterion = 2)
      temp_time_vec_simple <- c(temp_time_vec_simple,difftime(Sys.time(), start_time, units='mins'))
      start_time <- Sys.time()
      forest <- build_forest(train,val,treatment_list,response,control,n_trees = 20,n_features = 2,criterion = 1,
                             pruning = F,l = rep(1/n_treatments,n_treatments),
                             g = matrix(1/n_treatments^2,nrow = n_treatments, ncol = n_treatments))
      temp_time_vec_forest <- c(temp_time_vec_forest,difftime(Sys.time(), start_time, units='mins'))
      start_time <- Sys.time()
      for(t in treatment_list){
        tree <- causalForest(as.formula(paste(paste(response, " ~ ",sep = ""),
                                            paste(setdiff(remain_cols,c("treatment_Control")),collapse = "+"))),
                           data = train[,c(setdiff(remain_cols,c("treatment_Control")),t)],
                           treatment = train[,t], split.Rule = "CT", cv.option = "CT", split.Honest = T,
                           cv.Honest = T, split.Bucket = F, minsize = 20, propensity = 0.5, mtry = 2,
                           num.trees = 100, ncov_sample = 3, 
                           ncolx = (ncol(train[,c(setdiff(remain_cols,c("treatment_Control")))])-1))
      }
      temp_time_vec_causal_tree <- c(temp_time_vec_causal_tree,difftime(Sys.time(), start_time, units='mins'))
    }
    time_vector_rzp[[i-1]] <- temp_time_vec_rzp
    time_vector_simple[[i-1]] <- temp_time_vec_simple
    time_vector_forest[[i-1]] <- temp_time_vec_forest
    time_vector_causal_tree[[i-1]] <- temp_time_vec_causal_tree
  }
  rzp_df <- data.frame(cbind(time_vector_rzp[[1]],cbind(time_vector_rzp[[2]],time_vector_rzp[[3]])))
  simple_df <- data.frame(cbind(time_vector_simple[[1]],cbind(time_vector_simple[[2]],time_vector_simple[[3]])))
  forest_df <- data.frame(cbind(time_vector_forest[[1]],cbind(time_vector_forest[[2]],time_vector_forest[[3]])))
  causal_tree_df <- data.frame(cbind(time_vector_causal_tree[[1]],cbind(time_vector_causal_tree[[2]],
                                                                        time_vector_causal_tree[[3]])))
  colnames(rzp_df) <- colnames(simple_df) <- colnames(causal_tree_df) <- colnames(forest_df) <- c("2 Treatments","3 Treatments", "4 Treatments") 
  rownames(rzp_df) <- rownames(simple_df) <- rownames(causal_tree_df) <- rownames(forest_df) <- c("50k","150k", "300k")
  attach(mtcars)
  par(mfrow=c(1,3))
  for(c in c("2 Treatments","3 Treatments", "4 Treatments")){
    counts <- rbind(rbind(rbind(round(simple_df[,c],2),round(rzp_df[,c],2)),round(causal_tree_df[,c],2)),
                    round(forest_df[,c],2))
    x <- barplot(counts, main= c,
                 xlab="Number of Simples", col=c("darkblue","red","green","orange"),
                 beside=TRUE, ylab = "Time (mins)",
                 names.arg = c("40k","120k", "240k"),ylim = c(0,max(counts)+0.1*max(counts)))
    legend ("topleft", 
            c("Simple","Rzp","Causal Forest","Rzp-Forest"),
            fill = c("darkblue","red","green","orange"),
            cex = 0.7)
    y <- as.matrix(counts)
    text(x,y+0.07*max(counts),labels=as.character(y)) 
  }
}