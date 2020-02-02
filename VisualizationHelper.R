summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

visualize <- function(temp_data){
  values <- c()
  percentile <- c()
  model <- c()
  for(f in 1:nrow(temp_data)){
    if(length(values) == 0){
      values <- temp_data[f,1:11]
      percentile <- colnames(temp_data)[1:11]
      model <- rep(temp_data[f,12],11)
    } else{
      values <- c(values,temp_data[f,1:11])
      percentile <- c(percentile, colnames(temp_data)[1:11])
      model <- c(model,rep(temp_data[f,12],11))
    }
  }
  temp_df <- data.frame(cbind(values,percentile,model))
  rownames(temp_df) <- 1:nrow(temp_df)
  colnames(temp_df) <- c("values","percentile","model")
  for(c in 1:2){
    temp_df[,c] <- as.numeric(as.character(temp_df[,c]))
  }
  temp_df[,3] <- as.character(temp_df[,3])
  tgc <- summarySE(temp_df, measurevar="values", groupvars=c("percentile","model"))
  # new_tgc <- tgc[order(tgc$model),]
  # rownames(new_tgc) <- 1:nrow(new_tgc)
  pd <- position_dodge(0.1) # move them .05 to the left and right
  print(ggplot(tgc, aes(x=percentile, y=values,color=model)) + 
          geom_errorbar(aes(ymin=values-ci, ymax=values+ci), width=1) +
          geom_line() +
          geom_point() +
          xlab("Percent assigned according to Model Prediction") +
          ylab("Expected Outcome per Person") +
          ggtitle("Mean and Confidence Interval for Expected Outcome"))
}
