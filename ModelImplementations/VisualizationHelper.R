library("DiagrammeR")
library("gridExtra")
# This method visualizes the incremental expected outcome for one or more models.
# If there are multiple predictions for each model, the function aggregates those predictions and can create
# errorbars to show the confidence interval of the predictions.
visualize <- function(temp_data,n_treated = NULL,errorbars = TRUE,multiplot = FALSE){
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
  tgc <- summarySE(data=temp_df, measurevar="values", groupvars=c("percentile","model"))
  pd <- position_dodge(1) # move them .05 to the left and right
  p1 <- ggplot(tgc, aes(x=percentile, y=mean,color=model)) + 
  {if(errorbars && sum(is.na(tgc)) == 0) geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci), position = pd)} +
    geom_line() +
    geom_point() +
    xlab("Percent assigned according to Model Prediction") +
    ylab("Expected Outcome per Person") +
    {if(multiplot) facet_wrap(~model)} +
    {if(multiplot) theme(legend.position = "none")} +
    ggtitle("Mean and Confidence Interval for Expected Outcome")
  if(!is.null(n_treated)){
    agg_df<- aggregate(n_treated$PercTreated, by=list(n_treated$Treatment,n_treated$Model,n_treated$Decile), 
                       FUN=mean)
    colnames(agg_df) <- c("Treatment","Model","Decile","PercTreated")
    p2 <- ggplot(agg_df, aes(fill=Treatment, y=PercTreated, x=Decile)) + 
      geom_bar(position="stack", stat="identity") +
      facet_wrap(~Model) 
    print(p1)
    print(p2)
  } else{
    print(p1)
  }
}


visualize_qini_uplift <- function(temp_data,type,multiple_predictions = TRUE,errorbars = TRUE,multiplot=TRUE){
  values <- c()
  percentile <- c()
  model <- c()
  temp_df <- temp_data
  if(multiple_predictions){
    tgc <- summarySE(data=temp_df, measurevar="values", groupvars=c("percentile","model"))
    if(errorbars){
      print(ggplot(tgc, aes(x=percentile, y=mean,color=model)) + 
              geom_errorbar(aes(ymin=mean-ci, ymax=mean+ci)) +
              geom_line() +
              geom_point() +
              xlab("Percent assigned according to Model Prediction") +
              ylab(paste("Cummulated",type,sep = " ")) +
              {if(multiplot) facet_wrap(~model)} +
              {if(multiplot) theme(legend.position = "none")} +
              ggtitle(paste("Mean and Confidence Interval for",type,"score",sep=" ")))
    } else{
      print(ggplot(tgc, aes(x=percentile, y=mean,color=model)) + 
              geom_line() +
              geom_point() +
              xlab("Percent assigned according to Model Prediction") +
              ylab(paste("Cummulated",type,sep = " ")) +
              {if(multiplot) facet_wrap(~model)} +
              {if(multiplot) theme(legend.position = "none")} +
              ggtitle(paste("Mean",type,"score",sep=" ")))
    }
    
  } else{
    print(ggplot(temp_df, aes(x=percentile, y=values,color=model)) + 
            geom_line() +
            geom_point() +
            xlab("Percent assigned according to Model Prediction") +
            ylab(paste("Cummulated",type,sep = " ")) +
            {if(multiplot) facet_wrap(~model)} +
            {if(multiplot) theme(legend.position = "none")} +
            ggtitle(paste("Mean",type,"score",sep=" ")))
  }
}


# Function used to aggregate multiple predictions for the same model.
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
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}




# Methods used to plot a tree
# It is advisable to build a shallow tree (low value for max_depth) to make the graph somewhat readable.
visualize_tree <- function(tree, only_split = F){
  plot_list <- get_plot_params(tree)
  plot_string <- "digraph flowchart {  \n node [fontname = Helvetica, shape = rectangle]"
  
  label = ""
  split_counter <- 1
  if(only_split){
    for(x in 1:length(plot_list[[1]])){
      if(!(plot_list[[1]][x]) == "leaf"){
        label = paste(label, paste(as.character(plot_list[[5]][x]), "[label = '",
                                   "Split", names(plot_list[[2]][split_counter]), as.character(plot_list[[2]][[split_counter]]),
                                   "'] \n", sep = " "), sep = "")
        split_counter <- split_counter + 1
      } else{
        label = paste(label, paste(as.character(plot_list[[5]][x]), "[label = '", plot_list[[1]][x], "\n",
                                   "Number of Samples", as.character(plot_list[[3]][x]), "\n", 
                                   names(plot_list[[4]][((x-1) * 3 + 1)]), as.character(plot_list[[4]][[((x-1) * 3 + 1)]]), "\n",
                                   names(plot_list[[4]][((x-1) * 3 + 2)]), as.character(plot_list[[4]][[((x-1) * 3 + 2)]]), "\n",
                                   names(plot_list[[4]][((x-1) * 3 + 3)]), as.character(plot_list[[4]][[((x-1) * 3 + 3)]]), "\n",
                                   "'] \n", sep = " "), sep = "")
      }
    }
  } else{
    for(x in 1:length(plot_list[[1]])){
      if(!(plot_list[[1]][x]) == "leaf"){
        label = paste(label, paste(as.character(plot_list[[5]][x]), "[label = '", plot_list[[1]][x], "\n", 
                                   "Split", names(plot_list[[2]][split_counter]), as.character(plot_list[[2]][[split_counter]]), "\n", 
                                   "Number of Samples", as.character(plot_list[[3]][x]), "\n", 
                                   names(plot_list[[4]][((x-1) * 3 + 1)]), as.character(plot_list[[4]][[((x-1) * 3 + 1)]]), "\n",
                                   names(plot_list[[4]][((x-1) * 3 + 2)]), as.character(plot_list[[4]][[((x-1) * 3 + 2)]]), "\n",
                                   names(plot_list[[4]][((x-1) * 3 + 3)]), as.character(plot_list[[4]][[((x-1) * 3 + 3)]]), "\n",
                                   "'] \n", sep = " "), sep = "")
        split_counter <- split_counter + 1
      } else{
        label = paste(label, paste(as.character(plot_list[[5]][x]), "[label = '", plot_list[[1]][x], "\n",
                                   "Number of Samples", as.character(plot_list[[3]][x]), "\n", 
                                   names(plot_list[[4]][((x-1) * 3 + 1)]), as.character(plot_list[[4]][[((x-1) * 3 + 1)]]), "\n",
                                   names(plot_list[[4]][((x-1) * 3 + 2)]), as.character(plot_list[[4]][[((x-1) * 3 + 2)]]), "\n",
                                   names(plot_list[[4]][((x-1) * 3 + 3)]), as.character(plot_list[[4]][[((x-1) * 3 + 3)]]), "\n",
                                   "'] \n", sep = " "), sep = "")
      }
    }
  }
  links <- paste(plot_list[[6]], collapse = "; \n ")
  
  print(grViz(paste(plot_string,label,links,"}",sep = "\n")))
}


# Helperfunction for the tree visualization.
get_plot_params <- function(tree,counter = 0){
  type_list <- c(tree[["type"]])
  split_list <- c(tree[["split"]])
  n_sample_list <- c(tree[["n_samples"]])
  results_list <- c(tree[["results"]])
  counter_list <- c(counter)
  connection_list <- c()
  if(!is.null(tree[["left"]])){
    temp_list <- get_plot_params(tree[["left"]], counter = counter + 1)
    type_list <- c(type_list,temp_list[[1]])
    split_list <- c(split_list,temp_list[[2]])
    n_sample_list <- c(n_sample_list,temp_list[[3]])
    results_list <- c(results_list,temp_list[[4]])
    counter_list <- c(counter_list,temp_list[[5]])
    connection_list <- c(paste(as.character(counter), "->", as.character(temp_list[[5]][1]), sep = " "))
    if (length(temp_list) > 5) {
      connection_list <- c(connection_list, temp_list[[6]])
    }
    
  }
  if(!is.null(tree[["right"]])){
    temp_list <- get_plot_params(tree[["right"]], counter = max(counter_list) + 1)
    type_list <- c(type_list,temp_list[[1]])
    split_list <- c(split_list,temp_list[[2]])
    n_sample_list <- c(n_sample_list,temp_list[[3]])
    results_list <- c(results_list,temp_list[[4]])
    counter_list <- c(counter_list,temp_list[[5]])
    connection_list <- c(connection_list,paste(as.character(counter), "->", as.character(temp_list[[5]][1]), sep = " "))
    if (length(temp_list) > 5) {
      connection_list <- c(connection_list, temp_list[[6]])
    }
    
  }
  result_list <- list()
  result_list[[1]] <- type_list
  result_list[[2]] <- split_list
  result_list[[3]] <- n_sample_list
  result_list[[4]] <- round(results_list,2)
  result_list[[5]] <- counter_list
  result_list[[6]] <- connection_list
  return(result_list)
}