#This script contains the implementation of the decision tree following the paper 'Decision trees for uplift modeling
#with single and multiple treatments'

#Importing libraries
list.of.packages <- c("FNN",'LaplacesDemon','philentropy')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library('FNN')
library('LaplacesDemon')
library('philentropy')



#Transform the treatment in dummy variables
email$men_treatment <- ifelse(email$segment=='Mens E-Mail',1,0)
email$women_treatment <- ifelse(email$segment=='Womens E-Mail',1,0)
email$control <- ifelse(email$segment=='No E-Mail',1,0)



#Set up tests
set_up_tests <- function(x){
  type_list <- sapply(x, class)
  categorical_splits = list()
  numerical_splits = list()
  for(n in colnames(x)){
    if(type_list[[n]] == 'factor'){
      categorical_splits[[n]]<-levels(x[,n])
    }
    else{
      temp_list <- sort(unique(x[,n]))
      r <- 1
      s <- 2
      final_list <- c()
      while(r < length(temp_list)){
        final_list <- c(final_list,(temp_list[r]+temp_list[s])/2)
        r <- r+1
        s <- s+1
      }
      numerical_splits[[n]] <- final_list
    }
  }
  output = list()
  output[['categorical']] <- categorical_splits
  output[['numerical']] <- numerical_splits
  return(output)
}


select_split <- function(a,l,g,divergence,test_list,treatment,control,target,temp_data){
  gain_list <- c()
  name_list <- c()
  for(t in test_list$categorical){
    gain_list <- c(gain_list,gain(a,l,g,divergence,t[[1]],
                                  treatment,control,target,temp_data,'categorical',names(t)))
    name_list <- c(name_list,t)
  }
  for(t in test_list$numerical){
    gain_list <- c(gain_list,gain(a,l,g,divergence,number_of_treatments,t[[1]],
                                  treatment,control,target,temp_data,'numerical',names(t)))
    name_list <- c(name_list,t)
  }
  return(name_list[match(max(gain_list),gain_list)])
}

gain <- function(a,l,g,divergence, test_case,treatment,control,target,temp_data,test_type,test_col){
  conditional <- conditional_divergence(a,l,g,divergence,test_case,treatments,control,target,temp_data,
                                         test_type,test_col)
  multiple <- multiple_divergence(a,l,g,divergence,treatments,control,target,temp_data,test_col)
  return(conditional-multiple)
}

multiple_divergence <- function(a,l,g,divergence,treatments,control,target,temp_data){
  divergence_function <- match.fun(divergence)
  multiple <- 0
  for(t in length(treatments)){
    multiple <- multiple + a*l[t]*divergence_function(temp_data[temp_data[treatments[t]]==1,target],
                                                      temp_data[temp_data[control]==1,target])
    between_treatments <- 0
    for(s in length(treatments)){
      between_treatments <- between_treatments + g[t,s]*divergence_function(
        temp_data[temp_data[treatments[t]]==1,target],temp_data[temp_data[treatments[s]]==1,target])
    }
    multiple <- multiple + (1-a)*between_treatments
  }
  return(multiple)
}

conditional_divergence <- function(a,l,g,divergence,test_case,treatments,control,target,temp_data,test_type,
                                   test_col){
  div <- 0
  if(test_type == 'categorical'){
    for(t in test_case){
      div <- div + nrow(temp_data[temp_data[test_col]==t,])/nrow(temp_data)*
        multiple_divergence(a,l,g,divergence,treatments,control,target,
                            temp_data[temp_data[test_col]==t,])
    }
  }
  else{
    div <- div + nrow(temp_data[temp_data[test_col]<t,])/nrow(temp_data)*
      multiple_divergence(a,l,g,divergence,treatments,control,target,
                          temp_data[temp_data[test_col]<t,])
    div <- div + nrow(temp_data[temp_data[test_col]>t,])/nrow(temp_data)*
      multiple_divergence(a,l,g,divergence,treatments,control,target,
                          temp_data[temp_data[test_col]>t,])
  }
  return(div)
}

binary_KL_divergence <- function(x,y){
  p <- mean(x)
  q <- mean(y)
  return((p*log(p/q))+((1-p)*log((1-p)/1-q)))
}



####Test Area

divergence <-'binary_KL_divergence'
multiple <- 0
a <- 1
l<- c(0.5,0.5)
g<- matrix(0.25,nrow = 2,ncol = 2)
target <- 'visit'
treatment_list <- c('men_treatment','women_treatment')
test_list <- set_up_tests(email)
control <- 'control'
temp_test <- test_list$categorical[1]
test_col <- names(temp_test)

conditional_divergence(a,l,g,'binary_KL_divergence',temp_test,treatment_list,control,target,email,'categorical',
                    names(temp_test))




select_split(a,l,g,'binary_KL_divergence',test_list,treatment_list,'control','visit',email)

