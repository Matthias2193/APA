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
treatment_list <- c('men_treatment','women_treatment')


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


select_split <- function(divergence,treatments,control,target,temp_data){
  gain_list <- c()
  name_list <- c()
  for(t in test_list){
    gain_list <- c(gain_list,gain())
    name_list <- c(name_list,t)
  }
  return(name_list[match(max(gain_list),gain_list)])
}

gain <- function(a,l,g,divergence,number_of_treatments, test_case,treatment,control,target,temp_data){
  conditional <- conditional_convergance(a,l,g,divergence,test_case,treatments,control,target,temp_data)
  multiple <- multiple_divergence(a,l,g,divergence,treatments,control,target,temp_data)
  return(conditional-multiple)
}

multiple_divergence <- function(a,l,g,divergence,treatments,control,target,temp_data){
  divergence_function <- match.fun(divergence)
  multiple <- 0
  for(t in length(treatments)){
    multiple <- multiple + a*l[t]*divergence_function(temp_data[temp_data$treatments[t]==1,]$target,
                                                      temp_data[temp_data$control==1,]$target)
    between_treatments <- 0
    for(s in length(treatments)){
      between_treatments <- between_treatments + g[t,s]*divergence_function(
        temp_data[temp_data$treatments[t]==1,]$target,temp_data[temp_data$treatments[s]==1,]$target)
    }
    multiple <- multiple + (1-a)*between_treatments
  }
  return(multiple)
}

binary_KL_divergence <- function(x,y){
  p <- mean(x)
  q <- mean(y)
  return((p*log(p/q))+((1-p)*log((1-p)/1-q)))
}



####Test Area

multiple_divergence(1,c(0.5,0.5),matrix(0.25,nrow = 2,ncol = 2),'KL.divergence',treatment_list,'control',
                    'visit',email)

divergence_function <- match.fun('binary_KL_divergence')
multiple <- 0
a <- 1
l<- c(0.5,0.5)
g<- matrix(0.25,nrow = 2,ncol = 2)
target <- 'visit'
temp_data <- email
for(t in length(treatment_list)){
  multiple <- multiple + a*l[t]*divergence_function(temp_data[temp_data[,treatment_list[t]]==1,target],
                                                    temp_data[temp_data[,'control']==1,target])
  between_treatments <- 0
  for(s in length(treatments)){
    between_treatments <- between_treatments + g[t,s]*divergence_function(
      temp_data[temp_data$treatments[t]==1,]$target,temp_data[temp_data$treatments[s]==1,]$target)
  }
  multiple <- multiple + (1-a)*between_treatments
}
return(multiple)


email
X<- rexp(10000, rate=0.2)
Y<- rexp(10000, rate=0.4)
KL.divergence(X,Y)

conditional_divergence <- function(a,l,g,divergence,test_case,treatments,control,target,temp_data){
  divergence_function <- match.fun(divergence)
}

