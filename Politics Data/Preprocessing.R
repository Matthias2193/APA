#Preprocessing Script

library("mlr")
library("tidyr")
library("dplyr")

#####Preprocessing#####

setwd("/Users/tobychelton/Documents/MEMS/Predictive Analytics")

#read in individuals data
individuals <- read.csv("individuals.csv")

#convert past votes into dummy coding. Convert sex into factor. Create variable age.
individuals <- individuals%>%
  mutate(g2000=ifelse(g2000=="yes",1,0),
         g2002=ifelse(g2002=="yes",1,0),
         g2004=ifelse(g2004=="yes",1,0),
         p2000=ifelse(p2000=="yes",1,0),
         p2002=ifelse(p2002=="yes",1,0),
         p2004=ifelse(p2004=="Yes",1,0),
         female=ifelse(sex=="female",1,0),
         age=2006-yob
  )

#each level of "treatment" has a white space at the beginning. Remove it.
individuals$treatment <- as.factor(ifelse(individuals$treatment==" Control", "Control",
                                      ifelse(individuals$treatment==" Civic Duty", "Civic Duty",
                                          ifelse(individuals$treatment==" Hawthorne", "Hawthorne",
                                              ifelse(individuals$treatment==" Self", "Self", "Neighbors"
                                                  )))))


#drop unknown variables as well as yob
#also dropping g2004 because it's the same for every individual
individuals<-individuals[,c("female","age","g2000","g2002","p2000","p2002", "p2004", "treatment", "cluster", "voted", "hh_id", "hh_size")]

#convert treatment into 5 dummy columns and give these columns handy names
treatment_dummies <- createDummyFeatures(individuals[,"treatment"])
colnames(treatment_dummies) <- paste("treatment", colnames(treatment_dummies), sep="_")
colnames(treatment_dummies)[colnames(treatment_dummies)=="treatment_Civic Duty"] <- "treatment_CivicDuty"
individuals <- cbind(individuals, treatment_dummies)
#drop the original treatment variable
individuals <- individuals[,-which(colnames(individuals)=="treatment")]





#####Feature Engineering#####

#create new variable: sum of votes in the past
individuals$sum_votes <- individuals$g2000 + individuals$g2002 + individuals$p2000 + individuals$p2002 + individuals$p2004


#create new variables: mean by household

#mean female by household
mean_female_df <- aggregate(female~hh_id, data=individuals, mean)
colnames(mean_female_df)[colnames(mean_female_df)=="female"] <- "mean_female_hh"
individuals <- merge(individuals, mean_female_df, by="hh_id")

#mean g2000 by household
mean_g2000_df <- aggregate(g2000~hh_id, data=individuals, mean)
colnames(mean_g2000_df)[colnames(mean_g2000_df)=="g2000"] <- "mean_g2000_hh"
individuals <- merge(individuals, mean_g2000_df, by="hh_id")

#mean g2002 by household
mean_g2002_df <- aggregate(g2002~hh_id, data=individuals, mean)
colnames(mean_g2002_df)[colnames(mean_g2002_df)=="g2002"] <- "mean_g2002_hh"
individuals <- merge(individuals, mean_g2002_df, by="hh_id")

#mean p2000 by household
mean_p2000_df <- aggregate(p2000~hh_id, data=individuals, mean)
colnames(mean_p2000_df)[colnames(mean_p2000_df)=="p2000"] <- "mean_p2000_hh"
individuals <- merge(individuals, mean_p2000_df, by="hh_id")

#mean p2002 by household
mean_p2002_df <- aggregate(p2002~hh_id, data=individuals, mean)
colnames(mean_p2002_df)[colnames(mean_p2002_df)=="p2002"] <- "mean_p2002_hh"
individuals <- merge(individuals, mean_p2002_df, by="hh_id")

#mean p2004 by household
mean_p2004_df <- aggregate(p2004~hh_id, data=individuals, mean)
colnames(mean_p2004_df)[colnames(mean_p2004_df)=="p2004"] <- "mean_p2004_hh"
individuals <- merge(individuals, mean_p2004_df, by="hh_id")

#mean age by household
mean_age_df <- aggregate(age~hh_id, data=individuals, mean)
colnames(mean_age_df)[colnames(mean_age_df)=="age"] <- "mean_age_hh"
individuals <- merge(individuals, mean_age_df, by="hh_id")

#mean sum of votes by household
mean_sum_votes_df <- aggregate(sum_votes~hh_id, data=individuals, mean)
colnames(mean_sum_votes_df)[colnames(mean_sum_votes_df)=="sum_votes"] <- "mean_sum_votes_hh"
individuals <- merge(individuals, mean_sum_votes_df, by="hh_id")




#remove all objects but the individuals data frame
rm(list = ls()[ls() != "individuals"])