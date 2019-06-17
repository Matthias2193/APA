#Double Machine Learning for Multiple Treatments

install.packages('devtools')
library(devtools)
install_github(repo="MCKnaus/dmlmt")
library(dmlmt)

Y <- email$conversion
D <- c(2,0,1)[as.numeric(email$segment)]
X <- model.matrix(~ -1 + recency + history_segment + history + mens + womens + zip_code + newbie +
                    channel, data = email)

stand_pl_mult <- dmlmt(X,D,Y)


se_rules <- c(-1,-.5,.5,1)

# Binary
ext_pl_bin <- dmlmt(X,D,Y,se_rule=se_rules,w=TRUE)

# Example how to plot the results
df <- data.frame(SE_rule = factor(colnames(ext_pl_bin$SE_rule[[1]])
                                  ,levels = colnames(ext_pl_bin$SE_rule[[1]]))
                 ,coef = ext_pl_bin$SE_rule[[1]][1,],se = ext_pl_bin$SE_rule[[2]][1,])
ggplot(df, aes(SE_rule, coef, ymin = coef-se, ymax = coef+se)) +
  geom_errorbar() + geom_point()

# Example how to check balancing with the package of your choice, e.g. cobalt
library(cobalt)
balance <- bal.tab(as.data.frame(X), treat = D,weights=ext_pl_bin$weights,method = "weighting",
                   s.d.denom = "pooled", disp.v.ratio = TRUE, disp.ks = TRUE, un = TRUE)
love.plot(balance,abs = TRUE, line=TRUE, var.order="unadjusted")


library(grf)

# Initialize nuisance matrices
values <- sort(unique(D))
ps_mat <- t_mat <- y_mat <- matrix(NA,length(Y),length(values))

# Get nuisance parameter predictions
for (tr in 1:length(values)){
  t_mat[,tr] <- as.numeric(D == values[tr])
  rf_p <- regression_forest(X,t_mat[,tr])
  ps_mat[,tr] <- predict(rf_p, X)$predictions
  rf_y <- regression_forest(X[t_mat[,tr] == 1,],Y[t_mat[,tr] == 1])
  y_mat[,tr] <- predict(rf_y, X)$predictions
}

# Calculate generalized p-score and enforce common support
rf_gps <- gps_cs(ps_mat,t_mat)

# Potential outcomes
rf_PO <- PO_dmlmt(t_mat,Y,y_mat,rf_gps$p,cs_i=rf_gps$cs)
# ATE
rf_ATE <- TE_dmlmt(rf_PO$mu,rf_gps$cs)

