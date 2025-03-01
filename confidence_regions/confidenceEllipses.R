'############ RESULTS BUILD ACCORDING TO MY PYTHON IMPLAMANTATION

#rm(list = ls())
rct_no_het <- fread("rct_results_no_heterogeneity.csv")
rct_no_het$V1 <- NULL


rct_het <- fread("rct_results.csv")
rct_het$V1 <- NULL
#library(RColorBrewer)
#colors = brewer.pal(8, "Dark2")


library(ggplot2)

# Your data

# Calculate confidence intervals
rct_het$lower_ci <- rct_het$ate1 - 1.96 * rct_het$std_1
rct_het$upper_ci <- rct_het$ate1 + 1.96 * rct_het$std_1


ggplot(rct_het, aes(x = theta_11, y = theta_12)) +
  geom_point(alpha = 0.8)  + theme_classic()  +
  geom_point(aes(x = ate1, y = ate2), alpha = 0.8, color = "red") 




ggplot(df, aes(x = x, y = y, color = group)) +
  geom_point() +
  stat_ellipse()


###########


library(car)

# Set seed for reproducibility
set.seed(123)

# Example data
mean_vector <- c(1, 2)
covariance_matrix <- matrix(c(0.2, 0.3, 0.3, 0.5), nrow = 2, byrow = TRUE)
covariance_matrix_shape <- matrix(c(0.2, 0, 0, 0.5), nrow = 2, byrow = TRUE)

# Confidence level
confidence_level <- 0.95

# Create ellipse data
ellipse_data <- car::ellipse(
  shape = covariance_matrix_shape,
  center = mean_vector,
  level = confidence_level, 
  radius = sqrt(qchisq(confidence_level, df = 2)),
)

# Scatter plot of data points
data_points <- MASS::mvrnorm(n = 100, mu = mean_vector, Sigma = covariance_matrix)
plot(data_points, pch = 16, col = "blue", xlab = "X-axis", ylab = "Y-axis", main = "Scatter Plot with Confidence Ellipse")

# Draw confidence ellipse
lines(ellipse_data, col = "red", lty = 2)

# Add legend
legend("topright", legend = c("Data Points", "95% Confidence Ellipse"), col = c("blue", "red"), lty = c(0, 2), pch = c(16, NA))

# Show the plot
# Scatter plot with ggplot2
df <- as.data.frame(data_points)

ggplot(df, aes(x = V1, y = V2)) +
  geom_point(color = "blue") +
  geom_path(data = as.data.frame(ellipse_data), aes(x, y), color = "red", linetype = "dashed") +
  labs(x = "X-axis", y = "Y-axis", title = "Scatter Plot with Confidence Ellipse") +
  theme_classic() +
  theme(legend.position = "top") +
  guides(color = guide_legend(title = "Legend"))

'

library(data.table)
library(ggplot2)


############ RESULTS BUILD ACCORDING GRF Multi_armed_causal_forest 
############## VISUALIZE NETWORKS ######################

#install.packages('igraph')
library(igraph)
#install.packages('network')
#install.packages('sna')
library(sna)
#install.packages('ndtv')
library(ndtv)
#install.packages('visNetwork')
library(visNetwork)


net.bg <- sample_pa(80)
V(net.bg)$size <- 8
V(net.bg)$frame.color <- "white"
V(net.bg)$color <- "orange"
V(net.bg)$label <- ""
E(net.bg)$arrow.mode <- 0
plot(net.bg)


l <- layout_with_fr(net.bg)
plot(net.bg, layout=l)


plot(net.bg, layout=layout_randomly)
par(mfrow=c(2,2), mar=c(0,0,0,0)) # plot four figures - 2 rows, 2 columns
plot(net.bg, layout=layout_with_fr, vertex.color="#4863A0", vertex.size = 10)
plot(net.bg, layout=layout_with_fr, vertex.color="#307D7E", vertex.size = 10)
plot(net.bg, layout=layout_with_fr, vertex.color="#B3446C", vertex.size = 10)
plot(net.bg, layout=l, vertex.color = "#BAB86C", vertex.size = 10)





### Build a function that runs 1000 times for a single tree
### that becomes equivalent to the results of 1000 trees
simulate_data <- function(N, cov_dim){
  set.seed(0)
  seeds <- seq(10, 1000, 10)  
  
  
  p = cov_dim
  # generate a network from the BA model with alpha = 1, N = 1000+2, m = 1
  net <- generate_BA(N = 1000+2, multiple_node = 1, alpha = 1)
  str(net)
  #plot(net)
  
  
  ## basically we have a graph which contains from_id (id of the source) 
  # to to_id (id of the destination)
  #head(net$graph) 
  
  ## we now can calculate a distance matrix - i.e., how many nodes it 
  # takes to get from_id to to_id (so we can measure the closure)
  # the highest correlated ones are with the less distance.
  
  net_data <- as.data.frame(net$graph)
  net_data$dist <- net_data$V1 - net_data$V2
  
  eye <- diag(p)
  tri <- c(net_data$dist)
  length(tri)
  
  b = matrix(0, p, p) 
  b[lower.tri(b, diag=FALSE)] <- tri
  b[upper.tri(b, diag=FALSE)] <- tri
  b <- b + eye
  covar <- matrix(0, p, p)
  
  
  ## now convert distances to a covariance matrix (Yunchuan Kong & TianweiYu, 2018)
  for (i in 1:nrow(b)){
    
    for (j in 1:ncol(b)){
      
      covar[i, j] <- 0.5^b[i, j]
      
    }
  }
  
  print(covar)
  # Extract diagonal elements
  diagonal_elements <- diag(covar)
  
  #  Replace diagonal elements with zones
  diagonal_elements[] <- rep(1,  length(diagonal_elements))
  
  #  Assign modified diagonal elements back to the matrix
  covar <- `diag<-`(covar, diagonal_elements)
  
  ### Now simulate the data
  #N <- 500 #number of obs
  NX <- cov_dim # number of covariates
  
  ## instruments
  mean_X <- rep(0, len = NX)
  #var_X <- crossprod(covar, covar*(NX:1))
  var_X <- covar + diag(3)*0.5

  
  X <- as.matrix(MASS::mvrnorm(N, c(mean_X), var_X))
  
  
  
  #set.seed(30)
  P = rbinom(N, 1, 0.5)
  
  
  ## simulate outcome
  y1 <- 0.5 +  P*(X[, 2]+ 0.5*X[, 1]) + rnorm(N, 0, 1)
  y2 <- 0.5 + P*(X[, 3] + 0.5*X[, 1]) + rnorm(N, 0, 1)
  print(cov(X[, 2]+ 0.5*X[, 1], X[, 3] + 0.5*X[, 1]))
  
  data = as.data.frame(cbind(y1, y2, P, X))
  
  return(data)
  
}



sim_data <- simulate_data(N = 1000, cov_dim = 3)
head(sim_data)


########### Construct a multivariate causal forest that saves 1000 results for
#### a single tree grown on the same data
Y <- as.matrix(sim_data[, c(1, 2)])
W <- as.factor(sim_data[, 3])
X <- as.matrix(sim_data[, c(4:ncol(sim_data))])


predictions1 <- list()
predictions2 <- list()

library(grf)
## this is as if bootstrapped. O
for (i in 1:1000){
  model <- multi_arm_causal_forest(X, Y, W, num.trees = 1, seed = i)
  predictions1[[i]] <- model$predictions[, 1]
  predictions2[[i]] <- model$predictions[, 2]
  
}




# Assuming predictions1 and predictions2 are your bootstrapped lists of vectors
# Each list element is a vector of length N
predictions1_mean <- colMeans(do.call(rbind, predictions1), na.rm = T)
predictions2_mean <- colMeans(do.call(rbind, predictions2), na.rm = T)


temp_list1_var <- list()
temp_list2_var <- list()
temp_list3_cov <- list()


for (i in 1:length(predictions1)){
  
  temp_list1_var[[i]] <- (predictions1[[i]] - predictions1_mean)^2
  temp_list2_var[[i]] <- (predictions2[[i]] - predictions2_mean)^2
  temp_list3_cov[[i]] <- (predictions1[[i]] - predictions1_mean)*(predictions2[[i]] - predictions2_mean)
  
}



predictions1_var <- colMeans(do.call(rbind, temp_list1_var), na.rm = T)
predictions2_var <- colMeans(do.call(rbind, temp_list2_var), na.rm = T)
bootstrapped_cov <- colMeans(do.call(rbind, temp_list3_cov), na.rm = T)


mean(round(bootstrapped_cov, 3)) ### real covariance is 0.125, while estimated mean cov 0.139
round(bootstrapped_cov, 3) ### real covariance is 0.125, while estimated mean cov 0.139




library(ggplot2)
ggplot(data = sim_data, aes(x = X[, 2] + 0.5 * X[, 1], y = X[, 3] + 0.5 * X[, 1])) +
  geom_point(aes(color = "simulated"), alpha = 0.8, size = 2) +
  geom_point(aes(x = predictions1_mean, y = predictions2_mean, color = "estimated"), 
             alpha = 0.8, size = 2) +
  theme_classic() +
  scale_color_manual(values = c("simulated" = "#666362", "estimated" = "#87AFC7"), name = "policy effect") +
  labs(x = "policy effect (first outcome)", y = "policy effect (second outcome)") +
  theme(axis.text = element_text(size = 12))  # Adjust the size (14 is just an example, you can change it)




mean(predictions1_var)
mean(data_final$bootstrapped_cov)



data_final <- as.data.frame(cbind(ate1_true  = X[, 2]+ 0.5*X[, 1], 
                                  ate2_true =  X[, 3]+ 0.5*X[, 1], 
                                  predictions1_mean, predictions2_mean, 
                                  predictions1_var, predictions2_var, bootstrapped_cov))



data_final$cov_binary <- ifelse(data_final$bootstrapped_cov <= 0, "non-positive", "positive")


cov(data_final$ate1_true, data_final$ate2_true)
cov(data_final$predictions1_mean, data_final$predictions2_mean)


head(data_final)
table(data_final$cov_binary)


mean((data_final$ate1_true - data_final$predictions1_mean)^2)
mean((data_final$ate2_true - data_final$predictions2_mean)^2)


library(data.table)
data_final <- as.data.table(data_final)
data_final[, .(covariance = mean(bootstrapped_cov), 
               var_pred1 = mean(predictions1_var), 
               var_pred2 = mean(predictions2_var)), .(cov_binary)]



data_final$var_mean1 <- ifelse(data_final$cov_binary == "positive", 0.3640521, 0.2966197 )
data_final$var_mean2 <- ifelse(data_final$cov_binary == "positive", 0.3655752, 0.3226346 )
data_final$cov_mean <- ifelse(data_final$cov_binary == "positive", 0.08639157 , -0.01717121 )




############ CONFIDENCE ELLIPSES #######################################
# Example dataframe (replace this with your actual data)
library(car)

# Set seed for reproducibility
set.seed(123)
par(mfrow=c(1,1), mar=c(0,0,0,0)) # plot four figures - 2 rows, 2 colum


mean_vector <- c(mean(data_final$predictions1_mean), mean(data_final$predictions2_mean))


covariance_matrix <- matrix(c(1.709836, 0.7723078, 0.7723078, 1.086786), nrow = 2, byrow = TRUE)
covariance_matrix_shape <- matrix(c(1.709836, 0, 0, 1.086786), nrow = 2, byrow = TRUE)

# Confidence level
confidence_level <- 0.95


# Create ellipse data
ellipse_data <- car::ellipse(
  shape = covariance_matrix,
  center = mean_vector,
  level = confidence_level, 
  radius = sqrt(qchisq(confidence_level, df = 2)),
)

# Show the plot
# Scatter plot with ggplot2
#install.packages("sp")
library(sp)
data_final$inside_ellipse <- ifelse(point.in.polygon(data_final$predictions1_mean, 
                                                     data_final$predictions2_mean,
                                                     ellipse_data[, 1], ellipse_data[, 2]) == 1, "Inside", "Outside")
table(data_final$inside_ellipse)

# Define colors
inside_color <- "#666362"  # Color for points inside the ellipse
outside_color <- "#E55451"  # Color for points outside the ellipse


ggplot(data_final, aes(x = predictions1_mean, y = predictions2_mean, color = inside_ellipse)) +
  geom_point(alpha = 0.8) +
  geom_path(data = as.data.frame(ellipse_data), aes(x, y), 
            linetype = "dashed", size = 1.5, color = inside_color) +
  scale_color_manual(values = c("Inside" = inside_color, "Outside" = outside_color)) +
  labs(x = "policy effect (first outcome)", 
       y = "policy effect (second outcome)", 
       title = "") +
  theme_classic() +
  theme(legend.position = "top") +
  guides(color = guide_legend(title = "Estimates")) + 
  theme(axis.text = element_text(size = 15), 
        axis.title=element_text(size=15), 
        legend.text=element_text(size=15), 
        legend.title = element_text(size = 15))  # Adjust the size




'ggplot(data_final, aes(x = predictions1_mean, y = predictions2_mean)) +
  geom_point(color = "#666362") +
  geom_path(data = as.data.frame(ellipse_data), aes(x, y), color = "#728FCE", 
            linetype = "dashed", size = 1.5) +
  labs(x = "policy effect (first outcome)", 
       y = "policy effect (second outcome)", 
       title = "") +
  theme_classic() +
  theme(legend.position = "top") +
  guides(color = guide_legend(title = "Legend")) + 
  theme(axis.text = element_text(size = 12))  # Adjust the size (14 is just an example, you can change it)


'




########### NEGATIVE ####
# Assuming you have a dataframe named df with columns group, V1, and V2

# Compute means, variances, and covariances for each group
group_nonpositive <- data_final[data_final$cov_binary == "non-positive"]

mean_vector <- c(mean(group_nonpositive$predictions1_mean), mean(group_nonpositive$predictions2_mean))
covariance_matrix <- matrix(c(0.2966197, -0.01717121, -0.01717121, 1.086786), nrow = 2, byrow = TRUE)



confidence_level <- 0.95

# Create ellipse data
ellipse_data_1 <- car::ellipse(
  shape = covariance_matrix,
  center = mean_vector,
  radius = sqrt(qchisq(confidence_level, df = 2)),
)

# Show the plot
# Scatter plot with ggplot2
'
ggplot(group_nonpositive, aes(x = predictions1_mean, y = predictions2_mean)) +
  geom_point(color = "blue") +
  geom_path(data = as.data.frame(ellipse_data_1), aes(x, y), color = "red", linetype = "dashed") +
  labs(x = "X-axis", y = "Y-axis", title = "Scatter Plot with Confidence Ellipse") +
  theme_classic() +
  theme(legend.position = "top") +
  guides(color = guide_legend(title = "Legend"))

'


########## POSITIVE ################### 
group_positive <- data_final[data_final$cov_binary == "positive"]

mean_vector <- c(mean(group_positive$predictions1_mean), 
                 mean(group_positive$predictions2_mean))
covariance_matrix <- matrix(c(group_positive$var_mean1[1], 
                              group_positive$cov_mean[1], group_positive$cov_mean[1], 
                              group_positive$var_mean2[1]), nrow = 2, byrow = TRUE)



confidence_level <- 0.95

# Create ellipse data
ellipse_data_2 <- car::ellipse(
  shape = covariance_matrix,
  center = mean_vector,
  radius = sqrt(qchisq(confidence_level, df = 2)),
)




#########
# Combine dataframes for both groups
library(dplyr)
combined_data <- rbind(group_nonpositive, group_positive)
head(ellipse_data_1)
ellipse_data_1_new <- as.data.frame(cbind(ellipse_data_1, 
                                          rep("non-positive", 
                                              length(ellipse_data_1))))
ellipse_data_2_new <- as.data.frame(cbind(ellipse_data_2, 
                                          rep("positive", 
                                              length(ellipse_data_2))))

ellipse_combined <-rbind(ellipse_data_1_new, ellipse_data_2_new)
str(ellipse_combined)
ellipse_combined$x <- as.numeric(ellipse_combined$x )
ellipse_combined$y <- as.numeric(ellipse_combined$y )


# Create a new variable for grouping colors
combined_data$group_color <- as.character(combined_data$cov_binary, 
                                    levels = c("non-positive", "positive"), 
                                    labels = c("blue", "green"))

cov_binary_colors <- c("positive" = "#52595D", "non-positive" = "#990012")

# Plot the combined data with confidence ellipses
ggplot(combined_data, aes(x = predictions1_mean, 
                          y = predictions2_mean, color = group_color)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_path(data = ellipse_combined, aes(x, y, color = V3), 
            linetype = "dashed", size = 2) +
  labs(x = "policy effect (first outcome)", y = "policy effect (second outcome)", title = "") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_manual(name = "Legend", values = cov_binary_colors) +
  guides(color = guide_legend(title = "covariance")) +
  theme(legend.position = "right") +
  theme(axis.text = element_text(size = 12))  # Adjust the size (14 is just an example, you can change it)


  

#666362", "estimated" = "#728FCE




############ ASSUMING ZERO COVARIANCE
########### NEGATIVE ####
# Assuming you have a dataframe named df with columns group, V1, and V2

# Compute means, variances, and covariances for each group
group_nonpositive <- data_final[data_final$cov_binary == "non-positive"]

mean_vector <- c(mean(group_nonpositive$predictions1_mean), mean(group_nonpositive$predictions2_mean))
covariance_matrix <- matrix(c(0.2966197, 0, 0, 1.086786), nrow = 2, byrow = TRUE)



confidence_level <- 0.95

# Create ellipse data
ellipse_data_1 <- car::ellipse(
  shape = covariance_matrix,
  center = mean_vector,
  radius = sqrt(qchisq(confidence_level, df = 2)),
)

# Show the plot
# Scatter plot with ggplot2
'
ggplot(group_nonpositive, aes(x = predictions1_mean, y = predictions2_mean)) +
  geom_point(color = "blue") +
  geom_path(data = as.data.frame(ellipse_data_1), aes(x, y), color = "red", linetype = "dashed") +
  labs(x = "X-axis", y = "Y-axis", title = "Scatter Plot with Confidence Ellipse") +
  theme_classic() +
  theme(legend.position = "top") +
  guides(color = guide_legend(title = "Legend"))

'


########## POSITIVE ################### 
group_positive <- data_final[data_final$cov_binary == "positive"]

mean_vector <- c(mean(group_positive$predictions1_mean), 
                 mean(group_positive$predictions2_mean))
covariance_matrix <- matrix(c(group_positive$var_mean1[1], 
                              0, 0, 
                              group_positive$var_mean2[1]), nrow = 2, byrow = TRUE)



confidence_level <- 0.95

# Create ellipse data
ellipse_data_2 <- car::ellipse(
  shape = covariance_matrix,
  center = mean_vector,
  radius = sqrt(qchisq(confidence_level, df = 2)),
)




#########
# Combine dataframes for both groups
library(dplyr)
combined_data <- rbind(group_nonpositive, group_positive)
head(ellipse_data_1)
ellipse_data_1_new <- as.data.frame(cbind(ellipse_data_1, 
                                          rep("non-positive", 
                                              length(ellipse_data_1))))
ellipse_data_2_new <- as.data.frame(cbind(ellipse_data_2, 
                                          rep("positive", 
                                              length(ellipse_data_2))))

ellipse_combined <-rbind(ellipse_data_1_new, ellipse_data_2_new)
str(ellipse_combined)
ellipse_combined$x <- as.numeric(ellipse_combined$x )
ellipse_combined$y <- as.numeric(ellipse_combined$y )


# Create a new variable for grouping colors
combined_data$group_color <- as.character(combined_data$cov_binary, 
                                          levels = c("non-positive", "positive"), 
                                          labels = c("blue", "green"))

cov_binary_colors <- c("positive" = "#52595D", "non-positive" = "#990012")

# Plot the combined data with confidence ellipses
ggplot(combined_data, aes(x = predictions1_mean, 
                          y = predictions2_mean, color = group_color)) +
  geom_point(alpha = 0.4, size = 2) +
  geom_path(data = ellipse_combined, aes(x, y, color = V3), 
            linetype = "dashed", size = 2) +
  labs(x = "policy effect (first outcome)", y = "policy effect (second outcome)", title = "") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_manual(name = "Legend", values = cov_binary_colors) +
  guides(color = guide_legend(title = "covariance")) +
  theme(legend.position = "right") +
  theme(axis.text = element_text(size = 12))  # Adjust the size (14 is just an example, you can change it)




########### Hypothesis Testing for one particular observation:
library(car)

# Set seed for reproducibility
set.seed(123)
par(mfrow=c(1,1), mar=c(0,0,0,0)) # plot four figures - 2 rows, 2 colum


mean_vector <- c(1.64, -0.140)


covariance_matrix <- matrix(c(0.596, 0.445, 0.445, 0.839), nrow = 2, byrow = TRUE)
covariance_matrix_shape <- matrix(c(0.596, 0, 0, 0.839), nrow = 2, byrow = TRUE)

# Confidence level
confidence_level <- 0.90

library(MASS)
df <- as.data.frame(mvrnorm(n = nrow(data_final), 
                            mu = mean_vector, Sigma = covariance_matrix))

# Create ellipse data
ellipse_data_cov <- car::ellipse(
  shape = covariance_matrix,
  center = mean_vector,
  level = confidence_level, 
  radius = sqrt(qchisq(confidence_level, df = 2)),
)

ellipse_data_no_cov <- car::ellipse(
  shape = covariance_matrix_shape,
  center = mean_vector,
  level = confidence_level, 
  radius = sqrt(qchisq(confidence_level, df = 2)),
)
# Show the plot
# Scatter plot with ggplot2

total_ellipse <- as.data.frame(rbind(cbind(ellipse_data_cov, 
                             rep("with covariance", nrow(ellipse_data_cov))), 
                       cbind(ellipse_data_no_cov, 
                             rep("without covariance", nrow(ellipse_data_no_cov)))))

total_ellipse$x <- as.numeric(total_ellipse$x)
total_ellipse$y <- as.numeric(total_ellipse$y)


cov_binary_colors <- c("with covariance" = "#4863A0", "without covariance" = "#9E4638")


ggplot(df, aes(x = V1, y = V2)) +
  geom_point(color = "#666362", alpha = 0.5, size = 2) +
  geom_path(data = total_ellipse, aes(x, y, color = V3), 
            linetype = "dashed", size = 1.5) +
  labs(x = "policy effect (first outcome)", 
       y = "policy effect (second outcome)", 
       title = "") +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_manual(name = "Covariance", values = cov_binary_colors) +
  guides(color = guide_legend(title = "covariance"))  +
  theme(axis.text = element_text(size = 12))  # Adjust the size (14 is just an example, you can change it)



