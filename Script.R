# Load necessary libraries
library(tidyverse)
library(cluster)
library(factoextra)
library(caret)  


# Load the dataset
data <- read.csv("/home/joaomonteiro/Desktop/EAD/Trabalho2/water_pollution_disease.csv")

# View structure and summary
str(data)
summary(data)

X <- data[, -12]

dummies <- dummyVars("~ .", data = X)
data_encoded <- predict(dummies, newdata = X)
data_encoded <- as.data.frame(data_encoded)

data_scaled <- scale(data_encoded)



# Elbow Method 

wss <- function(k) {
  kmeans(data_scaled, centers = k, nstart = 10)$tot.withinss
}

k.values <- 1:20
wss_values <- sapply(k.values, wss)

plot(k.values, wss_values, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares")


# Silhouette

silhouette_score <- function(k){
  km <- kmeans(data_scaled, centers = k, nstart = 10)
  ss <- silhouette(km$cluster, dist(data_scaled))
  mean(ss[, 3])
}

sil_scores <- sapply(2:20, silhouette_score)

plot(1:19, sil_scores, type = "b", pch = 19,
     xlab = "Number of clusters K",
     ylab = "Average Silhouette Score")



best_k <- k.values[which.max(sil_scores)]
cat("Best k according to silhouette method is:", best_k, "\n")



# Gap Stat
gap_stat <- clusGap(data_scaled, FUN = kmeans, nstart = 25,
                    K.max = 15, B = 5)
print(gap_stat, method = "firstmax")
plot(gap_stat)


# Kmeans
k <- 10
km_result <- kmeans(data_scaled, centers = k, nstart = 25)



# Hierarchical clustering (complete linkage)
dist_matrix <- dist(data_scaled)

hc <- hclust(dist_matrix, method = "complete")

plot(hc, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")

data$hclust_cluster <- cutree(hc, k = 10)

rect.hclust(hc, k = 10, border = 2:11)  # adds color-coded cluster boxes






# LDA


library(MASS)      # For lda()
library(caret)     # For data splitting and evaluation
library(ggplot2)   # For visualization


# Read the dataset
data <- read.csv("/home/joaomonteiro/Desktop/EAD/Trabalho2/water_pollution_disease.csv")

# Check structure
str(data)


table(data$Water.Treatment.Method)

# Ensure 'Water Treatment Method' is a factor (target variable)
data$Water.Treatment.Method <- as.factor(data$Water.Treatment.Method)

# Remove non-predictive variables (e.g., Country, Region, Year)
data_clean <- subset(data, select = -c(Country, Region, Year))

# Remove rows with missing values
data_clean <- na.omit(data_clean)

# Split the data into training and testing sets
set.seed(123)
trainIndex <- createDataPartition(data_clean$Water.Treatment.Method, p = 0.8, list = FALSE)
trainData <- data_clean[trainIndex, ]
testData <- data_clean[-trainIndex, ]

preProc <- preProcess(trainData[, -which(names(trainData) == "Water.Treatment.Method")], method = c("center", "scale"))

train_scaled <- predict(preProc, trainData)
test_scaled <- predict(preProc, testData)

# Apply LDA
lda_model <- lda(Water.Treatment.Method ~ ., data = train_scaled)

# View model details
print(lda_model)

# Predict on test data
pred <- predict(lda_model, newdata = test_scaled)

# Confusion matrix
confusionMatrix(pred$class, test_scaled$Water.Treatment.Method)


# Plot LDA
lda_df <- data.frame(pred$x, Water.Treatment.Method = testData$Water.Treatment.Method)
ggplot(lda_df, aes(x = LD1, y = LD2, color = Water.Treatment.Method)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "LDA: Water Treatment Method Classification")


# QDA

qda_model <- qda(Water.Treatment.Method ~ ., data = train_scaled)

pred_qda <- predict(qda_model, newdata = test_scaled)

# Evaluate performance
confusionMatrix(pred_qda$class, test_scaled$Water.Treatment.Method)


# Random Forest

library(randomForest)

# Fit a Random Forest
rf_model <- randomForest(Water.Treatment.Method ~ ., data = trainData)

# Predict
rf_pred <- predict(rf_model, newdata = testData)

# Evaluate
confusionMatrix(rf_pred, testData$Water.Treatment.Method)











# Com nova target


library(MASS)      # For qda
library(caret)     # For data split and evaluation
library(dplyr)     # For data manipulation
library(ggplot2)

df <- read.csv("/home/joaomonteiro/Desktop/EAD/Trabalho2/water_pollution_disease.csv")

df$disease_burden <- df$Diarrheal.Cases.per.100.000.people +
  df$Cholera.Cases.per.100.000.people +
  df$Typhoid.Cases.per.100.000.people

ggplot(df, aes(x = disease_burden)) +
  geom_histogram(binwidth = 50, fill = "steelblue", color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Disease Burden",
       x = "Disease Burden (per 100,000 people)",
       y = "Frequency")


df$disease_burden_class <- cut(df$disease_burden,
                               breaks = c(0, 150, 300, 450, Inf),
                               labels = c("Very Low", "Low", "High", "Very High"),
                               include.lowest = TRUE)

df$disease_burden_class <- as.factor(df$disease_burden_class)


df_clean <- df[, !(names(df) %in% "disease_burden")]

df_clean <- na.omit(df_clean)

df_clean$Country <- as.factor(df_clean$Country)
df_clean$Region <- as.factor(df_clean$Region)
df_clean$Water.Source.Type <- as.factor(df_clean$Water.Source.Type)
df_clean$Water.Treatment.Method <- as.factor(df_clean$Water.Treatment.Method)

table(df_clean$disease_burden_class)

set.seed(123)
trainIndex <- createDataPartition(df_clean$disease_burden_class, p = 0.8, list = FALSE)
trainData <- df_clean[trainIndex, ]
testData <- df_clean[-trainIndex, ]


# LDA

lda_model <- lda(disease_burden_class ~ ., data = trainData)

# Step 9: Predict on test set
lda_pred <- predict(lda_model, newdata = testData)

# Step 10: Evaluate performance
confusionMatrix(lda_pred$class, testData$disease_burden_class)

# Optional: Visualize the result
lda_df <- data.frame(lda_pred$x, Class = testData$disease_burden_class)
ggplot(lda_df, aes(x = LD1, y = LD2, color = Class)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "LDA Projection: Disease Burden Classes")


# QDA


qda_model <- qda(disease_burden_class ~ ., data = trainData)

# Predict
qda_pred <- predict(qda_model, newdata = testData)

# Evaluate performance
confusionMatrix(qda_pred$class, testData$disease_burden_class)



