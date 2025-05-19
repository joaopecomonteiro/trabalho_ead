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



# KMeans


# Assuming df_clean is already prepared
numeric_features <- df_clean[, sapply(df_clean, is.numeric)]

# Remove the target (if still present)
numeric_features <- numeric_features[, !names(numeric_features) %in% "disease_burden_class"]

# Scale the features
scaled_data <- scale(numeric_features)

wss <- vector()
for (k in 1:20) {
  set.seed(123)
  km <- kmeans(scaled_data, centers = k, nstart = 10)
  wss[k] <- km$tot.withinss
}

# Plot the elbow
plot(1:20, wss, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters K",
     ylab = "Total Within-Cluster Sum of Squares (WSS)",
     main = "Elbow Method for Optimal K")


library(cluster)

avg_sil <- vector()
for (k in 2:10) {
  set.seed(123)
  km <- kmeans(scaled_data, centers = k, nstart = 10)
  sil <- silhouette(km$cluster, dist(scaled_data))
  avg_sil[k] <- mean(sil[, 3])
}

# Plot silhouette scores
plot(2:10, avg_sil[2:10], type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters K",
     ylab = "Average Silhouette Width",
     main = "Silhouette Method for Optimal K")


best_k <- which.max(avg_sil)
cat("Best number of clusters (k) by silhouette:", best_k, "\n")



# Kmeans


set.seed(123)
kmeans_result <- kmeans(scaled_data, centers = 2, nstart = 25)

# Add cluster labels
df_clean$KMeans_Cluster <- as.factor(kmeans_result$cluster)



# Hierarchical

dist_matrix <- dist(scaled_data)

# Perform hierarchical clustering using Ward's method
hc <- hclust(dist_matrix, method = "ward.D2")

# Cut the dendrogram to get 2 clusters
df_clean$HClust_Cluster <- as.factor(cutree(hc, k = 2))


km_sil <- silhouette(as.integer(df_clean$KMeans_Cluster), dist(scaled_data))
cat("K-Means avg silhouette:", mean(km_sil[, 3]), "\n")

# Hierarchical clustering silhouette
hc_sil <- silhouette(as.integer(df_clean$HClust_Cluster), dist(scaled_data))
cat("Hierarchical avg silhouette:", mean(hc_sil[, 3]), "\n")

