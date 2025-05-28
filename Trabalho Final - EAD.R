# Load necessary packages
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(forcats)
library(stringr)
library(data.table)
library(cluster)
library(factoextra)
library(mclust)
library(fpc)
library(MASS)
library(caret) 
library(tidyverse)

# Load the data
df <- read_csv("water_pollution_disease.csv")

set.seed(123)




#------------------------- EDA -------------------------#
glimpse(df)
summary(df)
sapply(df, function(x) sum(is.na(x))) 

numeric_cols <- df %>% 
  dplyr::select(where(is.numeric)) %>% 
  dplyr::select(-Year)

# Convert to long format and plot histograms
numeric_cols %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 100, fill = "skyblue", color = "black") +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(title = "Histograms of All Numeric Variables", x = NULL, y = "Frequency")



cor_matrix <- cor(numeric_cols, use = "complete.obs")
melted_cor <- as.data.frame(as.table(cor_matrix)) %>%
  dplyr::rename(Var1 = Var1, Var2 = Var2, Correlation = Freq)
ggplot(melted_cor, aes(Var1, Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix Heatmap", x = "", y = "")



cat_df <- df %>% dplyr::select(where(~ is.character(.x) || is.factor(.x)))


cat_df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = value)) +
  geom_bar(fill = "steelblue") +
  facet_wrap(~ variable, scales = "free_x", ncol = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Bar Plots of All Categorical Variables", x = "", y = "Count")










data <- fastDummies::dummy_cols(df, 
                                select_columns = c("Country", "Region", "Water Source Type"),
                                remove_first_dummy = TRUE,  
                                remove_selected_columns = TRUE) 


num_cols <- sapply(data, is.numeric)

data <- data
data[num_cols] <- scale(data[num_cols])


#------------------------- KMeans -------------------------#

df_clustering <- data[, !(names(data) %in% "Water Treatment Method")]

set.seed(123)

wss_plot <- fviz_nbclust(df_clustering, kmeans, method = "wss") +
  labs(title = "Elbow Method for Optimal k")
print(wss_plot)

silhouette_plot <- fviz_nbclust(df_clustering, kmeans, method = "silhouette") +
  labs(title = "Silhouette Method for Optimal k")
print(silhouette_plot)



optimal_k <- 8

kmeans_result <- kmeans(df_clustering, centers = optimal_k, nstart = 25)

kmeans_clusters <- kmeans_result$cluster

kmeans_res <- cluster.stats(dist(df_clustering), kmeans_clusters)

cat("Calinski-Harabasz Index:", kmeans_res$ch, "\n")
cat("Silhouette médio:", kmeans_res$avg.silwidth, "\n")
cat("Índice de Dunn:", kmeans_res$dunn, "\n")
cat("Within-cluster sum of squares:", kmeans_res$within.cluster.ss, "\n")
variance_explained <- kmeans_result$betweenss / (kmeans_result$tot.withinss + kmeans_result$betweenss)
cat("Variance explained by clusters:", round(variance_explained * 100, 2), "%\n")

# ---- Clustering Evaluation Summary (K = 8) ----
# - Calinski-Harabasz Index: 106.73
#   -> Indicates moderate separation between clusters, but not very strong.
#
# - Average Silhouette Width: 0.11
#   -> Very low; suggests clusters overlap and are not well-defined.
#
# - Dunn Index: 0.36
#   -> Moderate cluster compactness and separation; room for improvement.
#
# - Within-Cluster Sum of Squares (WSS): 91,192
#   -> High internal dispersion; clusters are not tight.
#
# - Variance Explained by Clusters: 19.98%
#   -> Only a small portion of total variance is captured by the clustering.
#
# => Overall, clustering quality is weak with k = 8.
#    Consider trying fewer clusters (e.g., k = 3–5) or other algorithms like GMM or hierarchical clustering.


pca_result <- prcomp(df_clustering, scale. = FALSE)  
pca_data <- as.data.frame(pca_result$x[, 1:2])
pca_data$Cluster <- factor(kmeans_clusters)

ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "K-Means Clustering (PCA Visualization)",
       x = "PC 1",
       y = "PC 2") +
  scale_color_brewer(palette = "Set1")


fviz_cluster(kmeans_result, data = df_clustering,
             geom = "point", ellipse.type = "convex",
             palette = "jco", ggtheme = theme_minimal()) +
  labs(title = "K-Means Clustering Visualized")


#------------------------- Hclust -------------------------#

d_matrix <- dist(df_clustering)
hc <- hclust(d_matrix, method = "ward.D2")
plot(hc, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")

hc_clusters <- cutree(hc, k = optimal_k)

rect.hclust(hc, k = optimal_k, border = 2:11) 


hc_res <- cluster.stats(dist(df_clustering), hc_clusters)

cat("Calinski-Harabasz Index:", hc_res$ch, "\n")
cat("Silhouette médio:", hc_res$avg.silwidth, "\n")
cat("Índice de Dunn:", hc_res$dunn, "\n")
cat("Within-cluster sum of squares:", hc_res$within.cluster.ss, "\n")

# ---- Clustering Evaluation Summary (Hierarchical Clustering, k = 8) ----
# - Calinski-Harabasz Index: 100.35
#   -> Indicates moderate separation between clusters, but lower than K-Means (106.73).
#
# - Average Silhouette Width: 0.10
#   -> Very low; suggests significant overlap and weak cluster cohesion.
#
# - Dunn Index: 0.37
#   -> Moderate compactness and separation; similar to K-Means (0.36).
#
# - Within-Cluster Sum of Squares (WSS): 92,294
#   -> Slightly higher internal dispersion than K-Means (91,192).
#
# => Overall, clustering quality is weak and slightly lower than K-Means.
#    Silhouette and CH Index both suggest limited structure in the data at k = 8.






d_matrix <- dist(df_clustering)
hc_c <- hclust(d_matrix, method = "complete")
plot(hc_c, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")

hc_c_clusters <- cutree(hc_c, k = optimal_k)

rect.hclust(hc_c, k = optimal_k, border = 2:11) 


hc_c_res <- cluster.stats(dist(df_clustering), hc_c_clusters)

cat("Calinski-Harabasz Index:", hc_c_res$ch, "\n")
cat("Silhouette médio:", hc_c_res$avg.silwidth, "\n")
cat("Índice de Dunn:", hc_c_res$dunn, "\n")
cat("Within-cluster sum of squares:", hc_c_res$within.cluster.ss, "\n")


# ---- Clustering Evaluation Summary (Hierarchical Clustering - Complete Linkage, k = 8) ----
# - Calinski-Harabasz Index: 82.34
#   -> Lower than both Ward's method and K-Means; indicates weaker separation between clusters.
#
# - Average Silhouette Width: 0.07
#   -> Very low; clusters are poorly defined with significant overlap.
#
# - Dunn Index: 0.34
#   -> Slightly weaker than Ward’s method (0.37); moderate compactness and separation.
#
# - Within-Cluster Sum of Squares (WSS): 95,554
#   -> Higher internal dispersion than both K-Means (91,192) and Ward (92,294), suggesting less compact clusters.
#
# => Overall, clustering quality using complete linkage is the weakest among the tested methods.
#    All evaluation metrics point to poor cohesion and separation.




#------------------------- GMM -------------------------#


gmm_model <- Mclust(df_clustering) # automatically choose k based on BIC

summary(gmm_model)

gmm_clusters <- gmm_model$classification

d_gmm <- dist(df_clustering)

# Evaluate GMM clustering
res_gmm <- cluster.stats(d_gmm, gmm_clusters)

# Print evaluation metrics
cat("Calinski-Harabasz Index:", res_gmm$ch, "\n")
cat("Silhouette médio:", res_gmm$avg.silwidth, "\n")
cat("Índice de Dunn:", res_gmm$dunn, "\n")
cat("Within-cluster sum of squares:", res_gmm$within.cluster.ss, "\n")

# ---- Clustering Evaluation Summary (Gaussian Mixture Model - GMM) ----
# - Calinski-Harabasz Index: 108.33
#   -> Highest among all tested methods; suggests relatively good separation between clusters.
#
# - Average Silhouette Width: 0.11
#   -> Low, but slightly better than Hierarchical (Ward: 0.10, Complete: 0.07); indicates weak-to-moderate cluster definition.
#
# - Dunn Index: 0.38
#   -> Highest (slightly) among all methods; indicates relatively better compactness and separation.
#
# - Within-Cluster Sum of Squares (WSS): 88,359
#   -> Lowest WSS so far, suggesting more compact clusters.
#
# => Overall, GMM outperforms K-Means and Hierarchical Clustering based on CH, Dunn, and WSS.
#    Despite low silhouette width, the other metrics indicate that GMM forms the most structured and well-separated clusters in this dataset.



pca_gmm <- prcomp(df_clustering, scale. = FALSE) 
pca_df_gmm <- as.data.frame(pca_gmm$x[, 1:2])
pca_df_gmm$Cluster <- factor(gmm_clusters)

ggplot(pca_df_gmm, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(alpha = 0.8, size = 2) +
  theme_minimal() +
  scale_color_brewer(palette = "Set3") +
  labs(title = "GMM Clustering (PCA Visualization)",
       x = "Principal Component 1",
       y = "Principal Component 2")


pca_gmm <- prcomp(df_clustering, scale. = FALSE)
pca_df <- as.data.frame(pca_gmm$x[, 1:2])
pca_df$Uncertainty <- gmm_model$uncertainty

ggplot(pca_df, aes(x = PC1, y = PC2, color = Uncertainty)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "GMM Classification Uncertainty (PCA projection)",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Uncertainty")








jaccard_index <- function(labels1, labels2) {
  n <- length(labels1)
  a <- b <- c <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      same1 <- labels1[i] == labels1[j]
      same2 <- labels2[i] == labels2[j]
      if (same1 && same2) a <- a + 1
      else if (same1 && !same2) b <- b + 1
      else if (!same1 && same2) c <- c + 1
    }
  }
  return(a / (a + b + c))
}



jaccard_km_vs_hc <- jaccard_index(kmeans_clusters, hc_clusters)
jaccard_km_vs_gmm <- jaccard_index(kmeans_clusters, gmm_clusters)
jaccard_hc_vs_gmm <- jaccard_index(hc_clusters, gmm_clusters)

jaccard_km_vs_hc_c <- jaccard_index(kmeans_clusters, hc_c_clusters)
jaccard_hc_vs_hc_c <- jaccard_index(hc_clusters, hc_c_clusters)
jaccard_gmm_vs_hc_c <- jaccard_index(gmm_clusters, hc_c_clusters)




jaccard_matrix <- matrix(c(
  1, jaccard_km_vs_hc, jaccard_km_vs_gmm, jaccard_km_vs_hc_c,
  jaccard_km_vs_hc, 1, jaccard_hc_vs_gmm, jaccard_hc_vs_hc_c,
  jaccard_km_vs_gmm, jaccard_hc_vs_gmm, 1, jaccard_gmm_vs_hc_c,
  jaccard_km_vs_hc_c, jaccard_hc_vs_hc_c, jaccard_gmm_vs_hc_c, 1
), nrow = 4, byrow = TRUE)

rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- c("K-Means", "HClust-Ward", "GMM", "HClust-Complete")

print(round(jaccard_matrix, 3))

# ---- Jaccard Index Comparison Between Clustering Algorithms ----
#                K-Means   HClust-Ward   GMM     HClust-Complete
# K-Means           1.000       0.522     0.756   0.366
# HClust-Ward       0.522       1.000     0.636   0.378
# GMM               0.756       0.636     1.000   0.430
# HClust-Complete   0.366       0.378     0.430   1.000
#
# - K-Means and GMM have the highest similarity (Jaccard = 0.756), suggesting reasonable alignment in cluster structure.
# - HClust-Ward is moderately similar to both K-Means (0.522) and GMM (0.636), which makes sense since Ward’s method tends to create compact, spherical clusters similar to those found by K-Means
# - HClust-Ward and HClust-Complete are quite different from each other (Jaccard index = 0.378), showing that the choice of linkage method in hierarchical clustering significantly impacts the resulting cluster structure.
#
# => GMM provides the most consistent clustering when compared to other methods.
#    HClust-Ward aligns better with the overall structure than HClust-Complete.

#-------------------------------------------------------#

data$`Water Treatment Method` <- as.factor(data$`Water Treatment Method`)
names(data) <- make.names(names(data))


set.seed(123)
trainIndex <- createDataPartition(data$Water.Treatment.Method, p = 0.8, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]
y_test <- testData$Water.Treatment.Method
X_test <- testData %>%
  dplyr::select(-Water.Treatment.Method)

#------------------------- LDA -------------------------#


lda_model <- lda(Water.Treatment.Method ~ ., data = trainData)

lda_pred <- predict(lda_model, newdata = X_test)

table(Predicted = lda_pred$class, Actual = y_test)

mean(lda_pred$class == y_test)



lda_values <- predict(lda_model)

lda_df <- data.frame(lda_values$x)
lda_df$Class <- trainData$Water.Treatment.Method

ggplot(lda_df, aes(x = LD1, y = LD2, color = Class)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "LDA Projection (Training Set)",
       x = "Linear Discriminant 1",
       y = "Linear Discriminant 2") +
  scale_color_brewer(palette = "Set1")



lda_test_df <- data.frame(lda_pred$x)
lda_test_df$Class <- y_test

ggplot(lda_test_df, aes(x = LD1, y = LD2, color = Class)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "LDA Projection (Test Set)",
       x = "Linear Discriminant 1",
       y = "Linear Discriminant 2") +
  scale_color_brewer(palette = "Set1")


#------------------------- QDA -------------------------#

qda_model <- qda(Water.Treatment.Method ~ ., data = trainData)

qda_pred <- predict(qda_model, newdata = X_test)

table(Predicted = qda_pred$class, Actual = y_test)

mean(qda_pred$class == y_test)


qda_train_pred <- predict(qda_model, newdata = trainData)

qda_train_df <- data.frame(qda_train_pred$posterior)  # Usar as probabilidades de pertença
qda_train_df$Class <- trainData$Water.Treatment.Method

pca_train <- prcomp(qda_train_df[ , -ncol(qda_train_df)], center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_train$x[, 1:2])
pca_df$Class <- qda_train_df$Class

ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "QDA Projection (Training Set via Posterior Probabilities)",
       x = "PC1", y = "PC2") +
  scale_color_brewer(palette = "Set1")



qda_test_df <- data.frame(qda_pred$posterior)
qda_test_df$Class <- y_test

pca_test <- prcomp(qda_test_df[ , -ncol(qda_test_df)], center = TRUE, scale. = TRUE)
pca_df_test <- as.data.frame(pca_test$x[, 1:2])
pca_df_test$Class <- qda_test_df$Class

ggplot(pca_df_test, aes(x = PC1, y = PC2, color = Class)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "QDA Projection (Test Set via Posterior Probabilities)",
       x = "PC1", y = "PC2") +
  scale_color_brewer(palette = "Set1")

