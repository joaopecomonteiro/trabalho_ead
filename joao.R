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
df <- read_csv("/home/joaomonteiro/Desktop/EAD/Trabalho2/water_pollution_disease.csv")






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
melted_cor <- melt(cor_matrix)
ggplot(melted_cor, aes(Var1, Var2, fill = value)) +
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
                                remove_first_dummy = TRUE,  # Avoid dummy variable trap
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

#df_kmeans$kmeans_cluster <- kmeans_result$cluster
kmeans_clusters <- kmeans_result$cluster

kmeans_res <- cluster.stats(dist(df_clustering), kmeans_clusters)

cat("Calinski-Harabasz Index:", kmeans_res$ch, "\n")
cat("Silhouette médio:", kmeans_res$avg.silwidth, "\n")
cat("Índice de Dunn:", kmeans_res$dunn, "\n")
cat("Within-cluster sum of squares:", kmeans_res$within.cluster.ss, "\n")
variance_explained <- kmeans_result$betweenss / (kmeans_result$tot.withinss + kmeans_result$betweenss)
cat("Variance explained by clusters:", round(variance_explained * 100, 2), "%\n")

# ---- Clustering Evaluation Summary (K = 8) ----
# - Calinski-Harabasz Index: 106.95
#   -> Indicates moderate separation between clusters, but not very strong.
#
# - Average Silhouette Width: 0.11
#   -> Very low; suggests clusters overlap and are not well-defined.
#
# - Dunn Index: 0.35
#   -> Moderate cluster compactness and separation; room for improvement.
#
# - Within-Cluster Sum of Squares (WSS): 91,153.1
#   -> High internal dispersion; clusters are not tight.
#
# - Variance Explained by Clusters: 20.01%
#   -> Only a small portion of total variance is captured by the clustering.
#
# => Overall, clustering quality is weak with k = 8.
#    Consider trying fewer clusters (e.g., k = 3–5) or other algorithms like GMM or hierarchical clustering.



#------------------------- Hclust -------------------------#

d_matrix <- dist(df_clustering)
hc <- hclust(d_matrix, method = "ward.D2")
plot(hc, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")

hc_clusters <- cutree(hc, k = optimal_k)

rect.hclust(hc, k = 8, border = 2:11) 


hc_res <- cluster.stats(dist(df_clustering), hc_clusters)

cat("Calinski-Harabasz Index:", hc_res$ch, "\n")
cat("Silhouette médio:", hc_res$avg.silwidth, "\n")
cat("Índice de Dunn:", hc_res$dunn, "\n")
cat("Within-cluster sum of squares:", hc_res$within.cluster.ss, "\n")

# ---- Clustering Evaluation Summary (Hierarchical Clustering, k = 8) ----
# - Calinski-Harabasz Index: 96.39
#   -> Indicates moderate separation between clusters, but lower than K-Means (106.95).
#
# - Average Silhouette Width: 0.097
#   -> Very low; suggests significant overlap and weak cluster cohesion.
#
# - Dunn Index: 0.34
#   -> Moderate compactness and separation; similar to K-Means (0.35).
#
# - Within-Cluster Sum of Squares (WSS): 92,991.24
#   -> Slightly higher internal dispersion than K-Means (91,153.1).
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
# - Calinski-Harabasz Index: 70.16
#   -> Lower than both Ward's method and K-Means; indicates weaker separation between clusters.
#
# - Average Silhouette Width: 0.061
#   -> Very low; clusters are poorly defined with significant overlap.
#
# - Dunn Index: 0.31
#   -> Slightly weaker than Ward’s method (0.34); moderate compactness and separation.
#
# - Within-Cluster Sum of Squares (WSS): 97,892.68
#   -> Higher internal dispersion than both K-Means (91,153.1) and Ward (92,991.2), suggesting less compact clusters.
#
# => Overall, clustering quality using complete linkage is the weakest among the tested methods.
#    All evaluation metrics point to poor cohesion and separation.




#------------------------- GMM -------------------------#


gmm_model <- Mclust(df_clustering) #Escolhe k automaticamente baseado no BIC

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
# - Calinski-Harabasz Index: 108.49
#   -> Highest among all tested methods; suggests relatively good separation between clusters.
#
# - Average Silhouette Width: 0.110
#   -> Low, but slightly better than Hierarchical (Ward: 0.097, Complete: 0.061); indicates weak-to-moderate cluster definition.
#
# - Dunn Index: 0.35
#   -> Highest among all methods; indicates relatively better compactness and separation.
#
# - Within-Cluster Sum of Squares (WSS): 88,329.62
#   -> Lowest WSS so far, suggesting more compact clusters.
#
# => Overall, GMM outperforms K-Means and Hierarchical Clustering based on CH, Dunn, and WSS.
#    Despite low silhouette width, the other metrics indicate that GMM forms the most structured and well-separated clusters in this dataset.



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
# K-Means           1.000       0.410     0.545   0.357
# HClust-Ward       0.410       1.000     0.494   0.309
# GMM               0.545       0.494     1.000   0.390
# HClust-Complete   0.357       0.309     0.390   1.000
#
# - K-Means and GMM have the highest similarity (Jaccard = 0.545), suggesting moderate alignment in cluster structure.
# - HClust-Ward is moderately similar to both K-Means (0.410) and GMM (0.494), indicating some shared structure.
# - HClust-Complete is the most divergent method overall, with the lowest similarity scores across the board.
#
# => GMM provides the most consistent clustering when compared to other methods.
#    HClust-Ward aligns better with the overall structure than HClust-Complete.


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





#------------------------- QDA -------------------------#

qda_model <- qda(Water.Treatment.Method ~ ., data = trainData)

qda_pred <- predict(qda_model, newdata = X_test)

table(Predicted = qda_pred$class, Actual = y_test)

mean(qda_pred$class == y_test)




















