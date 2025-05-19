library(MASS)      # For qda
library(caret)     # For data split and evaluation
library(dplyr)     # For data manipulation
library(ggplot2)
library(factoextra)
library(mclust)
library(cluster)
library(fpc)


df <- read.csv("water_pollution_disease.csv")

# Pre-Processing

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
kmeans_result <- kmeans(scaled_data, centers = best_k, nstart = 25)

# Add cluster labels
df_clean$KMeans_Cluster <- as.factor(kmeans_result$cluster)

# Visualização com PCA
pca <- prcomp(scaled_data, scale. = TRUE)
pca_data <- as.data.frame(pca$x[, 1:3])
pca_data$cluster <- as.factor(kmeans_result$cluster)

fviz_cluster(kmeans_result, data = scaled_data, geom = "point", ellipse.type = "convex") +
  theme_minimal() +
  ggtitle(paste("K-means Clustering com K =", best_k))

# -----------------------------
# Avaliação de qualidade do clustering
res <- cluster.stats(dist(scaled_data), kmeans_result$cluster)

cat("Calinski-Harabasz Index:", res$ch, "\n")
cat("Silhouette médio:", res$avg.silwidth, "\n")
# Um valor de 0.0385 é muito baixo, sugerindo que os clusters estão mal definidos, com muita sobreposição. Este é um sinal claro de que o agrupamento pode não estar a capturar estrutura real nos dados.
cat("Índice de Dunn:", res$dunn, "\n")
# Um valor abaixo de 1 (como 0.248) é baixo, indicando que os clusters estão muito próximos uns dos outros
cat("Within-cluster sum of squares:", res$within.cluster.ss, "\n")


set.seed(123)

gap_stat <- clusGap(scaled_data,
                    FUN = kmeans,
                    nstart = 25,
                    K.max = 10,
                    B = 5)  # número de amostragens bootstrap

# Ver os resultados
print(gap_stat, method = "firstmax")

# Plot do Gap Statistic
fviz_gap_stat(gap_stat) +
  ggtitle("Gap Statistic para Determinar o Número Ótimo de Clusters")


# Hierarchical

dist_matrix <- dist(scaled_data)

# Perform hierarchical clustering using Ward's method
hc <- hclust(dist_matrix, method = "ward.D2")

plot(hc, main = "Dendrograma do Clustering Hierárquico")

# Cut the dendrogram to get 2 clusters
df_clean$HClust_Cluster <- as.factor(cutree(hc, k = 2))


km_sil <- silhouette(as.integer(df_clean$KMeans_Cluster), dist(scaled_data))
cat("K-Means avg silhouette:", mean(km_sil[, 3]), "\n")

# Hierarchical clustering silhouette
hc_sil <- silhouette(as.integer(df_clean$HClust_Cluster), dist(scaled_data))
cat("Hierarchical avg silhouette:", mean(hc_sil[, 3]), "\n")


hc2 <- hclust(dist_matrix, method = "complete")

plot(hc2, labels = FALSE, hang = -1, main = "Hierarchical Clustering Dendrogram")

data$hclust_cluster <- cutree(hc2, k = 10)

rect.hclust(hc2, k = 10, border = 2:11) 




# Escolha do número de clusters (método da Silhueta)
sil_widths <- sapply(2:10, function(k) {
  clusters <- cutree(hc, k = k)
  mean(silhouette(clusters, dist_matrix)[, 3])
})

plot(2:10, sil_widths, type = "b", pch = 19,
     xlab = "Número de clusters",
     ylab = "Silhueta média",
     main = "Silhueta para escolher K")

# k=7
abline(v = 7, col = "red", lty = 2)

# Cortar a árvore em 5 clusters
clusters_hc <- cutree(hc, k = 7)

# Adicionar clusters ao dataset
df_clean$HC_Cluster <- as.factor(clusters_hc)

# Avaliação das métricas de qualidade do clustering
res <- cluster.stats(dist_matrix, clusters_hc)

cat("Calinski-Harabasz Index:", res$ch, "\n")
#  Um valor de 18.8 é baixo, sugerindo baixa separação entre grupos neste modelo hierárquico. No KMeans tinhas um valor 121.76, bem mais alto.
cat("Silhouette médio:", res$avg.silwidth, "\n")
# Um valor negativo indica que muitos pontos estão mais próximos de clusters vizinhos do que do seu próprio → estrutura de cluster fraca ou inexistente
cat("Índice de Dunn:", res$dunn, "\n")
# Este valor, embora ligeiramente melhor que no KMeans (0.248), continua fraco.
cat("Within-cluster sum of squares:", res$within.cluster.ss, "\n")
# Este valor é ligeiramente pior que o do KMeans (57,639.07), sugerindo menos coesão interna.

# Estabilidade dos clusters (Jaccard Index)
set.seed(123)
clusterboot_result <- clusterboot(
  dist_matrix, 
  clustermethod = disthclustCBI,
  method = "ward.D2", 
  k = 5,
  B = 100
)

cat("Estabilidade (Jaccard médio por cluster):\n")
print(clusterboot_result$bootmean)
# Estes valores (~0.15–0.20) são muito baixos, o que significa que:
#Os clusters não são consistentes entre diferentes amostras (bootstraps).
#A estrutura de agrupamento é provavelmente espúria ou sensível a ruído.


################### Gaussian Mixture Models (EM algorithm) ###############

# G = 7 clusters (escolhido previamente)
mod7 <- Mclust(scaled_data, G = 7)
summary(mod7)

# Visualização dos clusters em PC1 vs PC2
df7 <- as.data.frame(scaled_data)
df7$em_cluster <- as.factor(mod7$classification)

pc <- princomp(scaled_data)
scores <- as.data.frame(pc$scores)
scores$em_cluster <- df7$em_cluster

ggplot(scores, aes(x = Comp.1, y = Comp.2, color = em_cluster)) +
  geom_point() +
  theme_minimal() +
  labs(title = "EM Clustering com 7 Clusters", x = "PC1", y = "PC2", color = "Cluster")

################## Avaliação com métricas de cluster (para G = 7)

cl7 <- mod7$classification
dist_mat <- dist(scaled_data)
stats7 <- cluster.stats(dist_mat, cl7)
sil7 <- mean(silhouette(cl7, dist_mat)[, 3])

cat("Métricas para G = 7:\n")
cat("BIC:", mod7$bic, "\n")
cat("Calinski-Harabasz Index:", stats7$ch, "\n")
cat("Silhouette médio:", sil7, "\n")
cat("Dunn Index:", stats7$dunn, "\n")


#-------------------------
df <- read.csv("water_pollution_disease.csv")

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




