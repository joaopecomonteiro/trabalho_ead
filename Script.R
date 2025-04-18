library(tidyverse)
library(corrplot)
library(ggplot2)
library(factoextra)
library(cluster)
library(psych)
library(RColorBrewer)


ai_data <- read.csv("C:/Users/diogo/OneDrive/Ambiente de Trabalho/FCUP - Data Science/1º Ano/2º Semestre/Estatística e Análise de Dados/1st Assignment/AI_index_db.csv")


# Basic inspection of the data
head(ai_data)
str(ai_data)
summary(ai_data)

# Check for missing values
sum(is.na(ai_data))

# Convert categorical variables to factors
ai_data$Region <- as.factor(ai_data$Region)
ai_data$Cluster <- as.factor(ai_data$Cluster)
ai_data$Income.group <- as.factor(ai_data$Income.group)
ai_data$Political.regime <- as.factor(ai_data$Political.regime)

# 1. UNIVARIATE ANALYSIS

# 1.1 Descriptive statistics for numerical variables
num_vars <- ai_data %>% 
  select(Talent, Infrastructure, Operating.Environment, Research, 
         Development, Government.Strategy, Commercial, Total.score)

# Summary statistics
summary_stats <- data.frame(
  Variable = names(num_vars),
  Mean = sapply(num_vars, mean),
  Trimmed_Mean = sapply(num_vars, mean, trim = 0.1),
  Median = sapply(num_vars, median),
  Min = sapply(num_vars, min),
  Max = sapply(num_vars, max),
  Range = sapply(num_vars, function(x) max(x) - min(x)),
  IQR = sapply(num_vars, IQR),
  SD = sapply(num_vars, sd),
  CV = sapply(num_vars, function(x) sd(x)/mean(x)),
  Skewness = sapply(num_vars, function(x) psych::skew(x)),
  Kurtosis = sapply(num_vars, function(x) psych::kurtosi(x))
)

print(summary_stats)

# 1.2 Histograms for all numerical variables
par(mfrow = c(2, 4))
for (var in names(num_vars)) {
  hist(ai_data[[var]], breaks=10, main = paste("Histogram of", var), xlab = var, col = "lightblue")
}


# 1.3 Boxplots for all numerical variables
par(mfrow = c(2, 4))
for (var in names(num_vars)) {
  boxplot(ai_data[[var]], main = paste("Boxplot of", var), col = "lightblue")
}

# 1.4 Check for outliers (values beyond 1.5 * IQR)
outliers <- lapply(num_vars, function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  
  x[x < lower_bound | x > upper_bound]
})

# Print countries with outlier values
for (i in 1:length(outliers)) {
  if (length(outliers[[i]]) > 0) {
    cat("\nOutliers in", names(num_vars)[i], ":\n")
    outlier_indices <- which(num_vars[[i]] %in% outliers[[i]])
    print(ai_data$Country[outlier_indices])
    print(outliers[[i]])
  }
}


# 1.5 Frequency tables for categorical variables
cat_vars <- c("Region", "Cluster", "Income.group", "Political.regime")

for (var in cat_vars) {
  cat("\nFrequency table for", var, ":\n")
  print(table(ai_data[[var]]))
  cat("\nRelative frequency for", var, ":\n")
  print(prop.table(table(ai_data[[var]])))
}

# 1.6 Bar charts for categorical variables
pastel_colors <- brewer.pal(5, "Pastel1")

par(mfrow = c(2, 2))
for (var in cat_vars) {
  counts <- table(ai_data[[var]])
  n_colors <- length(counts)
  barplot(counts, 
          main = paste("Distribution of", var), 
          col = pastel_colors[1:n_colors], 
          las = 2, 
          cex.names = 0.7)
}

# 2. BIVARIATE ANALYSIS

# 2.1 Correlation analysis for numerical variables
correlation_matrix <- cor(num_vars)
print(correlation_matrix)

# Visualize the correlation matrix
corrplot(correlation_matrix, method = "circle", type = "upper", 
         tl.col = "black", tl.srt = 45)

# 2.2 Scatter plots for selected pairs
# Let's examine the relationship between Total score and some key variables
pairs(~ Total.score + Talent + Infrastructure + Research + Development, 
      data = ai_data, main = "Scatter Plot Matrix of Key Variables")

# 2.3 Contingency tables for categorical variables
for (i in 1:(length(cat_vars)-1)) {
  for (j in (i+1):length(cat_vars)) {
    cat("\nContingency table for", cat_vars[i], "and", cat_vars[j], ":\n")
    print(table(ai_data[[cat_vars[i]]], ai_data[[cat_vars[j]]]))
  }
}

# 2.4 Comparison of numerical variables across groups
# For example, comparing Total scores across different regions
boxplot(Total.score ~ Region, data = ai_data, 
        main = "Total AI Score by Region", col = rainbow(length(levels(ai_data$Region))),
        las = 2, cex.axis = 0.7)

# ANOVA to test if mean Total scores differ by region
anova_region <- aov(Total.score ~ Region, data = ai_data)
summary(anova_region)

# Compare Total scores across clusters
boxplot(Total.score ~ Cluster, data = ai_data, 
        main = "Total AI Score by Cluster", col = rainbow(length(levels(ai_data$Cluster))),
        las = 2)

# ANOVA for clusters
anova_cluster <- aov(Total.score ~ Cluster, data = ai_data)
summary(anova_cluster)

# 3. MULTIVARIATE ANALYSIS

# 3.1 Principal Component Analysis (PCA)
# Select only numerical variables for PCA
pca_data <- num_vars
pca_result <- prcomp(pca_data, scale = TRUE)

# Summary of PCA
summary(pca_result)

# Scree plot to visualize eigenvalues
fviz_eig(pca_result, addlabels = TRUE)

# Visualize variables on principal components
fviz_pca_var(pca_result, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Biplot: Visualize both variables and observatihttp://127.0.0.1:36483/graphics/plot_zoom_png?width=460&height=253ons
fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")  # Individuals color

# 3.2 Cluster Analysis
# Standardize the data for clustering
scaled_data <- scale(pca_data)

# Determine optimal number of clusters using silhouette method
fviz_nbclust(scaled_data, kmeans, method = "silhouette")

# K-means clustering with the suggested number of clusters (let's say k=3)
set.seed(123) # For reproducibility
k_means_result <- kmeans(scaled_data, centers = 3, nstart = 25)

# Add cluster assignments to the original data
ai_data$cluster <- as.factor(k_means_result$cluster)

# Visualize clusters on PCA plot
fviz_cluster(k_means_result, data = scaled_data,
             geom = "point",
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())

# Examine cluster characteristics
cluster_means <- aggregate(pca_data, by = list(Cluster = k_means_result$cluster), mean)
print(cluster_means)

# Compare clusters with original categorizations
table(ai_data$Cluster, ai_data$cluster)

# 3.3 Hierarchical Clustering
# Calculate distance matrix
dist_matrix <- dist(scaled_data)

# Perform hierarchical clustering
hc_result <- hclust(dist_matrix, method = "ward.D2")

# Plot dendrogram
plot(hc_result, main = "Hierarchical Clustering Dendrogram", 
     xlab = "", sub = "", cex = 0.6)

# Cut tree to form clusters
hc_clusters <- cutree(hc_result, k = 3)

# Add hierarchical cluster assignments to data
ai_data$hc_cluster <- as.factor(hc_clusters)

# Compare hierarchical clusters with K-means clusters
table(ai_data$cluster, ai_data$hc_cluster)

# 4. CONCLUSION
# This section will be filled in by you based on the results of the analysis

# Create a data frame with Country names and their cluster assignments for reference
country_clusters <- data.frame(
  Country = ai_data$Country,
  KMeans_Cluster = ai_data$cluster,
  HC_Cluster = ai_data$hc_cluster,
  Original_Group = ai_data$Cluster
)

print(country_clusters)

# Export results if needed
write.csv(summary_stats, "ai_index_summary_stats.csv", row.names = FALSE)
write.csv(country_clusters, "ai_index_cluster_assignments.csv", row.names = FALSE)









