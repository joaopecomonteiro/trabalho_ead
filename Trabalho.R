library(tidyverse)
library(corrplot)
library(ggplot2)
library(factoextra)
library(cluster)
library(psych)
library(RColorBrewer)


ai_data <- read.csv("/home/joaomonteiro/Desktop/EAD/Trabalho/AI_index_db.csv")


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
  ymax <- max(counts) + 7
  
  bp <- barplot(counts,
                main = paste("Distribution of", var),
                col = pastel_colors[1:n_colors],
                ylim = c(0, ymax),
                xaxt = "n")  
  
  text(x = bp,
       y = par("usr")[3] - 1.5,  
       labels = names(counts),
       srt = 45,
       adj = 1,
       xpd = TRUE,
       cex = 0.8)
}


# 2. BIVARIATE ANALYSIS

# 2.1 Correlation analysis for numerical variables

cor_matrix <- cor(num_vars, method = "pearson")
round(cor_matrix, 2)  

corrplot(cor_matrix, method = "circle", type = "upper",
         tl.cex = 0.8, number.cex = 0.7, addCoef.col = "black",
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("red", "white", "blue"))(200))

# 2.2 Scatter plots for selected pairs
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

boxplot(Total.score ~ Region, data = ai_data, 
        main = "Total AI Score by Region", col = rainbow(length(levels(ai_data$Region))),
        las = 2, cex.axis = 0.7)

anova_region <- aov(Total.score ~ Region, data = ai_data)
summary(anova_region)

boxplot(Total.score ~ Cluster, data = ai_data, 
        main = "Total AI Score by Cluster", col = rainbow(length(levels(ai_data$Cluster))),
        las = 2)

anova_cluster <- aov(Total.score ~ Cluster, data = ai_data)
summary(anova_cluster)


# 3. MULTIVARIATE ANALYSIS

# 3.1 Principal Component Analysis (PCA)

data_numeric <- ai_data %>%
  select(where(is.numeric)) %>%
  select(-Total.score) %>%
  na.omit()

pca_result <- prcomp(data_numeric, scale = TRUE)

summary(pca_result)

fviz_eig(pca_result, addlabels = TRUE) + ggtitle("")


fviz_pca_var(pca_result, col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)


fviz_pca_biplot(pca_result, repel = TRUE,
                col.var = "#2E9FDF", 
                col.ind = "#696969") 

fviz_contrib(pca_result, choice = "var", axes = 1, top = 10) + ggtitle("") # PC1
fviz_contrib(pca_result, choice = "var", axes = 2, top = 10) + ggtitle("") # PC2
fviz_contrib(pca_result, choice = "var", axes = 3, top = 10) + ggtitle("") # PC3

loadings <- pca_result$rotation[, 1:3]
contrib <- loadings^2
contrib <- sweep(contrib, 2, colSums(contrib), "/") * 100
round(contrib, 2)


# 3.2 Factor Analysis
fa.parallel(data_numeric, fa = "fa") + ggtitle("")

fa_result <- fa(data_numeric, nfactors = 2, rotate = "varimax", fm = "ml")

print(fa_result)


# 3.3 Multidimensional Scaling (MDS)


data_numeric <- ai_data %>%
  select(where(is.numeric)) %>%
  scale()


d <- dist(data_numeric, method = "euclidean")

mds_result <- mds(d, type = "interval")
nmds_result <- mds(d, type = "ordinal")

mds_result$stress
nmds_result$stress


mds_df <- as.data.frame(mds_result$conf)

# Non-metric MDS
nmds_df <- as.data.frame(nmds_result$conf)

# Add country names (assuming you have them)
mds_df$Label <- ai_data$Country
nmds_df$Label <- ai_data$Country

mds_df$Cluster <- ai_data$Cluster
nmds_df$Cluster <- ai_data$Cluster



mds_df$Total.score = ai_data$Total.score

mds_labels <- mds_df %>%
  group_by(Cluster) %>%
  arrange(desc(Total.score)) %>%
  slice(1:5)


ggplot(mds_df, aes(x = D1, y = D2, color = Cluster)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.8, linetype = "dashed") +
  geom_text_repel(data = mds_labels, aes(label = Label), color = "black", size = 3, max.overlaps = 20, box.padding = 0.5) +
  #ggtitle("NMDS Plot with Group Ellipses and 5 Black Labels per Group") +
  theme_minimal() +
  theme(legend.title = element_blank())



nmds_df$Total.score = ai_data$Total.score

nmds_labels <- nmds_df %>%
  group_by(Cluster) %>%
  arrange(desc(Total.score)) %>%
  slice(1:5)


ggplot(nmds_df, aes(x = D1, y = D2, color = Cluster)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.8, linetype = "dashed") +
  geom_text_repel(data = nmds_labels, aes(label = Label), color = "black", size = 3, max.overlaps = 50, box.padding = 0.5) +
  #ggtitle("NMDS Plot with Group Ellipses and 5 Black Labels per Group") +
  theme_minimal() +
  theme(legend.title = element_blank())



