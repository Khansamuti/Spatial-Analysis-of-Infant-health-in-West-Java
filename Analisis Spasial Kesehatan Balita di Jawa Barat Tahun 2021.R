library(readxl)
library(rgdal)

# Prepare Dataset
df <- read_excel("kesehatan balita di jawa barat 2021.xlsx")

X1 <- df$luas_wilayah
X2 <- df$presentase_wilayah
X3 <- df$jumlah_peserta_kb
X4 <- df$jumlah_posyandu_aktif
X5 <- df$presentase_pemberian_asi
X6 <- df$jumlah_balita_mendapat_vit_a
X7 <- df$jumlah_balita_gizi_buruk

y_cont <- df$presentase_bb_lahir_rendah
y_pois <- df$jumlah_balita_mendapat_imunisasi_lengkap

## Clustering
k <- 2
cluster_result <- kmeans(df$jumlah_balita_kurang_gizi, centers=k)
df$cluster_balita_kurang_gizi <- cluster_result[1]$cluster - 1
df$cluster_balita_kurang_gizi
y_binary <- df$cluster_balita_kurang_gizi

coordinates(df) <- c("longitude", "latitude")
class(df)
n <- length(df)
k <- 7

# ===========GWR===========
# global regression
glm <- glm(y_cont ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, data=df,
           family="gaussian")

bw_gwr <- bw.gwr(y_cont ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, df, 
                 kernel="bisquare", adaptive=T)
model_gwr <- gwr.basic(y_cont ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, df, 
                       bw=bw_gwr, kernel="bisquare", adaptive=T)

# Uji kesesuaian
# H0: bk (ui, vi) = bk
# H1: minimal ada bk (ui, vi) != bk

# Perform likelihood ratio test
lr_test <- anova(glm, model_gwr, test = "Chisq")

# Extract relevant information from the test result
chi2_statistic <- lr_test$Deviance[2]
dof <- lr_test$Df[2]
p_value <- pchisq(chi2_statistic, dof, lower.tail = FALSE)
decision <- ifelse(p_value < 0.05, "Reject H0", "Fail to reject H0")

# Print the results
cat("Test Statistic (Chi-Square):", chi2_statistic, "\n")
cat("Degrees of Freedom:", dof, "\n")
cat("P-value:", p_value, "\n")
cat("Decision:", decision, "\n")

## Uji parameter parsial pada lokasi ke-i
# H0: beta(ui, vi) = 0
# H1: beta(ui, vi) != 0
# Create a list of variable names
variable_names <- c("Intercept_TV", "X1_TV", "X2_TV", "X3_TV", "X4_TV", "X5_TV", "X6_TV", "X7_TV")

# Initialize empty data frames for t-values and p-values
t_values_df <- data.frame(matrix(nrow = length(model_gwr$SDF[[variable_names[1]]]), ncol = length(variable_names)))
p_values_df <- data.frame(matrix(nrow = length(model_gwr$SDF[[variable_names[1]]]), ncol = length(variable_names)))

# Assign t-values and p-values to data frames
for (i in 1:length(variable_names)) {
  t_values_df[, i] <- model_gwr$SDF[[variable_names[i]]]
  p_values_df[, i] <- 2 * pt(abs(model_gwr$SDF[[variable_names[i]]]), df = n - k - 1, lower.tail = FALSE)
}

# Set the column names
colnames(t_values_df) <- variable_names
colnames(p_values_df) <- variable_names

# Print the data frames
print(t_values_df)
print(p_values_df)

# Uji Serentak
# H0: b1 (ui, vi) = ... = bp (ui, vi) = 0
# H1: minimal ada satu bk (ui, vi) != 0
decision <- ifelse(glm$deviance > qchisq(0.05, n-k-1), "Reject H0", "Fail to reject H0")
decision


# ===========GWLR===========
DM <- gw.dist(dp.locat=coordinates(df), longlat=T)
bw_gwlr <- bw.ggwr(y_binary ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, df, 
                   family="binomial", kernel="bisquare", adaptive=T, 
                   dMat=DM, longlat=T, approach="AICc")
model_gwlr <- ggwr.basic(y_binary ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, df, 
                         bw=bw_gwlr, family="binomial", kernel="bisquare", 
                         dMat=DM, adaptive=T, longlat=T)

## Uji parameter parsial pada lokasi ke-i
# H0: beta(ui, vi) = 0
# H1: beta(ui, vi) != 0
# Create a list of variable names
variable_names <- c("Intercept_TV", "X1_TV", "X2_TV", "X3_TV", "X4_TV", "X5_TV", "X6_TV", "X7_TV")

# Initialize empty data frames for t-values and p-values
t_values_df <- data.frame(matrix(nrow = length(model_gwlr$SDF[[variable_names[1]]]), ncol = length(variable_names)))
p_values_df <- data.frame(matrix(nrow = length(model_gwlr$SDF[[variable_names[1]]]), ncol = length(variable_names)))

# Assign t-values and p-values to data frames
for (i in 1:length(variable_names)) {
  t_values_df[, i] <- model_gwlr$SDF[[variable_names[i]]]
  p_values_df[, i] <- 2 * pt(abs(model_gwlr$SDF[[variable_names[i]]]), df = n - k - 1, lower.tail = FALSE)
}

# Set the column names
colnames(t_values_df) <- variable_names
colnames(p_values_df) <- variable_names

# Print the data frames
print(t_values_df)
print(p_values_df)

# Uji Serentak
# H0: b1 (ui, vi) = ... = bp (ui, vi) = 0
# H1: minimal ada satu bk (ui, vi) != 0
decision <- ifelse(model_gwlr$glms$deviance > qchisq(0.05, n-k-1), "Reject H0", "Fail to reject H0")
decision


# ===========GWPR===========
# Poissonn Regression
reg_pois <- glm(y_cont ~ X1 + X2 + X3 +X4 + X5 + X6 + X7, df, family="poisson")

library(AER)
dispersiontest(reg_pois, alternative="greater") # equidispersion

# GWPR
DM <- gw.dist(dp.locat=coordinates(df), longlat=T)
bw_gwpr <- bw.ggwr(y_pois ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, df, 
                   family="poisson", kernel="gaussian", adaptive=T, 
                   dMat=DM, longlat=T, approach="AICc")
model_gwpr <- ggwr.basic(y_pois ~ X1 + X2 + X3 + X4 + X5 + X6 + X7, df, 
                         bw=bw_gwlr, family="poisson", kernel="gaussian", 
                         dMat=DM, adaptive=T, longlat=T)

# Uji kesesuaian
# H0: The Poisson regression model and the GWPR model are equivalent or have equal fit. There is no significant difference between the two models.
# H1: The Poisson regression model and the GWPR model are not equivalent or do not have equal fit. There is a significant difference between the two models.

# Perform likelihood ratio test
lr_test <- anova(reg_pois, model_gwpr, test = "Chisq")

# Extract relevant information from the test result
chi2_statistic <- lr_test$Deviance[2]
dof <- lr_test$Df[2]
p_value <- pchisq(chi2_statistic, dof, lower.tail = FALSE)
decision <- ifelse(p_value < 0.05, "Reject H0", "Fail to reject H0")

# Print the results
cat("Test Statistic (Chi-Square):", chi2_statistic, "\n")
cat("Degrees of Freedom:", dof, "\n")
cat("P-value:", p_value, "\n")
cat("Decision:", decision, "\n")

## Uji parameter parsial pada lokasi ke-i
variable_names <- c("Intercept_TV", "X1_TV", "X2_TV", "X3_TV", "X4_TV", "X5_TV", "X6_TV", "X7_TV")

### Initialize empty data frames for t-values and p-values
t_values_df <- data.frame(matrix(nrow = length(model_gwpr$SDF[[variable_names[1]]]), ncol = length(variable_names)))
p_values_df <- data.frame(matrix(nrow = length(model_gwpr$SDF[[variable_names[1]]]), ncol = length(variable_names)))

### Assign t-values and p-values to data frames
for (i in 1:length(variable_names)) {
  t_values_df[, i] <- model_gwpr$SDF[[variable_names[i]]]
  p_values_df[, i] <- 2 * pt(abs(model_gwpr$SDF[[variable_names[i]]]), df = n - k - 1, lower.tail = FALSE)
}

### Set the column names
colnames(t_values_df) <- variable_names
colnames(p_values_df) <- variable_names

### Print the data frames
print(t_values_df)
print(p_values_df)

# Uji Serentak
# H0: b1 (ui, vi) = ... = bp (ui, vi) = 0
# H1: minimal ada satu bk (ui, vi) != 0
decision <- ifelse(model_gwpr$glms$deviance > qchisq(0.05, n-k-1), "Reject H0", "Fail to reject H0")
decision

