source("clustering_helper.R")

data_path <- "./data"
result_path <- "./results"
pathway <- "hallmark"

D <- 100
K <- 15

X <- log2(all_mrna + 1)
common_patients <- intersect(rownames(X), rownames(all_tissue))
X <- X[common_patients,]
sigma_N <- 1000

pathways <- read_pathways(pathway)
gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
X <- X[, which(colnames(X) %in% gene_names)]
N_pathway <- length(pathways)

valid_features <- as.numeric(which(apply(X, 2, sd) != 0))
X <- X[, valid_features]
X <- scale(X)
set.seed(1606)
distance_indices <- sample(1:nrow(X), sigma_N)
Z <- list()
K_train_silhouette <- list()
for(m in 1:length(pathways)) {
  feature_indices <- which(colnames(X) %in% pathways[[m]]$symbols)
  D_train <- pdist(X[distance_indices, feature_indices], X[distance_indices, feature_indices])
  sigma <- mean(D_train)
  K_train_silhouette[[m]] <- exp(-D_train^2 / (2 * sigma^2))
  
  gamma <- 1 / (2 * sigma^2)
  sigma_p <- sqrt(2 * gamma)
  
  W <- matrix(rnorm(n = length(feature_indices) * D, mean = 0, sd = sigma_p), nrow = length(feature_indices), ncol = D)
  b <- runif(D, min = 0, max = 2 * pi)
  Z[[m]] <- sqrt(1 / D) * cbind(cos(X[, feature_indices] %*% W + matrix(b, nrow = nrow(X), ncol = D, byrow = TRUE)), sin(X[, feature_indices] %*% W + matrix(b, nrow = nrow(X), ncol = D, byrow = TRUE)))
}

indices <- 1:N_pathway
npathways_desired <- 5

set.seed(1606)
selected <- numeric()
value <- numeric()
final_cluster_objects <- list()
while(length(selected) < npathways_desired) {
  print(length(selected))
  candidates <- setdiff(indices, selected)
  measure <- numeric()
  cluster_object <- list()
  for(candidate in candidates) {
    Z_union <- do.call(cbind, Z[union(selected, candidate)])
    Ks_to_use <- K_train_silhouette[union(selected, candidate)]
    K_avg <- Reduce('+', Ks_to_use) / length(Ks_to_use)
    cluster_object[[candidate]] <- kmeans(Z_union, centers = K, nstart = 5, iter.max = 100, algorithm = "Hartigan-Wong")
    si <- cluster::silhouette(cluster_object[[candidate]]$cluster[distance_indices], dmatrix = 2 - 2 * K_avg)
    measure[candidate] <- mean(si[,3])
  }
  to_add <- which.max(measure)
  value[length(selected) + 1] <- measure[to_add]
  final_cluster_objects[[length(selected) + 1]] <- cluster_object[[to_add]]
  selected <- c(selected, to_add)
}

best_cluster_object <- final_cluster_objects[[which.max(value)]]
Z_used <- do.call(cbind, Z[selected[1:which.max(value)]])

result <- list()
result$best_cluster_object <- best_cluster_object
result$selected <- selected
result$value <- value
result$pathway <- pathway
result$D <- D
result$Z_used <- Z_used
result$K <- K

save("result", file = sprintf("%s/forward_makkmeans_D_%d_pathway_%s_K_%d_result.RData", result_path, D, pathway, K))
