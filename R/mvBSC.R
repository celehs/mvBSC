#' @import stats
#' @importFrom Matrix sparseMatrix
NULL

#' @title ...
#' @param cosK ...
#' @param h ...
#' @param k ...
#' @param R ...
#' @export
get_U <- function(cosK, h, k, R) {
  cosK_band <- cosK * (R <= h)
  eigen_cosK_band <- eigen(cosK_band)
  idx <- order(abs(eigen_cosK_band$values), decreasing = TRUE)
  U <- eigen_cosK_band$vectors[, idx[1:k]]
  rownames(U) <- rownames(cosK)
  return(U)
}

#' @title ...
#' @param codes ...
#' @param labels ...
#' @export
get_Z <- function(codes, labels) {
  uniqueLabels <- sort(unique(labels))
  Z <- 1 * outer(labels, uniqueLabels, "==")  
  rownames(Z) <- codes
  colnames(Z) <- uniqueLabels
  return(Z)
}

#' @title ...
#' @param U ...
#' @param Zvec group membership vector
#' @export
ratio <- function(U, Zvec) {
  n <- length(Zvec)
  mu0 <- colMeans(U) # grand mean of all codes
  sstot <- sum(U^2) - n * sum(mu0^2)
  ssbtw <- 0
  g <- unique(Zvec)
  nk <- table(Zvec)[g]
  for (i in 1:length(g)) {
    muk <- colMeans(subset(U, Zvec == g[i])) # group mean
    ssbtw <- ssbtw + nk[i] * sum((muk - mu0)^2)
  }
  names(ssbtw) <- NULL
  ratio <- ssbtw / sstot
  return(ratio)
}

#' @title ...
#' @param Z1 ...
#' @param Z2 ...
#' @param digits ...
#' @export
NMI <- function(Z1, Z2, digits = 3) {
  Z1 <- as.matrix(Z1) 
  Z2 <- as.matrix(Z2)
  n <- nrow(Z1)
  A <- t(Z1) %*% Z2
  probA = A / sum(A) #joint prob
  probZ1 = apply(probA, 1, sum)
  probZ2 = apply(probA, 2, sum)
  probB = log(probA / (probZ1 %*% t(probZ2)))
  probB[is.infinite(probB)] <- 0
  top <- sum(probA * probB)
  bottom <- sqrt(t(probZ1) %*% log(probZ1) * t(probZ2) %*% log(probZ2))
  nmi <- as.numeric(top / bottom)
  return(round(nmi, digits = digits))
}

#' @title ...
#' @param Z1 ...
#' @param Z2 ...
#' @param digits ...
#' @export
ARI <- function(Z1, Z2, digits = 3) {
  n <- nrow(Z1)
  contingencyTb <- t(Z1) %*% Z2
  one <- rep(1, n)
  a <- t(Z1) %*% one
  b <- t(Z2) %*% one
  A <- 0.5 * sum(apply(contingencyTb, c(1, 2), function(x) x^2 - x)) 
  B <- 0.25 * sum(a * (a - 1)) * sum(b * (b - 1))
  C <- 0.5 * n * (n - 1)
  D <- 0.25 * (sum(a * (a - 1)) + sum(b * (b - 1)))
  ARI <- (A - B / C) / (D - B / C)
  return(round(ARI, digits = digits))
}

#' @title ...
#' @param Z ...
#' @param Z0 ...
#' @param digits ...
#' @export
F_measure <- function(Z, Z0, digits = 3) {
  PZ <- Z %*% t(Z)
  PZ0 <- Z0 %*% t(Z0)
  precision <- round(sum(PZ * PZ0) / sum(PZ), digits = digits)
  recall <- round(sum(PZ * PZ0) / sum(PZ0), digits = digits)
  f_measure <- round(2 * precision * recall / (precision + recall), digits = digits)
  return(c("precision" = precision, "recall" = recall, "f" = f_measure))
}

#' @title ...
#' @param codes ...
#' @param distance ...
#' @param similarity ...
#' @param k ...
#' @param ncluster ...
#' @param delta ...
#' @param h ...
#' @param wt ...
#' @param seed ...
#' @export
mvbsc_fit0 <- function(codes, distance, similarity, k, ncluster, delta, h, wt, seed) {
  R <- distance[codes, codes]
  W <- 0
  for (i in 1:length(wt)) {
    if (wt[i] > 0) {
      cosK <- similarity[[i]][codes, codes]
      tmp <- get_U(cosK = cosK, h = h, k = k, R = R)
      W <- W + wt[i] * tcrossprod(tmp)
    }
  }
  W_eigen <- eigen(W)
  idx <- order(abs(W_eigen$values), decreasing = TRUE)
  U <- W_eigen$vectors[, idx[1:k]]
  rownames(U) <- rownames(W)
  set.seed(seed)
  fit <- kmeans(U, ncluster, iter.max = 100, nstart = 25)
  tbl <- data.frame(cluster = 1:ncluster, size = NA, max_dist = NA)
  for (i in 1:ncluster) {
    v <- names(fit$cluster)[fit$cluster == i]
    tbl$size[i] <- length(v)
    tbl$max_dist[i] <- max(R[v, v])
  }
  tbl2 <- tbl[order(tbl$max_dist), ]
  rownames(tbl2) <- 1:ncluster
  list(cluster = fit$cluster, cluster_info = tbl2, U = U)
}

#' @title ...
#' @param codes ...
#' @param distance ...
#' @param similarity ...
#' @param k ...
#' @param delta ...
#' @param h ...
#' @param wt ...
#' @param seed ...
#' @export
mvbsc_fit <- function(codes, distance, similarity, k, delta, h, wt, seed) {
  initial <- mvbsc_fit0(
    codes = codes, 
    distance = distance,
    similarity = similarity, 
    k = k, 
    ncluster = k,
    delta = delta,
    h = h,
    wt = wt, 
    seed = seed)
  cluster0 <- subset(initial$cluster_info, initial$cluster_info$max_dist > delta)$cluster
  regroup <- vector("list", length(cluster0))
  DF <- data.frame(level1 = initial$cluster, level2 = 0)
  for (i in 1:length(cluster0)) {
    codes0 <- names(initial$cluster[initial$cluster == cluster0[i]])
    for (nc in 2:(length(codes0) - 1)) {
      fit0 <- mvbsc_fit0(
        codes = codes0,
        distance = distance, 
        similarity = similarity, 
        k = nc,
        ncluster = nc,
        delta = delta,           
        h = h, 
        wt = wt,   
        seed = seed)
      if (all(fit0$cluster_info$max_dist <= delta)) break
    }
    regroup[[i]] <- fit0
    DF[names(regroup[[i]]$cluster), "level2"] <- regroup[[i]]$cluster
  }
  cluster <- paste0("C", DF$level1, ".", DF$level2)
  cluster <- factor(cluster, labels = 1:length(unique(cluster)))
  names(cluster) <- rownames(DF)
  ratio <- ratio(initial$U[names(cluster), ], cluster)
  # summary <- data.frame(delta = delta, h = h, ratio = ratio)
  list(ratio = ratio,
       k = k,
       delta = delta,
       h = h,
       wt = wt, 
       size = table(cluster),
       cluster = cluster)  
}

#' @title ...
#' @param codes ...
#' @param distance ...
#' @param similarity ...
#' @param k ...
#' @param delta ...
#' @param h ...
#' @param wt ...
#' @param seed ...
#' @export
mvbsc <- function(codes, distance, similarity, k, 
                  delta = NULL, h = NULL, wt = NULL, seed = 123) {
  if (is.null(delta)) delta <- min(apply(distance, 1, max))
  if (is.null(h)) h <- seq(0, max(delta), length.out = 11)[-1]
  if (is.null(wt)) wt <- rep(1, length(similarity)) / length(similarity)
  DF <- expand.grid(k = k, delta = delta, h = h, ratio = NA)
  N <- NROW(DF)
  fit <- vector("list", N)
  for (i in 1:N) {
    fit[[i]] <- mvbsc_fit(codes = codes,
                          distance = distance, 
                          similarity = similarity,
                          k = DF$k[i],
                          delta = DF$delta[i], 
                          h = DF$h[i],
                          wt = wt,
                          seed = seed)  
    DF$ratio[i] <- fit[[i]]$ratio
  }
  idx <- which.max(DF$ratio)
  list(summary = DF[order(-DF$ratio), ],
       optimal = fit[[idx]])
}

#' @title ...
#' @param W UU'
#' @param k ...
#' @param constraint_set ...
#' @param alpha ...
#' @param beta ...
#' @param seed ...
#' @export
semi_kmeans <- function(W, k, constraint_set, alpha, beta, seed) {
  n <- nrow(W)
  A <- Matrix::sparseMatrix(
    i = constraint_set$i,
    j = constraint_set$j,
    x = constraint_set$link,
    dims = c(n, n))
  A <- A + t(A)
  diag(A) <- 1
  K <- W + alpha * (A > 0) - beta * (A < 0)
  eigK <- eigen(K)
  large <- sort.int(eigK$values, decreasing = TRUE, index.return = TRUE)
  # newU <- eigK$vectors[,large$ix[1:k]] # %*%diag(sqrt(large$x[1:k])) # should not use eigenvectors only
  newk <- min(k, sum(eigK$values > 0))
  newU <- eigK$vectors[, large$ix[1:newk]] %*% diag(sqrt(large$x[1:newk]))
  rownames(newU) <- rownames(W)
  set.seed(seed)
  res <- kmeans(newU, k, iter.max = 1000, nstart = 25) 
  res
}
