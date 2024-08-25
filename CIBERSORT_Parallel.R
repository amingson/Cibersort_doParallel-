CIBERSORT_Parallel <- function (sig_matrix = lm22, mixture_file, perm, QN = TRUE, 
                                absolute, abs_method = "sig.score",num_cores=50,parallel_doperm=parallel_doperm) 
{
  if (length(intersect(rownames(mixture_file), rownames(sig_matrix))) == 
      0) {
    stop("None identical gene between eset and reference had been found.\n         Check your eset using: intersect(rownames(eset), rownames(reference))")
  }
  if (absolute && abs_method != "no.sumto1" && abs_method != 
      "sig.score") 
    stop("abs_method must be set to either 'sig.score' or 'no.sumto1'")
  X <- sig_matrix
  Y <- rownames_to_column(mixture_file, var = "symbol")
  dups <- dim(Y)[1] - length(unique(Y[, 1]))
  if (dups > 0) {
    warning(paste(dups, " duplicated gene symbol(s) found in mixture file!", 
                  sep = ""))
    rownames(Y) <- make.names(Y[, 1], unique = TRUE)
  }
  else {
    rownames(Y) <- Y[, 1]
  }
  Y <- Y[, -1]
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  X <- X[order(rownames(X)), ]
  Y <- Y[order(rownames(Y)), ]
  P <- perm
  if (max(Y) < 50) {
    Y <- 2^Y
  }
  if (QN == TRUE) {
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  Yorig <- Y
  Ymedian <- max(median(Yorig), 1)
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, ]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, ]
  X <- (X - mean(X))/sd(as.vector(X))
  cat('Run perm')
  # absolute <- absolute
  # abs_method <-abs_method
  # num_cores <- num_cores
  if (P > 0) {
    # nulldist <- sort(doPerm(P, X, Y, absolute, abs_method)$dist)
    kk <- parallel_doperm(perm1 = P, X1=X, Y1=Y, absolute1=absolute, abs_method1=abs_method,num_cores1 = num_cores)
    nulldist <-  sort(kk$dist)
  }
  header <- c("Mixture", colnames(X), "P-value", "Correlation", 
              "RMSE")
  if (absolute) 
    header <- c(header, paste("Absolute score (", abs_method, 
                              ")", sep = ""))
  library(foreach)
  library(doParallel)
  
  # 设置并行计算的核心数
  num_cores <- num_cores
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  mixtures <- dim(Y)[2]
  pval <- 9999
  cat('Run foreach')
  # 并行计算
  results <- foreach(itor = 1:mixtures, .combine = rbind, .packages = c('IOBR')) %dopar% {
    # 每个工作线程的本地变量
    y <- Y[, itor]
    y <- (y - mean(y)) / sd(y)
    result <- CoreAlg(X, y, absolute, abs_method)
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse
    
    if (absolute && abs_method == "sig.score") {
      w <- w * median(Y[, itor]) / Ymedian
    }
    
    pval <- 9999
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }
    
    out <- c(colnames(Y)[itor], w, pval, mix_r, mix_rmse)
    if (absolute) 
      out <- c(out, sum(w))
    
    out
  }
  
  # 关闭并行计算
  stopCluster(cl)
  closeAllConnections()
  # 输出结果
  output <- results
  
  obj <- rbind(header, output)
  obj <- obj[, -1]
  obj <- obj[-1, ]
  obj <- matrix(as.numeric(unlist(obj)), nrow = nrow(obj))
  rownames(obj) <- colnames(Y)
  if (!absolute) {
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", 
                       "RMSE")
  }
  else {
    colnames(obj) <- c(colnames(X), "P-value", "Correlation", 
                       "RMSE", paste("Absolute score (", abs_method, ")", 
                                     sep = ""))
  }
  obj
}
