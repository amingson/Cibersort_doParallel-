# 定义并行版本的函数
parallel_doperm <- function(perm1, X1, Y1, absolute1, abs_method1,num_cores1=num_cores) {
  library(foreach)
  library(doParallel)
  # 设置并行计算的核心数 
  # 使用所有可用核心中的一个减少负担
  cl <- makeCluster(num_cores1)
  registerDoParallel(cl)
  
  # 将Y转换为列表格式
  Ylist <- as.list(data.matrix(Y1))
  
  # 并行计算
  dist <- foreach(itor = 1:perm1, .combine = rbind, .packages = c('IOBR')) %dopar% {
    # 每个工作线程的本地变量
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X1)[1])])
    yr <- (yr - mean(yr)) / sd(yr)
    result <- CoreAlg(X1, yr, absolute1, abs_method1)
    mix_r <- result$mix_r
    mix_r
  }
  
  # 关闭并行计算
  stopCluster(cl)
  closeAllConnections()
  # 返回结果
  newList <- list(dist = dist)
  return(newList)
}

