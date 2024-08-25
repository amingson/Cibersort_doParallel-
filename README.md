#从单细胞矩阵开始运行 对TCGA进行预测
'
library(data.table)
library(tidyverse)
library(tibble)
library(Seurat)
sce <- readRDS("../03.TAMsub/results/sce.rds")
sce.norm <- sce
DefaultAssay(sce.norm) <- 'RNA'
sce.norm <- NormalizeData(sce.norm)
sce.norm <-FindVariableFeatures(sce.norm,selection.method = "vst",nfeatures = 2000)
sce.norm <- ScaleData(sce.norm,features = rownames(sce.norm))

#制作卷积参考矩阵 
Idents(sce.norm) <- "cell_type3"
#制作细胞类型的系数表
X <- AverageExpression(sce.norm)[[1]]    
write.table(X,"./results/cibsort_ref.txt",sep = "\t",col.names = T,row.names = T)
ref_cib <- as.data.frame(fread("./results/cibsort_ref.txt")) %>% column_to_rownames('V1') %>% as.data.frame()
#===========================================================对TCGA数据进行解卷积
library(IOBR)
tcga_dat <- as.data.frame(fread("../01.datapre/TCGA/results/tcga_dat_T.txt")) %>% column_to_rownames('V1') %>% as.data.frame()
source("./parallel_doperm.R")
source("./CIBERSORT_Parallel.R")
library(preprocessCore)
closeAllConnections()
#自定义运行函数  批量运行 加快速度
tcga_dat_cib <- CIBERSORT_Parallel(sig_matrix = ref_cib,mixture_file = tcga_dat,perm = 100,absolute = F,num_cores = 50,parallel_doperm = parallel_doperm)
'
