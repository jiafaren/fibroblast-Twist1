

######################################downstream analysis#######01 create_object
##create seurat file 为了接下来的数据整合（使用step2的 /DimReduction/scOpen/scOpen_barcodes.txt） 简单质控 data integration label transfer 聚类 下游分析  这里和step1不同
#step1 只是单纯的为了step2的降维操作，因为scopen是作者辛苦开发出来的降维工具
#整好之后保存为/data/uuo.Rds


```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
```

## Load data
```{r load_data}
counts <- readRDS("../PeakMatrix.Rds")
dim(counts)

metadata <- read.csv(
  file = "../meta_data.csv",
  header = TRUE,
  row.names = 1)


chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("_", "_"),
  genome = 'mm10'
)




########3以上文件在step1完成创建#####################3

obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata,
  names.field = 1, 
  names.delim = "#")

DefaultAssay(obj)
obj
Idents(obj) %>%table()

obj$nFrags_log10 <- log10(obj$nFrags)
obj$sample <- Idents(obj)


```


## Computing QC Metrics
```{r qc, fig.height=4, fig.width=12}
p=VlnPlot(
  object = obj,
  features = c('TSSEnrichment', 'nFrags_log10',
               'FRIP', 'DoubletScore'),
  pt.size = 0,
  ncol = 4
)
pdf("Computing QC Metrics.pdf",width=20,height=20)
print(p)
dev.off()



```



## save data
```{r}
saveRDS(obj, "../uuo.Rds")
```







############title:############################# 02  data_integration####################
###使用step2 scopen数据来进行降维  并通过harmony去除批次效应
#保存为/data/uuo.integrated.scOpen.Rds"


/home/data/t040413/scope_abc/second

#####安装archr包##别处复制
.libPaths(c("/home/data/t040413/R/yll/usr/local/lib/R/site-library",  "/home/data/t040413/R/x86_64-pc-linux-gnu-library/4.2", "/usr/local/lib/R/library"))

library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(harmony)


## Load data

cols.samples <- c("D0_1" = "#a6cee3",
                  "D0_2" = "#1f78b4",
                  "D2_1" = "#b2df8a",
                  "D2_2" = "#33a02c",
                  "D10_1" = "#fb9a99",
                  "D10_2" = "#e31a1c")

.libPaths(c("/home/data/t040413/R/yll/usr/local/lib/R/site-library",  "/home/data/t040413/R/x86_64-pc-linux-gnu-library/4.2", "/usr/local/lib/R/library"))

## dimension reduction scOpen and batch correction with Harmony

obj <- readRDS("../uuo.Rds")



.libPaths(c("/home/data/t040413/R/yll/usr/local/lib/R/site-library",  "/home/data/t040413/R/x86_64-pc-linux-gnu-library/4.2", "/usr/local/lib/R/library"))

library(Seurat)
Idents(obj) %>%table()
obj$sample <- Idents(obj)





##读取step2降维数据
df <- read.table("../scopen__barcodes.txt",#step创建的
                 header = TRUE, check.names = FALSE, sep = "\t",
                 comment.char="")

head(df)[,1:8]


df <- df[, colnames(obj)]
dim(df)
#[1]    20 30129

seq_len(nrow(df))
rownames(df) <- paste0("PC_", seq_len(nrow(df)))
n.features <- nrow(df)
head(df)[,1:8]



'''
333
df <- read.csv("../uuo.integrated.annotated.scOpen.csv")
rownames(df)=df$X
table(rownames(df) %in% colnames(obj))
table(colnames(obj) %in%  rownames(df))
colnames(obj) %in%  rownames(df)
obj=obj[,colnames(obj) %in%  rownames(df)]
#> table(colnames(obj) %in%  rownames(df))
 TRUE 
22419 

uuo.atac=readRDS("../uuo.integrated.annotated.scOpen.Rds")
uuo.atac@reductions$scOpen %>% dim()
#[1] 22472    20

table(colnames(obj) %in% rownames(uuo.atac@reductions$scOpen ))
obj[,colnames(obj) %in% rownames(uuo.atac@reductions$scOpen )]

333
'''

obj[['scOpen']] <- CreateDimReducObject(embeddings = as.matrix(t(df)),
                                        key = "PC_", 
                                        assay = "peaks")
DepthCor(obj, assay = "peaks", 
         reduction = "scOpen", n = n.features)

obj@reductions

library(harmony)
obj <- RunHarmony(
  object = obj,
  group.by.vars = 'sample',
  reduction = 'scOpen',
  project.dim = FALSE,
  assay.use = 'peaks',
  plot_convergence = FALSE,
  verbose = TRUE
)

library(reticulate)
use_python("/home/data/t040413/.virtualenvs/Renv/bin/python")


library(Seurat)

memory.limit(size = 999999)
'''
.libPaths(c("/home/data/t040413/R/yll/usr/local/lib/R/site-library",  "/home/data/t040413/R/x86_64-pc-linux-gnu-library/4.2", "/usr/local/lib/R/library"))
obj=readRDS("../uuo.Fibroblast.scOpen.Rds")
obj <- RunUMAP(object = obj, reduction = 'scOpen', 
               dims = 1:3, n.neighbors = 3L,
               metric = "correlation",
               umap.method = "uwot")
'''




memory.limit(size = 999999999999999999)

##conda activate five #pip install umap-learn
obj <- RunUMAP(object = obj, reduction = 'scOpen', 
               dims = 1:n.features, 
               reduction.name = "umap_harmony"
               metric = "correlation",
               umap.method = "uwot") #umap-learn

obj <- RunUMAP(obj, 
               dims = 1:n.features, 
               reduction = 'harmony',
               reduction.name = "umap_harmony",
               metric = "correlation",
               umap.method = "uwot")##umap-learn
obj@reductions

p1 <- DimPlot(object = obj, reduction = "umap",
              group.by = "sample") +
    scale_color_manual(values = cols.samples) +
    xlab("UMAP1") + ylab("UMAP2")

p2 <- DimPlot(object = obj, reduction = "umap_harmony", 
              group.by = "sample") +
    scale_color_manual(values = cols.samples) +
    xlab("UMAP1") + ylab("UMAP2")

p3 <- FeaturePlot(object = obj, reduction = "umap_harmony", 
              feature = "nFrags_log10") +
    xlab("UMAP1") + ylab("UMAP2")

p4 <- FeaturePlot(object = obj, reduction = "umap_harmony", 
              feature = "DoubletScore") +
    xlab("UMAP1") + ylab("UMAP2")

pdf("p1.pdf",width=20,height=20)
print(p1)
dev.off()

pdf("p2.pdf",width=20,height=20)
print(p2)
dev.off()

pdf("p3.pdf",width=20,height=20)
print(p3)
dev.off()

pdf("p4.pdf",width=20,height=20)
print(p4)
dev.off()

ggAlignPlots(p1, p2, type = "h")
ggAlignPlots(p3, p4, type = "h")
saveRDS(obj, "../uuo.integrated.scOpen.Rds")

```

##############################3333######03  label transfer##############
```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
```

## Load data
```{r load_data}
cols.samples <- c("D0_1" = "#a6cee3",
                  "D0_2" = "#1f78b4",
                  "D2_1" = "#b2df8a",
                  "D2_2" = "#33a02c",
                  "D10_1" = "#fb9a99",
                  "D10_2" = "#e31a1c")
cols.celltypes <- c("Pod" = "#2e3091", 
                "EC" = "#45c0ba",
                "PT(S1)" = "#8ec73e",
                "PT(S2)" = "#0eaf4f",
                "PT(S3)" = "#3ab44a",
                "Dediff. PT" = "#007d3b",
                "Prolif. PT" = "#00491b",
                "DL and tAL" = "#e0e332",
                "TAL" = "#cba214",
                "DCT" = "#F79123",
                "CNT" = "#f26221",
                "CD-PC" = "#f05325",
                "IC" = "#8f0b13",
                "Fib.1" = "#ed98b0",
                "Fib.2" = "#ca4d9c",
                "JGA" = "#d776ab",
                "Ma" = "#00ff00")

```




"/home/data/t040413/scope_abc/second"
## Load snRNA-seq####直接从官网下载rds 文件 不需要这一步
```{r load_rna}

gene.expression <- read.csv(file = "./GSE119531_UUO.dge.txt",
                            header = TRUE, row.names = 1, sep="\t")

metadata <- read.csv(file = "./GSE119531_UUO.cell.annotation.txt", 
                     header = TRUE, row.names = 1,sep="\t",encoding = "UTF-8")
metadata[metadata$CellType=="M\xcc\xfc",]

table(str_replace_all(string =metadata$CellType,pattern = "M\xcc\xfc",replacement = "Mac"))

#myrownames=rownames(metadata)

metadata$CellType=str_replace_all(string =metadata$CellType,pattern = "M\xcc\xfc",replacement = "Mac")
head(metadata)
table(metadata$CellType)
#rownames(metadata)=myrownames

head(gene.expression)[,1:7]
head(metadata)
#identical(rownames(metadata),colnames(gene.expression))

uuo.rna <- CreateSeuratObject(counts = gene.expression, 
                              assay = 'gene.expression',
                              project = 'UUO', min.cells = 1, 
                              names.field = 2,
                              names.delim = "_")

uuo.rna$CellType <- metadata$CellType
uuo.rna <- NormalizeData(uuo.rna) %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA()
```
> dim(uuo.rna)
[1] 19492  6147



## Label transfer
```{r, fig.height=6, fig.width=8}
gene.activity <- readRDS("../GeneScoreMatrix.Rds")
head(gene.activity)[,1:9]


for (dr_method in c( "scOpen")) { #"LSI", "cisTopic", , "snapATAC"
    dr_method="scOpen"
    uuo.atac <- readRDS(sprintf("../uuo.integrated.%s.Rds",
                                dr_method))
    
    uuo.atac[["RNA"]] <- CreateAssayObject(counts = gene.activity)
    
    DefaultAssay(uuo.atac) <- 'RNA'
    dim(uuo.atac)
 > dim(uuo.atac)
[1] 24333 30129   


    uuo.atac <- NormalizeData(object = uuo.atac)
    
    transfer.anchors <- FindTransferAnchors(reference = uuo.rna, 
                                        query = uuo.atac,
                                        features = VariableFeatures(uuo.rna),
                                        reference.assay = "gene.expression",
                                        query.assay = "RNA", 
                                        reduction = "cca")
    
    predicted.labels <- TransferData(anchorset = transfer.anchors,
                                     refdata = uuo.rna$CellType, #注意编码utf-8
                                     weight.reduction = uuo.atac[["harmony"]],
                                     dims = 1:ncol(uuo.atac[["harmony"]]))
    
    uuo.atac <- AddMetaData(object = uuo.atac, 
                            metadata = predicted.labels)
    
    p1 <-DimPlot(uuo.atac, group.by = "predicted.id", 
             reduction = "umap_harmony",
             label = TRUE, pt.size = 0.1, label.size = 4,
             cols = cols.celltypes) +
    xlab("UMAP1") + ylab("UMAP2") +
        ggtitle(sprintf("predicted labels: %s", dr_method))
    
    pdf("umap1_umap2.pdf",width=20,height=20)
  	print(p1)
	dev.off()

    
    saveRDS(uuo.atac, file = sprintf("../uuo.integrated.%s.Rds", dr_method))
}



```








#######################################33########04  clustering_k_medoids
#使用##########################################scopen数据来聚类
```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
```

## Load data
```{r load_data}
cols.celltypes <- c("Pod" = "#2e3091", 
                "EC" = "#45c0ba",
                "PT(S1)" = "#8ec73e",
                "PT(S2)" = "#0eaf4f",
                "PT(S3)" = "#3ab44a",
                "Dediff. PT" = "#007d3b",
                "Prolif. PT" = "#00491b",
                "DL and tAL" = "#e0e332",
                "TAL" = "#cba214",
                "DCT" = "#F79123",
                "CNT" = "#f26221",
                "CD-PC" = "#f05325",
                "IC" = "#8f0b13",
                "Fib.1" = "#ed98b0",
                "Fib.2" = "#ca4d9c",
                "JGA" = "#d776ab",
                "Ma" = "#00ff00")
```


## Clustering with K-medoids
```{r, fig.width=12, fig.height=5}
library(parallel)
library(doParallel)
library(foreach)
cl <- parallel::makeForkCluster(40)
doParallel::registerDoParallel(cl)

foreach(dr_method = c( "scOpen")) %dopar% {#"LSI", "cisTopic", , "snapATAC"
    
    dr_method = c( "scOpen")
    uuo.atac <- readRDS(sprintf("../uuo.integrated.%s.Rds",
                                dr_method))
    
    # here we use 1 - cor as distance metric
    
    df.dist <- 1 - cor(as.matrix(t(uuo.atac@reductions$harmony@cell.embeddings)))
   

    pa <- pam(df.dist, k = 17, diss = TRUE, cluster.only = TRUE)
    df_cluster <- as.data.frame(pa)
    




    colnames(df_cluster) <- c("pam_clusters")
    uuo.atac <- AddMetaData(uuo.atac, metadata = df_cluster)
   
    p1 <- DimPlot(uuo.atac, reduction = "umap_harmony",
              label = TRUE, label.size = 4,
              group.by = "pam_clusters")
    
    p2 <-DimPlot(uuo.atac, group.by = "predicted.id", 
             reduction = "umap_harmony",
             label = TRUE, pt.size = 0.1, label.size = 4,
             cols = cols.celltypes) +
    xlab("UMAP1") + ylab("UMAP2") +
        ggtitle(sprintf("predicted labels: %s", dr_method))
    
    
 	pdf("04pam_clusters1.pdf",width=20,height=20)
	print(p1)
	dev.off()

 	pdf("04pam_clusters2.pdf",width=20,height=20)
	print(p2)
	dev.off()

 
    ggAlignPlots(p1, p2, type = "h")
    saveRDS(uuo.atac, file = sprintf("../uuo.integrated.%s.Rds", dr_method))

}
parallel::stopCluster(cl)
```


## visualize
```{r, fig.height=6, fig.width=8}
for(dr_method in c("LSI", "cisTopic", "scOpen", "snapATAC")) {
    

    dr_method="scOpen"
    uuo.atac <- readRDS(sprintf("../data/uuo.integrated.%s.Rds",
                                dr_method))
    
    p1 <- DimPlot(uuo.atac, reduction = "umap_harmony",
              label = TRUE, label.size = 4,
              group.by = "pam_clusters")
    
 	pdf("04pam_clustersviz.pdf",width=20,height=20)
	print(p1)
	dev.off()



}
```
#######################################################33## 05 cluster_k_medoids
###########################33#降维之后 画图 umap 和tsne

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
```

## Load data
```{r load_data}
cols.celltypes <- c("Pod" = "#2e3091", 
                "EC" = "#45c0ba",
                "PT(S1)" = "#8ec73e",
                "PT(S2)" = "#0eaf4f",
                "PT(S3)" = "#3ab44a",
                "Dediff. PT" = "#007d3b",
                "Prolif. PT" = "#00491b",
                "DL and tAL" = "#e0e332",
                "TAL" = "#cba214",
                "DCT" = "#F79123",
                "CNT" = "#f26221",
                "CD-PC" = "#f05325",
                "IC" = "#8f0b13",
                "Fib.1" = "#ed98b0",
                "Fib.2" = "#ca4d9c",
                "JGA" = "#d776ab",
                "Ma" = "#00ff00")
```


## Clustering with K-medoids
```{r, fig.width=12, fig.height=5}
library(parallel)
library(doParallel)
library(foreach)
cl <- parallel::makeForkCluster(5)
doParallel::registerDoParallel(cl)

foreach(dr_method = c("scOpen")) %dopar% {#"LSI", "cisTopic", , "snapATAC"

    uuo.atac <- readRDS(sprintf("../uuo.integrated.%s.Rds",
                                dr_method))
    
    # here we use 1 - cor as distance metric
    df.dist <- 1 - cor(as.matrix(t(uuo.atac@reductions$harmony@cell.embeddings)))
   
    pa <- pam(df.dist, k = 17, diss = TRUE, cluster.only = TRUE) #计算距离
    
    df_cluster <- as.data.frame(pa)
    
    colnames(df_cluster) <- c("pam_clusters")
    
    uuo.atac <- AddMetaData(uuo.atac, metadata = df_cluster)###最好运行下
    
    p1 <- DimPlot(uuo.atac, reduction = "umap_harmony",
              label = TRUE, label.size = 4,
              group.by = "pam_clusters")
    
    p2 <-DimPlot(uuo.atac, group.by = "predicted.id", 
             reduction = "umap_harmony",
             label = TRUE, pt.size = 0.1, label.size = 4,
             cols = cols.celltypes) +
    
    xlab("UMAP1") + ylab("UMAP2") +
        ggtitle(sprintf("predicted labels: %s", dr_method))
    
pdf(" 05 cluster_k_medoids1.pdf",width=20,height=20)
print(p1)
dev.off()

pdf(" 05 cluster_k_medoids2.pdf",width=20,height=20)
print(p2)
dev.off()


    ggAlignPlots(p1, p2, type = "h")
    
    saveRDS(uuo.atac, file = sprintf("../uuo.integrated.%s.Rds", dr_method))
}
parallel::stopCluster(cl)








```


## visualize
```{r, fig.height=6, fig.width=8}
for(dr_method in c("scOpen")) {#"LSI", "cisTopic", , "snapATAC"
    uuo.atac <- readRDS(sprintf("../data/uuo.integrated.%s.Rds",
                                dr_method))
    
    p1 <- DimPlot(uuo.atac, reduction = "umap_harmony",
              label = TRUE, label.size = 4,
              group.by = "pam_clusters")
    print(p1)

    pdf(" 05 cluster_vis.pdf",width=20,height=20)
	print(p1)
	dev.off()
}
```


## ARI
##可以不做
```{r}
cols <- c("scOpen" = "#e31a1c",
          "cisTopic" = "#33a02c",
          "SnapATAC" = "#386cb0",
          "Cusanovich2018" = "#984ea3")
df_list <- list()
idx <- 1
for (dr_method in c("LSI", "cisTopic", "scOpen", "snapATAC")) {
    uuo.atac <- readRDS(sprintf("../data/uuo.integrated.%s.Rds",
                                dr_method))
    meta.data <- as.data.frame(uuo.atac@meta.data)
    meta.data$sample2 <- stringr::str_replace_all(meta.data$sample,
                                                  c("D0_1" = "Day 0",
                                                    "D0_2" = "Day 0",
                                                    "D2_1" = "Day 2",
                                                    "D2_2" = "Day 2",
                                                    "D10_1" = "Day 10",
                                                    "D10_2" = "Day 10"))
    
    
    ARI <- adjustedRandIndex(meta.data$pam_clusters,
                             meta.data$predicted.id)
    
    df_list[[idx]] <- data.frame(ARI = ARI,
                                 dr_method = dr_method,
                                 data = "Integrated")
    idx <- idx + 1
    for (sample in c("Day 0", "Day 2", "Day 10")) {
        meta.data.sub <- meta.data[meta.data$sample2 == sample, ]
        ARI <- adjustedRandIndex(meta.data.sub$pam_clusters,
                                 meta.data.sub$predicted.id)
        df_list[[idx]] <- data.frame(ARI=ARI,
                                        dr_method=dr_method,
                                        data = sample)
        idx <- idx + 1
    }
    
}
df <- Reduce(rbind, df_list)
df$dr_method <- stringr::str_replace_all(df$dr_method, c("LSI" = "Cusanovich2018",
                                                         "snapATAC" = "SnapATAC"))
p <- df %>%
    ggplot(aes(x = data, y = ARI)) + 
    geom_bar(aes(fill = dr_method), position = "dodge", stat = "identity") +
    scale_fill_manual(values = cols) +
    xlab("") + ylab("ARI") +
    theme_cowplot() +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1))
p
```

## Ranking
```{r, fig.height=4, fig.width=3}
df <- group_by(df, data) %>%
    mutate(rank = rank(ARI))
p <- ggplot(data = df, aes(x = reorder(dr_method, -rank, FUN = median), 
                               y = rank, fill = dr_method)) +
        geom_boxplot() +
    scale_fill_manual(values = cols) +
        xlab("") + ylab("") +
        theme_cowplot() +
        theme(legend.title = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(angle = 60, hjust = 1),
              legend.position = "none")
p
```

########################################06 silhousette_score
#计算house得分 距离得分  评价聚类好坏

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
```

## Load data
```{r load_data}
cols.celltypes <- c("Pod" = "#2e3091", 
                "EC" = "#45c0ba",
                "PT(S1)" = "#8ec73e",
                "PT(S2)" = "#0eaf4f",
                "PT(S3)" = "#3ab44a",
                "Dediff. PT" = "#007d3b",
                "Prolif. PT" = "#00491b",
                "DL and tAL" = "#e0e332",
                "TAL" = "#cba214",
                "DCT" = "#F79123",
                "CNT" = "#f26221",
                "CD-PC" = "#f05325",
                "IC" = "#8f0b13",
                "Fib.1" = "#ed98b0",
                "Fib.2" = "#ca4d9c",
                "JGA" = "#d776ab",
                "Ma" = "#00ff00")
cols <- c("scOpen" = "#e31a1c"
          )#,"cisTopic" = "#33a02c","SnapATAC" = "#386cb0","Cusanovich2018" = "#984ea3"
```


## Compute silhouette score
```{r, fig.width=12, fig.height=5}

for (dr_method in c("scOpen")) {#"LSI", "cisTopic", , "snapATAC"
    
    uuo.atac <- readRDS(sprintf("../data/uuo.integrated.%s.Rds",
                                dr_method))
    
    df.dist <- 1 - cor(as.matrix(t(uuo.atac@reductions$harmony@cell.embeddings)))
    meta.data <- uuo.atac@meta.data
    
    meta.data$sample2 <- stringr::str_replace_all(meta.data$sample,
                                                  c("D0_1" = "Day_0",
                                                    "D0_2" = "Day_0",
                                                    "D2_1" = "Day_2",
                                                    "D2_2" = "Day_2",
                                                    "D10_1" = "Day_10",
                                                    "D10_2" = "Day_10"))
    
    si <- silhouette(x = as.numeric(factor(meta.data$predicted.id)), 
                     dmatrix = df.dist)
    
    saveRDS(si, file = sprintf("../data/si.%s.Rds",
                                dr_method))
    
    for (sample in c("Day_0", "Day_2", "Day_10")) {
        meta.data.sub <- meta.data[meta.data$sample2 == sample, ]
        cell.embedding <- uuo.atac@reductions$harmony@cell.embeddings[rownames(meta.data.sub), ]
        
        df.dist <- 1 - cor(t(cell.embedding))
        si <- silhouette(x = as.numeric(factor(meta.data.sub$predicted.id)), 
                     dmatrix = df.dist)
        
        saveRDS(si, file = sprintf("../data/si.%s.%s.Rds",
                                dr_method, sample))
    }
}
```

## plot
```{r}
df_list <- list()
idx <- 1

for (dr_method in c("scOpen")) {#, "snapATAC"  , "LSI", "cisTopic" 
    obj <- readRDS(sprintf("../data/si.%s.Rds", dr_method))
    si <- mean(obj[, 'sil_width'])
    
    df_list[[idx]] <- data.frame(si = si,
                                 dr_method = dr_method,
                                 data = "Integrated")
    idx <- idx + 1
    
    for (sample in c("Day_0", "Day_2", "Day_10")) {
        obj <- readRDS(sprintf("../data/si.%s.%s.Rds", dr_method, sample))
        si <- mean(obj[, 'sil_width'])
        
        df_list[[idx]] <- data.frame(si=si,
                                        dr_method=dr_method,
                                        data = sample)
        idx <- idx + 1
    }
}

df <- Reduce(rbind, df_list)
#df$dr_method <- stringr::str_replace_all(df$dr_method, c("LSI" = "Cusanovich2018",
#                                                        "snapATAC" = "SnapATAC"))

p <- df %>%
    ggplot(aes(x = data, y = si)) + 
    geom_bar(aes(fill = dr_method), position = "dodge", stat = "identity") +
    scale_fill_manual(values = cols) +
    xlab("") + ylab("Silhouette score") +
    theme_cowplot() +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1))
p
```

## Ranking
```{r, fig.height=4, fig.width=3}
df <- group_by(df, data) %>%
    mutate(rank = rank(si))
df$dr_method <- factor(df$dr_method, levels = c("scOpen"
                                                ))# , "cisTopic","SnapATAC", "Cusanovich2018"
p <- ggplot(data = df, aes(x = reorder(dr_method, -rank, FUN = median), 
                               y = rank, fill = dr_method)) +
        geom_boxplot() +
    scale_fill_manual(values = cols) +
        xlab("") + ylab("") +
        theme_cowplot() +
        theme(legend.title = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(angle = 60, hjust = 1),
              legend.position = "none")
p
```

#############################3##############07 viz_markers
#画图  可以不画
if (F){```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
```


## Viz markers
```{r, fig.width=20, fig.height=4}
library(factoextra)
library(cluster)
for (dr_method in c("scOpen")) { #,"LSI", "cisTopic", "snapATAC" 
 

 dr_method="scOpen"
    uuo.atac <- readRDS(sprintf("../uuo.integrated.%s.Rds",
                                dr_method))
    
    DefaultAssay(uuo.atac) 
    dim(uuo.atac)
    DefaultAssay(uuo.atac) <- "RNA"
    
    plot_list <- lapply(c("Scara5", "Fbln2", "Dcn", "Fn1", "Col1a1"), function(x){
        p <- FeaturePlot(uuo.atac, features = x, 
                     reduction = "umap_harmony",
                     pt.size = 0.1,
                     min.cutoff = "q1",
                     max.cutoff = "q99") +
            xlab("") + ylab("") +
         ggtitle(sprintf("%s", x)) +
         NoLegend() +
            theme(axis.text = element_blank(),
                  axis.ticks = element_blank())
        
        
        pdf(paste0("07_vizmarkers",x, ".pdf"),width=15,height=15)
		print(p)
		dev.off()

    }  )
    
    ArchR::ggAlignPlots(plotList = plot_list, type = "h")
}

```}


##########################3333333 08  find_markers#####33333
#总群的marker基因寻找 需要archr对象

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 10)
addArchRGenome("mm10")
```


```{r}
proj <- loadArchRProject(path = "../UUO/")

uuo.atac <- readRDS("../uuo.integrated.scOpen.Rds")
uuo.integrated.annotated.scOpen=readRDS("../uuo.integrated.annotated.scOpen_forminternet.Rds")
dim(uuo.atac)
dim(uuo.integrated.annotated.scOpen)



 uuo.atac@meta.data %>%head()
 uuo.integrated.annotated.scOpen@meta.data %>%head()

table(colnames(uuo.atac) %in% colnames(uuo.integrated.annotated.scOpen))
uuo.atac=uuo.atac[,colnames(uuo.atac) %in% colnames(uuo.integrated.annotated.scOpen)]



from.inte.meta=uuo.integrated.annotated.scOpen@meta.data 

dim(from.inte.meta[rownames(from.inte.meta) %in% colnames(uuo.atac),])

from.inte.meta=from.inte.meta[rownames(from.inte.meta) %in% colnames(uuo.atac),]
dim(from.inte.meta)
identical(rownames(from.inte.meta),colnames(uuo.atac))
dim(uuo.atac)



uuo.atac <- AddMetaData(object = uuo.atac, 
                            metadata = from.inte.meta)
uuo.atac@meta.data %>%head()
uuo.integrated.annotated.scOpen@meta.data %>%head()



df <- uuo.atac@meta.data
#df <- df[rownames(proj), ]####??????????

dim(df)
getCellColData(proj)
dim(getCellColData(proj))


proj=proj[rownames(uuo.atac@meta.data),]

dim(proj[rownames(uuo.atac@meta.data),])
 proj@cellColData %>%dim()

identical(rownames(proj@cellColData),rownames(df))

proj <- addCellColData(proj, data = as.character(df$pam_clusters),
                       cells = rownames(proj@cellColData),
                       name = "pam_clusters",
                       force = TRUE)

proj <- addCellColData(proj, data = as.character(df$predicted.id),
                       cells = rownames(proj@cellColData),
                       name = "predicted.id",
                       force = TRUE)
```

## confusion matrix
```{r}
library(tidyr)
library(circlize)
library(ComplexHeatmap)

df.plot <- df %>%
    group_by(pam_clusters, predicted.id) %>%
   # summarise() %>%
    mutate(counts = n(),frac=round(counts/sum(counts), 4)) %>%
    subset(., select = c("pam_clusters", "predicted.id", "frac")) %>%
    spread(., predicted.id, 
           value = frac, fill = 0)
rownames(df.plot) <- paste0("C", df.plot$pam_clusters) 

df.plot$pam_clusters <- NULL
mat <- as.matrix(t(df.plot))
colnames(mat) <- rownames(df.plot)
cols <- colorRamp2(c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), 
                   c("white", "#ffffd9", '#edf8b1', "#c7e9b4", "#7fcdbb", "#41b6c4",
                     "#1d91c0", "#225ea8", "#253494", "#081d58"))
p <- Heatmap(mat,
             name = "Perc.",
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             show_column_names = TRUE,
             show_row_names = TRUE,
             row_order = c("PT(S2)", "TAL", "PT(S3)", "DCT", "PT(S1)",
                           "CD-PC", "EC", "Fib.1", "IC", "CNT",
                           "DL and tAL", "Ma", "Pod", "JGA", "Fib.2",
                           "Dediff. PT"),
             row_names_side = "left",
             column_names_side = "bot",
             col = cols)
draw(p)
```


# find marker
```{r, fig.height=8, fig.width=6}
known.markers <- c("Flt1", "Erg", # EC
                   "Tmem117", # IC
                   "Slc14a2", # DL & TAL
                   "C1qa", "C1qb", # Mac
                   "Cd1d", "Lta", "Ltb", # Lymphoid
                   "Aqp4", # CD-PC
                   "Slc12a1", "Umod", # TAL
                   "Slc12a3", # DCT
                   "Scnn1g", # CNT
                   "Pdgfrb", # Fibroblast
                   "Fbln2", "Dcn", # MyoFibroblast
                   "Wt1", "Nphs1", "Nphs2", # Pod
                   "Vcam1", "Havcr1", # Injured PT
                   "Lrp2", # PT (S1)
                   "Slc5a12", # PT (S2)
                   "Slc27a2"# PT (S3)
                   )
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "pam_clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    verbose = FALSE
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = known.markers
)

ComplexHeatmap::draw(heatmapGS, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot",
                     show_column_dend = FALSE)

markerList <- lapply(markerList, as.data.frame)
saveRDS(markersGS, file = "../MarkerGenesForClusters.Rds")

for(i in 1:length(markerList)){
    markerList[[i]] <- markerList[[i]][order(-markerList[[i]]$Log2FC), ]
}

WriteXLS(markerList,
         ExcelFileName = "../MarkerGenesForClusters.xlsx",
         SheetNames = names(markerList))
```


## save data
```{r}
saveArchRProject(ArchRProj = proj, 
                 load = FALSE)
```


#
##############################33####09 add_annatiaion##################################
#给总群添加注释  去双细胞 去不要的群 得到 改名字


#write.csv(df, file = "../data/uuo.integrated.annotated.scOpen.csv")
#saveRDS(uuo.atac, "../data/uuo.integrated.annotated.scOpen.Rds")
```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 10)
addArchRGenome("mm10")
```


```{r}
uuo.atac <- readRDS("../uuo.integrated.scOpen.Rds")

table(Idents(uuo.atac))
Idents(uuo.atac) <- "pam_clusters"
table(Idents(uuo.atac))
dim(uuo.atac)

```


## qc per cluster
```{r qc_cluster, fig.width=8, fig.height=6}
p1 <- VlnPlot(uuo.atac, features = "TSSEnrichment", pt.size = 0) +
    NoLegend()
    
p2 <- VlnPlot(uuo.atac, features = "nFrags_log10", pt.size = 0) +
    NoLegend()

p3 <- VlnPlot(uuo.atac, features = "DoubletEnrichment", pt.size = 0) +
    NoLegend()
p1
p2
p3

df <- as.data.frame(uuo.atac@meta.data)
p <- ggplot(data = df, aes(x = DoubletEnrichment)) +
    geom_histogram(bins = 100) +
    geom_vline(xintercept = 2.5, color = "red")
#p <- hist(uuo.atac$DoubletEnrichment, breaks = 200)
print(p)


```


## subset
```{r subset}
uuo.atac <- subset(uuo.atac, idents = '4', invert = TRUE)
# remove cells with high doublet score
uuo.atac <- subset(uuo.atac, DoubletEnrichment < 2.5)
```

## add new umap after cell filtering
```{r}
uuo.atac@meta.data %>%head()
uuo.atac <- RunUMAP(uuo.atac, 
               dims = 1:20, 
               reduction = 'harmony',
               reduction.name = "umap_harmony_v2",
               metric = "correlation",
               umap.method = "uwot") ####
```

dim(uuo.atac)

table(Idents(uuo.atac))
## add annotation
```{r add_annotation, fig.height=6, fig.width=6}
newLabels <- c("1" = "PT(S2)",
               "2" = "TAL",
               "3" = "PT(S3)",
               "5" = "DCT",
               "6" = "PT(S1)",
               "7" = "CD-PC",
               "8" = "EC",
               "9" = "CNT",
               "10" = "IC",
               "11" = "Fibroblast",
               "12" = "DL & TAL",
               "13" = "Mac",
               "14" = "Injured PT",
               "15" = "Lymphoid",
               "16" = "Pod",
               "17" = "TAL")
uuo.atac <- RenameIdents(uuo.atac, newLabels)

uuo.atac$CellType <- as.character(Idents(uuo.atac))

pal <- ArchR::paletteDiscrete(values = uuo.atac$CellType)

p1 <- DimPlot(uuo.atac, reduction = "umap_harmony_v2",
              label = TRUE,
              group.by = "CellType", pt.size = 0.1) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    scale_color_manual(values = pal) +
    theme(axis.ticks = element_blank(),
          axis.text = element_blank()) +
    ggtitle("")


pdf("_all_umap.pdf",width=10,height=10)
print(p1)
dev.off()
```


## heatmap
```{r heatmap, fig.height=6, fig.width=12}
known.markers <- c("Aqp2", "Aqp4", # CD-PC,
                   "Scnn1g", # CNT
                   "Slc12a3", # DCT
                   "Flt1", "Erg", # EC
                   "Pdgfrb", "Col3a1",  # Fibroblast
                   "Tmem117", # IC
                   "Vcam1", "Havcr1", # Injured PT
                   "Slc14a2", # DL & TAL
                   "Cd3d", "Lta", "Ltb", # Lymphoid
                   "C1qa", "C1qb", # Mac
                   "Wt1", "Nphs1", "Nphs2", # Pod
                   "Slc5a2", # PT (S1)
                   "Slc17a3", # PT (S2)
                   "Slc27a2", # PT (S3),
                   "Slc12a1", "Umod" # TAL
                   )
uuo.atac <- uuo.atac %>%
    ScaleData(features = known.markers)
p=DoHeatmap(uuo.atac, 
          features = known.markers,
          group.by = "CellType",
          disp.min = -1,
          disp.max = 1.5,
          group.colors = pal,
          size = 4)


pdf("heatmap_all.pdf",width=6,height=12)
print(p)
dev.off()


```

## save data
```{r}
df <- uuo.atac@meta.data
write.csv(df, file = "../uuo.integrated.annotated.scOpen.csv")
saveRDS(uuo.atac, "../uuo.integrated.annotated.scOpen.Rds")
```



######################################10 cell_propotion

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
```


```{r}
df <- read.csv(file="../uuo.integrated.annotated.scOpen.csv")
pal <- ArchR::paletteDiscrete(values = df$CellType)
df$Sample <- stringr::str_replace_all(df$Sample, c("D0_1" = "Day 0",
                                                   "D0_2" = "Day 0",
                                                   "D2_1" = "Day 2",
                                                   "D2_2" = "Day 2",
                                                   "D10_1" = "Day 10",
                                                   "D10_2" = "Day 10"))
```

```{r, fig.height=6, fig.width=10}
df.plot <- df %>%
    group_by(Sample, CellType) %>%
    summarise(counts = n()) %>%
    mutate(frac = counts / sum(counts))

df.plot <- df %>% ########我改动的
    group_by(Sample, CellType,counts = n()) %>%
    mutate(frac = counts / sum(counts))



df.plot$Sample <- factor(df.plot$Sample, levels = c("Day 0", "Day 2", "Day 10"))

p <- ggplot(data = df.plot, aes(x = Sample, y = frac, fill = CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    facet_wrap(~CellType, nrow = 3, scales = "free_y") +
    xlab("") + ylab("Fraction of cells") +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
pdf("fraction.pdf",width=20,height=20)
print(p)
dev.off()

```

##############################################333#11 fibroblast_subset
#从总群取出成纤维子集  并把从seurat取出的成纤维放入ArchRproject对象
#把seurat对象转换为archrproject对象转换  取子集

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
```


library(Seurat)
library(dplyr)
```{r load_data}
uuo.atac <- readRDS("../uuo.integrated.annotated.scOpen.Rds")


DefaultAssay(uuo.atac)
assay(uuo.atac)

Idents(uuo.atac) %>%table()
Idents(uuo.atac) <- "CellType"

GetAssayData(uuo.atac, slot = "data")


uuo.atac.fib <- subset(uuo.atac, idents = 'Fibroblast')
saveRDS(uuo.atac.fib, "../uuo.Fibroblast.scOpen.Rds")



addArchRThreads(threads = 40)
addArchRGenome("mm10")



proj <- loadArchRProject(path = "../UUO")
# select all fibroblasts

proj <- subsetArchRProject(proj, 
                           cells = colnames(uuo.atac.fib),
                           outputDirectory = "../Fibroblast",
                           force = TRUE)


proj@cellColData$Sample2 <- stringr::str_replace_all(proj@cellColData$Sample,
                                                    c("D0_1" = "Day 0",
                                                   "D0_2" = "Day 0",
                                                   "D2_1" = "Day 2",
                                                   "D2_2" = "Day 2",
                                                   "D10_1" = "Day 10",
                                                   "D10_2" = "Day 10"))

proj@cellColData$Sample2 <- factor(proj@cellColData$Sample2, levels = c("Day 0",
                                                      "Day 2",
                                                      "Day 10"))
cols.samples <- c("Day 0" = "#e41a1c",
                  "Day 2" = "#377eb8",
                  "Day 10" = "#4daf4a")
```


## save data
```{r}
saveArchRProject(ArchRProj = proj, 
                 load = FALSE)
```


######################################12 fibroblast clustering#####33####################3
#子集聚类 成纤维聚类
##转换 格式


```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
library(harmony)
```


```{r load_data}
#obj <- readRDS("../uuo.Fibroblast.scOpen.Rds")

obj=readRDS("../uuo.Fibroblast.scOpen.Rds") #uuo.atac.fib

proj <- loadArchRProject(path = "../Fibroblast")
```


## add umap and dim. reduction to ArchR project
```{r, fig.width=10, fig.height=4}


proj@reducedDims[['scOpen']] <- SimpleList(matDR = obj@reductions$scOpen@cell.embeddings,
                                      params = NULL, date = Sys.time(), scaleDims = NA, 
        corToDepth = NA)
proj@reducedDims[['harmony']] <- SimpleList(matDR = obj@reductions$harmony@cell.embeddings,
                                      params = NULL, date = Sys.time(), scaleDims = NA, 
        corToDepth = NA)




proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "harmony", 
    name = "UMAP_Fib", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
)

embedding <- as.data.frame(obj@reductions$umap_harmony_v2@cell.embeddings)

colnames(embedding) <- c("Harmony#UMAP_Dimension_1",
                         "Harmony#UMAP_Dimension_2")

proj@embeddings[['umap_harmony']] <- SimpleList(df = embedding,
                                      params = NULL)
p1 <- plotEmbedding(ArchRProj = proj,
                   embedding = "umap_harmony",
                   colorBy = "cellColData",
                   name = "Sample2",
                   plotAs = "points",
                   labelAsFactors = FALSE) +
    theme_cowplot() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle("")

p2 <- plotEmbedding(ArchRProj = proj,
                   embedding = "UMAP_Fib",
                   colorBy = "cellColData",
                   name = "Sample2",
                   plotAs = "points",
                   labelAsFactors = FALSE) +
    theme_cowplot() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle("")

plotPDF(p1,p2, 
	name = "umap_fib.pdf", ArchRProj = proj, 
	addDOC = FALSE, width = 8, height = 8) ##成纤维画图

ggAlignPlots(p1, p2, type = "h")
```


```{r}
proj <- addImputeWeights(proj, reducedDims = "harmony") #赋予权重 更好看

p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Col1a1", 
    embedding = "UMAP_Fib",
    quantCut = c(0.01, 0.95),
) +theme_cowplot()

p2 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Postn", 
    embedding = "UMAP_Fib",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()

p3 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Scara5", 
    embedding = "UMAP_Fib",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()
p4 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Pdgfrb", 
    embedding = "UMAP_Fib",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()

p5 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Notch3", 
    embedding = "UMAP_Fib",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()

p6 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Rgs5", 
    embedding = "UMAP_Fib",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()

p7 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Cspg4", 
    embedding = "UMAP_Fib",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()

plotPDF(p4,p5,p6,p7,p3,p1,p2,
	name = "umap_fib_imputed.pdf", ArchRProj = proj, 
	addDOC = FALSE, width = 48, height = 48) ##成纤维画图




```
## clustering
```{r}
proj <- addClusters(
    input = proj,
    reducedDims = "harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8,
    force = TRUE,
    annoy.metric = "cosine",
    maxClusters = 3
)
p1 <- plotEmbedding(ArchRProj = proj,
                   embedding = "UMAP_Fib",
                   colorBy = "cellColData",
                   name = "Clusters",
                   plotAs = "points",
                   labelAsFactors = FALSE) +
    theme_cowplot() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle("")
p2 <- plotGroups(ArchRProj = proj, groupBy = "Clusters", colorBy = "cellColData",
                 name = "DoubletEnrichment", plotAs = "violin")

plotPDF(p1,p2,
	name = "umap_fib_imputed_cluster.pdf", ArchRProj = proj, 
	addDOC = FALSE, width = 18, height = 18) ##成纤维画图



```



## Diffusiion map
```{r, fig.width=10, fig.height=4}
library(destiny)
#matDR <- proj@embeddings$UMAP_Fib$df
matDR <- proj@reducedDims$harmony$matDR
dm <- DiffusionMap(as.matrix(matDR),
                   verbose = TRUE)

plot(dm)

embedding <- as.data.frame(dm)[, c("DC1", "DC2")]
colnames(embedding) <- c("Harmony#DC_Dimension_1",
                         "Harmony#DC_Dimension_2")
proj@embeddings[["dm"]] <- SimpleList(df = as.data.frame(embedding),
                                      params = NULL)


p1 <- plotEmbedding(ArchRProj = proj,
                   embedding = "dm",
                   colorBy = "cellColData",
                   name = "Sample2",
                   plotAs = "points",
                   labelAsFactors = FALSE) +
    theme_cowplot() +
    xlab("DC 1") + ylab("DC 2") +
    ggtitle("")



p2 <- plotEmbedding(ArchRProj = proj,
                   embedding = "dm",
                   colorBy = "cellColData",
                   name = "Clusters",
                   plotAs = "points",
                   labelAsFactors = FALSE) +
    theme_cowplot() +
    xlab("DC 1") + ylab("DC 2") +
    ggtitle("")
ggAlignPlots(p1, p2, type = "h")

```

```{r}

p1 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Twist1", 
    embedding = "dm",
    quantCut = c(0.01, 0.95)
) +theme_cowplot()

pdf("Twist1_DIFUSION_MAP.pdf",width=20,height=20)
print(p1)
dev.off()





p2 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Postn", 
    embedding = "dm",
    quantCut = c(0.01, 0.95)
) +
    theme_cowplot()



p3 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Scara5", 
    embedding = "dm",
    quantCut = c(0.01, 0.95)
) +
    theme_cowplot()
p4 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Pdgfrb", 
    embedding = "dm",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()
p5 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Notch3", 
    embedding = "dm",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()


p6 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Rgs5", 
    embedding = "dm",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()

p7 <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "GeneScoreMatrix", 
    name = "Cspg4", 
    embedding = "dm",
    quantCut = c(0.01, 0.95),
) +
    theme_cowplot()



```


## find marker genes
```{r, fig.height=8, fig.width=6}

markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon",
    verbose = FALSE
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)

ComplexHeatmap::draw(heatmapGS, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot",
                     show_column_dend = FALSE)

markerList <- lapply(markerList, as.data.frame)
saveRDS(markersGS, file = "../Fibroblast/MarkerGenesForClusters.Rds")




for(i in 1:length(markerList)){
    markerList[[i]] <- markerList[[i]][order(-markerList[[i]]$Log2FC), ]
}
WriteXLS(markerList,
         ExcelFileName = "../Fibroblast/MarkerGenesForClusters.xlsx",
         SheetNames = names(markerList))

```


## add annotation
```{r, fig.height=4, fig.width=6}


proj@cellColData$CellType <- stringr::str_replace_all(proj@cellColData$Clusters,
                                                      c("C1" = "Pericytes",
                                                        "C2" = "Myofibroblast",
                                                        "C3" = "Fib_Scara5"))
df <- as.data.frame(proj@cellColData)
pal <- ArchR::paletteDiscrete(values = df$CellType)


df.plot <- df %>%
    group_by(Sample2, CellType) %>%
    summarise(counts = n()) %>%
    mutate(frac = counts / sum(counts))
df.plot$Sample2 <- factor(df.plot$Sample2, levels = c("Day 0", "Day 2", "Day 10"))

p <- ggplot(data = df.plot, aes(x = Sample2, y = frac, fill = CellType)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal) +
    facet_wrap(~CellType, nrow = 1, scales = "free_y") +
    xlab("") + ylab("Fraction of cells") +
    theme_cowplot() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))

pdf("Annotation_fib_map.pdf",width=10,height=10)
print(p)
dev.off()




```

## Peak calling
```{r peak_calling}

proj <- addGroupCoverages(ArchRProj = proj, 
                          groupBy = "pam_clusters",
                          force = TRUE)

#conda activate my_masc2
#pathToMacs2 <- findMacs2()
#
pathToMacs2="/home/data/t040413/anaconda3/envs/my_masc2/bin/macs2"

proj <- addReproduciblePeakSet( #######ZUIHAOO YIBUBU YUNXING
    ArchRProj = proj, 
    groupBy = "pam_clusters", 
    pathToMacs2 = pathToMacs2
)

getPeakSet(proj)

proj <- addPeakMatrix(proj)

getAvailableMatrices(proj)

peakMatrix <- getMatrixFromProject(proj,
                                   useMatrix = "PeakMatrix")

counts <- peakMatrix@assays@data$PeakMatrix

df_rangers <- as.data.frame(peakMatrix@rowRanges@ranges)

rownames(counts) <- paste(peakMatrix@rowRanges@seqnames,
                          df_rangers$start,
                          df_rangers$end,
                          sep = "_") 

saveRDS(counts, file = "../Fibroblast/PeakMatrix.Rds")


if(!dir.exists("../Fibroblast/filtered_peak_bc_matrix")){
    dir.create("../Fibroblast/filtered_peak_bc_matrix")
}

writeMM(counts, file = "../Fibroblast/filtered_peak_bc_matrix/matrix.mtx")
barcodes <- as.data.frame(colnames(counts))

peaks <- as.data.frame(stringr::str_split_fixed(rownames(counts), "_", 3))

write.table(barcodes, file = "../Fibroblast/filtered_peak_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(peaks, file = "../Fibroblast/filtered_peak_bc_matrix/peaks.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


```


## add motif matrix
```{r}
#devtools::install_github("GreenleafLab/chromVARmotifs")

#Error: peakSet is NULL. You need a peakset to run addMotifAnnotations! See addReproduciblePeakSet!

proj <- addMotifAnnotations(ArchRProj = proj, 
                                 motifSet = "cisbp", 
                                 name = "Motif",
                            force = TRUE)
proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)


```

## save data
```{r}
saveArchRProject(ArchRProj = proj, 
                 load = FALSE)
```



#######################################################13轨迹分析###############



if(F){
	
```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 80)
addArchRGenome("mm10")
```


```{r load_data}
proj <- loadArchRProject(path = "../Fibroblast")
```


## add trajectory
```{r, fig.height=6, fig.width=6}
trajectory <- c("C3", "C2")

proj <- addTrajectory(
    ArchRProj = proj, 
    name = "Fib", 
    groupBy = "Clusters",
    trajectory = trajectory, 
    preFilterQuantile = 0.9,
    postFilterQuantile = 0.95,
    useAll = TRUE,
    embedding = "dm", 
    force = TRUE
)
p <- plotTrajectory(proj, trajectory = "Fib", 
                    colorBy = "cellColData",
                    continuousSet = "blueYellow",
                    name = "Fib",
                    embedding = "dm",
                    plotAs = "points",
                    rastr = FALSE,
                    smoothWindow = 35,
                    size = 1)


p1 <- p[[1]] + 
    theme_cowplot() +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle("")

pdf("tracy.pdf",width=10,height=10)
print(p1)
dev.off()

p1
```

##
```{r, fig.height=10, fig.width=6}
trajMM  <- getTrajectory(ArchRProj = proj,
                         name = "Fib",
                         useMatrix = "MotifMatrix",
                         log2Norm = FALSE)

trajMM <- trajMM[!grepl("deviations", rownames(trajMM)), ]

p1 <- plotTrajectoryHeatmap(trajMM, 
                            varCutOff = 0.3,
                            pal = paletteContinuous(set = "solarExtra"),
                            limits = c(-2, 2))

pdf("tracy_heatmap.pdf",width=10,height=10)
print(p1)
dev.off()

```


```{r, fig.height=10, fig.width=6}

trajGSM <- getTrajectory(ArchRProj = proj,
                         name = "Fib",
                         useMatrix = "GeneScoreMatrix",
                         log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2))

pdf("tracy_heatmap_trajectoryHeatmap.pdf",width=10,height=10)
print(p2)
dev.off()



```


##
```{r}


corGSM_MM <- correlateTrajectories(trajGSM, trajMM,
                                   varCutOff1 = 0.3,
                                   varCutOff2 = 0.3,
                                   corCutOff = 0)
#install.packages("magick")

df <- as.data.frame(corGSM_MM[[1]])
df <- df[!grepl("-AS", df$name1), ]
df <- df[df$FDR < 0.1, ]
df <- df[!duplicated(df$matchname1), ]
trajGSM2 <- trajGSM[df$name1, ]
trajMM2 <- trajMM[df$name2, ]

trajCombined <- trajGSM2

assay(trajCombined, withDimnames=FALSE) <- t(apply(assay(trajGSM2), 1, scale)) + 
    t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined,
                                     returnMat = TRUE,
                                     varCutOff = 0,
                                     force = TRUE)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,
                             pal = paletteContinuous(set = "horizonExtra"),
                             varCutOff = 0,
                             rowOrder = rowOrder,
                             limits = c(-2, 2),
                             labelRows = TRUE)
ht2 <- plotTrajectoryHeatmap(trajMM2,
                             pal = paletteContinuous(set = "solarExtra"),
                             varCutOff = 0,
                             rowOrder = rowOrder,
                             limits = c(-2, 2),
                             labelRows = TRUE)
ht1 + ht2

pdf("ht1.pdf",width=20,height=20)
print(ht1)
dev.off()

pdf("ht2.pdf",width=20,height=20)
print(ht2)
dev.off()


df <- df[rowOrder, ]
write.table(df, file = "../Fibroblast/sel_tf_motif_gex.txt",
             quote = FALSE, sep = "\t",
             row.names = FALSE)


```

}


##轨迹分析  相关性散点图

## plot correlation for each TF
```{r, fig.width=6, fig.height=4}

library(ggrepel)
df <- read.table("../Fibroblast/sel_tf_motif_gex.txt",
                 header = TRUE)
p <- ggplot(df, aes(x = reorder(matchname1, -Correlation), 
                    y = Correlation)) +
    geom_text_repel(aes(label = matchname2)) +
    geom_point() +
    xlab("TFs") + ylab("Correlation") +
    theme_cowplot() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank())

pdf("p_plot correlation for each TF.pdf",width=20,height=20)
print(p)
dev.off()





```



## save data
```{r}
saveArchRProject(ArchRProj = proj, 
                 load = FALSE)
```

#####################################14 fibroblast p2g peak to gene links##########3###
#
#
#
#
```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
source("./HiddenUtils.R")



```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 90)
addArchRGenome("mm10")
```

# load data
```{r load_data}


proj <- loadArchRProject(path = "../Fibroblast")


```

## save peaks
```{r}


peaks <- getPeakSet(proj)
rgs <- as.data.frame(peaks@ranges)
df_peaks <- data.frame(chr = peaks@seqnames,
                       start = rgs$start,
                       end = rgs$end)


write.table(df_peaks, file = "../Fibroblast/peaks.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


```


## define function
```{r define_func}
library(Rcpp)
Rcpp::sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  using namespace std;
  // Adapted from https://github.com/AEBilgrau/correlateR/blob/master/src/auxiliary_functions.cpp
  // [[Rcpp::export]]
  Rcpp::NumericVector rowCorCpp(IntegerVector idxX, IntegerVector idxY, Rcpp::NumericMatrix X, Rcpp::NumericMatrix Y) {
    
    if(X.ncol() != Y.ncol()){
      stop("Columns of Matrix X and Y must be equal length!");
    }
    if(max(idxX)-1 > X.nrow()){
      stop("Idx X greater than nrow of Matrix X");
    }
    if(max(idxY)-1 > Y.nrow()){
      stop("Idx Y greater than nrow of Matrix Y");
    }
    // Transpose Matrices
    X = transpose(X);
    Y = transpose(Y);
    
    const int nx = X.ncol();
    const int ny = Y.ncol();
    // Centering the matrices
    for (int j = 0; j < nx; ++j) {
      X(Rcpp::_, j) = X(Rcpp::_, j) - Rcpp::mean(X(Rcpp::_, j));
    }
    for (int j = 0; j < ny; ++j) {
      Y(Rcpp::_, j) = Y(Rcpp::_, j) - Rcpp::mean(Y(Rcpp::_, j));
    }
    // Compute 1 over the sample standard deviation
    Rcpp::NumericVector inv_sqrt_ss_X(nx);
    for (int i = 0; i < nx; ++i) {
      inv_sqrt_ss_X(i) = 1/sqrt(Rcpp::sum( X(Rcpp::_, i) * X(Rcpp::_, i) ));
    }
    Rcpp::NumericVector inv_sqrt_ss_Y(ny);
    for (int i = 0; i < ny; ++i) {
      inv_sqrt_ss_Y(i) = 1/sqrt(Rcpp::sum( Y(Rcpp::_, i) * Y(Rcpp::_, i) ));
    }
    //Calculate Correlations
    const int n = idxX.size();
    Rcpp::NumericVector cor(n);
    for(int k = 0; k < n; k++){
      cor[k] = Rcpp::sum( X(Rcpp::_, idxX[k] - 1) * Y(Rcpp::_, idxY[k] - 1) ) * inv_sqrt_ss_X(idxX[k] - 1) * inv_sqrt_ss_Y(idxY[k] - 1);    } 
    return(cor);
  }'
)
Rcpp::sourceCpp(code='
  #include <Rcpp.h>
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
Rcpp::IntegerVector determineOverlapCpp(IntegerMatrix m, int overlapCut){
  int k2 = 2 * m.ncol();
  int nr = m.nrow();
  int nUnion;
  int maxOverlap;
  IntegerVector unionVector;
  IntegerVector testVector = IntegerVector(nr);
  IntegerVector nOverlap = IntegerVector(nr);
  NumericVector maxOverlapVector = NumericVector(nr);
  IntegerVector vi;
  IntegerVector vj;
  for (int i = 1; i < nr; i++){
   
    if (i % 500 == 0) Rcpp::Rcout << "Completed Computing KNN Overlap " << i << " of " << nr << endl;
    
    for(int j = 0; j < i; j++){
      
      if(testVector(j) == 0){
        vi = m(i, _);
        vj = m(j, _);
        unionVector = union_( vi , vj );
        nUnion = unionVector.size();
        nOverlap(j) = k2 - nUnion;
      }else{
        nOverlap(j) = 0;
      }
    }
    maxOverlap = max( nOverlap );
    maxOverlapVector(i) = maxOverlap;
    if(maxOverlap > overlapCut){
      testVector(i) = -1;
    }
  }
  return testVector;
}'
)



```
## aggregate data along the trajectory 12d ######
```{r}


trajectory <- getCellColData(proj, "Fib")
trajectory

trajectory <- trajectory[!is.na(trajectory[, 1]), , drop = FALSE]

breaks <- seq(0, 100, 1)
breaks
seq_along(breaks)

groupList <- lapply(seq_along(breaks), function(x) {
        if (x == 1) {
            NULL
        }
        else {
            rownames(trajectory)[which(trajectory[, 1] > breaks[x - 
                1] & trajectory[, 1] <= breaks[x])]
        }
    })[-1]
groupList


names(groupList) <- paste0("T.", breaks[-length(breaks)], 
        "_", breaks[-1])
knnObj <- groupList

type(knnObj)
length(knnObj)

```




## create aggregated gene matrix
```{r agg_gex}
geneMatrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
geneMatrix
assay(geneMatrix) %>%dim()

geneSet <- geneMatrix@elementMetadata
geneSet

geneStart <- GRanges(geneSet$seqnames, 
                   IRanges(geneSet$start, 
        width = 1), name = geneSet$name, 
        idx = geneSet$idx)
geneStart

mat <- geneMatrix@assays@data$GeneScoreMatrix
dim(mat)

groupMatRNA <- lapply(seq_along(knnObj), function(x){
    cellNames <- knnObj[[x]]
    z <- Matrix::rowSums(mat[, cellNames])
    z
}) %>% Reduce(cbind, .)

rawMatRNA <- groupMatRNA
groupMatRNA <- t(t(groupMatRNA)/colSums(groupMatRNA)) * 10000


saveRDS(groupMatRNA, file = "../Fibroblast/aggr_rna.rds")




seRNA <- SummarizedExperiment(assays = SimpleList(RNA = groupMatRNA, 
                                                  RawRNA = rawMatRNA), 
                              rowRanges = geneStart)




```


## create aggregated atac-seq with raw data
```{r agg_raw}


PeakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
peakSet <- getPeakSet(ArchRProj = proj)
mat <- PeakMatrix@assays@data$PeakMatrix
groupMatATAC <- lapply(seq_along(knnObj), function(x){
    cellNames <- knnObj[[x]]
    z <- Matrix::rowSums(mat[, cellNames])
    z
}) %>% Reduce(cbind, .)
rawMatATAC <- groupMatATAC
groupMatATAC <- t(t(groupMatATAC)/colSums(groupMatATAC)) * 10000

dim(groupMatATAC)
dim(rawMatATAC)
dim(peakSet)



seATACRaw <- SummarizedExperiment(assays = SimpleList(ATAC = groupMatATAC, 
                                                      RawATAC = rawMatATAC),
                                  rowRanges = peakSet)

saveRDS(groupMatATAC, file = "../Fibroblast/aggr_raw.rds")

```

## find putative peak-to-gene p2g
```{r overlap}


maxDist = 250000

o <- DataFrame(findOverlaps(resize(seRNA, 2 * maxDist + 1, "center"),
                            resize(rowRanges(seATACRaw), 1, "center"),
                            ignore.strand = TRUE))

#install.packages("philentropy")
#library(philentropy)

#destiny::distance
o$distance <- GenomicFeatures::distance(rowRanges(seRNA)[o[, 1]], 
                       rowRanges(seATACRaw)[o[, 2]])


colnames(o) <- c("B", "A", "distance")
df <- rowRanges(seATACRaw)[o$A, ]
o$gene <- rowData(seRNA)[o$B, ]$name
o$peak <- paste0(df@seqnames, "_",
                 as.data.frame(df@ranges)$start, "_",
                 as.data.frame(df@ranges)$end)




```


## compute correlation using raw data
```{r cor_raw}

o$Correlation <- rowCorCpp(as.integer(o$A), 
                           as.integer(o$B), 
                           assay(seATACRaw), assay(seRNA))
o$VarAssayA <- .getQuantiles(matrixStats::rowVars(assay(seATACRaw)))[o$A]
o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
o$TStat <- (o$Correlation/sqrt((pmax(1 - o$Correlation^2, 
        1e-17, na.rm = TRUE))/(ncol(seATACRaw) - 2)))
o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATACRaw) - 2)
o$FDR <- p.adjust(o$Pval, method = "fdr")
out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", 
        "VarAssayB", "distance")]
colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", 
        "VarQATAC", "VarQRNA", "Distance")
  out$gene <- o$gene
  out$peak <- o$peak
  
out <- out[!is.na(out$FDR), ]
write.csv(out, file = "../Fibroblast/p2g_raw.csv")



```



##这个地方要小心 找到相应的文件 scopen.txt  应该读取peaks 读取fibro的降维矩阵
## create aggregated atac-seq with imputed data
```{r agg_imputed}

for(method in c("scOpen")) {#, "MAGIC", "cisTopic", "SCALE"
  
  method="scOpen"

'''#读取大文件需要改写
  mat <- read.table(glue::glue("../mytmp/{method}.txt"), 
                    header = TRUE, 
                    check.names = FALSE, 
                    sep = "\t", 
                    comment.char="")
'''

library(data.table)
mat <- fread(file = glue::glue("../mytmp/{method}.txt"), 
             header = TRUE, 
             sep = "\t",  
             verbose = FALSE)

####





 dim(mat)
 head(mat)[,1:7]
#[1]    152299     21  


seq_along(knnObj)

  groupMatATAC <- lapply(seq_along(knnObj), function(x){
      cellNames <- knnObj[[x]]
      z <- Matrix::rowSums(mat[, cellNames])
      z
  }) %>% Reduce(cbind, .)
  


  rawMatATAC <- groupMatATAC
  groupMatATAC <- t(t(groupMatATAC)/colSums(groupMatATAC)) * 10000
#####################################这里失败了 
#
dim(groupMatATAC)
dim(rawMatATAC)
dim(peakSet)

head(groupMatATAC)
 
 seATAC <- SummarizedExperiment(assays = SimpleList(ATAC = groupMatATAC, 
          RawATAC = rawMatATAC), rowRanges = peakSet)
    

  saveRDS(groupMatATAC, 
          file = glue::glue("../Fibroblast/aggr_{method}.rds"))




  
  o$Correlation <- rowCorCpp(as.integer(o$A), as.integer(o$B), 
        assay(seATAC), assay(seRNA))
  o$VarAssayA <- .getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
  
  o$VarAssayB <- .getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
  
  o$TStat <- (o$Correlation/sqrt((pmax(1 - o$Correlation^2, 
          1e-17, na.rm = TRUE))/(ncol(seATAC) - 2)))
  o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) - 2)
  o$FDR <- p.adjust(o$Pval, method = "fdr")
  out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", 
          "VarAssayB", "distance")]
  colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", 
          "VarQATAC", "VarQRNA", "Distance")
  out$gene <- o$gene
  out$peak <- o$peak
  
  out <- out[!is.na(out$FDR), ]
 
  write.csv(out, file = glue::glue("../Fibroblast/p2g_{method}.csv"))



}





```


### ################################333###15 fibroblat_Twist1#######
#主要是评价预测结果 可以跳过


```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
source("HiddenUtils.R")



```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 90)
addArchRGenome("mm10")
```

```{r load_data}
proj <- loadArchRProject(path = "../Fibroblast")


```


## loadd p2g links
```{r load_p2g}


df_p2g_raw <- read.csv("../Fibroblast/p2g_raw.csv", row.names = 1)

df_p2g_scopen <- read.csv("../Fibroblast/p2g_scOpen.csv", row.names = 1)

'''
df_p2g_magic <- read.csv("../Fibroblast/p2g_MAGIC.csv", row.names = 1)
df_p2g_cisTopic <- read.csv("../Fibroblast/p2g_cisTopic.csv", row.names = 1)
df_p2g_SCALE <- read.csv("../Fibroblast/p2g_SCALE.csv", row.names = 1)
'''


df_p2g_raw <- subset(df_p2g_raw, Distance > 0)
df_p2g_scopen <- subset(df_p2g_scopen, Distance > 0)


df_p2g_magic <- subset(df_p2g_magic, Distance > 0)
df_p2g_cisTopic <- subset(df_p2g_cisTopic, Distance > 0)
df_p2g_SCALE <- subset(df_p2g_SCALE, Distance > 0)



df_p2g_raw$data <- "Raw"
df_p2g_scopen$data <- "scOpen"


df_p2g_magic$data <- "MAGIC"
df_p2g_cisTopic$data <- "cisTopic"
df_p2g_SCALE$data <- "SCALE"



df <- rbind(df_p2g_raw, df_p2g_scopen)#df_p2g_magic,df_p2g_cisTopic, df_p2g_SCALE
                                      


p1 <- ggplot(data = df, aes(x = reorder(data, -abs(Correlation), FUN = mean), 
                            y = Correlation, fill = data)) +
  geom_violin() +
  theme_cowplot()

pdf("p2g_raw scopen.pdf",width=20,height=20)
print(p1)
dev.off()




df.plot <- df %>%
  group_by(data) %>%
  summarise(AvgAbsCor = mean(abs(Correlation)))

p2 <- ggplot(data = df.plot, aes(x = reorder(data, -AvgAbsCor), 
                                y = AvgAbsCor, fill = data)) +
  geom_bar(stat = "identity") +
  theme_cowplot()


pdf("p2g_raw scopen_avg.pdf",width=10,height=10)
print(p2)
dev.off()

AvgAbsCor
1 0.1022388

```



## get Twist1 binding sites
```{r}


motifInPeaks <- readRDS("../Fibroblast/Annotations/Motif-Matches-In-Peaks.rds")
###################
#getPositions(ArchRProj = proj, annoName = Twist1_785)

Twist1InPeaks <- motifInPeaks[, grep("Twist1_785", colnames(motifInPeaks))]
#Twist1_785
motifInPeaks
Twist1InPeaks


peaksWithTwist1 <- which(as.data.frame(Twist1InPeaks@assays@data$matches)$Twist1_785)

motifPositions <- getPositions(proj)#motifInPeaks



Twist1BindingSites <- read.table("../HINT/FootprintingFibroblast/Twist1.bed")
motifPositions <- getPositions(proj)

Twist1BindingSites=motifPositions[["Twist1_785"]] #改动代码
dir.create("../HINT/FootprintingFibroblast/",recursive=T)
#save(Twist1BindingSites,file="../HINT/FootprintingFibroblast/Twist1.bed")




gr <- GRanges(seqnames = seqnames(Twist1BindingSites),  ###=gr=改写
              IRanges(start = as.numeric(Twist1BindingSites$V2),
                      end = as.numeric(Twist1BindingSites$V3)))

gr=Twist1BindingSites
```


## load DE genes
```{r}


library(DESeq2)
load("../../../Twist1_RNA/results/star_salmon/deseq2_qc/deseq2.dds.RData")
de_res <- results(dds) %>%
  as.data.frame() %>%
    subset(., !is.na(padj))
#de_res$DE <- ifelse(de_res$padj < 1e-05 & abs(de_res$log2FoldChange) > 1.5, 1, 0)
de_res$DE <- ifelse(de_res$padj < 1e-05 & abs(de_res$log2FoldChange) > 1, 1, 0)
de_res$gene <- rownames(de_res)
```


## define function
```{r}
#install.packages("PRROC")

library(PRROC)
get_evaluate <- function(p2g_list, final_prediction){
    pr_list <- lapply(seq_along(p2g_list), function(x){
    df_p2g <- p2g_list[[x]]
    df_p2g$gene <- toupper(df_p2g$gene)
    
    method <- names(p2g_list)[[x]]
    
    df_p2g$peak_chr <- stringr::str_split_fixed(df_p2g$peak, "_", 3)[, 1]
    df_p2g$peak_start <- stringr::str_split_fixed(df_p2g$peak, "_", 3)[, 2]
    df_p2g$peak_end <- stringr::str_split_fixed(df_p2g$peak, "_", 3)[, 3]
    peaks <- GRanges(seqnames = df_p2g$peak_chr,
                         IRanges(start = as.numeric(df_p2g$peak_start),
                                 end = as.numeric(df_p2g$peak_end)))
    o <- findOverlaps(query = peaks, subject = gr,
                          ignore.strand = TRUE)
    Twist1_p2g <- df_p2g[as.integer(o@from), ]
    
    Twist1_p2g <- Twist1_p2g %>%
        group_by(gene) %>%
        summarise(Correlation = abs(max(Correlation))) %>%
        as.data.frame()
    
    rownames(Twist1_p2g) <- Twist1_p2g$gene
    commonGene <- intersect(Twist1_p2g$gene, de_res$gene)
    
    de_res2 <- de_res[commonGene, ]
    Twist1_p2g <- Twist1_p2g[commonGene, ]
    df <- merge.data.frame(de_res2, Twist1_p2g)
    
    fg <- df[df$DE == 1, ]
    bg <- df[df$DE == 0, ]
    pr <- pr.curve(scores.class0 = fg$Correlation, 
                   scores.class1 = bg$Correlation, curve = T)
    
    pr$backgroud <- nrow(fg) / (nrow(fg) + nrow(bg))
    
    pr
})
    
    df <- lapply(seq_along(pr_list), function(x){
        pr <- pr_list[[x]]
        method <- names(p2g_list)[[x]]
      
          curve <- as.data.frame(pr$curve[, 1:2])
          colnames(curve) <- c("Recall", "Precision")
          curve$Method <- method
        
        curve
    })%>% Reduce(rbind, .)
    df.aupr <- lapply(seq_along(pr_list), function(x){
      pr <- pr_list[[x]]
        method <- names(p2g_list)[[x]]
        
        aupr <- pr$auc.integral
        
        df <- data.frame(aupr = aupr, method = method)
    })%>% Reduce(rbind, .)
    background_pre <- pr_list[[1]]$backgroud
    
    return(list(df, df.aupr, background_pre))
}


```


## evaluation of peak-to-gene links using DE genes
```{r, fig.height=5, fig.width=10}


res_list <- get_evaluate(p2g_list, final_prediction = "abs_max")
df <- res_list[[1]]
df.aupr2 <- res_list[[2]]
background_pre <- res_list[[3]]
p1 <- ggplot(data = df, aes(x = Recall, y = Precision)) +
  geom_line(aes(color = Method)) +
  geom_hline(yintercept = background_pre) +
  theme_cowplot()+
    ggtitle("abs(max)")
p2 <- ggplot(data = df.aupr2, aes(x =reorder(method, -aupr), 
                                 y = aupr, fill = method)) +
  geom_bar(stat = "identity")+
    theme_cowplot() +
        xlab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ggtitle("abs(max)")
p1 + p2
saveRDS(res_list, file = "../data/Fibroblast/abs_max.Rds")


```

## Session information
```{r}
sessionInfo()
```

################################################16_fibroblast_Twist1

```{r setup, include=FALSE}


library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
source("~/scope_abc/third/HiddenUtils.R")


```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 80)
addArchRGenome("mm10")
```


```{r load_data}


proj <- loadArchRProject(path = "../Fibroblast")
geneMatrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
peakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")


```

## find marker peaks
```{r}
# markersPeaks <- getMarkerFeatures(
#     ArchRProj = proj, 
#     useMatrix = "PeakMatrix", 
#     groupBy = "CellType",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# heatmapPeaks <- plotMarkerHeatmap(
#   seMarker = markersPeaks, 
#   cutOff = "FDR <= 0.05 & Log2FC >= 0.5")
# 
# draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
# 
# markerList <- getMarkers(markersPeaks, 
#                          cutOff = "FDR <= 0.05 & Log2FC >= 0.5", returnGR = TRUE)
# markerList$Myofibroblast
# 
# peaks_Myofibroblast <- as.data.frame(markerList$Myofibroblast)
# 
# peaks_Myofibroblast$peaks <- paste0(peaks_Myofibroblast$seqnames, "_",
#                                     peaks_Myofibroblast$start, "_",
#                                     peaks_Myofibroblast$end)
```

## loadd p2g links
```{r load_p2g}


df_p2g_raw <- read.csv("../Fibroblast/p2g_raw.csv", row.names = 1)
df_p2g_scopen <- read.csv("../Fibroblast/p2g_scOpen.csv", row.names = 1)

'''
df_p2g_magic <- read.csv("../data/Fibroblast/p2g_MAGIC.csv", row.names = 1)
df_p2g_cisTopic <- read.csv("../data/Fibroblast/p2g_cisTopic.csv", row.names = 1)
df_p2g_SCALE <- read.csv("../data/Fibroblast/p2g_SCALE.csv", row.names = 1)
'''

df_p2g_raw <- subset(df_p2g_raw, Distance > 0)
df_p2g_scopen <- subset(df_p2g_scopen, Distance > 0)


df_p2g_magic <- subset(df_p2g_magic, Distance > 0)
df_p2g_cisTopic <- subset(df_p2g_cisTopic, Distance > 0)
df_p2g_SCALE <- subset(df_p2g_SCALE, Distance > 0)


df_p2g_raw$data <- "Raw"
df_p2g_scopen$data <- "scOpen"
df_p2g_raw %>%head()
df_p2g_scopen %>%head()


'''
df_p2g_magic$data <- "MAGIC"
df_p2g_cisTopic$data <- "cisTopic"
df_p2g_SCALE$data <- "SCALE"
'''

df <- rbind(df_p2g_raw, df_p2g_scopen)# df_p2g_magic, df_p2g_cisTopic, df_p2g_SCALE
           


p1 <- ggplot(data = df, aes(x = data, y = Correlation)) +
  geom_boxplot()


pdf("Correlation1_raw_scopen.pdf",width=10,height=10)
print(p1)
dev.off()



p2 <- ggplot(data = df, aes(x = data, y = Correlation)) +
  geom_violin()

pdf("Correlation2Correlation1_raw_scopen.pdf",width=10,height=10)
print(p2)
dev.off()






```


## get Twist1 binding sites
```{r}

load("../HINT/FootprintingFibroblast/Twist1.bed") #Twist1BindingSites <- 
gr=Twist1BindingSites


gr <- GRanges(seqnames = Twist1BindingSites$V1,
              IRanges(start = as.numeric(Twist1BindingSites$V2),
                      end = as.numeric(Twist1BindingSites$V3)))





```


## define function
```{r}

getCellColData(proj, "Fib")


get_trajectory <- function(proj, mat){
    trajectory <- getCellColData(proj, "Fib")
    trajectory <- trajectory[!is.na(trajectory[, 1]), , drop = FALSE]
    breaks <- seq(0, 100, 1)

groupList <- lapply(seq_along(breaks), function(x) {
        if (x == 1) {
            NULL
        }
        else {
            rownames(trajectory)[which(trajectory[, 1] > breaks[x - 
                1] & trajectory[, 1] <= breaks[x])]
        }
    }  ) [-1]


names(groupList) <- paste0("T.", breaks[-length(breaks)], 
        "_", breaks[-1])
    
    groupMat <- lapply(seq_along(groupList), function(x){
    	#x=2
        cellNames <- groupList[[x]]  #############代码我改动了一下
   
   #############改动位置
   #https://blog.csdn.net/yijiaobani/article/details/83784733
   head(mat)[,1:6]
   mat[, c("D0_2#ACTGTCCTCGGCAATT-1","D2_1#GCAGCCAAGCCGCTGT-1")]
   Matrix::rowSums( mat[, c("D0_2#ACTGTCCTCGGCAATT-1","D2_1#GCAGCCAAGCCGCTGT-1")] )
   cellNames=cellNames[cellNames %in% colnames(mat)]
   table(cellNames %in% colnames(mat))
   print(paste0("这是第=====",x,"--=====个循环"))
  mat[,c("D10_2#ACTTCCGGTAGAGAGA-1","D0_2#TGCTATTTCCAACGCG-1")]
   ###############################

        z <- Matrix::rowSums( mat[, cellNames,with=F]  )
        z
    }) %>% Reduce(cbind, .)
    


    groupMat <- t(t(groupMat)/colSums(groupMat)) * 10000
    
    
    message("Smoothing...")
    smoothGroupMat <- as.matrix(t(apply(groupMat, 1, function(x) centerRollMean(x, k = 3))))
    colnames(smoothGroupMat) <- paste0(colnames(groupMat))
        colnames(groupMat) <- paste0(colnames(groupMat))
        seTrajectory <- SummarizedExperiment(assays = SimpleList(smoothMat = as.matrix(smoothGroupMat), 
            mat = as.matrix(groupMat)), rowData = NULL)
    rownames(seTrajectory) <- rownames(mat)
    metadata(seTrajectory)$Params <- list(useMatrix = "peaks", 
        matrixClass = "Sparse.Integer.Matrix", scaleTo = 10000, log2Norm = TRUE, 
        smoothWindow = 11, date = Sys.Date())
        
    seTrajectory

}

#https://rdrr.io/github/CostaLab/scMEGA/src/R/utils.R
```

centerRollMean <- function(v = NULL, k = NULL){
  o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
  if(k%%2==0){
    o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else if(k%%2==1){
    o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
  }else{
    stop("Error!")
  }
  o2
}





## heatmap of all peaks and genes
```{r}


peaks <- as.data.frame(peakMatrix@rowRanges)
peaks

mat <- peakMatrix@assays@data$PeakMatrix
rownames(mat) <- paste0(peaks$seqnames, "_", 
                                 peaks$start, "_", 
                                 peaks$end)
peaks %>%head()
mat %>%head()


trajPM  <- get_trajectory(proj, mat)




p <- plotTrajectoryHeatmap(trajPM, 
                             pal = paletteContinuous(set = "solarExtra"))
 

pdf("solarExtra_trajectory_second.pdf",width=10,height=10)
print(p)
dev.off()



trajGSM  <- getTrajectory(proj, name = "Fib",
                          useMatrix = "GeneScoreMatrix",
                          smoothWindow = 11)



rownames(trajGSM) <- stringr::str_split_fixed(rownames(trajGSM), ":", 2)[, 2]
rownames(trajGSM) <- toupper(rownames(trajGSM))



p <- plotTrajectoryHeatmap(trajGSM, 
                             pal = paletteContinuous(set = "solarExtra"))
 
pdf("solarExtra_trajectorytrajGSM.pdf",width=10,height=10)
print(p)
dev.off()




print(nrow(trajPM))
print(nrow(trajGSM))


geneSet <- geneMatrix@elementMetadata
groupMatRNA <- readRDS(file = "../Fibroblast/aggr_rna.rds")



rownames(groupMatRNA) <- toupper(geneSet$name)

groupMatRNA %>%head()


```
## add peak-to-gene links using raw count matrix
```{r raw, fig.height=10, fig.width=14}

df_p2g_raw <- subset(df_p2g_raw, FDR < 0.05)
df_p2g_raw$peak_chr <- stringr::str_split_fixed(df_p2g_raw$peak, "_", 3)[, 1]
df_p2g_raw$peak_start <- stringr::str_split_fixed(df_p2g_raw$peak, "_", 3)[, 2]
df_p2g_raw$peak_end <- stringr::str_split_fixed(df_p2g_raw$peak, "_", 3)[, 3]


peaks <- GRanges(seqnames = df_p2g_raw$peak_chr,
                         IRanges(start = as.numeric(df_p2g_raw$peak_start),
                                 end = as.numeric(df_p2g_raw$peak_end)))

peaks



o <- findOverlaps(query = peaks, subject = gr,
                          ignore.strand = TRUE)






Twist1_2_g_raw <- df_p2g_raw[as.integer(o@from), ]
peaks <- as.data.frame(peakMatrix@rowRanges)
peaks %>%head()



groupMatATAC <- readRDS(file = "../Fibroblast/aggr_raw.rds")
groupMatATAC %>%head()
dim(groupMatATAC)
dim(peaks)

rownames(groupMatATAC) <- paste0(peaks$seqnames, "_", 
                                 peaks$start, "_", 
                                 peaks$end)
Twist1_2_g_raw$gene <- toupper(Twist1_2_g_raw$gene)
dim(Twist1_2_g_raw)
Twist1_2_g_raw %>%head() #下游靶点



dir.out <- "../Fibroblast/Raw"

if(!dir.exists(dir.out)){
    dir.create(dir.out)
}







for (i in 1:nrow(Twist1_2_g_raw)) {
	#i=5
    gene <- Twist1_2_g_raw$gene[[i]]
    peak <- Twist1_2_g_raw$peak[[i]]
    

    gene
    peak
    correlation <- Twist1_2_g_raw$Correlation[[i]]
    
    rna <- groupMatRNA[gene, ]
    atac <- groupMatATAC[peak, ]
    rna
    atac
groupMatATAC %>%head()
groupMatRNA%>%head()

    df <- data.frame(rna = rna,
                     atac = atac)
    df
    df$PseudoTime <- 1:nrow(df)
    
    pal <- rev(ArchR::paletteContinuous(set = "blueYellow"))
    
    p <- ggplot(data = df, aes(x = rna, y = atac, 
                               color = PseudoTime)) +
        geom_point() +
        scale_color_gradientn(colours = pal) +
        xlab("RNA") + ylab("ATAC") +
        ggtitle(glue::glue("Cor: {correlation}")) +
        theme_cowplot()
    
    peak <- stringr::str_replace_all(peak, ":", "_")
   

    peak
    pdf(file = glue::glue("{dir.out}/{gene}_{peak}.pdf"),
        width = 5, height = 4)
    print(p)
    dev.off()


}








trajGSM_sub <- trajGSM[Twist1_2_g_raw$gene, ]
trajPM_sub <- trajPM[Twist1_2_g_raw$peak, ]
print(nrow(trajPM_sub))
print(nrow(trajGSM_sub))


ht1 <- plotTrajectoryHeatmap(trajPM_sub, 
                           varCutOff = 0,
                            pal = paletteContinuous(set = "solarExtra"))
ht2 <- plotTrajectoryHeatmap(trajGSM_sub, 
                           varCutOff = 0,
                            pal = paletteContinuous(set = "horizonExtra"))
ht1 + ht2


pdf("ht1_solarExtra.pdf",width=10,height=10)
print(h1)
dev.off()

pdf("ht2_horizonExtra.pdf",width=10,height=10)
print(h2)
dev.off()

plotPDF(ht1,ht2, name = "ht1_solarExtra_ht2_horizonExtra.pdf", 
	ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)



  
write.csv(Twist1_2_g_raw, 
          file = "../Fibroblast/Twist1_target_genes_raw.csv")


```

## add peak-to-gene links using scopen matrix
```{r scopen, fig.height=10, fig.width=14}

df_p2g_scopen <- subset(df_p2g_scopen, FDR < 0.05)
df_p2g_scopen$peak_chr <- stringr::str_split_fixed(df_p2g_scopen$peak, "_", 3)[, 1]
df_p2g_scopen$peak_start <- stringr::str_split_fixed(df_p2g_scopen$peak, "_", 3)[, 2]
df_p2g_scopen$peak_end <- stringr::str_split_fixed(df_p2g_scopen$peak, "_", 3)[, 3]
peaks <- GRanges(seqnames = df_p2g_scopen$peak_chr,
                         IRanges(start = as.numeric(df_p2g_scopen$peak_start),
                                 end = as.numeric(df_p2g_scopen$peak_end)))


o <- findOverlaps(query = peaks, subject = gr,
                          ignore.strand = TRUE)

Twist1_2_g_scopen <- df_p2g_scopen[as.integer(o@from), ]
peaks <- as.data.frame(peakMatrix@rowRanges)


groupMatATAC <- readRDS(file = "../Fibroblast/aggr_scOpen.rds")
rownames(groupMatATAC) <- paste0(peaks$seqnames, "_", 
                                 peaks$start, "_", 
                                 peaks$end)
Twist1_2_g_scopen$gene <- toupper(Twist1_2_g_scopen$gene)
dir.out <- "../Fibroblast/scOpen"

if(!dir.exists(dir.out)){
    dir.create(dir.out)
}



for (i in 1:nrow(Twist1_2_g_scopen)) {
    gene <- Twist1_2_g_scopen$gene[[i]]
    peak <- Twist1_2_g_scopen$peak[[i]]
    
    correlation <- Twist1_2_g_scopen$Correlation[[i]]
    
    rna <- groupMatRNA[gene, ]
    atac <- groupMatATAC[peak, ]
    
    df <- data.frame(rna = rna,
                     atac = atac)
    df$PseudoTime <- 1:nrow(df)
    
    pal <- ArchR::paletteContinuous(set = "blueYellow")
    
    p <- ggplot(data = df, aes(x = rna, y = atac, 
                               color = PseudoTime)) +
        geom_point() +
        scale_color_gradientn(colours = pal) +
        xlab("RNA") + ylab("ATAC") +
        ggtitle(glue::glue("Cor: {correlation}")) +
        theme_cowplot()
    
    peak <- stringr::str_replace_all(peak, ":", "_")
    
    pdf(file = glue::glue("{dir.out}/{gene}_{peak}.pdf"),
        width = 5, height = 4)
    print(p)
    dev.off()
}

##############peudo-bulks-seq????????
mat <- read.table(glue::glue("../mytmp/scOpen.txt"),
                   header = TRUE, check.names = FALSE, sep = "\t",
                 comment.char="")



library(data.table)
mat <- fread(file = glue::glue("../mytmp/scOpen.txt"), 
             header = TRUE, 
             sep = "\t",  
             verbose = FALSE)



rownames(mat) <- paste0(peaks$seqnames, "_", 
                                 peaks$start, "_", 
                                 peaks$end)
dim(mat)
head(mat)


trajPM  <- get_trajectory(proj, mat)



trajPM_sub <- trajPM[Twist1_2_g_scopen$peak, ]
trajGSM_sub <- trajGSM[Twist1_2_g_scopen$gene, ]
print(nrow(trajPM_sub))
print(nrow(trajGSM_sub))


ht1 <- plotTrajectoryHeatmap(trajPM_sub,
                           varCutOff = 0,
                           limits = c(-2, 2),
                            pal = paletteContinuous(set = "solarExtra"))
ht2 <- plotTrajectoryHeatmap(trajGSM_sub,
                           varCutOff = 0,
                            pal = paletteContinuous(set = "horizonExtra"))
ht1 + ht2

plotPDF(ht1,ht2, name = "ht1_solarExtra_ht2_horizonExtra_scopen.pdf", 
  ArchRProj = proj, addDOC = FALSE, width = 8, height = 8)


write.csv(Twist1_2_g_scopen, file = "../Fibroblast/Twist1_target_genes_scopen.csv")



```

##不用做            fc
#########################17_fibroblast_Twist1###################################
################17_fibroblast_Twist1################################
#主要就是加载

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
source("HiddenUtils.R")
```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 80)
addArchRGenome("mm10")
```


```{r load_data}


proj <- loadArchRProject(path = "../data/Fibroblast")
geneMatrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
peakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")
peaks <- as.data.frame(peakMatrix@rowRanges)
peaks$peaks <- paste0(peaks$seqnames, ":", peaks$start, "_", peaks$end)




```


## loadd p2g links
```{r load_p2g}


df_p2g_raw <- read.csv("../Fibroblast/p2g_raw.csv", row.names = 1)
df_p2g_scopen <- read.csv("../Fibroblast/p2g_scOpen.csv", row.names = 1)



df_p2g_magic <- read.csv("../data/Fibroblast/p2g_MAGIC.csv", row.names = 1)
df_p2g_cisTopic <- read.csv("../data/Fibroblast/p2g_cisTopic.csv", row.names = 1)
df_p2g_SCALE <- read.csv("../data/Fibroblast/p2g_SCALE.csv", row.names = 1)



df_p2g_raw <- subset(df_p2g_raw, Distance > 0)
df_p2g_scopen <- subset(df_p2g_scopen, Distance > 0)




df_p2g_magic <- subset(df_p2g_magic, Distance > 0)
df_p2g_cisTopic <- subset(df_p2g_cisTopic, Distance > 0)
df_p2g_SCALE <- subset(df_p2g_SCALE, Distance > 0)
```


## get Twist1 binding sites
```{r}
load("../HINT/FootprintingFibroblast/Twist1.bed") #Twist1BindingSites <- 
gr=Twist1BindingSites

gr <- GRanges(seqnames = Twist1BindingSites$V1, 
              IRanges(start = as.numeric(Twist1BindingSites$V2),
                      end = as.numeric(Twist1BindingSites$V3)))



df_p2g_raw$peak_chr <- stringr::str_split_fixed(df_p2g_raw$peak, "_", 3)[, 1]
df_p2g_raw$peak_start <- stringr::str_split_fixed(df_p2g_raw$peak, "_", 3)[, 2]
df_p2g_raw$peak_end <- stringr::str_split_fixed(df_p2g_raw$peak, "_", 3)[, 3]
peaks <- GRanges(seqnames = df_p2g_raw$peak_chr,
                         IRanges(start = as.numeric(df_p2g_raw$peak_start),
                                 end = as.numeric(df_p2g_raw$peak_end)))
o <- findOverlaps(query = peaks, subject = gr,
                          ignore.strand = TRUE)



```


```{r}
Twist1_2_g_raw <- df_p2g_raw[as.integer(o@from), ]
Twist1_2_g_scopen <- df_p2g_scopen[as.integer(o@from), ]


Twist1_2_g_magic <- df_p2g_magic[as.integer(o@from), ]
Twist1_2_g_cisTopic <- df_p2g_cisTopic[as.integer(o@from), ]
Twist1_2_g_SCALE <- df_p2g_SCALE[as.integer(o@from), ]


Twist1_2_g_raw$gene <- toupper(Twist1_2_g_raw$gene)
Twist1_2_g_scopen$gene <- toupper(Twist1_2_g_scopen$gene)


Twist1_2_g_magic$gene <- toupper(Twist1_2_g_magic$gene)
Twist1_2_g_cisTopic$gene <- toupper(Twist1_2_g_cisTopic$gene)
Twist1_2_g_SCALE$gene <- toupper(Twist1_2_g_SCALE$gene)
```


## correlation to FC
```{r, fig.height=12, fig.width=12}
library(DESeq2)
library(ggrepel)
load("../../../Twist1_RNA/results/star_salmon/deseq2_qc/deseq2.dds.RData")
de_res <- results(dds) %>%
  as.data.frame() %>%
  subset(., !is.na(log2FoldChange) &!is.na(padj))
de_res$gene <- rownames(de_res)
df_raw <- merge.data.frame(de_res, Twist1_2_g_raw)
df_scopen <- merge.data.frame(de_res, Twist1_2_g_scopen)
df_magic <- merge.data.frame(de_res, Twist1_2_g_magic)
df_cisTopic <- merge.data.frame(de_res, Twist1_2_g_cisTopic)
df_SCALE <- merge.data.frame(de_res, Twist1_2_g_SCALE)
cor_raw <- cor(df_raw$log2FoldChange, df_raw$Correlation)
cor_scopen <- cor(df_scopen$log2FoldChange, df_scopen$Correlation)
cor_magic <- cor(df_magic$log2FoldChange, df_magic$Correlation)
cor_cisTopic <- cor(df_cisTopic$log2FoldChange, df_cisTopic$Correlation)
cor_SCALE <- cor(df_SCALE$log2FoldChange, df_SCALE$Correlation)
horizontal.lins <- c(-0.4, 0.4)
vertical.lins <- c(-2.5, 2.5)
df_text <- subset(df_raw, abs(log2FoldChange) > 2.5 & 
                      abs(Correlation) > 0.4)
p1 <- ggplot(data = df_raw, aes(x = log2FoldChange, y = Correlation)) +
  geom_point() + ggtitle("Raw") +
    geom_text_repel(data = df_text, aes(x = log2FoldChange, y = Correlation,
                                        label = gene)) +
    theme_cowplot() +
    geom_hline(yintercept = horizontal.lins) +
    geom_vline(xintercept = vertical.lins)
df_text <- subset(df_scopen, abs(log2FoldChange) > 2.5 & 
                      abs(Correlation) > 0.4)
p2 <- ggplot(data = df_scopen, aes(x = log2FoldChange, y = Correlation)) +
  geom_point() + ggtitle("scOpen")+
        geom_text_repel(data = df_text, aes(x = log2FoldChange, y = Correlation,
                                        label = gene)) +
    theme_cowplot() +
    geom_hline(yintercept = horizontal.lins) +
    geom_vline(xintercept = vertical.lins)
df_text <- subset(df_magic, abs(log2FoldChange) > 2.5 & 
                      abs(Correlation) > 0.4)
p3 <- ggplot(data = df_magic, aes(x = log2FoldChange, y = Correlation)) +
  geom_point() + ggtitle("MAGIC") +
        geom_text_repel(data = df_text, aes(x = log2FoldChange, y = Correlation,
                                        label = gene)) +
    theme_cowplot() +
    geom_hline(yintercept = horizontal.lins) +
    geom_vline(xintercept = vertical.lins)
df_text <- subset(df_cisTopic, abs(log2FoldChange) > 2.5 & 
                      abs(Correlation) > 0.4)
p4 <- ggplot(data = df_cisTopic, aes(x = log2FoldChange,
                                     y = Correlation)) +
    geom_point() + ggtitle("cisTopic") +
        geom_text_repel(data = df_text, aes(x = log2FoldChange, y = Correlation,
                                        label = gene)) +
    theme_cowplot() +
    geom_hline(yintercept = horizontal.lins) +
    geom_vline(xintercept = vertical.lins)
df_text <- subset(df_SCALE, abs(log2FoldChange) > 2.5 & 
                      abs(Correlation) > 0.4)
p5 <- ggplot(data = df_SCALE, aes(x = log2FoldChange,
                                  y = Correlation)) +
  geom_point() + ggtitle("SCALE") +
        geom_text_repel(data = df_text, aes(x = log2FoldChange, y = Correlation,
                                        label = gene)) +
    theme_cowplot() +
    geom_hline(yintercept = horizontal.lins) +
    geom_vline(xintercept = vertical.lins)
p1
p2
p3
p4
p5
```


####################################18_fibroblast_Twist1#################3
#

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
```


## loadd p2g links
```{r load_p2g}

df_p2g <- read.csv("../Fibroblast/p2g_scOpen.csv", row.names = 1) %>%
    subset(., Distance > 0 & FDR < 0.05) %>%
    subset(., select = c("gene", "peak", "Correlation", "Distance"))
df_p2g$gene <- toupper(df_p2g$gene)





```

## load DE results
```{r}
library(DESeq2)
load("../../../Twist1_RNA/results/star_salmon/deseq2_qc/deseq2.dds.RData")

res <- results(dds) %>%
  as.data.frame() %>%
  subset(., !is.na(log2FoldChange) &!is.na(padj)) %>%
    subset(., select = c("log2FoldChange", "padj")) %>%
    subset(., padj < 0.05)
res$gene <- rownames(res)
res$regulation_Twist1overexpression <- ifelse(res$log2FoldChange > 0, "up", "down")

# gene  log2FoldChange  padj  regulation_Twist1overexpression
```
merge.data.frame("1","3")
merge.data.frame("1",c("3","1"))

res=data.frame()
df_p2g %>%head()

#########################zichuang #######

res=df_p2g
colnames(res)=c("gene",  "log2FoldChange" , "padj", "regulation_Twist1overexpression")
res %>%head()
## merge data frames
```{r}

df <- merge.data.frame(res, df_p2g)
df %>%head()

```


## get Twist1 binding sites
```{r}

Twist1BindingSites <- read.table("../../HINT/FootprintingFibroblast/Twist1.bed")
gr <- GRanges(seqnames = Twist1BindingSites$V1, 
              IRanges(start = as.numeric(Twist1BindingSites$V2),
                      end = as.numeric(Twist1BindingSites$V3)))


peak_chr <- stringr::str_split_fixed(df$peak, "_", 3)[, 1]
peak_start <- stringr::str_split_fixed(df$peak, "_", 3)[, 2]
peak_end <- stringr::str_split_fixed(df$peak, "_", 3)[, 3]
peaks <- GRanges(seqnames = peak_chr,
                         IRanges(start = as.numeric(peak_start),
                                 end = as.numeric(peak_end)))
o <- findOverlaps(query = peaks, subject = gr,
                          ignore.strand = TRUE)
df$Twist1Binding <- 0
df$Twist1Binding[o@from] <- 1


```


## load ArchR project
```{r}

proj <- loadArchRProject(path = "../data/Fibroblast")

markers <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "CellType",#########
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Myofibroblast",
  bgdGroups = "Fib_Scara5"
)



pv <- plotMarkers(seMarker = markers, 
                  name = "Myofibroblast", 
                  cutOff = "FDR <= 0.1", 
                  plotAs = "Volcano")
pv #差异marker基因


df_de <- assays(markers) %>% Reduce(cbind, .)
colnames(df_de) <- names(assays(markers))


df_de$gene <- toupper(rowData(markers)$name)
df_de <- subset(df_de, select = c("gene", "Log2FC", "FDR"))
colnames(df_de) <- c("gene", "Log2FC_UUO", "FDR_UUO")
df_de$regulation_UUO <- ifelse(df_de$Log2FC > 0, "up", "down")


```

## merge
```{r}
head(df)
head(df_p2g)
head(df_de)


df <- merge.data.frame(df, df_de)
head(df)

colnames(df) <- c("gene",
                  "log2FC_Twist1overexpression",
                  "padj_Twist1overexpression",
                  "regulation_Twist1overexpression",
                  
                  "peak",
                  "Correlation",
                  "Distance",
                  
                  "Twist1Binding",
                  
                  "log2FC_UUO",
                  "padj_UUO",
                  "regulation_UUO")

write.csv(df, file = "../Fibroblast/scOpen_p2g_DE_Twist1.csv",
          row.names = FALSE, quote = FALSE)


```

#############################19_gene_along_trajectory######################################
##不用运行

```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)


```

## load data
```{r}
proj <- loadArchRProject("../data/Fibroblast")
proj <- addImputeWeights(proj, reducedDims = "harmony")


```

## plot
```{r, fig.width=4, fig.height=4}

##########想要展示的基因 选择的基因
for (gene in c("Scara5", "Fbln2", "Dcn", "Fn1", "Col1a1", "Twist2",
               "Tgfbr1", "Twist1")) {
    
    
        p <- plotEmbedding(
            ArchRProj = proj,
            colorBy = "GeneScoreMatrix",
            name = gene,
            embedding = "dm",
            quantCut = c(0.01, 0.99),
            plotAs = "points",
            labelAsFactors = FALSE,
            labelSize = 0,
            size = 1.5
        ) + xlab("") + ylab("") +
            theme_cowplot() +
            theme(legend.position = "none",
		 axis.line = element_blank(),
                  axis.text = element_blank(),
                  axis.ticks = element_blank()) +
            ggtitle(gene)
        

pdf(paste0(gene,"_plotEmbedding.pdf"),width=5,height=5)
print(p)
dev.off()

}




```


########################20_Gviz#################
#
#
```{r setup, include=FALSE}
library(ggplot2)
library(stringr)
library(magrittr)
library(cowplot)
library(ArchR)
library(Gviz)
library(cicero) #there is no package called 'cicero'
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(genomation) #


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("genomation")

```

## load reference genome
```{r}


nohup wget -c ftp://ftp.ensembl.org/pub/release-84/gtf/mus_musculus/Mus_musculus.GRCm38.84.gtf.gz &
gunzip Mus_musculus.GRCm38.84.gtf.gz


#gene_model <- readGFF("../../../Reference/refdata-cellranger-atac-mm10-1.2.0/genes/genes.gtf")

gene_model <- readGFF("../Reference/genes.gtf")

table(gene_model$seqid)

gene_model$chromosome <- gene_model$seqid
gene_model$gene <- gene_model$gene_id
gene_model$transcript <- gene_model$transcript_id
gene_model$symbol <- gene_model$gene_name
gene_model$exon <- gene_model$exon_id
gene_model$width <- gene_model$end - gene_model$start + 1
gene_model$feature <- gene_model$transcript_type

gene_model <- subset(gene_model, !is.na(transcript) & !is.na(exon))

gene_model %>%head()

```


## load data
```{r}


df_links <- read.csv("../Fibroblast/p2g_scOpen.csv", row.names = 1)
df_links <- subset(df_links, Distance > 0 & FDR < 0.05 & Correlation > 0)

df_links %>%head()
## load gene information
#
#


proj <- loadArchRProject("../Fibroblast")
geneMatrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")




df_gene <- as.data.frame(rowData(geneMatrix))
df_gene$gene_pos <- paste(df_gene$seqnames, 
                        df_gene$start, 
                        df_gene$start + 1,
                        sep = "_")

df_gene <- subset(df_gene, select = c("name", "gene_pos")) 
colnames(df_gene) <- c("gene", "gene_pos")
df <- merge.data.frame(df_links, df_gene)
df %>%head()



```



## load bigwig files
```{r}

data_track_ylim <- c(0, 4)



dt_t1 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T1.bw", 
                   genome = "mm10", type = "h", name = "T1", col = "#F8FA0D",
                   ylim = data_track_ylim)

dt_t2 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T2.bw", 
                   genome = "mm10", type = "h", name = "T2", col = "#F8BA43",
                   ylim = data_track_ylim)
dt_t3 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T3.bw", 
                   genome = "mm10", type = "h", name = "T3", col = "#2DB7A3",
                   ylim = data_track_ylim)
dt_t4 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T4.bw", 
                   genome = "mm10", type = "h", name = "T4", col = "#0262E0",
                   ylim = data_track_ylim)

dt_t5 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T5.bw", 
                   genome = "mm10", type = "h", name = "T5", col = "#352A86",
                   ylim = data_track_ylim)





load("../HINT/FootprintingFibroblast/Twist1.bed") #df_tf <- ,header = FALSE
 Twist1BindingSites                 
gr=Twist1BindingSites

gr <- GRanges(seqnames = seqnames(Twist1BindingSites), 
              IRanges(start = as.numeric(Twist1BindingSites$V2),
                      end = as.numeric(Twist1BindingSites$V3)))
###########################wojiade
granges(df_tf)
IRanges(df_tf)
mcols(gr)
granges(gr)
granges(gr) %>% as.data.frame() %>%head() %>%rownames()
df_tf1=granges(gr) %>% as.data.frame()
dim(df_tf1)
dim(mcols(gr))
df_tf=cbind(df_tf1,mcols(gr))
##################################
head(df_tf)
####列名有变化

colnames(df_tf) <- c("chromosome", "start", "end", "name",  "strand","score")
gr_tf <- GRanges(seqnames = df_tf$chromosome,
                     ranges = IRanges(start = df_tf$start,
                                      end = df_tf$end))


```

## plot 

# , "Twist2", "Col15a1"
# 
Twist1_2_g_raw$gene %>%head()

## plot 
```{r, fig.height=6, fig.width=6}
# , "Twist2", "Col15a1"
# 



for (sel_gene in Twist1_2_g_raw$gene ) {
    message(sprintf("plotting links for gene %s", sel_gene))
    
    df.plot <- subset(df, gene == sel_gene)
    
    peak_chr <- stringr::str_split_fixed(df.plot$peak, "_", 3)[, 1]
    peak_start <- stringr::str_split_fixed(df.plot$peak, "_", 3)[, 2]
    peak_end <- stringr::str_split_fixed(df.plot$peak, "_", 3)[, 3]
    
    ## highlight Twist1 binding sites
    gr_peaks <- GRanges(seqnames = peak_chr,
                     ranges = IRanges(start = as.numeric(peak_start),
                                      end = as.numeric(peak_end)))
    gr_overlap <- findOverlaps(query = gr_tf, 
                               subject = gr_peaks, 
                                type = "any",
                               ignore.strand = TRUE)
    df_tf_sub <- df_tf[gr_overlap@from, ]
    gr_tf_sub <- gr_tf[gr_overlap@from, ]
    
    print(df_tf_sub)
    
    df_tf_sub$symbol <- ""
    Twist1_track <- GeneRegionTrack(range = df_tf_sub,
                            genome = "mm10",
                            name = "Twist1",
                            shape = "box",
                            col = "blue",
                            fill = "blue",
                            fontsize = 12,
                            fontcolor = "black")
        
    gr_view <- resize(gr_tf_sub, width = 2000, fix = "center")
    peak_center <- (as.numeric(peak_start) + as.numeric(peak_end)) / 2
    
    df.plot$peak_pos <- paste(peak_chr, 
                            peak_center, 
                            peak_center + 1,
                            sep = "_")
    
    df.plot <- subset(df.plot, 
                      select = c("peak_pos", "gene_pos", "FDR"))
    colnames(df.plot) <- c("Peak1", "Peak2", "coaccess")
    df.plot$coaccess <- -log10(df.plot$coaccess)
    
    coaccess_cutoff <- 3
    
    peaks <- as.data.frame(str_split_fixed(df.plot$Peak2, "_", 3))
    
    chr <- unique(peak_chr)
    if(sel_gene == "Col15a1"){
        minbp <- min(as.numeric(as.character(peaks$V2))) - 2000
    } else{
        minbp <- min(as.numeric(df_tf_sub$start)) - 10000
        maxbp <- max(as.numeric(df_tf_sub$end)) + 10000        
    }
    
    connection_ymax <- max(df.plot$coaccess)
    
    link_tracks <- plot_connections(connection_df = df.plot, 
                                    gene_model = gene_model,
                                    gene_model_shape = "smallArrow",
                                    collapseTranscripts = "longest",
                                    chr = chr, 
                     minbp = minbp, 
                     maxbp = maxbp,
                     alpha_by_coaccess = FALSE,
                     connection_ymax = connection_ymax,
                     coaccess_cutoff = coaccess_cutoff,
                     connection_width = 1,
                     comparison_connection_color = "#9A6324",
                     include_axis_track = FALSE,
                     return_as_list = TRUE)
    link_track <- link_tracks[[1]]
    peak_track <- link_tracks[[2]]
    gene_model_track <- link_tracks[[3]]
    gene_axis_track <- GenomeAxisTrack(fontsize = 4)
    link_track@dp@pars$fontsize <- 12
    dt_t1@dp@pars$fontsize <- 12
    dt_t2@dp@pars$fontsize <- 12
    dt_t3@dp@pars$fontsize <- 12
    dt_t4@dp@pars$fontsize <- 12
    dt_t5@dp@pars$fontsize <- 12
    gene_model_track@dp@pars$fontsize <- 12
    trackList <- list(link_track,
                      peak_track,
                      dt_t1,
                      dt_t2,
                      dt_t3,
                      dt_t4,
                      dt_t5,
                      gene_model_track,
                      Twist1_track)
    
    trackList <- HighlightTrack(trackList = trackList, 
                                range = gr_view,
                                chromosome = chr,
                                col = "#F0544F", 
                                fill = "white",
                                #fill = "#EFD8D7", 
                                inBackground = FALSE, 
                                alpha = 0.3)
#grid.newpage()
plotTracks(trackList, title.width = 0.5, showTitle = TRUE, 
               from = minbp, to = maxbp, chromosome = chr, 
               sizes = c(0.2, 0.05, 
                         0.1, 0.1, 0.1, 0.1, 0.1, 
                         0.1, 0.05), 
               transcriptAnnotation = "symbol", 
               background.title = "white", 
               col.border.title = "transparent", 
               lwd.border.title = "transparent", 
               col.axis = "black", fontcolor.legend = "black",
               innerMargin = 3)
           }



```



##########################21 21_viz_p2g####################
#


```{r setup, include=FALSE}
library(ggplot2)
library(Seurat)
library(stringr)
library(magrittr)
library(readr)
library(Matrix)
library(tidyr)
library(dplyr)
library(plotly)
library(Signac)
library(cluster)
library(mclust)
library(cowplot)
library(gridExtra)
library(ArchR)
library(foreach)
library(WriteXLS)
library(caret)
```

```{r set_parameters, echo=FALSE}
## set parameters
set.seed(42)
addArchRThreads(threads = 100)
addArchRGenome("mm10")
```


```{r load_data}


proj <- loadArchRProject(path = "../data/Fibroblast")
geneMatrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
peakMatrix <- getMatrixFromProject(proj, useMatrix = "PeakMatrix")

peakMatrix
geneMatrix

```


## loadd p2g links
```{r load_p2g}



df_p2g <- read.csv("../Fibroblast/p2g_scOpen.csv", row.names = 1)
df_p2g <- subset(df_p2g, Distance > 0 & FDR < 0.05 & Correlation > 0)

df_p2g$peak_chr <- stringr::str_split_fixed(df_p2g$peak, "_", 3)[, 1]
df_p2g$peak_start <- stringr::str_split_fixed(df_p2g$peak, "_", 3)[, 2]
df_p2g$peak_end <- stringr::str_split_fixed(df_p2g$peak, "_", 3)[, 3]

gr_peaks <- GRanges(seqnames = df_p2g$peak_chr,
                         IRanges(start = as.numeric(df_p2g$peak_start),
                                 end = as.numeric(df_p2g$peak_end)))




```


## get Twist1 binding sites
```{r}



load("../HINT/FootprintingFibroblast/Twist1.bed") #Twist1BindingSites <- 
gr_Twist1=Twist1BindingSites

gr_Twist1 <- GRanges(seqnames = Twist1BindingSites$V1,
              IRanges(start = as.numeric(Twist1BindingSites$V2),
                      end = as.numeric(Twist1BindingSites$V3)))
```

## find overlapping
```{r}

o <- findOverlaps(query = gr_peaks, subject = gr_Twist1,
                  ignore.strand = TRUE)
o



## it's possible that one peak can have multiple Twist1 binding sites


df_p2g <- df_p2g[unique(as.integer(o@from)), ]
df_p2g$gene <- toupper(df_p2g$gene)
rownames(df_p2g) <- paste0(df_p2g$peak, "_", df_p2g$gene)




```

## load p2g links generated by using raw data
```{r}

df_p2g_raw <- read.csv("../Fibroblast/p2g_raw.csv", row.names = 1)
df_p2g_raw$gene <- toupper(df_p2g_raw$gene)
rownames(df_p2g_raw) <- paste0(df_p2g_raw$peak, "_", df_p2g_raw$gene)
df_p2g_raw <- df_p2g_raw[rownames(df_p2g), ]
df_p2g_raw %>%head()



```


## load data
```{r}

peaks <- as.data.frame(peakMatrix@rowRanges)
geneSet <- geneMatrix@elementMetadata
groupMatATAC_raw <- readRDS(file = "../Fibroblast/aggr_raw.rds")
rownames(groupMatATAC_raw) <- paste0(peaks$seqnames, "_", 
                                 peaks$start, "_", 
                                 peaks$end)



groupMatATAC_scOpen <- readRDS(file = "../Fibroblast/aggr_scOpen.rds")

rownames(groupMatATAC_scOpen) <- paste0(peaks$seqnames, "_", 
                                 peaks$start, "_", 
                                 peaks$end)

groupMatRNA <- readRDS(file = "../Fibroblast/aggr_rna.rds")

rownames(groupMatRNA) <- toupper(geneSet$name)
stopifnot(nrow(groupMatATAC_scOpen) == nrow(groupMatATAC_raw))


```


## plot individual link
```{r}


if(!dir.exists("../Fibroblast/Raw")){
    dir.create("../Fibroblast/Raw")
}
if(!dir.exists("../Fibroblast/scOpen")){
    dir.create("../Fibroblast/scOpen")
}


nrow(df_p2g)

df_p2g %>%head()

for (i in 1:nrow(df_p2g)) {
    
    gene <- df_p2g$gene[[i]]
    peak <- df_p2g$peak[[i]]
    
    correlation1 <- df_p2g$Correlation[[i]]
    correlation2 <- df_p2g_raw$Correlation[[i]]
    
    rna <- scale(groupMatRNA[gene, ])
    atac_raw <- scale(groupMatATAC_raw[peak, ])
    atac_scopen <- scale(groupMatATAC_scOpen[peak, ])
    df <- data.frame(rna = rna,
                     atac_raw = atac_raw,
                     atac_scopen = atac_scopen)
    
    df$PseudoTime <- 1:nrow(df)
    pal <- rev(ArchR::paletteContinuous(set = "blueYellow"))
    
    p1 <- ggplot(data = df, aes(x = rna, y = atac_scopen,
                               color = PseudoTime)) +
        geom_point() +
        scale_color_gradientn(colours = pal) +
        xlab("RNA") + ylab("ATAC") +
        ggtitle(glue::glue("Cor: {correlation1}")) +
        theme_cowplot()
    
    p2 <- ggplot(data = df, aes(x = rna, y = atac_raw,
                               color = PseudoTime)) +
        geom_point() +
        scale_color_gradientn(colours = pal) +
        xlab("RNA") + ylab("ATAC") +
        ggtitle(glue::glue("Cor: {correlation2}")) +
        theme_cowplot()

        
    write.table(df, file = glue::glue("../Fibroblast/scOpen/{gene}_{peak}.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)

    
    pdf(file = glue::glue("../Fibroblast/scOpen/{gene}_{peak}.pdf"),
        width = 5, height = 4)
    print(p1)
    dev.off()
    
    pdf(file = glue::glue("../Fibroblast/Raw/{gene}_{peak}.pdf"),
        width = 5, height = 4)
    print(p2)
    dev.off()


}





```

#############################22
#
```{r setup, include=FALSE}
library(ggplot2)
library(stringr)
library(magrittr)
library(cowplot)
library(ArchR)
library(Gviz)
library(cicero)
library(rtracklayer)
library(stringr)
library(GenomicRanges)
library(genomation)




```

## define function
```{r}
    gene_model <- readGFF("../../../Reference/refdata-cellranger-atac-mm10-1.2.0/genes/genes.gtf")
    gene_model$chromosome <- gene_model$seqid
    gene_model$gene <- gene_model$gene_id
    gene_model$transcript <- gene_model$transcript_id
    gene_model$symbol <- gene_model$gene_name
    gene_model$exon <- gene_model$exon_id
    gene_model$width <- gene_model$end - gene_model$start + 1
    gene_model$feature <- gene_model$transcript_type
    gene_model <- subset(gene_model, !is.na(transcript) & !is.na(exon))
    
```

## load bigwig files
```{r}
data_track_ylim <- c(0, 4)
dt_t1 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T1.bw", 
                   genome = "mm10", type = "h", name = "T1", col = "#F8FA0D",
                   ylim = data_track_ylim)
dt_t2 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T2.bw", 
                   genome = "mm10", type = "h", name = "T2", col = "#F8BA43",
                   ylim = data_track_ylim)
dt_t3 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T3.bw", 
                   genome = "mm10", type = "h", name = "T3", col = "#2DB7A3",
                   ylim = data_track_ylim)
dt_t4 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T4.bw", 
                   genome = "mm10", type = "h", name = "T4", col = "#0262E0",
                   ylim = data_track_ylim)
dt_t5 <- DataTrack(range = "../../PseudoBulkTrajectory/BigWig/T5.bw", 
                   genome = "mm10", type = "h", name = "T5", col = "#352A86",
                   ylim = data_track_ylim)




df_tf <- read.table("../HINT/FootprintingFibroblast/Twist1.bed")#, header = FALSE

colnames(df_tf) <- c("chromosome", "start", "end", "name", "score", "strand")

gr_tf <- GRanges(seqnames = df_tf$chromosome,
                     ranges = IRanges(start = df_tf$start,
                                      end = df_tf$end))
```

## load gene
```{r}
## load gene information

proj <- loadArchRProject("../data/Fibroblast")
geneMatrix <- getMatrixFromProject(proj, useMatrix = "GeneScoreMatrix")
df_gene <- as.data.frame(rowData(geneMatrix))
df_gene$gene_pos <- paste(df_gene$seqnames, 
                        df_gene$start, 
                        df_gene$start + 1,
                        sep = "_")

df_gene <- subset(df_gene, select = c("name", "gene_pos")) 
colnames(df_gene) <- c("gene", "gene_pos")



```


## define function
```{r}

plot_links <- function(df_links_p2g,
                       coaccess_cutoff_p2g,
                       df_links_p2p_1,
                       coaccess_cutoff_p2p_1,
                       df_links_p2p_2,
                       coaccess_cutoff_p2p_2,
                       sel_gene) {
    df_links_p2g <- subset(df_links_p2g, gene == sel_gene)
    
    peak_chr <-
        stringr::str_split_fixed(df_links_p2g$peak, "_", 3)[, 1]
    peak_start <-
        stringr::str_split_fixed(df_links_p2g$peak, "_", 3)[, 2]
    peak_end <-
        stringr::str_split_fixed(df_links_p2g$peak, "_", 3)[, 3]
    
    ## highlight binding sites
    gr_peaks <- GRanges(seqnames = peak_chr,
                        ranges = IRanges(start = as.numeric(peak_start),
                                         end = as.numeric(peak_end)))
    gr_overlap <- findOverlaps(
        query = gr_tf,
        subject = gr_peaks,
        type = "any",
        ignore.strand = TRUE
    )
    df_tf_sub <- df_tf[gr_overlap@from,]
    gr_tf_sub <- gr_tf[gr_overlap@from,]
    df_tf_sub$symbol <- ""
    Twist1_track <- GeneRegionTrack(
        range = df_tf_sub,
        genome = "mm10",
        name = "Twist1",
        shape = "box",
        col = "blue",
        fill = "blue",
        fontsize = 12,
        fontcolor = "black"
    )
    
    # genomic coordinates
    chr <- unique(peak_chr)
    minbp <- min(as.numeric(df_tf_sub$start)) - 10000
    maxbp <- max(as.numeric(df_tf_sub$end)) + 10000
    
    
    # subset peak-to-peak links
    gene_pos <- unique(df_links_p2g$gene_pos)
    
    o <- find_overlapping_coordinates(c(df_links_p2p_1$Peak1,
                                        df_links_p2p_1$Peak2),
                                      gene_pos,
                                      maxgap = 1000)
    df_links_p2p_1 <- df_links_p2p_1[df_links_p2p_1$Peak1 %in% o |
                                         df_links_p2p_1$Peak2 %in% o,]
    df_links_p2p_1 <-
        subset(df_links_p2p_1, coaccess > coaccess_cutoff_p2p_1)
    
    o <- find_overlapping_coordinates(c(df_links_p2p_2$Peak1,
                                        df_links_p2p_2$Peak2),
                                      gene_pos,
                                      maxgap = 1000)
    df_links_p2p_2 <- df_links_p2p_2[df_links_p2p_2$Peak1 %in% o |
                                         df_links_p2p_2$Peak2 %in% o,]
    df_links_p2p_2 <-
        subset(df_links_p2p_2, coaccess > coaccess_cutoff_p2p_1)
    
    # resize peaks for peak-to-gene links
    peak_center <-
        (as.numeric(peak_start) + as.numeric(peak_end)) / 2
    df_links_p2g$peak_pos <- paste(peak_chr,
                                   peak_center,
                                   peak_center + 1,
                                   sep = "_")
    df_links_p2g <-
        subset(df_links_p2g,  select = c("peak_pos", "gene_pos", "FDR"))
    colnames(df_links_p2g) <- c("Peak1", "Peak2", "coaccess")
    df_links_p2g$coaccess <- -log10(df_links_p2g$coaccess)
    
    p2g_track <- plot_connections(
        connection_df = df_links_p2g,
        gene_model = gene_model,
        gene_model_shape = "smallArrow",
        collapseTranscripts = "longest",
        chr = chr,
        minbp = minbp,
        maxbp = maxbp,
        alpha_by_coaccess = FALSE,
        connection_ymax = max(df_links_p2g$coaccess),
        coaccess_cutoff = coaccess_cutoff_p2g,
        connection_width = 1,
        include_axis_track = FALSE,
        return_as_list = TRUE
    )
    
    p2p_track1 <- plot_connections(
        connection_df = df_links_p2p_1,
        gene_model = gene_model,
        gene_model_shape = "smallArrow",
        collapseTranscripts = "longest",
        chr = chr,
        minbp = minbp,
        maxbp = maxbp,
        alpha_by_coaccess = FALSE,
        connection_ymax = max(df_links_p2p_1$coaccess),
        coaccess_cutoff = coaccess_cutoff_p2p_1,
        connection_width = 1,
        include_axis_track = FALSE,
        return_as_list = TRUE
    )
    
     p2p_track2 <- plot_connections(
        connection_df = df_links_p2p_2,
        gene_model = gene_model,
        gene_model_shape = "smallArrow",
        collapseTranscripts = "longest",
        chr = chr,
        minbp = minbp,
        maxbp = maxbp,
        alpha_by_coaccess = FALSE,
        connection_ymax = max(df_links_p2p_1$coaccess),
        coaccess_cutoff = coaccess_cutoff_p2p_1,
        connection_width = 1,
        include_axis_track = FALSE,
        return_as_list = TRUE
    )   
    
    link_track1 <- p2g_track[[1]]
    peak_track1 <- p2g_track[[2]]
    gene_model_track <- p2g_track[[3]]
    
    link_track2 <- p2p_track1[[1]]
    peak_track2 <- p2p_track1[[2]]
    
    link_track3 <- p2p_track2[[1]]
    peak_track3 <- p2p_track2[[2]]
    
    link_track1@dp@pars$fontsize <- 12
    link_track2@dp@pars$fontsize <- 12
    link_track3@dp@pars$fontsize <- 12
    
    dt_t1@dp@pars$fontsize <- 12
    dt_t2@dp@pars$fontsize <- 12
    dt_t3@dp@pars$fontsize <- 12
    dt_t4@dp@pars$fontsize <- 12
    dt_t5@dp@pars$fontsize <- 12
    gene_model_track@dp@pars$fontsize <- 12
    
    trackList <- list(
        link_track1,
        peak_track1,
        link_track2,
        peak_track2,
        link_track3,
        peak_track3,
        dt_t1,
        dt_t2,
        dt_t3,
        dt_t4,
        dt_t5,
        gene_model_track,
        Twist1_track
    )
    
    plotTracks(
        trackList,
        title.width = 0.5,
        showTitle = TRUE,
        from = minbp,
        to = maxbp,
        chromosome = chr,
        sizes = c(0.2, 0.05, 0.2, 0.05, 0.2, 0.05,
                  0.1, 0.1, 0.1, 0.1, 0.1,
                  0.1, 0.05),
        transcriptAnnotation = "symbol",
        background.title = "white",
        col.border.title = "transparent",
        lwd.border.title = "transparent",
        col.axis = "black",
        fontcolor.legend = "black",
        innerMargin = 3
    )
}




```

## load data
```{r}


df_links_p2g <- read.csv("../Fibroblast/p2g_scOpen.csv", row.names = 1)
df_links_p2g <- subset(df_links_p2g, Distance > 0 & FDR < 0.05 & Correlation > 0)
df_links_p2g <- merge.data.frame(df_links_p2g, df_gene)



df_links_p2p_scopen <- read.table("../../Cicero/scOpen/Cicero.txt", 
                                  header = TRUE) %>%
    subset(., coaccess > 0)
df_links_p2p_cistopic <- read.table("../../Cicero/cisTopic/Cicero.txt", 
                                  header = TRUE) %>%
    subset(., coaccess > 0)
df_links_p2p_raw <- read.table("../../Cicero/Raw/Cicero.txt", 
                                  header = TRUE) %>%
    subset(., coaccess > 0)
```



## plot for Tgfbr1
```{r, fig.height=8, fig.width=6}

plot_links(df_links_p2g = df_links_p2g, 
           coaccess_cutoff_p2g = 3,
           df_links_p2p_1 = df_links_p2p_scopen, 
           coaccess_cutoff_p2p_1 = 0.2,
           df_links_p2p_2 = df_links_p2p_raw, 
           coaccess_cutoff_p2p_2 = 0.0,
           sel_gene = "Tgfbr1")

plot_links(df_links_p2g = df_links_p2g, 
           coaccess_cutoff_p2g = 3,
           df_links_p2p_1 = df_links_p2p_scopen, 
           coaccess_cutoff_p2p_1 = 0.2,
                      df_links_p2p_2 = df_links_p2p_raw, 
           coaccess_cutoff_p2p_2 = 0.0,
           sel_gene = "Twist2")


```

