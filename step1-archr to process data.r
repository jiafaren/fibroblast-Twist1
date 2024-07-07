conda create -n scopen python=3.8
conda activate scopen
pip install scopen -i  https://pypi.douban.com/simple/
pip install MACS2 -i  https://pypi.douban.com/simple/ 




#####install package archr
.libPaths(c("/home/data/t040413/R/yll/usr/local/lib/R/site-library",  "/home/data/t040413/R/x86_64-pc-linux-gnu-library/4.2", "/usr/local/lib/R/library"))




nohup wget -c ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119531/suppl/GSE119531_Healthy.combined.dge.txt.gz  &
 
###################01 create archRproject from compressed file##
 tss.enrichemnt 

#tar file
tar -xvf GSE139950_RAW

#setwd("/home/yll/yll")

#setwd("/home/abc/scope_abc")

fs=list.files('./mydownloads/','^GSM') 
fs


#
library(stringr)
#samples=str_split(fs,'_',simplify = T)[,2] 

samples=paste0(str_split(fs,'_',simplify = T)[,2],"_",str_split(fs,'_',simplify = T)[,3])
samples
getwd()#"/home/data/t040413/scope_abc"

lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  #y=c(TRUE  ,TRUE  ,TRUE  ,TRUE, FALSE ,FALSE ,FALSE, FALSE ,FALSE ,FALSE, FALSE)
  #####
  folder=paste0("mydownloads/", paste0(str_split(y[1],'_',simplify = T)[,2], 
                                       str_split(y[1],'_',simplify = T)[,3])     )
  dir.create(folder,recursive = T)
  
  
  file.rename(paste0("mydownloads/",y[1]),file.path(folder,"barcodes.tsv.gz"))
 
  file.rename(paste0("mydownloads/",y[2]),file.path(folder,"fragments.tsv.gz"))
  file.rename(paste0("mydownloads/",y[3]),file.path(folder,"matrix.mtx.gz"))
  file.rename(paste0("mydownloads/",y[4]),file.path(folder,"peaks.bed.gz"))
  
  })
getwd()
#samples=list.files("mydownloads//")
samples



# "/home/abc/scope_abc"
## Creating Arrow Files

inputFiles <- c("D0_1" = "./mydownloads/Day01/fragments.tsv.gz",
                "D0_2" = "./mydownloads/Day02/fragments.tsv.gz",
                "D2_1" = "./mydownloads/Day21/fragments.tsv.gz",
                "D2_2" = "./mydownloads/Day22/fragments.tsv.gz",
                "D10_1" = "./mydownloads/Day101/fragments.tsv.gz",
                "D10_2" = "./mydownloads/Day102/fragments.tsv.gz")
filterTSS <- 8
filterFrags <- 1000


library(ggplot2)
library(stringr)
library(magrittr)
library(WriteXLS)
library(tidyr)
library(dplyr)
library(plotly)
library(cluster)
library(cowplot)
library(gridExtra)
library(viridis)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)
library(ArchR)


## set parameters
set.seed(42)
addArchRThreads(threads = parallel::detectCores() - 12)
addArchRGenome("mm10")
```

## Creating Arrow Files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  outputNames = names(inputFiles),
  minTSS = filterTSS, 
  minFrags = filterFrags, 
  QCDir = "QualityControl",
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles
```



## Plotting
```{r, fig.width=6, fig.height=6}
for(sample in unique(names(inputFiles))){
    input_filename <- sprintf("./QualityControl/%s/%s-Pre-Filter-Metadata.rds", sample, sample)
    
    if(file.exists(input_filename)){
        Metadata <- readRDS(input_filename)
    
        ggtitle <- sprintf("%s\n%s\n%s",
            paste0(sample, "\nnCells Pass Filter = ", sum(Metadata$Keep)),
            paste0("Median Frags = ", median(Metadata$nFrags[Metadata$Keep==1])),
            paste0("Median TSS Enrichment = ", median(Metadata$TSSEnrichment[Metadata$Keep==1]))
          )
    
        gg <- ggPoint(
          x = pmin(log10(Metadata$nFrags), 5) + rnorm(length(Metadata$nFrags), sd = 0.00001),
          y = Metadata$TSSEnrichment + rnorm(length(Metadata$nFrags), sd = 0.00001), 
          colorDensity = TRUE,
          xlim = c(2.5, 5),
          ylim = c(0, max(Metadata$TSSEnrichment) * 1.05),
          baseSize = 6,
          continuousSet = "sambaNight",
          xlabel = "Log 10 (Unique Fragments)",
          ylabel = "TSS Enrichment",
          title = ggtitle,
          rastr = TRUE) + 
          geom_hline(yintercept=filterTSS, lty = "dashed", size = 0.25) +
          geom_vline(xintercept=log10(filterFrags), lty = "dashed", size = 0.25)
        
        print(gg)
    }
}
```


## Inferring Doublets
```{r, fig.width=6, fig.height=6}
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)
```


## Creating an ArchRProject
```{r, fig.width=6, fig.height=6}
# With our Arrow files in hand, we are now ready to create an ArchRProject. An ArchRProject is associated with a set of Arrow files and is the backbone of nearly all ArchR analyses.
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "UUO",
  showLogo = FALSE,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

p1 <- plotGroups(ArchRProj = proj, 
                 groupBy = "Sample", 
                 colorBy = "cellColData", 
                 name = "TSSEnrichment",
                 alpha = 0.4,
                 plotAs = "violin",
                 addBoxPlot = TRUE)
p2 <- plotGroups(ArchRProj = proj, 
                 groupBy = "Sample", 
                 colorBy = "cellColData", 
                 name = "log10(nFrags)",
                 plotAs = "violin",
                 alpha = 0.4,
                 addBoxPlot = TRUE)
print(p1)
print(p2)
```

## save data
```{r}
# Now we can filter putative doublets based on the previously determined doublet scores using the filterDoublets() function. This doesn’t physically remove data from the Arrow files but rather tells the ArchRProject to ignore these cells for downstream analysis.
proj <- filterDoublets(ArchRProj = proj)
saveArchRProject(ArchRProj = proj, 
                 load = FALSE)

```






########################################——##02_create_mat.Rmd########
#1.peak calling（pseudo，addGroupCoverages），mscs2
#2 addPeakMatrix, save as the file named "./PeakMatrix.Rds"
#

#3save  "./filtered_peak_bc_matrix/barcodes.tsv" "./filtered_peak_bc_matrix/peaks.bed"   "./filtered_peak_bc_matrix/matrix.mtx"    "GeneScoreMatrix.Rds")


#####install package archr##
.libPaths(c("/home/data/t040413/R/yll/usr/local/lib/R/site-library",  "/home/data/t040413/R/x86_64-pc-linux-gnu-library/4.2", "/usr/local/lib/R/library"))

```{r setup, include=FALSE}
library(ggplot2)
library(stringr)
library(magrittr)
library(WriteXLS)
library(tidyr)
library(dplyr)
library(plotly)
library(cluster)
library(cowplot)
library(gridExtra)
library(viridis)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)
library(ArchR)
```


## set parameters
set.seed(42)
addArchRThreads(threads = 50)
addArchRGenome("mm10")

proj2 <- loadArchRProject(path = "./Save-ProjHeme2")

## Load project
```{r}
proj <- loadArchRProject(path = "./UUO")
```

## Peak calling
```{r}
proj <- addGroupCoverages(ArchRProj = proj, 
                          groupBy = "Sample",
                          maxCells = 3000,
                          force = TRUE)

pathToMacs2 <- "/home/data/t040413/anaconda3/envs/six/bin/macs2"#findMacs2() #/home/data/t040413/anaconda3/envs/six/bin/macs3  /home/data/t040413/anaconda3/envs/six/bin/macs2
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Sample", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(proj)
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)
```


## save data
```{r}
peakMatrix <- getMatrixFromProject(proj,
                                   useMatrix = "PeakMatrix")
counts <- peakMatrix@assays@data$PeakMatrix
df_rangers <- as.data.frame(peakMatrix@rowRanges@ranges)
rownames(counts) <- paste(peakMatrix@rowRanges@seqnames,
                          df_rangers$start,
                          df_rangers$end,
                          sep = "_") 
saveRDS(counts, file = "./PeakMatrix.Rds")
if(!dir.exists("./filtered_peak_bc_matrix")){
    dir.create("./filtered_peak_bc_matrix")
}

writeMM(counts, file = "./filtered_peak_bc_matrix/matrix.mtx")

barcodes <- as.data.frame(colnames(counts))

peaks <- as.data.frame(stringr::str_split_fixed(rownames(counts), "_", 3))

write.table(barcodes, file = "./filtered_peak_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(peaks, file = "./filtered_peak_bc_matrix/peaks.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
atac <- getMatrixFromProject(ArchRProj = proj,
                             useMatrix = "GeneScoreMatrix")
counts <- atac@assays@data$GeneScoreMatrix

rownames(counts) <- atac@elementMetadata$name

saveRDS(counts, file = "GeneScoreMatrix.Rds")

df <- as.data.frame(proj@cellColData)
write.csv(df, file = "meta_data.csv", quote = FALSE, row.names = TRUE)
saveArchRProject(ArchRProj = proj, 
                 load = FALSE)
```


###########
##############################3
#################################################################333
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-ProjHeme2",
                 load = FALSE)

setwd("/home/data/t040413/scope_abc/")
#####install package archr##
.libPaths(c("/home/data/t040413/R/yll/usr/local/lib/R/site-library",  "/home/data/t040413/R/x86_64-pc-linux-gnu-library/4.2", "/usr/local/lib/R/library"))


library(ggplot2)
library(stringr)
library(magrittr)
library(WriteXLS)
library(tidyr)
library(dplyr)
library(plotly)
library(cluster)
library(cowplot)
library(gridExtra)
library(viridis)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)
library(ArchR)


## set parameters
set.seed(42)
addArchRThreads(threads = 50)
addArchRGenome("mm10")

proj <- loadArchRProject(path = "./UUO")



genome=proj@genomeAnnotation$genome 
genome
library(BSgenome)
library(BSgenome.Mmusculus.UCSC.mm10)
BSgenome <- eval(parse(text = genome)) 
BSgenome <- validBSgenome(BSgenome) 

proj <- addGroupCoverages(ArchRProj = proj, 
                          groupBy = "Sample",
                          maxCells = 3000,
                          force = TRUE)

pathToMacs2 <- "/home/data/t040413/anaconda3/envs/six/bin/macs2"

proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "Sample", 
  pathToMacs2 = pathToMacs2
)
getPeakSet(proj)
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

#################save data
peakMatrix <- getMatrixFromProject(proj,
                                   useMatrix = "PeakMatrix")
peakMatrix
assay(peakMatrix)
rowRanges(peakMatrix)@seqnames


counts <- peakMatrix@assays@data$PeakMatrix
counts[1:4,1:5]
df_rangers <- as.data.frame(peakMatrix@rowRanges@ranges)
rownames(counts) <- paste(peakMatrix@rowRanges@seqnames,
                          df_rangers$start,
                          df_rangers$end,
                          sep = "_") 
saveRDS(counts, file = "./PeakMatrix.Rds")
dim(counts) #152299  30129
##

if(!dir.exists("./filtered_peak_bc_matrix")){
  dir.create("./filtered_peak_bc_matrix")
}

writeMM(counts, file = "./filtered_peak_bc_matrix/matrix.mtx")

barcodes <- as.data.frame(colnames(counts))

peaks <- as.data.frame(stringr::str_split_fixed(rownames(counts), "_", 3))

write.table(barcodes, file = "./filtered_peak_bc_matrix/barcodes.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(peaks, file = "./filtered_peak_bc_matrix/peaks.bed", 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

atac <- getMatrixFromProject(ArchRProj = proj,
                             useMatrix = "GeneScoreMatrix")


##gsm
counts <- atac@assays@data$GeneScoreMatrix

rownames(counts) <- atac@elementMetadata$name

saveRDS(counts, file = "GeneScoreMatrix.Rds")

df <- as.data.frame(proj@cellColData)
write.csv(df, file = "meta_data.csv", quote = FALSE, row.names = TRUE)
head(df)
dim(df)
saveArchRProject(ArchRProj = proj, outputDirectory = "Save-ProjHeme2",
                 load = FALSE)
