library(dplyr)
library(Seurat)
install.packages("devtools")
devtools::install_github("jlmelville/uwot")
pkgbuild::check_build_tools(debug=TRUE)
library(uwot)
library(patchwork)
library(Matrix)

sample.data <- Read10X(data.dir = "C:\\Users\\USER\\Downloads\\filtered_gene_bc_matrices\\hg19\\")
sample <- CreateSeuratObject(counts=sample.data,project ="pbmc3k",min.cells=3)
sample
sample.data[0:3,0:3]
rawData <- read.csv(file = "", sep = ",", header = FALSE )
rawData[c("100","10"),1:3]
rawData <- CreateSeuratObject(counts=rawData)
head(rawData@meta.data, 5)
#create sce object
rawSce <- SingleCellExperiment(list(counts = data.matrix(rawData)))
rawSce
dataz <- CreateSeuratObject(counts = rawData)
dataz[["percent.mt"]] <- PercentageFeatureSet(dataz, pattern = "^MT")
VlnPlot(dataz, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


rawData <- transform(rawData, geneId=paste("MT",gene_id,sep=""))
gene <- as.matrix(subset(rawData, select = -c(gene_id, geneId)))
gene <- cbind(subset(rawData, select = c(geneId)),gene)



geneData <- as.matrix(subset(rawData, select = -c(gene_id)))

temp <- subset(rawData,select=c(gene_id))
temp <- transform(temp, geneId=paste("MT",gene_id,sep=""))
temp[3,]
geneData <- cbind(subset(temp,select=c(geneId)),geneData)
geneData[3,]

geneData <- cbind(subset(rawData, select = c(gene_id)),geneData)

cancer <- CreateSeuratObject(counts = geneData)

#calculate QC metrics
cancer[["percent.mt"]] <- PercentageFeatureSet(cancer, pattern = "^MT-")
VlnPlot(cancer, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(cancer, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cancer, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#feature selection
cancer <- FindVariableFeatures(object = cancer, nfeatures = 10000)
plot1 <- VariableFeaturePlot(cancer)
plot1

#pca
all.genes <- rownames(cancer)
cancer <- ScaleData(cancer, features = all.genes)
cancer <- RunPCA(cancer, features = VariableFeatures(object=cancer))
print(cancer[["pca"]],dims=1:3,nfeatures=3)

#pca 시각화
VizDimLoadings(cancer, dims = 1:2, reduction = "pca")
DimPlot(cancer, reduction = "pca")
DimHeatmap(cancer, dims = 1, cells = 70, balanced = TRUE)

#choose pc
cancer <- JackStraw(cancer, num.replicate = 100)
cancer <- ScoreJackStraw(cancer, dims = 1:20)
JackStrawPlot(cancer, dims = 1:15)

ElbowPlot(cancer)
