library(Seurat)
library(hdf5r)
library(zellkonverter)
library(SingleCellExperiment)

# Replace 'path/to/your/file.h5ad' with the actual file path
file_path <- "GSE180665_hb_integrated_normalized_annotated_harmony.h5ad"

# Read the .h5ad file into an AnnData object
sce <- readH5AD("cairo.h2ad")

scee <- sce

names(assays(scee))=c("counts")

assay(scee,"logcounts")<-assay(scee,"counts")

data<-as.seurat(scee)


load("cairo.rda")
str(data)
x<-data[[]]

meta<-read.csv("meta.csv",h=T,row.names=1)

all(row.names(meta)==row.names(x))

data <- AddMetaData(object = data, metadata = meta)

data<-FindVariableFeatures(data)
data<-ScaleData(data)
data<-RunPCA(data)
ElbowPlot(data,ndims = 50)
data<-RunUMAP(data,dims=1:30)



Idents(object = data) <- 'Cell.class_concise'
Idents(object = data)

library(pals)
DimPlot(data,reduction="umap",group.by ="sample",pt.size=1,cols=cols25(),label=T)



### integration
data[["originalexp"]] <- split(data[["originalexp"]], f = data$sample)

data <- IntegrateLayers(
  object = data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


DimPlot(data, reduction = "harmony",group.by="sample",cols=cols25(),label=F)
DimPlot(data, reduction = "harmony",group.by="Cell.class_concise",cols=cols25(),label=F)
save(all,file="harmony.rda")



all$class<-x$class
DimPlot(data, reduction = "umap",group.by="Cell.class_concise")


data <- JoinLayers(data)
data
data<-RunPCA(data)
ElbowPlot(data,ndims = 50)


data<-RunUMAP(data,dims=1:30, reduction = "harmony")

ElbowPlot(data, ndims = 50, reduction = "harmony")
#obj <- RunUMAP(obj, dims = 1:30, reduction = "harmony")
save(data,file="all.rda")

DimPlot(data,reduction="umap",group.by ="sample_group",pt.size=0.5,cols=cols25(),label=F)

DimPlot(data,reduction="umap",group.by ="RNA_snn_res.0.5",pt.size=0.5,cols=cols25(),label=T)
save(data,file="harmony.rda")


## ferroptosis score down

library(Seurat)

list.files()

load("harmony.rda")

cross<-read.table("crossresults.tsv",h=T,sep="\t")

library(dplyr)
cross%>%filter(group=="down.down")%>%distinct()%>%pull(GSE104766_hs.gene)->vector
cross%>%distinct()%>%pull(GSE104766_hs.gene)->vector
equation<-paste(vector,collapse="','")






ferroptosis.down <- list(c('GLS2','SOCS2','GOT1','ENO3','DHODH','HPX','TFR2','DECR1','GCH1','ANO1','ESR1','SCP2','CP','ERN1',
'STEAP3','AQP3','RND1','MPC1','SATB1','HMGCL','FNDC5','AR','TLR4','ACSL1','CBS','MAP3K14','ANGPTL4',
'TTPA','IL1B','AADAC','PTGS2'))
data<- AddModuleScore(
  object = data,
  features = ferroptosis.down, name = 'ferroptosis.down',
ctrl=31
)

library(pals)
VlnPlot(data, features = c("ferroptosis.down1"), slot = "data", log = TRUE,pt.size=1,split.by="Cell.class_concise",col=cols25())

colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)
FeaturePlot(data,reduction="umap",pt.size=0.1,features="ferroptosis.down1",min.cutoff = "q9",split.by="sample_group",col=c("grey90","darkred"))+
	plotTheme+coord_fixed()

DotPlot(data,group.by="Cell.class",features=vector)+
coord_flip() + scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -45, hjust = 0))

DotPlot(data,group.by="sample_group",features=vector)+
 scale_color_gradientn(colors = colgex) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0))




library(viridis)
library(RColorBrewer)
library(ggplot2)
FeaturePlot(data,"ferroptosis.score1",reduction = "umap",pt.size=0.1,min.cutoff = "q9") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) + 
ggtitle("ferroptosis.score1")

DimPlot(data,reduction="umap",group.by ="sample_group",pt.size=0.1,cols=cols25(),label=F)


## neuronal cell subset analyses
library(Seurat)

load("harmony.rda")

neu<-subset(data,idents="Neuro")

rm(data)
gc()

neu<-NormalizeData(neu)

neu<-FindVariableFeatures(neu)

neu<-ScaleData(neu)

neu<-RunPCA(neu)

ElbowPlot(neu)

neu<-RunTSNE(neu,dims=1:5)

DimPlot(neu,reduction="tsne",group.by ="sample_group")

vector<-c("NQO1","TGFB2","PRKAA2","ACSL4","IGF2BP3","SRC","SLC7A11","BEX1")

colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=12)
library(ggplot2)
FeaturePlot(neu,reduction="tsne",pt.size=2,features="TRIB2",min.cutoff = "q9",split.by="sample_group",col=c("grey85","navy"))+plotTheme


gene_vector <- vector[vector %in% rownames(neu)]

expression <- FetchData(neu, vars = gene_vector)

library(dplyr)

expression%>%mutate(NQO1_cat=ifelse(NQO1>0,"POSITIVE","negative"))->expression
expression%>%mutate(TGFB2_cat=ifelse(TGFB2>0,"POSITIVE","negative"))->expression
expression%>%mutate(ACSL4_cat=ifelse(ACSL4>0,"POSITIVE","negative"))->expression
expression%>%mutate(PRKAA2_cat=ifelse(PRKAA2>0,"POSITIVE","negative"))->expression
expression%>%mutate(IGF2BP3_cat=ifelse(IGF2BP3>0,"POSITIVE","negative"))->expression
expression%>%mutate(SRC_cat=ifelse(SRC>0,"POSITIVE","negative"))->expression
expression%>%mutate(SLC7A11_cat=ifelse(SLC7A11>0,"POSITIVE","negative"))->expression
expression%>%mutate(BEX1_cat=ifelse(BEX1>0,"POSITIVE","negative"))->expression

df<-preprocess(expression)

results<-cltest(df)

nlpplot(results,titleplot="",txt=18)

chinet(results,fold=2,cex=1.5,distance=3,family="sans",layout=layout_nicely)


neu <- AddMetaData(neu, metadata = expression)

save(neu,file="neuro.rda")


## single cell trajectory

load("neuro.rda")
library(Seurat)

##cell cyle phase prediction
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

neu <- CellCycleScoring(neu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(neu[[]])

##red score
red.score <- list(c('TGFB2','PRKAA2','SLC7A11'))
neu<- AddModuleScore(
  object = neu,
  features = red.score, name = 'red.score',
ctrl=3
)


library(viridis)
library(RColorBrewer)
library(ggplot2)
colgex=c("grey90",brewer.pal(7,"Reds"))
plotTheme=theme_classic(base_size=16)
FeaturePlot(neu,reduction="tsne",pt.size=2,features="red.score1",min.cutoff = "q9",split.by="sample_group",col=c("grey90","darkred"))+plotTheme


##blue score
blue.score <- list(c('BEX1','NQO1','ACSL4','SRC','IGFBP3'))
neu<- AddModuleScore(
  object = neu,
  features = blue.score, name = 'blue.score',
ctrl=5
)
FeaturePlot(neu,reduction="tsne",pt.size=2,features="blue.score1",min.cutoff = "q9",split.by="sample_group",col=c("grey90","darkred"))+plotTheme

library(pals)
DimPlot(neu,reduction="tsne",group.by ="Phase",pt.size=2,cols=cols25(),label=F)

VlnPlot(neu, features = c("red.score1"), slot = "data", log = TRUE,pt.size=1,split.by="Phase",col=cols25())


## single cell experiment
library(SingleCellExperiment)
neural.sce <- as.SingleCellExperiment(neu)

## cell entropy
library(TSCAN)
entropy<-perCellEntropy(neural.sce)


ent.data<-data.frame(cluster=neu$Phase,entropy=entropy)
head(ent.data)
library(dplyr)


ent.data$code<-as.numeric(ent.data$code)
ent.data$cluster<-as.factor(ent.data$Phase)
library(ggplot2)

ggplot(ent.data,aes(reorder(cluster,cluster),entropy))+
geom_violin(aes(fill = cluster), trim = FALSE) + 
geom_boxplot(width = 0.4,outlier.shape = NA)+
scale_fill_manual(values = c("blue", "red","green"))+
	#scale_fill_brewer(palette="Set1")+
  #geom_point(aes(fill=factor(cluster),size=0.5),shape = 21, alpha = .8, position = position_dodge2(width = .5))+
  theme_classic(base_size=18) +
  theme(legend.position = "none")+xlab("HB neural cells")+ggtitle("Cell entropy")+ 
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

neu$entropy<-ent.data$entropy

FeaturePlot(neu,reduction="tsne",pt.size=2,features="entropy",min.cutoff = "q9",split.by="sample_group",col=c("grey90","darkred"))+plotTheme

library(pals)
DimPlot(neu,reduction="tsne",group.by ="Phase",pt.size=2,cols=cols25(),label=F)

VlnPlot(neu, features = c("red.score1"), slot = "data", log = TRUE,pt.size=1,split.by="Phase",col=cols25())


library(slingshot)
library(Seurat)


sce <- slingshot(neural.sce, reducedDim = 'TSNE')  # no clusters

# Plot PC1 vs PC2 colored by Slingshot pseudotime.
colors <- rainbow(50, alpha = 1)
plot(reducedDims(sce)$TSNE, col = colors[cut(sce$slingPseudotime_1,breaks=50)], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2)

library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(dplyr)
library(cowplot)
library(RColorBrewer)
library(knitr)

library(ggbeeswarm)
library(ggthemes)
# Plot Slingshot pseudotime vs Phase. 
ggplot(as.data.frame(colData(sce)), aes(x = sce$slingPseudotime_1, y = Phase, 
                              colour = Phase)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau() + theme_classic() +
    xlab("Slingshot pseudotime") + ylab("Timepoint") +
    ggtitle("Cells ordered by Slingshot pseudotime")

ggsave(paste0("pseudotime_slingshot.svg"))

library(gam)

# Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
Y <- log2(counts(sce) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:1000]
Y <- Y[var1K, ]  # only counts for variable genes

# Fit GAM for each gene using pseudotime as independent variable.
t <- sce$slingPseudotime_1
gam.pval <- apply(Y, 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(gam.pval, decreasing = FALSE))[00:70]  
top70 <- as.data.frame(sort(gam.pval, decreasing = FALSE)[00:70]) 

colnames(top70)<-"gam_pvalues"
top70$gene<-row.names(top70)


# Identify genes with the most significant time-dependent model fit.
topgenes <- names(sort(sce$slingPseudotime_1, decreasing = FALSE))[00:70]  

topgenes<-rownames(sorted)[1:70]
# Prepare and plot a heatmap of the top genes that vary their expression over pseudotime.

## best gam genes
library(scran)
library(scater)
sce$entropy<-entropy
on.first.path <- !is.na(sce$slingPseudotime_1)

plotHeatmap(sce[,on.first.path], order_columns_by=c("slingPseudotime_1","red.score1","entropy","blue.score1"), 
    colour_columns_by=c("Phase"), features=topgenes,
    center=TRUE)

plotHeatmap(sce[,on.first.path], order_columns_by=c("slingPseudotime_1","red.score1"), 
    colour_columns_by=c("Phase"), features=topgenes,
    center=TRUE)


write.table(top70,file="top70.tsv",row.names=T,sep="\t")

write.table(gam.pval,file="gam.tsv",row.names=T,sep="\t")


plotExpression(sce, "NQO1", x = "slingPseudotime_1", 
               colour_by = "Phase", show_violin = FALSE,
               show_smooth = TRUE)+scale_color_tableau()


save(sce,file="scegood.rda")

