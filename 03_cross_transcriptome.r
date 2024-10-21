list.files()

GSE104766<-read.table("GSE104766.tsv",h=T,sep="\t")

GSE131329<-read.table("GSE131329.tsv",h=T,sep="\t")


head(GSE104766)

library(dplyr)

GSE104766%>%select(hs.gene,logFC,P.Value,regulation)%>%distinct()->GSE104766
colnames(GSE104766) <- paste0("GSE104766_", colnames(GSE104766))

GSE131329%>%select(hs.gene,logFC,P.Value,regulation)%>%distinct()->GSE131329
colnames(GSE131329) <- paste0("GSE131329_", colnames(GSE131329))

GSE104766%>%inner_join(GSE131329,by=c("GSE104766_hs.gene"="GSE131329_hs.gene"))->cross 

cross$group<-paste(cross$GSE104766_regulation,cross$GSE104766_regulation,sep=".")
library(ggplot2)
library(ggbeeswarm)
ggplot(cross,aes(GSE104766_logFC,GSE131329_logFC,color=group))+geom_point(aes(fill=group))+
scale_fill_brewer(palette="Pastel1")



library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
ggplot(cross,aes(GSE104766_logFC,GSE131329_logFC))+ scale_fill_brewer(palette="Pastel1")+
  geom_point(aes(fill=factor(group),size=1),shape = 21, alpha = 1, position = position_dodge2(width = 0.1))+
	geom_text_repel(data=cross, aes(label=GSE104766_hs.gene),size=6)+
  theme_classic(base_size=24) +
  theme(legend.position = "none")+ggtitle("")+
geom_vline(xintercept = 0, color = "navy", linetype = "dashed", size = 1)+
geom_hline(yintercept = 0, color = "navy", linetype = "dashed", size = 1)

library(readr)
write_tsv(cross,file="crossresults.tsv")

cross<-read_tsv("crossresults.tsv")
library(dplyr)

cross$fc= (cross$GSE104766_logFC+cross$GSE131329_logFC)/2

library(clusterProfiler)
library(org.Hs.eg.db)


ego <- enrichGO(gene= cross$GSE104766_hs.gene, OrgDb= org.Hs.eg.db,keyType= 'SYMBOL',ont= "BP",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
                
                
           
                


goplot(ego)

dotplot(ego, showCategory=10)
library(DOSE)
#BiocManager::install("enrichplot")
library(enrichplot)

genelist<-split(cross$fc,cross$GSE104766_hs.gene)

cnetplot(ego,foldChange=genelist,colorEdge=TRUE,showCategory=10)
,showCategory=10)

ego2 <- pairwise_termsim(ego)
treeplot(ego2)

upsetplot(ego)
















