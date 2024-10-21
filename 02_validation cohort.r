library(GEOquery)
# download GEO GSE104766
gse <- getGEO("GSE104766", GSEMatrix = TRUE)

# dataset access
set <- gse[[1]]



# phenotypes
annot <- pData(set)



write.csv(annot,file="pheno.csv",row.names=T)
## phenotype extraction GEOquery but RNAseq matrix was downloaded on NCBI GEO website
data<-read.csv("matrix.csv",h=T)
list.files()



library(transpipe14)


all(row.names(pheno)==colnames(data))

library(ferroviz)



res<-deg(data,pheno$tissue,control="Normal liver")

sig<-filtresig(res)

vollimma(res,nb=500,fc=1,p=0.02560168,size=4,alpha=1)

hsvol(res,ferrdbhs,color="grey",size=16,x=1,label=1.5)


process<-reducedf(sig,data,n=1933)

pcatransellipses(process,pheno,group="tissue",pal="Set1",alpha=1,names=F,x=1,y=2,level=0.75)

pheno%>%select(tissue)->annot
bestheat(process,annot,font=10,rownames=F)


barploths(res,ferrdbhs,fc=1,size=16)

df<-lisths(res,ferrdbhs,fc=1)
df

write.table(df,file="ferroptosis.tsv",row.names=F,sep="\t")


## GO annotation


load("human_annotation.rda")

df<-lisths(resulths,ferrdbhs,fc=0.125)
df

annotation%>%mutate(id=replace(id,id=="tf","transcription_factor"))->annotation
annotation%>%select(hs.gene,id)%>%distinct()%>%inner_join(df,by="hs.gene",relationship="many-to-many")->all




  
## annotation

all$conca<-paste(all$regulation,all$ferroptosis,sep="_")
all %>%
  count(id, sort = TRUE) %>%  # Count and sort by frequency
  slice_head(n=20)%>%pull(id)->vector

all%>%filter(id%in%vector)->all
library(pals)
ggplot(all,aes(x=factor(regulation),fill=factor(id)))+
    geom_bar(position="fill")+
    geom_text(aes(label=..count..),stat='count',position=position_fill(vjust=0.5),size=6,colour="white")+
    scale_fill_manual(values=cols25())+
    labs(fill = "Ontology")+  
    ggtitle("") +
    xlab("regulation_ferroptosis") + ylab("proportions/numbers")+theme_classic(base_size=16)+ coord_flip()

write.table(all,file="annotationFerro.tsv",row.names=F,sep="\t")

chisq.test(all$regulation,all$id)


vector<-as.vector(unique(df$hs.gene))
small<-data[row.names(data)%in%vector,]

bestheat(small,annot,font=8,rownames=T)

save(vector,file="GSE131329.rda")
pcatransellipses(small,pheno,group="tissue",pal="Set1",alpha=1,names=F,x=1,y=2,level=0.75)



