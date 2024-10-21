library(GEOquery)
# download GEO GSE131329
gse <- getGEO("GSE131329", GSEMatrix = TRUE)

# dataset access
set <- gse[[1]]

# expression data and features annotations
data<-exprs(set)
features<-fData(set)
# phenotypes
annot <- pData(set)
tail(features,n=100)


all<-merge(features,data,by="row.names")

write.csv(all,file="matrix.csv",row.names=F)
write.csv(annot,file="pheno.csv",row.names=T)




## ferroptosis analysis

list.files()

data<-read.csv("matrix.csv",h=T)

data<-data[complete.cases(data$gene),]
library(transpipe14)
ok<-filtermatrix(data)
pheno<-read.csv("phenoshort.csv",row.names=1)

all(row.names(pheno)==colnames(ok))

library(ferroviz)

res<-deg(ok,pheno$tissue,control="noncancerous_liver_tissue")

sig<-filtresig(res)

vollimma(res,nb=500,fc=1,p=1.668178e-02,size=4,alpha=1)

hsvol(res,ferrdbhs,color="grey",size=16,x=1,label=1.5)


process<-reducedf(sig,ok,n=1384)

pcatransellipses(process,pheno,group="tissue",pal="Set1",alpha=1,names=F,x=1,y=2,level=0.75)


bestheat(process,pheno,font=10,rownames=F)


barploths(res,ferrdbhs,fc=1,size=16)

df<-lisths(res,ferrdbhs,fc=1)
df

write.table(df,file="ferroptosis.tsv",row.names=F,sep="\t")


## ferroptosis GO annotations

load("human_annotation.rda")



annotation%>%mutate(id=replace(id,id=="tf","transcription_factor"))->annotation
annotation%>%select(hs.gene,id)%>%distinct()%>%inner_join(df,by="hs.gene",relationship="many-to-many")->all


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

chisq.test(all$conca,all$id)


vector<-as.vector(unique(df$hs.gene))
small<-ok[row.names(ok)%in%vector,]

bestheat(small,scale="none",pheno,font=8,rownames=T)

save(vector,file="GSE131329.rda")
pcatransellipses(small,pheno,group="tissue",pal="Set1",alpha=1,names=F,x=1,y=2,level=0.75)



library(vcd)

struct <- structable(~regulation+id, data = all)
mosaic(struct,  direction = "h", pop = FALSE,colorize = T, shade = TRUE,
       labeling_args = list(gp_labels = gpar(fontsize = 10, fontface = "bold"), # Taille et style des labels
                            rot_labels = c(left = 90, top = 90))) # Rotation des étiquettes

chisq.test(struct)
