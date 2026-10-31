


#### Needed functions and packages 
####
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(clusterProfiler)
require(ggplot2)
require(ggseqlogo)
library(igraph)
library(qgraph)

###calculate jaccerd index
jaccard.IB<-function(set1,set2){
  jac=sum(set1%in%set2)/sum(!duplicated(c(set1,set2)))
  return(jac)
}

# Function to count edges between two vertex sets
edges_between_groups <- function(g, A, B) {
  edges <- ends(g, E(g))
  sum(
    (edges[,1] %in% A & edges[,2] %in% B) |
      (edges[,2] %in% A & edges[,1] %in% B)
  )
}





####
#### Loading required tables
####

tab_4mir=readRDS("output/tab_4mir.rds")



####
#### Figure S4A (heatmap and hierarchical clustering)
####

###select columns for heatmap
temp_fig=tab_4mir[,grep("GUY11|DCDS1",colnames(tab_4mir))]
temp_fig=apply(temp_fig,2,as.numeric)
rownames(temp_fig)=tab_4mir[,1]
temp_fig_moc=temp_fig[,grep("moc",colnames(temp_fig))]
temp_fig_moc=na.omit(temp_fig_moc)

temp_fig=tab_4mir[,grep("GUY11|DCDS1",colnames(tab_4mir))]
temp_fig=apply(temp_fig,2,as.numeric)
rownames(temp_fig)=tab_4mir[,1]
temp_fig_hu=temp_fig[,grep("HU",colnames(temp_fig))]
temp_fig_hu=na.omit(temp_fig_hu)


###keep diff expressed 

reg=rep(0,nrow(temp_fig_moc))

for (i in 1:length(reg)){
  
  x=temp_fig_moc[i,grep("diff",colnames(temp_fig_moc))]
  y=temp_fig_moc[i,grep("padj",colnames(temp_fig_moc))]
  reg[i]=sum(abs(x)>=1 & y<=0.05)
  
}

temp_fig_moc=temp_fig_moc[reg>0,grep("diff",colnames(temp_fig_moc))]

#######
reg=rep(0,nrow(temp_fig_hu))

for (i in 1:length(reg)){
  
  x=temp_fig_hu[i,grep("diff",colnames(temp_fig_hu))]
  y=temp_fig_hu[i,grep("padj",colnames(temp_fig_hu))]
  reg[i]=sum(abs(x)>=1 & y<=0.05)
  
}

temp_fig_hu=temp_fig_hu[reg>0,grep("diff",colnames(temp_fig_hu))]

###Scale matrix
scaled_heat_moc=t(scale(t(temp_fig_moc)))
scaled_heat_hu=t(scale(t(temp_fig_hu)))

###calculate distances
d_moc <- dist(scaled_heat_moc)
hc_moc <- hclust(d_moc, method = "ward.D2")
clusters_moc <- cutree(hc_moc, k = 8)

d_hu <- dist(scaled_heat_hu)
hc_hu <- hclust(d_hu, method = "ward.D2")
clusters_hu <- cutree(hc_hu, k = 8)

#Color for the clusters
annotation_row_moc <- data.frame(Cluster = factor(clusters_moc))
rownames(annotation_row_moc) <- rownames(scaled_heat_moc)

ann_colors_moc <- list(
  Cluster = setNames(RColorBrewer::brewer.pal(8, "Pastel1"), levels(annotation_row_moc$Cluster))
)

annotation_row_hu <- data.frame(Cluster = factor(clusters_hu))
rownames(annotation_row_hu) <- rownames(scaled_heat_hu)

ann_colors_hu <- list(
  Cluster = setNames(RColorBrewer::brewer.pal(8, "Pastel1"), levels(annotation_row_hu$Cluster))
)

###plot heatmaps
pdf("output/FIGS4A_heatmap_MOC.pdf",height = 10,width = 10)
pheatmap(scaled_heat_moc,
         cluster_rows = hc_moc,
         cluster_cols =F,         
         annotation_row = annotation_row_moc,
         show_rownames = FALSE,      
         show_colnames = TRUE,  
         col=viridis(1000),
         main = "Hierarchical clusters MOC (k = 8)",
         cutree_rows = 8,
         annotation_colors = ann_colors_moc,cellwidth = 20)
dev.off()

pdf("output/FIGS4A_heatmap_HU.pdf",height = 10,width = 10)
pheatmap(scaled_heat_hu,
         cluster_rows = hc_hu,
         cluster_cols =F,         
         annotation_row = annotation_row_hu,
         show_rownames = FALSE,      
         show_colnames = TRUE,  
         col=viridis(1000),
         main = "Hierarchical clusters HU (k = 8)",
         cutree_rows = 8,
         annotation_colors = ann_colors_hu,cellwidth = 20)
dev.off()

####
#### Figure 2B Gene Ontology enrichment using hierarchichal clustering results
####

### Calculating enrichment and writing tables

####HU

all_prot=unique(gsub("T0|T1","",tab_4mir[,"proteinID"]))

cluster_tab=cbind(names(clusters_hu),"",clusters_hu)
colnames(cluster_tab)=c("site","prot","clust")
rownames(cluster_tab)=NULL
for (i in 1:nrow(cluster_tab)){cluster_tab[i,"prot"]=unlist(strsplit(cluster_tab[i,"site"],"-"))[1]}

fore=cluster_tab[,c("prot","clust")]
fore[,"prot"]=gsub("T0|T1","",fore[,"prot"])
fore=fore[!duplicated(fore),]

annotations <- read.delim("input/mart_export_GO_MG8.txt", sep="\t", header=TRUE)

####GOBP
##generate dataframe for mapping and naming
gene2go <- annotations[annotations[,"GO.domain"]=="biological_process", c("GO.term.accession", "Gene.stable.ID")]
go2name <- annotations[annotations[,"GO.domain"]=="biological_process", c("GO.term.accession", "GO.term.name")]

##run enrichment
ego <- enricher(
  gene = fore[fore[,"clust"]==1,"prot"],     
  TERM2GENE = gene2go,
  TERM2NAME = go2name,
  pvalueCutoff = 0.05)

gobp=as.matrix(cbind("1",ego@result))

for (i in 2:8){
  
  ego <- enricher(
    gene = fore[fore[,"clust"]==i,"prot"],     
    TERM2GENE = gene2go,
    TERM2NAME = go2name,
    pvalueCutoff = 0.05)
  
  gobp=rbind(gobp,as.matrix(cbind(i,ego@result)))
  
}
colnames(gobp)[1]="cluster"

####GOMF

##generate dataframe for mapping and naming
gene2go <- annotations[annotations[,"GO.domain"]=="molecular_function", c("GO.term.accession", "Gene.stable.ID")]
go2name <- annotations[annotations[,"GO.domain"]=="molecular_function", c("GO.term.accession", "GO.term.name")]

##run enrichment
ego <- enricher(
  gene = fore[fore[,"clust"]==1,"prot"],    
  TERM2GENE = gene2go,
  TERM2NAME = go2name,
  pvalueCutoff = 0.05)

gomf=as.matrix(cbind("1",ego@result))

for (i in 2:8){
  
  ego <- enricher(
    gene = fore[fore[,"clust"]==i,"prot"],     
    TERM2GENE = gene2go,
    TERM2NAME = go2name,
    pvalueCutoff = 0.05)
  
  gomf=rbind(gomf,as.matrix(cbind(i,ego@result)))
  
}
colnames(gomf)[1]="cluster"

####GOCC
##generate dataframe for mapping and naming
gene2go <- annotations[annotations[,"GO.domain"]=="cellular_component", c("GO.term.accession", "Gene.stable.ID")]
go2name <- annotations[annotations[,"GO.domain"]=="cellular_component", c("GO.term.accession", "GO.term.name")]

##run enrichment
ego <- enricher(
  gene = fore[fore[,"clust"]==1,"prot"],    
  TERM2GENE = gene2go,
  TERM2NAME = go2name,
  pvalueCutoff = 0.05)

gocc=as.matrix(cbind("1",ego@result))

for (i in 2:8){
  
  ego <- enricher(
    gene = fore[fore[,"clust"]==i,"prot"],    
    TERM2GENE = gene2go,
    TERM2NAME = go2name,
    pvalueCutoff = 0.05)
  
  gocc=rbind(gocc,as.matrix(cbind(i,ego@result)))
  
}
colnames(gocc)[1]="cluster"

gobp=gobp[as.numeric(gobp[,"pvalue"])<=0.05,]
gomf=gomf[as.numeric(gomf[,"pvalue"])<=0.05,]
gocc=gocc[as.numeric(gocc[,"pvalue"])<=0.05,]

enrich_HU=rbind(cbind("GOBP",gobp),
                cbind("GOMF",gomf),
                cbind("GOCC",gocc))
colnames(enrich_HU)[1]="GeneOntology"
rownames(enrich_HU)=NULL

write.table(cluster_tab,"output/hierar_tab_HU.txt",quote = F, row.names=F,sep="\t")
write.table(enrich_HU,"output/GO_enrich_HU.txt",quote = F, row.names=F,sep="\t")
saveRDS(enrich_HU,"output/GO_enrich_HU.rds")

####not HU

all_prot=unique(gsub("T0|T1","",tab_4mir[,"proteinID"]))

cluster_tab=cbind(names(clusters_moc),"",clusters_moc)
colnames(cluster_tab)=c("site","prot","clust")
rownames(cluster_tab)=NULL
for (i in 1:nrow(cluster_tab)){cluster_tab[i,"prot"]=unlist(strsplit(cluster_tab[i,"site"],"-"))[1]}

fore=cluster_tab[,c("prot","clust")]
fore[,"prot"]=gsub("T0|T1","",fore[,"prot"])
fore=fore[!duplicated(fore),]


####GOBP
##generate dataframe for mapping and naming
gene2go <- annotations[annotations[,"GO.domain"]=="biological_process", c("GO.term.accession", "Gene.stable.ID")]
go2name <- annotations[annotations[,"GO.domain"]=="biological_process", c("GO.term.accession", "GO.term.name")]

##run enrichment
ego <- enricher(
  gene = fore[fore[,"clust"]==1,"prot"],     # your DE genes, e.g. c("MGG_00123", "MGG_04567", ...)
  TERM2GENE = gene2go,
  TERM2NAME = go2name,
  pvalueCutoff = 0.05)

gobp=as.matrix(cbind("1",ego@result))

for (i in 2:8){
  
  ego <- enricher(
    gene = fore[fore[,"clust"]==i,"prot"],     # your DE genes, e.g. c("MGG_00123", "MGG_04567", ...)
    TERM2GENE = gene2go,
    TERM2NAME = go2name,
    pvalueCutoff = 0.05)
  
  gobp=rbind(gobp,as.matrix(cbind(i,ego@result)))
  
}
colnames(gobp)[1]="cluster"

####GOMF

##generate dataframe for mapping and naming
gene2go <- annotations[annotations[,"GO.domain"]=="molecular_function", c("GO.term.accession", "Gene.stable.ID")]
go2name <- annotations[annotations[,"GO.domain"]=="molecular_function", c("GO.term.accession", "GO.term.name")]

##run enrichment
ego <- enricher(
  gene = fore[fore[,"clust"]==1,"prot"],     # your DE genes, e.g. c("MGG_00123", "MGG_04567", ...)
  TERM2GENE = gene2go,
  TERM2NAME = go2name,
  pvalueCutoff = 0.05)

gomf=as.matrix(cbind("1",ego@result))

for (i in 2:8){
  
  ego <- enricher(
    gene = fore[fore[,"clust"]==i,"prot"],     # your DE genes, e.g. c("MGG_00123", "MGG_04567", ...)
    TERM2GENE = gene2go,
    TERM2NAME = go2name,
    pvalueCutoff = 0.05)
  
  gomf=rbind(gomf,as.matrix(cbind(i,ego@result)))
  
}
colnames(gomf)[1]="cluster"

####GOCC
##generate dataframe for mapping and naming
gene2go <- annotations[annotations[,"GO.domain"]=="cellular_component", c("GO.term.accession", "Gene.stable.ID")]
go2name <- annotations[annotations[,"GO.domain"]=="cellular_component", c("GO.term.accession", "GO.term.name")]

##run enrichment
ego <- enricher(
  gene = fore[fore[,"clust"]==1,"prot"],     # your DE genes, e.g. c("MGG_00123", "MGG_04567", ...)
  TERM2GENE = gene2go,
  TERM2NAME = go2name,
  pvalueCutoff = 0.05)

gocc=as.matrix(cbind("1",ego@result))

for (i in 2:8){
  
  ego <- enricher(
    gene = fore[fore[,"clust"]==i,"prot"],     # your DE genes, e.g. c("MGG_00123", "MGG_04567", ...)
    TERM2GENE = gene2go,
    TERM2NAME = go2name,
    pvalueCutoff = 0.05)
  
  gocc=rbind(gocc,as.matrix(cbind(i,ego@result)))
  
}
colnames(gocc)[1]="cluster"

gobp=gobp[as.numeric(gobp[,"pvalue"])<=0.05,]
gomf=gomf[as.numeric(gomf[,"pvalue"])<=0.05,]
gocc=gocc[as.numeric(gocc[,"pvalue"])<=0.05,]

enrich_MOC=rbind(cbind("GOBP",gobp),
                 cbind("GOMF",gomf),
                 cbind("GOCC",gocc))
colnames(enrich_MOC)[1]="GeneOntology"
rownames(enrich_MOC)=NULL

write.table(cluster_tab,"output/hierar_tab_MOC.txt",quote = F, row.names=F,sep="\t")
write.table(enrich_MOC,"output/GO_enrich_MOC.txt",quote = F, row.names=F,sep="\t")
saveRDS(enrich_MOC,"output/GO_enrich_MOC.rds")

############
#### heatmaps for the enrichments

temp_GO=enrich_HU[as.numeric(enrich_HU[,"qvalue"])<=0.05,]

####GOBP
heat=matrix(0,length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","Description"])),length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","cluster"])))
colnames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","cluster"])
rownames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","Description"])

for (i in 1:ncol(heat)){
  
  temp=as.numeric(temp_GO[temp_GO[,"cluster"]==colnames(heat)[i]  ,"FoldEnrichment"])
  names(temp)=temp_GO[temp_GO[,"cluster"]==colnames(heat)[i],"Description"]
  temp=temp[names(temp)%in%rownames(heat)]
  heat[names(temp),i]=temp
}
pdf("output/FIGS4B_GO_enrich_HU_GOBP.pdf",height = 6,width = 6)
pheatmap(heat,
         cluster_cols =F,  
         col=c("white",plasma(1000)),
         main = "GOBP enrichment_HU",cellwidth = 6, cellheight = 6, fontsize = 6,border_color = "black")
dev.off()

####GOMF
heat=matrix(0,length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","Description"])),length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","cluster"])))
colnames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","cluster"])
rownames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","Description"])

for (i in 1:ncol(heat)){
  
  temp=as.numeric(temp_GO[temp_GO[,"cluster"]==colnames(heat)[i]  ,"FoldEnrichment"])
  names(temp)=temp_GO[temp_GO[,"cluster"]==colnames(heat)[i],"Description"]
  temp=temp[names(temp)%in%rownames(heat)]
  heat[names(temp),i]=temp
}
pdf("output/FIGS4B_GO_enrich_HU_GOMF.pdf",height = 6,width = 6)
pheatmap(heat,
         cluster_cols =F,  
         col=c("white",plasma(1000)),
         main = "GOMF enrichment_HU",cellwidth = 6, cellheight = 6, fontsize = 6,border_color = "black")
dev.off()


####GOCC
heat=matrix(0,length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","Description"])),
            length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","cluster"])))
colnames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","cluster"])
rownames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","Description"])

for (i in 1:ncol(heat)){
  
  temp=as.numeric(temp_GO[temp_GO[,"cluster"]==colnames(heat)[i]  ,"FoldEnrichment"])
  names(temp)=temp_GO[temp_GO[,"cluster"]==colnames(heat)[i],"Description"]
  temp=temp[names(temp)%in%rownames(heat)]
  heat[names(temp),i]=temp
}
pdf("output/FIGS4B_GO_enrich_HU_GOCC.pdf",height = 6,width = 6)
pheatmap(heat,
         cluster_cols =F,  
         col=c("white",plasma(1000)),
         main = "GOCC enrichment_HU",cellwidth = 6, cellheight = 6, fontsize = 6,border_color = "black")
dev.off()

#####
#####
#####

temp_GO=enrich_MOC[as.numeric(enrich_MOC[,"qvalue"])<=0.05,]

####GOBP
heat=matrix(0,length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","Description"])),length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","cluster"])))
colnames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","cluster"])
rownames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOBP","Description"])

for (i in 1:ncol(heat)){
  
  temp=as.numeric(temp_GO[temp_GO[,"cluster"]==colnames(heat)[i]  ,"FoldEnrichment"])
  names(temp)=temp_GO[temp_GO[,"cluster"]==colnames(heat)[i],"Description"]
  temp=temp[names(temp)%in%rownames(heat)]
  heat[names(temp),i]=temp
}
pdf("output/FIGS4B_GOBP_enrich_MOC_GOBP.pdf",height = 6,width = 6)
pheatmap(heat,
         cluster_cols =F,  
         col=c("white",plasma(1000)),
         main = "GOBP enrichment_MOC",cellwidth = 6, cellheight = 6, fontsize = 6,border_color = "black")
dev.off()

####GOMF
heat=matrix(0,length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","Description"])),length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","cluster"])))
colnames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","cluster"])
rownames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOMF","Description"])

for (i in 1:ncol(heat)){
  
  temp=as.numeric(temp_GO[temp_GO[,"cluster"]==colnames(heat)[i]  ,"FoldEnrichment"])
  names(temp)=temp_GO[temp_GO[,"cluster"]==colnames(heat)[i],"Description"]
  temp=temp[names(temp)%in%rownames(heat)]
  heat[names(temp),i]=temp
}
pdf("output/FIGS4B_GO_enrich_MOC_GOMF.pdf",height = 6,width = 6)
pheatmap(heat,
         cluster_cols =F,  
         col=c("white",plasma(1000)),
         main = "GOMF enrichment_MOC",cellwidth = 6, cellheight = 6, fontsize = 6,border_color = "black")
dev.off()


####GOCC
heat=matrix(0,length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","Description"])),length(unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","cluster"])))
colnames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","cluster"])
rownames(heat)=unique(temp_GO[temp_GO[,"GeneOntology"]=="GOCC","Description"])

for (i in 1:ncol(heat)){
  
  temp=as.numeric(temp_GO[temp_GO[,"cluster"]==colnames(heat)[i]  ,"FoldEnrichment"])
  names(temp)=temp_GO[temp_GO[,"cluster"]==colnames(heat)[i],"Description"]
  temp=temp[names(temp)%in%rownames(heat)]
  heat[names(temp),i]=temp
}
pdf("output/FIGS4B_GO_enrich_MOC_GOCC.pdf",height = 6,width = 6)
pheatmap(heat,
         cluster_cols =F,  
         col=c("white",plasma(1000)),
         main = "GOCC enrichment_MOC",cellwidth = 6, cellheight = 6, fontsize = 6,border_color = "black")
dev.off()

#### JACCARD and PERMANOVA (Figure supplementary 4)


heat_jacc=matrix(0,8,8)
colnames(heat_jacc)=1:8
rownames(heat_jacc)=1:8

for (i in 1:nrow(heat_jacc)){
  
  A=names(clusters_hu)[clusters_hu==i]
  
  for (j in 1:ncol(heat_jacc)){
    
    B=names(clusters_moc)[clusters_moc==j]
    
    heat_jacc[i,j]=jaccard.IB(A,B)
    
    
  }
  
}

colnames(heat_jacc)=paste("MOC_",1:8,sep="")
rownames(heat_jacc)=paste("HU_",1:8,sep="")

A1=names(clusters_moc)[clusters_moc==1]
B1=names(clusters_hu)[clusters_hu==1]

A2=names(clusters_moc)[clusters_moc==1]
B2=names(clusters_hu)[clusters_hu==3]

A3=names(clusters_moc)[clusters_moc==3]
B3=names(clusters_hu)[clusters_hu==1]

A4=names(clusters_moc)[clusters_moc==5]
B4=names(clusters_hu)[clusters_hu==4]

A5=names(clusters_moc)[clusters_moc==7]
B5=names(clusters_hu)[clusters_hu==4]

A6=names(clusters_moc)[clusters_moc==8]
B6=names(clusters_hu)[clusters_hu==2]



mat=as.data.frame(cbind(scaled_heat_moc[A1[A1%in%B1],],scaled_heat_hu[B1[B1%in%A1],]))


sample_info <- data.frame(
  Sample = colnames(mat),
  Genotype = c("GUY11","GUY11","GUY11","DCDS1","DCDS1","DCDS1",
               "GUY11","GUY11","GUY11","DCDS1","DCDS1","DCDS1"),
  Treatment = c("moc","moc","moc","moc","moc","moc",
                "HU","HU","HU","HU","HU","HU"),
  Time = c("1.5h","2h","4h",
           "1.5h","2h","4h",
           "1.5h","2h","4h",
           "1.5h","2h","4h")
)


library(vegan)

#####PERMANOVA


perm1=adonis2(
  t(as.data.frame(cbind(scaled_heat_moc[A1[A1%in%B1],],scaled_heat_hu[B1[B1%in%A1],]))) ~ Genotype + Treatment + Time,
  data = sample_info,
  method = "euclidean",
  by = "margin"
)

perm2=adonis2(
  t(as.data.frame(cbind(scaled_heat_moc[A2[A2%in%B2],],scaled_heat_hu[B2[B2%in%A2],]))) ~ Genotype + Treatment + Time,
  data = sample_info,
  method = "euclidean",
  by = "margin"
)

perm3=adonis2(
  t(as.data.frame(cbind(scaled_heat_moc[A3[A3%in%B3],],scaled_heat_hu[B3[B3%in%A3],]))) ~ Genotype + Treatment + Time,
  data = sample_info,
  method = "euclidean",
  by = "margin"
)

perm4=adonis2(
  t(as.data.frame(cbind(scaled_heat_moc[A4[A4%in%B4],],scaled_heat_hu[B4[B4%in%A4],]))) ~ Genotype + Treatment + Time,
  data = sample_info,
  method = "euclidean",
  by = "margin"
)

perm5=adonis2(
  t(as.data.frame(cbind(scaled_heat_moc[A5[A5%in%B5],],scaled_heat_hu[B5[B5%in%A5],]))) ~ Genotype + Treatment + Time,
  data = sample_info,
  method = "euclidean",
  by = "margin"
)

perm6=adonis2(
  t(as.data.frame(cbind(scaled_heat_moc[A6[A6%in%B6],],scaled_heat_hu[B6[B6%in%A6],]))) ~ Genotype + Treatment + Time,
  data = sample_info,
  method = "euclidean",
  by = "margin"
)


pca1=prcomp(t(as.data.frame(cbind(scaled_heat_moc[A1[A1%in%B1],],scaled_heat_hu[B1[B1%in%A1],]))), scale. = TRUE)
pca2=prcomp(t(as.data.frame(cbind(scaled_heat_moc[A2[A2%in%B2],],scaled_heat_hu[B2[B2%in%A2],]))), scale. = TRUE)
pca3=prcomp(t(as.data.frame(cbind(scaled_heat_moc[A3[A3%in%B3],],scaled_heat_hu[B3[B3%in%A3],]))), scale. = TRUE)
pca4=prcomp(t(as.data.frame(cbind(scaled_heat_moc[A4[A4%in%B4],],scaled_heat_hu[B4[B4%in%A4],]))), scale. = TRUE)
pca5=prcomp(t(as.data.frame(cbind(scaled_heat_moc[A5[A5%in%B5],],scaled_heat_hu[B5[B5%in%A5],]))), scale. = TRUE)
pca6=prcomp(t(as.data.frame(cbind(scaled_heat_moc[A6[A6%in%B6],],scaled_heat_hu[B6[B6%in%A6],]))), scale. = TRUE)


###Representacion

pch_vect=c(17,17,17,16,16,16,17,17,17,16,16,16)
col_vect=c(rep("azure3",6),rep("gold",6))




pdf("output/FIGS4C_PERMANOVA_PCA.pdf",height = 6,width = 8)

par(mfrow=c(2,3))

plot(pca1$x[,1],pca1$x[,2],
     xlab = paste0("PC1 (", round(100*summary(pca1)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca1)$importance[2,2],1), "%)"),
     pch = pch_vect,col = col_vect,cex=1.5,main="G1")

plot(pca2$x[,1],pca2$x[,2],
     xlab = paste0("PC1 (", round(100*summary(pca2)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca2)$importance[2,2],1), "%)"),
     pch = pch_vect,col = col_vect,cex=1.5,main="G2")

plot(pca3$x[,1],pca3$x[,2],
     xlab = paste0("PC1 (", round(100*summary(pca3)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca3)$importance[2,2],1), "%)"),
     pch = pch_vect,col = col_vect,cex=1.5,main="G3")

plot(pca4$x[,1],pca4$x[,2],
     xlab = paste0("PC1 (", round(100*summary(pca4)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca4)$importance[2,2],1), "%)"),
     pch = pch_vect,col = col_vect,cex=1.5,main="G4")

plot(pca5$x[,1],pca5$x[,2],
     xlab = paste0("PC1 (", round(100*summary(pca5)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca5)$importance[2,2],1), "%)"),
     pch = pch_vect,col = col_vect,cex=1.5,main="G5")

plot(pca6$x[,1],pca6$x[,2],
     xlab = paste0("PC1 (", round(100*summary(pca6)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca6)$importance[2,2],1), "%)"),
     pch = pch_vect,col = col_vect,cex=1.5,main="G6")

dev.off()









