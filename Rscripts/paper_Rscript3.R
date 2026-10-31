
####The following code depends on rds tables calculated in Rscript1 & 2 (already provided in folder output)

####Here we provide the code for figures XXXXX . This code result in "raw" figures that were modified for aesthetic reasons using Inkscape.
####
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
#### Figure 2B
####

temp_normInt=apply(tab_4mir[,grep("guy11",colnames(tab_4mir))],2,as.numeric)
temp_vect=rep(0,nrow(temp_normInt))
for (i in 1:length(temp_vect)){temp_vect[i]=sum(!is.na(temp_normInt[i,]))}

counts=matrix(0,2,9)
colnames(counts)=c("Total identified","Total analyzed",unique(gsub("rep1|rep3|guy11_","",colnames(temp_normInt))))
rownames(counts)=c("Sites","Proteins")

counts[1,1]=sum(temp_vect>0)
counts[2,1]=length(unique(tab_4mir[temp_vect>0,"proteinID"]))
counts[1,2]=sum(!is.na(as.numeric(tab_4mir[,"GUY11_HU_2h_diff"])))
counts[2,2]=length(unique(tab_4mir[!is.na(as.numeric(tab_4mir[,"GUY11_HU_2h_diff"])),"proteinID"]))

for (i in 3:ncol(counts)){
  
  temp_vect=rep(0,nrow(temp_normInt))
  
  for (j in 1:length(temp_vect)){temp_vect[j]=sum(!is.na(temp_normInt[j,grep(colnames(counts)[i],colnames(temp_normInt))]))}
  
  counts[1,i]=sum(temp_vect>0)
  counts[2,i]=length(unique(tab_4mir[temp_vect>0,"proteinID"]))
}

pdf("output/FIG2B_barplot_counts.pdf",height = 6,width = 8)
x=barplot(counts,beside = T,las=3,ylim=c(0,7000))
text(x[1,],counts[1,],counts[1,])
text(x[2,],counts[2,],counts[2,])
dev.off()

####
#### Figure 2C (table for the pie chart and barplot)
####
count_pep=table(unlist(strsplit(tab_4mir[,"Aa"],";")))
pdf("output/FIG2C_barplot_residues.pdf",height = 6, width = 3)
x=barplot(count_pep)
text(x[,1],count_pep,count_pep)
dev.off()


####
#### Figure 2A (heatmap and hierarchical clustering)
####

###select columns for heatmap
temp_fig=tab_4mir[,grep("GUY11|DCDS1",colnames(tab_4mir))]
temp_fig=apply(temp_fig,2,as.numeric)
rownames(temp_fig)=tab_4mir[,1]
temp_fig_moc=temp_fig[,grep("GUY11",colnames(temp_fig))]
temp_fig_moc=na.omit(temp_fig_moc)

temp_fig=tab_4mir[,grep("GUY11|DCDS1",colnames(tab_4mir))]
temp_fig=apply(temp_fig,2,as.numeric)
rownames(temp_fig)=tab_4mir[,1]
temp_fig_hu=temp_fig[,grep("DCDS1",colnames(temp_fig))]
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

scaled_heat_moc=scaled_heat_moc[,c("GUY11_moc_1h30min_diff","GUY11_moc_2h_diff","GUY11_moc_4h_diff",
                                   "GUY11_HU_1h30min_diff","GUY11_HU_2h_diff","GUY11_HU_4h_diff")]


scaled_heat_hu=scaled_heat_hu[,c("DCDS1_moc_1h30min_diff","DCDS1_moc_2h_diff","DCDS1_moc_4h_diff",
                                 "DCDS1_HU_1h30min_diff","DCDS1_HU_2h_diff","DCDS1_HU_4h_diff")]


###plot heatmaps
pdf("output/FIG3A_heatmap_GUY11.pdf",height = 10,width = 10)
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

pdf("output/FIG3A_heatmap_deltaCDS1.pdf",height = 10,width = 10)
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
#### Figure 3B Gene Ontology enrichment using hierarchichal clustering results
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

write.table(cluster_tab,"output/hierar_tab_deltaCDS1.txt",quote = F, row.names=F,sep="\t")
write.table(enrich_HU,"output/GO_enrich_deltaCDS1.txt",quote = F, row.names=F,sep="\t")
saveRDS(enrich_HU,"output/GO_enrich_deltaCDS1.rds")

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

write.table(cluster_tab,"output/hierar_tab_GUY11.txt",quote = F, row.names=F,sep="\t")
write.table(enrich_MOC,"output/GO_enrich_GUY11.txt",quote = F, row.names=F,sep="\t")
saveRDS(enrich_MOC,"output/GO_enrich_GUY11.rds")

#########
######### barplots figure 3B
#########

enrich_guy=readRDS("output/GO_enrich_GUY11.rds")
enrich_cds=readRDS("output/GO_enrich_deltaCDS1.rds")

##selecting columns and rows (significant hits)

enrich_guy=enrich_guy[as.numeric(enrich_guy[,"qvalue"])<=0.05,c("GeneOntology","cluster","Description","qvalue","FoldEnrichment")]
enrich_cds=enrich_cds[as.numeric(enrich_cds[,"qvalue"])<=0.05,c("GeneOntology","cluster","Description","qvalue","FoldEnrichment")]

enrich_guy=enrich_guy[order(as.numeric(enrich_guy[,"cluster"]),
                            enrich_guy[,"GeneOntology"],
                            -1*as.numeric(enrich_guy[,"FoldEnrichment"])),]

enrich_cds=enrich_cds[order(as.numeric(enrich_cds[,"cluster"]),
                            enrich_cds[,"GeneOntology"],
                            -1*as.numeric(enrich_cds[,"FoldEnrichment"])),]
col_vect=RColorBrewer::brewer.pal(8,"Pastel1")
col_vect_GO=c("coral1","cornflowerblue","darkolivegreen3")
names(col_vect_GO)=c("GOBP","GOMF","GOCC")

clust_index=table(enrich_guy[,"cluster"])
target=max(clust_index)+1

for(i in 1:length(clust_index)){
  
  path_save=paste("output/FIG3B_barplot_GO_cluster_",
                  names(clust_index)[i],"_GUY11.pdf",sep="")
  
  sel=enrich_guy[,"cluster"]==names(clust_index)[i]
  
  pdf(path_save,height = 5,width = 10)
  
  names_plot=c(paste(enrich_guy[sel,"Description"],
                     enrich_guy[sel,"GeneOntology"],sep="_"),rep("",target-sum(sel)))
  number_plot=c(as.numeric(enrich_guy[sel,"FoldEnrichment"]),rep(0,target-sum(sel)))
  
  barplot(number_plot,
          border=NA,
          col=c(col_vect[as.numeric(enrich_guy[sel,"cluster"])],rep("grey",target-sum(sel))),
          ylim=c(0,max(as.numeric(enrich_guy[sel,"FoldEnrichment"]))+2),las=2,
          names=names_plot,
          ylab=paste("FoldEnrichment (Guy11)","cluster",names(clust_index)[i]))
  dev.off()
  
  
}


clust_index=table(enrich_cds[,"cluster"])
target=max(clust_index)+1

for(i in 1:length(clust_index)){
  
  path_save=paste("output/FIG3B_barplot_GO_cluster_",
                  names(clust_index)[i],"_deltaCDS1.pdf",sep="")
  
  sel=enrich_cds[,"cluster"]==names(clust_index)[i]
  
  pdf(path_save,height = 5,width = 10)
  
  names_plot=c(paste(enrich_cds[sel,"Description"],
                     enrich_cds[sel,"GeneOntology"],sep="_"),rep("",target-sum(sel)))
  number_plot=c(as.numeric(enrich_cds[sel,"FoldEnrichment"]),rep(0,target-sum(sel)))
  
  barplot(number_plot,
          border=NA,
          col=c(col_vect[as.numeric(enrich_cds[sel,"cluster"])],rep("grey",target-sum(sel))),
          ylim=c(0,max(as.numeric(enrich_cds[sel,"FoldEnrichment"]))+2),las=2,
          names=names_plot,
          ylab=paste("FoldEnrichment (deltaCDS1)","cluster",names(clust_index)[i]))
  dev.off()
  
}


####
#### Figure 3C JAccard index to measure site overlap between the two series for clustsered phosphosites
####

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

###Selecting pairs

jacc_v2=cbind(rep(1:8,8),rep(1:8,each=8),0)
colnames(jacc_v2)=c("moc","hu","jacc")
for (i in 1:nrow(jacc_v2)){
  
  A=names(clusters_hu)[clusters_hu==jacc_v2[i,"hu"]]
  B=names(clusters_moc)[clusters_moc==jacc_v2[i,"moc"]]
  
  jacc_v2[i,"jacc"]=jaccard.IB(A,B)
}

jacc_v2=jacc_v2[order(jacc_v2[,3],decreasing=T),]


pdf("output/FIG3C_jaccard_GUY11_deltaCDS1.pdf",height = 6,width = 6)
pheatmap(heat_jacc,  
         col=magma(1000),
         main = "jaccard_peptide",
         cellwidth = 7, cellheight = 7, fontsize = 7)
dev.off()


####Preparing table with both hieararchical clustering index together

jacc_mat=matrix("",length(unique(c(names(clusters_hu),names(clusters_moc)))),4)
colnames(jacc_mat)=c("Position","hier_clust_MOC","hier_clust_HU","jaccard_selection")
rownames(jacc_mat)=unique(c(names(clusters_hu),names(clusters_moc)))
jacc_mat[,"Position"]=unique(c(names(clusters_hu),names(clusters_moc)))
jacc_mat[names(clusters_moc),"hier_clust_MOC"]=clusters_moc
jacc_mat[names(clusters_hu),"hier_clust_HU"]=clusters_hu
jacc_mat[jacc_mat[,"hier_clust_MOC"]=="4" & 
           jacc_mat[,"hier_clust_HU"]=="3","jaccard_selection"]="4_3"
jacc_mat[jacc_mat[,"hier_clust_MOC"]=="1" & 
           jacc_mat[,"hier_clust_HU"]=="2","jaccard_selection"]="1_2"
jacc_mat[jacc_mat[,"hier_clust_MOC"]=="5" & 
           jacc_mat[,"hier_clust_HU"]=="2","jaccard_selection"]="5_2"
jacc_mat[jacc_mat[,"hier_clust_MOC"]=="8" & 
           jacc_mat[,"hier_clust_HU"]=="1","jaccard_selection"]="8_1"
jacc_mat[jacc_mat[,"hier_clust_MOC"]=="2" & 
           jacc_mat[,"hier_clust_HU"]=="3","jaccard_selection"]="2_3"
jacc_mat[jacc_mat[,"hier_clust_MOC"]=="1" & 
           jacc_mat[,"hier_clust_HU"]=="4","jaccard_selection"]="1_4"

colnames(jacc_mat)=c("Position","hier_clust_GUY11","hier_clust_deltaCDS1","jaccard_selection")


write.table(jacc_mat,"output/jacc_mat.txt",quote = F, row.names=F,sep="\t")

###Integrating to main tab

rownames(tab_4mir)=tab_4mir[,"IDpos"]

tab_4mir=cbind(tab_4mir,"","","")
colnames(tab_4mir)[colnames(tab_4mir)==""]=c("hier_clust_GUY11","hier_clust_deltaCDS1","jaccard_selection")
tab_4mir[jacc_mat[,"Position"],c("hier_clust_GUY11","hier_clust_deltaCDS1","jaccard_selection")]=jacc_mat[,c("hier_clust_GUY11","hier_clust_deltaCDS1","jaccard_selection")]

write.table(tab_4mir,"output/mainTab_mod.txt",quote = F, row.names=F,sep="\t")


####
#### Figure 3D Boxplot showing scaled abundances for top jaccard scoring cluster pairs
####

A1=names(clusters_moc)[clusters_moc==4]
B1=names(clusters_hu)[clusters_hu==3]

A2=names(clusters_moc)[clusters_moc==1]
B2=names(clusters_hu)[clusters_hu==2]

A3=names(clusters_moc)[clusters_moc==5]
B3=names(clusters_hu)[clusters_hu==2]

A4=names(clusters_moc)[clusters_moc==8]
B4=names(clusters_hu)[clusters_hu==1]

A5=names(clusters_moc)[clusters_moc==2]
B5=names(clusters_hu)[clusters_hu==3]

A6=names(clusters_moc)[clusters_moc==1]
B6=names(clusters_hu)[clusters_hu==4]

pdf("output/FIG3D_top6_jaccard_guy11_deltaCDS1.pdf",height = 5,width = 8)

par(mfrow=c(2,3))
boxplot(cbind(scaled_heat_moc[A1[A1%in%B1],],scaled_heat_hu[B1[B1%in%A1],]),
        col=c(rep(rgb(0.7,0.2,0.2),6),
              rep(rgb(0.2,0.5,0.8),6)),las=3,ylim=c(-2,2),main=paste("n=",sum(A1%in%B1), "GUY11_4 deltaCDS1_3"))
abline(h=0,lty=2,col=rgb(0.2,0.4,0.6))
boxplot(cbind(scaled_heat_moc[A2[A2%in%B2],],scaled_heat_hu[B2[B2%in%A2],]),
        col=c(rep(rgb(0.7,0.2,0.2),6),
              rep(rgb(0.2,0.5,0.8),6)),las=3,ylim=c(-2,2),main=paste("n=",sum(A2%in%B2), "GUY11_1 deltaCDS1_2"))
abline(h=0,lty=2,col=rgb(0.2,0.4,0.6))
boxplot(cbind(scaled_heat_moc[A3[A3%in%B3],],scaled_heat_hu[B3[B3%in%A3],]),
        col=c(rep(rgb(0.7,0.2,0.2),6),
              rep(rgb(0.2,0.5,0.8),6)),las=3,ylim=c(-2,2),main=paste("n=",sum(A3%in%B3), "GUY11_5 deltaCDS1_2"))
abline(h=0,lty=2,col=rgb(0.2,0.4,0.6))
boxplot(cbind(scaled_heat_moc[A4[A4%in%B4],],scaled_heat_hu[B4[B4%in%A4],]),
        col=c(rep(rgb(0.7,0.2,0.2),6),
              rep(rgb(0.2,0.5,0.8),6)),las=3,ylim=c(-2,2),main=paste("n=",sum(A4%in%B4), "GUY11_8 deltaCDS1_1"))
abline(h=0,lty=2,col=rgb(0.2,0.4,0.6))
boxplot(cbind(scaled_heat_moc[A5[A5%in%B5],],scaled_heat_hu[B5[B5%in%A5],]),
        col=c(rep(rgb(0.7,0.2,0.2),6),
              rep(rgb(0.2,0.5,0.8),6)),las=3,ylim=c(-2,2),main=paste("n=",sum(A5%in%B5), "GUY11_2 deltaCDS1_3"))
abline(h=0,lty=2,col=rgb(0.2,0.4,0.6))
boxplot(cbind(scaled_heat_moc[A6[A6%in%B6],],scaled_heat_hu[B6[B6%in%A6],]),
        col=c(rep(rgb(0.7,0.2,0.2),6),
              rep(rgb(0.2,0.5,0.8),6)),las=3,ylim=c(-2,2),main=paste("n=",sum(A6%in%B6), "GUY11_1 deltaCDS1_4"))
abline(h=0,lty=2,col=rgb(0.2,0.4,0.6))
dev.off()


####
#### Figure 4A heatmap showing cadidate sites 
####

####Annotation using the kinase prediction and nuclear localization (GO annotation from ensembl using biomart)

annotations <- as.matrix(read.delim("input/mart_export_GO_MG8.txt", sep="\t", header=TRUE))
nuclear=annotations[annotations[,"GO.domain"]=="cellular_component",c("GO.term.accession","GO.term.name")]
nuclear=nuclear[!duplicated(nuclear),]
nuc_prot=unique(annotations[annotations[,2]%in%c("GO:0005634","GO:0005730","GO:0000786","GO:0005654"),1])


annot=cbind(apply(tab_4mir[,c("canonicalPattern","NetPhospan","KinLib(CHK1-CHK2)")],2,as.numeric),0)
colnames(annot)=c("canonicalPattern","NetPhospan","KinLib(CHK1-CHK2)","nuclear")
rownames(annot)=tab_4mir[,"alt_name"]


for (i in 1:nrow(annot)){
 annot[i,"nuclear"]=sum(gsub("T0|T1","",tab_4mir[i,"proteinID"])%in%nuc_prot)
}

annot[annot[,1]>1,1]=1
annot[annot[,2]>1,2]=1
annot[annot[,3]>1,3]=1
annot[annot[,4]>1,4]=1

######################### Candidate sites 

sel1= tab_4mir[,"candidate_HU_1h30min"]=="2" | tab_4mir[,"candidate_HU_2h"]=="2" | tab_4mir[,"candidate_HU_4h"]=="2" |
      tab_4mir[,"candidate_moc_1h30min"]=="2" | tab_4mir[,"candidate_moc_2h"]=="2" | tab_4mir[,"candidate_moc_4h"]=="2"

new=apply(tab_4mir[sel1,grep("diff",colnames(tab_4mir))],2,as.numeric)
rownames(new)=tab_4mir[sel1,"alt_name"]

tab_prot=tab_4mir[sel1,c("alt_name","proteinID")]
tab_prot[,2]=gsub("T0|T1","",tab_prot[,2])

#####Ahora solo los valores relevantes como guia de distancia

mat_heat=matrix(0,nrow(new),ncol(new))
colnames(mat_heat)=colnames(new)
rownames(mat_heat)=rownames(new)

temp=new[tab_4mir[tab_4mir[,"candidate_HU_1h30min"]=="2","alt_name"],c("GUY11_HU_1h30min_diff","DCDS1_HU_1h30min_diff")]
mat_heat[rownames(temp),colnames(temp)]=temp
temp=new[tab_4mir[tab_4mir[,"candidate_HU_2h"]=="2","alt_name"],c("GUY11_HU_2h_diff","DCDS1_HU_2h_diff")]
mat_heat[rownames(temp),colnames(temp)]=temp
temp=new[tab_4mir[tab_4mir[,"candidate_HU_4h"]=="2","alt_name"],c("GUY11_HU_4h_diff","DCDS1_HU_4h_diff")]
mat_heat[rownames(temp),colnames(temp)]=temp

temp=new[tab_4mir[tab_4mir[,"candidate_moc_1h30min"]=="2","alt_name"],c("GUY11_moc_1h30min_diff","DCDS1_moc_1h30min_diff")]
mat_heat[rownames(temp),colnames(temp)]=temp
temp=new[tab_4mir[tab_4mir[,"candidate_moc_2h"]=="2","alt_name"],c("GUY11_moc_2h_diff","DCDS1_moc_2h_diff")]
mat_heat[rownames(temp),colnames(temp)]=temp
temp=new[tab_4mir[tab_4mir[,"candidate_moc_4h"]=="2","alt_name"],c("GUY11_moc_4h_diff","DCDS1_moc_4h_diff")]
mat_heat[rownames(temp),colnames(temp)]=temp

###annotating the heatmap
annot_frame=as.data.frame(annot[rownames(mat_heat),])

my_colors <-  rev(colorRampPalette(brewer.pal(n = 5, name = "RdBu"))(100))
# Set breaks so that 0 is at the center
max_val <- max(mat_heat, na.rm = TRUE)
min_val <- min(mat_heat, na.rm = TRUE)
lim <- max(abs(c(min_val, max_val))) # symmetric limits
breaks <- seq(-6, 6, length.out = length(my_colors) + 1)





mat_heat=mat_heat[,c("GUY11_moc_1h30min_diff","GUY11_moc_2h_diff","GUY11_moc_4h_diff",
                     "DCDS1_moc_1h30min_diff","DCDS1_moc_2h_diff","DCDS1_moc_4h_diff",
                     "GUY11_HU_1h30min_diff","GUY11_HU_2h_diff","GUY11_HU_4h_diff",
                     "DCDS1_HU_1h30min_diff","DCDS1_HU_2h_diff","DCDS1_HU_4h_diff")]

x=pheatmap(mat_heat,col=my_colors,cluster_cols = F, 
           breaks = breaks)
dev.off()

new=new[x$tree_row$labels[x$tree_row$order],c("GUY11_moc_1h30min_diff","GUY11_moc_2h_diff","GUY11_moc_4h_diff",
                                              "DCDS1_moc_1h30min_diff","DCDS1_moc_2h_diff","DCDS1_moc_4h_diff",
                                              "GUY11_HU_1h30min_diff","GUY11_HU_2h_diff","GUY11_HU_4h_diff",
                                              "DCDS1_HU_1h30min_diff","DCDS1_HU_2h_diff","DCDS1_HU_4h_diff")]

mat_heat=mat_heat[x$tree_row$labels[x$tree_row$order],]
mark_display=matrix("",nrow(mat_heat),ncol(mat_heat))
rownames(mark_display)=rownames(mat_heat)
colnames(mark_display)=colnames(mat_heat)

for (i in 1:ncol(mark_display)){
  mark_display[mat_heat[,i]!=0,i]="*"
}

mark_display=mark_display[,colnames(new)]
annot_frame=annot_frame[x$tree_row$labels[x$tree_row$order],]

###Colors for annotation 
ann_colors <- list(
  nuclear = c(
    "0" = "white",
    "1" = "gold"
  ),
  "KinLib(CHK1-CHK2)" = c(
    "0" = "white",
    "1" = "dodgerblue"
  ),
  "NetPhospan" = c(
    "0" = "white",
    "1" = "blueviolet"
  ),
  canonicalPattern = c(
    "0" = "white",
    "1" = "forestgreen"
  )
)



pdf("output/FIG4A_heat_candidate_filter_allQ.pdf",height = 14, width = 6)
pheatmap(new,col=my_colors,cluster_rows = F,cluster_cols = F, 
         cellwidth =6, cellheight = 6 , fontsize = 6,
         annotation_row = annot_frame, 
         display_numbers = mark_display,
         breaks = breaks,
         annotation_colors = ann_colors,border_color = "white")
dev.off()


####
#### Figure 4B sequence windows
####

####Selection of sequence windows for logo using all pool of candidates or the ones affected by deltaCDS1 and HU treatment only

#sites affected by CDS1 and

###all candidates
sel1=tab_4mir[,"candidate_HU_1h30min"]=="2" | tab_4mir[,"candidate_HU_2h"]=="2" | tab_4mir[,"candidate_HU_4h"]=="2" |
  tab_4mir[,"candidate_moc_1h30min"]=="2" | tab_4mir[,"candidate_moc_2h"]=="2" | tab_4mir[,"candidate_moc_4h"]=="2"
###candidates affected by CDS1 and HU treatment
selct_sites=tab_4mir[(tab_4mir[,"candidate_HU_1h30min"]=="2" | tab_4mir[,"candidate_HU_2h"]=="2" | tab_4mir[,"candidate_HU_4h"]=="2") &
                       tab_4mir[,"candidate_moc_1h30min"]!="2" & tab_4mir[,"candidate_moc_2h"]!="2" & tab_4mir[,"candidate_moc_4h"]!="2","Sequence"]
sel2=sel1 & tab_4mir[,"Sequence"]%in%selct_sites

wind1=unlist(strsplit(tab_4mir[sel1,"window7Aa"],";"))
wind2=unlist(strsplit(tab_4mir[sel2,"window7Aa"],";"))

####Vamos a probar seleccionando el site más probable con el pattern
pattern=".....R.[ST].......|....R..[ST].......|..R....[ST]......."

cand_sites=tab_4mir[sel1,c("alt_name","Sequence","window7Aa")]
cand_sites=cbind(cand_sites,"")
colnames(cand_sites)=c("alt_name","Sequence","window7Aa","best_wind")

for (i in c(1:nrow(cand_sites))[grep(";",cand_sites[,"window7Aa"])]){
  
  temp=unlist(strsplit(cand_sites[i,"window7Aa"],";"))
  
  if(length(grep(pattern,temp))>0){
    
    cand_sites[i,"best_wind"]=paste(temp[grep(pattern,temp)],collapse=";")
  }else{
    
    cand_sites[i,"best_wind"]=paste(temp,collapse=";")
    
  }
}

cand_sites[cand_sites[,"best_wind"]=="","best_wind"]=cand_sites[cand_sites[,"best_wind"]=="","window7Aa"]

wind1=unlist(strsplit(cand_sites[,"best_wind"],";"))
wind2=unlist(strsplit(cand_sites[cand_sites[,"Sequence"]%in%selct_sites,"best_wind"],";"))

test=list(wind1,wind2)
names(test)=c("all (n=106 sites, 98 peptides)","selection (n=35 sites, 32 peptides)")

pdf("output/FIG4B_motif_sites.pdf",height = 3, width = 8)
ggseqlogo( test , method = 'prob',ncol=2)+
  scale_x_continuous(breaks = 1:15, labels = -7:7)
dev.off()

####
#### Figure 4C boxplot for pathway annotation and network study
####
#### For this figure and the following one we need two extra tables: manual annotation of pathways in 
#### Magnaporthe provided by the Oses-Ruiz group and PPI from STRING 12.0 database (https://string-db.org/)

####paths
pathway=as.matrix(read.delim("input/MOR_pathways.txt"))

####string
string_id=as.matrix(read.delim("input/242507.protein.aliases.v12.0.txt"))
stringnet=as.matrix(read.delim("input/242507.protein.links.detailed.v12.0.txt",sep=" "))
colnames(string_id)=gsub("X.","",colnames(string_id))
stringnet=stringnet[as.numeric(stringnet[,"combined_score"])>=400,]
####Selecting regulated sites


###all candidates
sel1=tab_4mir[,"candidate_HU_1h30min"]=="2" | tab_4mir[,"candidate_HU_2h"]=="2" | tab_4mir[,"candidate_HU_4h"]=="2" |
  tab_4mir[,"candidate_moc_1h30min"]=="2" | tab_4mir[,"candidate_moc_2h"]=="2" | tab_4mir[,"candidate_moc_4h"]=="2"
###candidates affected by CDS1 and HU treatment
selct_sites=tab_4mir[(tab_4mir[,"candidate_HU_1h30min"]=="2" | tab_4mir[,"candidate_HU_2h"]=="2" | tab_4mir[,"candidate_HU_4h"]=="2") &
                       tab_4mir[,"candidate_moc_1h30min"]!="2" & tab_4mir[,"candidate_moc_2h"]!="2" & tab_4mir[,"candidate_moc_4h"]!="2","Sequence"]
sel2=sel1 & tab_4mir[,"Sequence"]%in%selct_sites

####Providing node table (proteins with regulated sites plus proteins assigned to a pathway using Oses-Ruiz annotation)
nodes=unique(c(gsub("T0","",tab_4mir[sel1,"proteinID"]),pathway[,"geneID"]))
nodes=cbind(nodes,"")
colnames(nodes)=c("proteinID","STRING_ID")

for (i in 1:nrow(nodes)){
  
  if(sum(nodes[i,"proteinID"]%in%string_id[,"alias"])>0){
    nodes[i,"STRING_ID"]=unique(string_id[string_id[,"alias"]%in%nodes[i,"proteinID"],"string_protein_id"])
    
  }else{
    
    nodes[i,"STRING_ID"]=""    
    
  }}

nodes=cbind(nodes,"0")
colnames(nodes)=c("proteinID","STRING_ID","cand")
nodes=nodes[nodes[,"STRING_ID"]!="",]
nodes[nodes[,"proteinID"]%in%gsub("T0","",tab_4mir[sel1,"proteinID"]),"cand"]="1"
table(nodes[,"cand"])

####Selecting nodes that are part of an annotated pathway

pathway_tab=matrix(0,length(unique(pathway[,1])),3)
colnames(pathway_tab)=c("total","node","perc")
temp=table(pathway[,1])

rownames(pathway_tab)=names(temp)
pathway_tab[,"total"]=as.numeric(temp)

temp=table(pathway[pathway[,2]%in%nodes[,1],1])
pathway_tab[names(temp),"node"]=as.numeric(temp)
pathway_tab=pathway_tab[names(temp),]

pathway_tab[,"perc"]=pathway_tab[,"node"]/pathway_tab[,"total"]
pathway_tab=pathway_tab[order(pathway_tab[,"perc"],decreasing = T),]

###Definning network and avoiding selfloops and duplicated edges
net_back= simplify(graph_from_data_frame(d=stringnet[,c("protein1","protein2")],directed=F), remove.multiple = TRUE, remove.loops = TRUE)

##definning the candidate group (group_A) and counting the edges between the candidates and each pathway    
group_A <- unique(nodes[nodes[,"cand"]=="1","STRING_ID"])
obs_edges=rep(0,nrow(pathway_tab))

for (j in 1:nrow(pathway_tab)){
  
  group_B <- unique(nodes[nodes[,"proteinID"]%in%pathway[pathway[,"term"]%in%rownames(pathway_tab)[j],"geneID"],"STRING_ID"])
  obs_edges[j] <- edges_between_groups(net_back, group_A, group_B)
  
}

names(obs_edges)=rownames(pathway_tab)

####1000 permutaciones to calculate pvalues after reshuffleing the network
#### This last part of the code could take long time
permut_mat=matrix(0,1000,nrow(pathway_tab))
colnames(permut_mat)=rownames(pathway_tab)

for (i in 1:nrow(permut_mat)) {
  # Rewire edges while preserving degree sequence
  g_random <- rewire(net_back, keeping_degseq(niter = ecount(net_back) * 10))
  
  for (j in 1:nrow(pathway_tab)){
    
    group_B <- unique(nodes[nodes[,"proteinID"]%in%pathway[pathway[,"term"]%in%rownames(pathway_tab)[j],"geneID"],"STRING_ID"])
    permut_mat[i,j] <- edges_between_groups(g_random, group_A, group_B)
    
  }
}
p_value=rep(0,length(obs_edges))
for (i in 1:length(obs_edges)){
  p_value[i] <- (sum(permut_mat[,i] >= obs_edges[i]) + 1) / (nrow(permut_mat) + 1)
}
names(p_value)=names(obs_edges)
padj=p.adjust(p_value,method="BH")

## Getting the table with the final results 

net_mat_98=cbind(obs_edges,p_value,padj)
colnames(net_mat_98)=c("edge_count","pval","padj")
permut_mat_98=permut_mat

saveRDS(net_mat_98,"output/net_mat_98.rds")
saveRDS(permut_mat_98,"output/permut_mat_98.rds")


order_98=order(net_mat_98[,"pval"])

###ORder of significants ptahs usubjected to 98 candidates
col_box_98=rep("dodgerblue3",nrow(net_mat_98))
col_box_98[1:sum(net_mat_98[,"padj"]<=0.05)]="brown2"

pdf("output/FIG4C_boxplot_edgeBetween.pdf",height = 6, width = 6)

boxplot(permut_mat_98[,order_98],las=3,ylab="Edge count",col="white",ylim=c(0,120),main="98 phospho")
points(net_mat_98[order_98,"edge_count"],col=col_box_98,pch=16)
legend("topright",c("adjusted pval<0.05","adjusted pval>0.05"),fill=c("brown2","dodgerblue3"))

dev.off()

####
#### Figure 4D Graphic representation of significand results
####

nodes=unique(c(gsub("T0","",tab_4mir[sel1,"proteinID"]),pathway[,"geneID"]))
nodes=cbind(nodes,"")
colnames(nodes)=c("proteinID","STRING_ID")

for (i in 1:nrow(nodes)){
  
  if(sum(nodes[i,"proteinID"]%in%string_id[,"alias"])>0){
    nodes[i,"STRING_ID"]=unique(string_id[string_id[,"alias"]%in%nodes[i,"proteinID"],"string_protein_id"])
    
  }else{
    
    nodes[i,"STRING_ID"]=""    
    
  }}

nodes=cbind(nodes,"0")
colnames(nodes)=c("proteinID","STRING_ID","cand")
nodes=nodes[nodes[,"STRING_ID"]!="",]
nodes[nodes[,"proteinID"]%in%gsub("T0","",tab_4mir[sel1,"proteinID"]),"cand"]="1"
table(nodes[,"cand"])

#### Selection of nodes coming from significand pathways
nodes_sel=c(nodes[nodes[,"cand"]=="1","STRING_ID"],
            nodes[nodes[,"proteinID"]%in%pathway[pathway[,"term"]%in%rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05],"geneID"],"STRING_ID"])
nodes_sel=unique(nodes_sel)
sel_sig=stringnet[,"protein1"]%in%nodes_sel & stringnet[,"protein2"]%in%nodes_sel
###net and filtering for important connexion
net_edge=graph_from_data_frame(d=stringnet[sel_sig,c("protein1","protein2")],directed=F)
net_edge<- simplify(net_edge, remove.multiple = TRUE, remove.loops = TRUE)
# Kamada-Kawai layout
e <- get.edgelist(net_edge,names = F)
layout3 <- qgraph.layout.fruchtermanreingold(e, vcount = vcount(net_edge),
                                             area = 4*(vcount(net_edge)^2), repulse.rad = vcount(net_edge)^3)
net_node=nodes[nodes[,"cand"]=="1" |
                 nodes[,"proteinID"]%in%pathway[pathway[,"term"]%in%rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05],"geneID"],]
net_node=net_node[net_node[,"STRING_ID"]%in%V(net_edge)$name,]
nombres=net_node[,"proteinID"]
names(nombres)=net_node[,"STRING_ID"]
col_nod=rep("white",length(nombres))
names(col_nod)=V(net_edge)$name
col_nam=rep("black",length(nombres))
names(col_nam)=V(net_edge)$name

##Coloring nodes
colores=c("springgreen3","azure2","bisque2","brown2","dodgerblue2","gold2","tan2")

col_nod[names(col_nod)%in%nodes[nodes[,"proteinID"]%in%pathway[pathway[,1]==rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05][1],2],"STRING_ID"]]=colores[1]
col_nod[names(col_nod)%in%nodes[nodes[,"proteinID"]%in%pathway[pathway[,1]==rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05][2],2],"STRING_ID"]]=colores[2]
col_nod[names(col_nod)%in%nodes[nodes[,"proteinID"]%in%pathway[pathway[,1]==rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05][3],2],"STRING_ID"]]=colores[3]
col_nod[names(col_nod)%in%nodes[nodes[,"proteinID"]%in%pathway[pathway[,1]==rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05][4],2],"STRING_ID"]]=colores[4]
col_nod[names(col_nod)%in%nodes[nodes[,"proteinID"]%in%pathway[pathway[,1]==rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05][5],2],"STRING_ID"]]=colores[5]
col_nod[names(col_nod)%in%nodes[nodes[,"proteinID"]%in%pathway[pathway[,1]==rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05][6],2],"STRING_ID"]]=colores[6]
col_nod[names(col_nod)%in%nodes[nodes[,"proteinID"]%in%pathway[pathway[,1]==rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05][7],2],"STRING_ID"]]=colores[7]

col_nam[names(col_nam)%in%nodes[nodes[,"cand"]=="1","STRING_ID"]]=rgb(0.9,0.2,0.2)

V(net_edge)$color= col_nod
V(net_edge)$frame.color=col_nam

pdf("output/FIG4D_network_sig.pdf",height = 10, width = 10)
plot(net_edge,
     layout=layout3,
     vertex.size=5,
     vertex.label=nombres[V(net_edge)$name],
     vertex.label.cex=0.5,
     vertex.label.color="black")

legend("topright",rownames(net_mat_98)[net_mat_98[,"padj"]<=0.05],fill=colores)
dev.off()

####Writing the tables for the network
edge_list <- get.edgelist(net_edge)

colnames(edge_list)=c("A","B")

comb_sco=stringnet[sel_sig,c("combined_score")]
names(comb_sco)=paste(stringnet[sel_sig,c("protein1")],
                      stringnet[sel_sig,c("protein2")],
                      sep="_")

edge_list=cbind(edge_list,comb_sco[paste(edge_list[,"A"],
                                         edge_list[,"B"],
                                         sep="_")])
colnames(edge_list)=c("A","B","score")

net_node=cbind(net_node,"")
colnames(net_node)=c("proteinID","STRING_ID","cand","pathway")

padj=net_mat_98[,"padj"]
names(padj)=rownames(net_mat_98)
net_node[net_node[,"proteinID"]%in%pathway[pathway[,1]==names(padj)[padj<=0.05][1],2],"pathway"]=names(padj)[padj<=0.05][1]
net_node[net_node[,"proteinID"]%in%pathway[pathway[,1]==names(padj)[padj<=0.05][2],2],"pathway"]=names(padj)[padj<=0.05][2]
net_node[net_node[,"proteinID"]%in%pathway[pathway[,1]==names(padj)[padj<=0.05][3],2],"pathway"]=names(padj)[padj<=0.05][3]
net_node[net_node[,"proteinID"]%in%pathway[pathway[,1]==names(padj)[padj<=0.05][4],2],"pathway"]=names(padj)[padj<=0.05][4]
net_node[net_node[,"proteinID"]%in%pathway[pathway[,1]==names(padj)[padj<=0.05][5],2],"pathway"]=names(padj)[padj<=0.05][5]
net_node[net_node[,"proteinID"]%in%pathway[pathway[,1]==names(padj)[padj<=0.05][6],2],"pathway"]=names(padj)[padj<=0.05][6]
net_node[net_node[,"proteinID"]%in%pathway[pathway[,1]==names(padj)[padj<=0.05][7],2],"pathway"]=names(padj)[padj<=0.05][67]

write.table(net_node,"output/node_tab.txt",quote = F, row.names=F,sep="\t")
write.table(edge_list,"output/edge_tab.txt",quote = F, row.names=F,sep="\t")











