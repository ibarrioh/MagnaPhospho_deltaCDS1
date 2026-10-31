

####The following code depends on rds tables calculated in Rscript1 (already provided in folder output)

####Here we assemble the final table containing all results for further analysis at phosphosite level

####
#### Loading required tables
####

sites=readRDS("output/sites.rds")
input_dcds1=readRDS("output/input_dcds1.rds")
input_guy11=readRDS("output/input_guy11.rds")
tab_norm=readRDS("output/tab_norm.rds")
final=readRDS("output/limma_res.rds")
pep=readRDS("input/pep_tab_exp1.rds")

####Renaming the sites

sites=cbind(sites,"")
colnames(sites)[ncol(sites)]="alt_name"
sites[,"alt_name"]=sites[,"IDpos"]

for (i in grep(";",sites[,"IDpos"])){
  temp=unlist(strsplit(sites[i,"IDpos"],";"))
  temp1=gsub(unlist(strsplit(sites[i,"IDpos"],"-"))[1],"",temp[2:length(temp)])
  sites[i,"alt_name"]=paste(c(temp[1],gsub("-",",",temp1)),collapse="")
}

####
#### 1.-Adding information concerning kinase prediction for the sites
####
##Possible kinase pattern (Arginine residue at -2,-3 or -5 position)
pattern=".....R.[ST].......|....R..[ST].......|..R....[ST]......."

###Predictions calculated using Netphospan using the protein kinse loop of CDS1

phos=as.matrix(read.delim("input/Ex1_10Aawindow_NetPhospan.txt"))

###All predictions from the kinase library. It is important to notice that this is human kinase prediction, limitations must be considered. https://kinase-library.mit.edu/home

pred_tab1=as.matrix(read.delim("input/score-sites-bulk-serine-threonine.tsv"))
pred_chr=pred_tab1[,1:6]
pred_num=apply(pred_tab1[,7:ncol(pred_tab1)],2,as.numeric)
num_perc=pred_num[,seq(from=3, to=ncol(pred_num),by=4)]
colnames(num_perc)=gsub("_percentile","",colnames(num_perc))

##keep top prediction for CHK1 & CHK2
chk1_chk2=pred_chr[num_perc[,"CHK1"]>=90 |
                     num_perc[,"CHK2"]>=90  ,1]

annot=matrix(0,nrow(sites),3)
colnames(annot)=c("canonicalPattern","NetPhospan","KinLib(CHK1-CHK2)")
rownames(annot)=sites[,"topPep"]

for (i in 1:nrow(annot)){
  temp=unlist(strsplit(sites[i,"window7Aa"],";"))
  annot[i,"canonicalPattern"]=length(grep(pattern,temp))
  annot[i,"KinLib(CHK1-CHK2)"]=sum(temp%in%chk1_chk2)
  temp=unlist(strsplit(sites[i,"window7Aa"],";"))
  annot[i,"NetPhospan"]=sum(temp%in%substr(phos[as.numeric(phos[,"Pred"])>=0.5,"Peptide"],start=4,stop=18))
}

annot[annot[,1]>1,1]=1
annot[annot[,2]>1,2]=1
annot[annot[,3]>1,3]=1

#########
######### Creating the template and adding the kinase prediction data

tab_4mir=cbind(sites,annot)
rownames(tab_4mir)=tab_4mir[,"topPep"]
tab_4mir=tab_4mir[tab_4mir[,"topPep"]%in%unique(final[,"ID"]),]

#########
######### Adding limma results

qtab=matrix(0,nrow(tab_4mir),36)
colnames(qtab)=c(paste(rep(paste("GUY11_",unique(final[,"treat"]),sep=""),each=3),c("diff","pval","padj"),sep="_"),
                 paste(rep(paste("DCDS1_",unique(final[,"treat"]),sep=""),each=3),c("diff","pval","padj"),sep="_"))
rownames(qtab)=rownames(tab_4mir)

index=unique(final[,"treat"])

for (i in 1:length(index)){
  
  temp=final[final[,"treat"]==index[i],"diff_guy11"]
  names(temp)=final[final[,"treat"]==index[i],"ID"]
  qtab[names(temp),paste("GUY11",index[i],"diff",sep="_")]=temp
  
  temp=final[final[,"treat"]==index[i],"pval_guy11"]
  names(temp)=final[final[,"treat"]==index[i],"ID"]
  qtab[names(temp),paste("GUY11",index[i],"pval",sep="_")]=temp
  
  temp=final[final[,"treat"]==index[i],"padj_guy11"]
  names(temp)=final[final[,"treat"]==index[i],"ID"]
  qtab[names(temp),paste("GUY11",index[i],"padj",sep="_")]=temp
  
  temp=final[final[,"treat"]==index[i],"diff_dcds1"]
  names(temp)=final[final[,"treat"]==index[i],"ID"]
  qtab[names(temp),paste("DCDS1",index[i],"diff",sep="_")]=temp
  
  temp=final[final[,"treat"]==index[i],"pval_dcds1"]
  names(temp)=final[final[,"treat"]==index[i],"ID"]
  qtab[names(temp),paste("DCDS1",index[i],"pval",sep="_")]=temp
  
  temp=final[final[,"treat"]==index[i],"padj_dcds1"]
  names(temp)=final[final[,"treat"]==index[i],"ID"]
  qtab[names(temp),paste("DCDS1",index[i],"padj",sep="_")]=temp
}

####ADding the information for the peptide table
pep=pep[,1:17]
rownames(pep)=pep[,"ID"]
pep=pep[rownames(tab_4mir),]

####Categorical columns for candidate filtering
candi=matrix(0,nrow(tab_4mir),6)
colnames(candi)=paste("candidate_",unique(final[,"treat"]),sep="")

candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="1" & final[,"treat"]=="HU_1h30min","ID"],1]="1"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="1" & final[,"treat"]=="HU_2h","ID"],2]="1"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="1" & final[,"treat"]=="HU_4h","ID"],3]="1"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="1" & final[,"treat"]=="moc_1h30min","ID"],4]="1"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="1" & final[,"treat"]=="moc_2h","ID"],5]="1"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="1" & final[,"treat"]=="moc_4h","ID"],6]="1"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="2" & final[,"treat"]=="HU_1h30min","ID"],1]="2"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="2" & final[,"treat"]=="HU_2h","ID"],2]="2"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="2" & final[,"treat"]=="HU_4h","ID"],3]="2"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="2" & final[,"treat"]=="moc_1h30min","ID"],4]="2"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="2" & final[,"treat"]=="moc_2h","ID"],5]="2"
candi[tab_4mir[,"topPep"]%in%final[final[,"candidate_CDS1"]=="2" & final[,"treat"]=="moc_4h","ID"],6]="2"



tab_4mir=cbind(tab_4mir[,c("IDpos","alt_name","proteinID","Aa","position","loc_score","window7Aa","canonicalPattern","NetPhospan","KinLib(CHK1-CHK2)","topPep")],
               pep[rownames(tab_4mir),c("Sequence","Modifications.in.Master.Proteins")],
               candi,
               qtab[rownames(tab_4mir),],
               tab_norm[rownames(tab_4mir),])

saveRDS(tab_4mir,"output/tab_4mir.rds")
write.table(tab_4mir,"output/PhosphoPeptides_limma_280825.txt",quote = F, row.names=F,sep="\t")











