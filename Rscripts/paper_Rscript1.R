
#### Data analysis from modified peptide table ()
#### As input we are using the excel table "" converted to .txt and renamed ""
#### Peptides intensities are already normalized

####
#### Needed functions and packages 
####

library(limma)

####Function ofor inpoutition of NA based on distribution shifted towards less intense peptides (adapted from Perseus PMID:27348712 DOI:10.1038/nmeth.3901)
perseus_impute <- function(x, down_shift = 1.8, width = 0.3) {
  # x: numeric vector with NAs (e.g., log2 intensity values for one sample or feature)
  x_obs <- x[!is.na(x)]  # non-missing values
  
  mu <- mean(x_obs)      # mean of observed values
  sigma <- sd(x_obs)     # standard deviation of observed values
  
  # parameters for the imputation distribution
  imp_mean <- mu - down_shift * sigma
  imp_sd <- width * sigma
  
  # how many values to impute?
  n_miss <- sum(is.na(x))
  
  # impute missing values
  imputed_values <- rnorm(n_miss, mean = imp_mean, sd = imp_sd)
  
  # replace NAs
  x_imputed <- x
  x_imputed[is.na(x)] <- imputed_values
  
  return(x_imputed)
}

####
#### 0.- Generation of table to select the most intense phosphopeptide per phosphosite to avoid redundancy 
#### (misscleavages or modifications other than phospho, multiphosphorylated peptides are considered as "sites" of their own)
####
setwd()
#### Generating the sites table from pep table (peptides)
#### You can skip this step and use the provided "sites" table located in the output folder or you can generatedfrom scratch using the following lines 

##KEeping relevant columns for the site conversion 
pep=readRDS("input/pep_tab_exp1.rds")
pep=pep[,c("ID",
           "Annotated.Sequence",
           "Sequence",
           "Modifications",
           "Qvality.q.value",
           "X..Proteins",
           "Master.Protein.Accessions",
           "Positions.in.Master.Proteins",
           "Modifications.in.Master.Proteins")]
###Making sure all rows contain phosphorylated peptides
pep=pep[grep("Phospho",pep[,"Modifications"]),]

###template with the master accession
sites=cbind(pep[,c("Master.Protein.Accessions","Modifications.in.Master.Proteins","ID")],"")
colnames(sites)=c("proteinID","position","ID","phospho")
patterns=c("S/T/Y","T/S/Y","S/Y/T","; S]","; T]","; Y]",
           "S/T","S/Y","T/Y","T/S","Y/T","Y/S")

####### we are going to exclude the multi phosphos if not all sites have been identified with confidence
####### This action excludes 279 multiphospho peptides 
sites=sites[grep(paste(patterns,collapse="|"),sites[,"position"],invert=T),]
###We keep the first proteinID when more than one are present
ref_cut=grep(";",sites[,"proteinID"])

for (i in 1:length(ref_cut)){
  
  temp=unlist(strsplit(sites[ref_cut[i],"proteinID"],";"))[1]
  sites[ref_cut[i],"proteinID"]=temp
  
  temp=unlist(strsplit(sites[ref_cut[i],"position"],"; MGG_"))[1]
  sites[ref_cut[i],"position"]=temp
}
##clean the position column, erase protein name 
for (i in 1:nrow(sites)){
  
  sites[i,"position"]=gsub(sites[i,"proteinID"],"",sites[i,"position"])
  sites[i,"position"]=gsub(" ","",sites[i,"position"])
  
}

sites=sites[grep("Phospho",sites[,"position"]),]

# Extract only the 'NxPhospho[...]' segments
pattern <- "\\d+xPhospho\\[[^\\]]+\\]"

for (i in 1:nrow(sites)){
  temp=unlist(regmatches(sites[i,"position"], gregexpr(pattern, sites[i,"position"], perl = TRUE)))
  sites[i,"phospho"]=temp[1]
}

####Final ID for sites
sites=cbind(sites,"","")
colnames(sites)=c("proteinID","position","pepID","phospho","IDpos","loc_score")

for (i in c(1:nrow(sites))[grep("1xPhospho",sites[,"phospho"])]){
  
  temp=sites[i,"phospho"]
  temp=unlist(strsplit(gsub("1xPhospho\\[|\\]","",temp),"\\("))
  
  sites[i,"IDpos"]=paste(sites[i,"proteinID"],
                         temp[1],sep="-")
  sites[i,"loc_score"]=gsub(")","",temp[2])
}

for (i in c(1:nrow(sites))[grep("2xPhospho",sites[,"phospho"])]){
  
  temp=sites[i,"phospho"]
  temp=gsub("2xPhospho\\[|\\]","",temp)
  temp=unlist(strsplit(temp,";"))
  temp1=paste(sites[i,"proteinID"],
              unlist(strsplit(temp,"\\("))[c(1,3)],sep="-")
  temp2=gsub(")","",unlist(strsplit(temp,"\\("))[c(2,4)])
  
  sites[i,"IDpos"]=paste(temp1,collapse=";")
  sites[i,"loc_score"]=paste(temp2,collapse=";")
}

for (i in c(1:nrow(sites))[grep("3xPhospho",sites[,"phospho"])]){
  
  temp=sites[i,"phospho"]
  temp=gsub("3xPhospho\\[|\\]","",temp)
  temp=unlist(strsplit(temp,";"))
  temp1=paste(sites[i,"proteinID"],
              unlist(strsplit(temp,"\\("))[c(1,3,5)],sep="-")
  temp2=gsub(")","",unlist(strsplit(temp,"\\("))[c(2,4,6)])
  
  sites[i,"IDpos"]=paste(temp1,collapse=";")
  sites[i,"loc_score"]=paste(temp2,collapse=";")
  
}

###Ready for assembly

final=cbind(unique(sites[,"IDpos"]),matrix("",length(unique(sites[,"IDpos"])),7))
colnames(final)=c("IDpos","proteinID","Aa","position","loc_score","window7Aa","pepIDs","Pepseqs")

#site information single
sel=grep(";",final[,"IDpos"],invert=T)
temp=unlist(strsplit(final[sel,"IDpos"],"-"))
final[sel,"proteinID"]=temp[c(T,F)]
final[sel,"Aa"]=substr(temp[c(F,T)], 1, 1)
final[sel,"position"]=gsub("S|T|Y","",temp[c(F,T)])

#site information multiple
sel=grep(";",final[,"IDpos"])

for (i in sel){
  temp=unlist(strsplit(final[i,"IDpos"],";"))
  temp=unlist(strsplit(temp,"-"))
  final[i,"proteinID"]=temp[1]
  final[i,"Aa"]=paste(substr(temp[seq(from=2,to=length(temp),by=2)], 1, 1),collapse=";")
  final[i,"position"]=paste(gsub("S|T|Y","",temp[seq(from=2,to=length(temp),by=2)]),collapse=";")
}

#score
for (i in 1:nrow(final)){
  final[i,"loc_score"]=max(as.numeric(unlist(strsplit(sites[sites[,"IDpos"]%in%final[i,"IDpos"],"loc_score"],";"))))
}

##Lets go for the pepIDs 

for (i in 1:nrow(final)){
  final[i,"pepIDs"]=paste(sites[sites[,"IDpos"]%in%final[i,"IDpos"],"pepID"],collapse=";")
}
for (i in 1:nrow(final)){
  temp_pep=unique(pep[pep[,"ID"]%in%unlist(strsplit(final[i,"pepIDs"],";")),"Sequence"])
  final[i,"Pepseqs"]=paste(temp_pep,collapse=";")
}

######### Sequence windows 
###proteome of magnaporthe that was used for the mass spec 
fasta=readLines("input/MagnaportheOryzaeMG8_20200715.fasta")
##generating table with seqs
ID_pos=grep(">",fasta)
fasta_ids=fasta[ID_pos]

start=ID_pos+1
stop=c(ID_pos[2:length(ID_pos)]-1,length(fasta))

fasta_seq=rep("",length(fasta_ids))

for (i in 1:length(ID_pos)){
  fasta_seq[i]=paste(fasta[start[i]:stop[i]],collapse="")
}
###Modifying for seq windows in protein ends (avoiding errors when cutting, careful when mapping)
fasta_seq=paste("__________",fasta_seq,"__________",sep="")
for (i in 1:length(fasta_ids)){fasta_ids[i]=gsub(">","",unlist(strsplit(fasta_ids[i], " "))[1])}
fasta_seq=fasta_seq[fasta_ids%in%final[,"proteinID"]]
fasta_ids=fasta_ids[fasta_ids%in%final[,"proteinID"]]

####Cutting fasta

for (i in grep(";",final[,"IDpos"],invert=T)){
  
  pos=as.numeric(final[i,"position"])+10
  seq_temp=fasta_seq[fasta_ids==final[i,"proteinID"]]
  final[i,"window7Aa"]=substr(seq_temp,start=pos-7,stop=pos+7)
}

for (i in grep(";",final[,"IDpos"])){
  pos=as.numeric(unlist(strsplit(final[i,"position"],";")))+10
  seq_temp=fasta_seq[fasta_ids==final[i,"proteinID"]]
  temp=unlist(mapply(substr, MoreArgs = list(x = seq_temp), start = pos-7, stop = pos+7))
  final[i,"window7Aa"]=paste(temp,collapse=";")
}


########### Renaming the final site table, saving the results and erasing all the intermediates

sites=final
saveRDS(sites,"output/phosphoSTY_tab_exp1.rds")
rm(final,temp,temp1,temp2,fasta,fasta_seq,fasta_ids,ID_pos,pep)
gc()

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


######You can start from here using the precalculated sites tables provided in the output folder

pep=readRDS("input/pep_tab_exp1.rds")
sites=readRDS("output/phosphoSTY_tab_exp1.rds")


#### Renaming tabs
tab_chr=sites
tab_chr_pep=pep[,"ID"]
##selecting numeric columns with intensities
tab_num_pep=pep[,151:192]
##renaming columns
colnames(tab_num_pep)=c("dcds1_HU_1h30min_rep1","dcds1_HU_1h30min_rep2","dcds1_HU_1h30min_rep3",
                        "dcds1_HU_2h_rep1","dcds1_HU_2h_rep2","dcds1_HU_2h_rep3",
                        "dcds1_HU_4h_rep1","dcds1_HU_4h_rep2","dcds1_HU_4h_rep3",
                        "dcds1_moc_1h_rep1","dcds1_moc_1h_rep2","dcds1_moc_1h_rep3",
                        "dcds1_moc_1h30min_rep1","dcds1_moc_1h30min_rep2","dcds1_moc_1h30min_rep3",
                        "dcds1_moc_2h_rep1","dcds1_moc_2h_rep2","dcds1_moc_2h_rep3",
                        "dcds1_moc_4h_rep1","dcds1_moc_4h_rep2","dcds1_moc_4h_rep3",
                        "guy11_HU_1h30min_rep1","guy11_HU_1h30min_rep2","guy11_HU_1h30min_rep3",
                        "guy11_HU_2h_rep1","guy11_HU_2h_rep2","guy11_HU_2h_rep3",
                        "guy11_HU_4h_rep1","guy11_HU_4h_rep2","guy11_HU_4h_rep3",
                        "guy11_moc_1h_rep1","guy11_moc_1h_rep2","guy11_moc_1h_rep3",
                        "guy11_moc_1h30min_rep1","guy11_moc_1h30min_rep2","guy11_moc_1h30min_rep3",
                        "guy11_moc_2h_rep1","guy11_moc_2h_rep2","guy11_moc_2h_rep3",
                        "guy11_moc_4h_rep1","guy11_moc_4h_rep2","guy11_moc_4h_rep3")

##Correcting a swap in two of the sample name tags
colnames(tab_num_pep)[colnames(tab_num_pep)%in%c("guy11_moc_2h_rep1","guy11_moc_2h_rep2")]=c("guy11_moc_2h_rep2","guy11_moc_2h_rep1")
##keeping replicas 1 & 3 (replica 2 is less abundant and shows a different dynamic)
tab_num_pep=tab_num_pep[,grep("rep2",invert=T,colnames(tab_num_pep))]
##chaging decimal comma pojnts and NA to 0
for (i in 1:ncol(tab_num_pep)){tab_num_pep[,i]=gsub(",",".",tab_num_pep[,i])}
tab_num_pep=apply(tab_num_pep,2,as.numeric)
for (i in 1:ncol(tab_num_pep)){tab_num_pep[is.na(tab_num_pep[,i]),i]=0}
####Vector with sum intensities to choose the best peptides per site
all_int=apply(tab_num_pep,1,sum)
names(all_int)=tab_chr_pep
####Selecting the ID for the top peptide
sites=cbind(sites,"")
colnames(sites)[ncol(sites)]="topPep"
sites[grep(";",sites[,"pepIDs"],invert=T),"topPep"]=sites[grep(";",sites[,"pepIDs"],invert=T),"pepIDs"]
###Preparing the final table with the top intensities
tab_num=matrix(0,nrow(sites),ncol(tab_num_pep))
colnames(tab_num)=colnames(tab_num_pep)

for (i in grep(";",sites[,"pepIDs"])){
  
  temp=all_int[unlist(strsplit(sites[i,"pepIDs"],";"))]
  temp=temp[order(temp,decreasing=T)]
  tab_num[i,]=tab_num_pep[tab_chr_pep==names(temp)[1],]
  
  sites[i,"topPep"]=names(temp)[1]
}

for (i in grep(";",sites[,"pepIDs"],invert=T)){
  tab_num[i,]=tab_num_pep[tab_chr_pep%in%unlist(strsplit(sites[i,"pepIDs"],";")),]
}

rownames(tab_num)=sites[,"topPep"]

####missing values asigned 0 back to NAs for the inputation step and log transform
for (i in 1:ncol(tab_num)){tab_num[tab_num[,i]==0,i]=NA}

tab_norm=log2(tab_num)

####cleaning intermediate tables and nam 
rm(pep, tab_num_pep, tab_chr_pep)

####
#### 1.- Filtering and inputation of missing values  
#### Peptides that are present in both replicas from at least one condition are considered for analysis
#### Wildtype (GUY11) and deltaCDS1 are segregated for filtering and limma test (step2)

#### input matrixes for each experiment

input_dcds1=tab_norm[,grep("dcds1",colnames(tab_norm))]
input_guy11=tab_norm[,grep("guy11",colnames(tab_norm))]

#### Setting up the filters

filt_dcds1=rep(0,nrow(tab_norm))
for (i in 1:length(filt_dcds1)){ 
  
  temp=input_dcds1[i,]
  
  if(sum(is.na(temp[1:2]))==0 | sum(is.na(temp[3:4]))==0 | 
     sum(is.na(temp[5:6]))==0 | sum(is.na(temp[7:8]))==0 | 
     sum(is.na(temp[9:10]))==0 | sum(is.na(temp[11:12]))==0 | 
     sum(is.na(temp[13:14]))==0){
    
    filt_dcds1[i]=2
    
  }else{
    filt_dcds1[i]=0
    
  }}


filt_guy11=rep(0,nrow(tab_norm))

for (i in 1:length(filt_guy11)){ 
  
  temp=input_guy11[i,]
  
  if(sum(is.na(temp[1:2]))==0 | sum(is.na(temp[3:4]))==0 | 
     sum(is.na(temp[5:6]))==0 | sum(is.na(temp[7:8]))==0 | 
     sum(is.na(temp[9:10]))==0 | sum(is.na(temp[11:12]))==0 | 
     sum(is.na(temp[13:14]))==0){
    
    filt_guy11[i]=2
    
  }else{
    
    filt_guy11[i]=0
    
  }}


filt_dcds1=filt_dcds1==2
filt_guy11=filt_guy11==2

######Inputation

input_dcds1=input_dcds1[filt_dcds1,]
input_guy11=input_guy11[filt_guy11,]

for (i in 1:ncol(input_dcds1)){ input_dcds1[,i]=perseus_impute(input_dcds1[,i]) }
for (i in 1:ncol(input_guy11)){ input_guy11[,i]=perseus_impute(input_guy11[,i]) }


####
#### 2.- Limma test
####

####  GUY11

ctrl="guy11_moc_1h"
treat=colnames(input_guy11)[grep(paste(ctrl,c("rep1","rep3"),sep="_",collapse="|"),colnames(input_guy11),invert=T)]
treat=unique(gsub("_rep1|_rep2|_rep3","",treat))

res_boot=matrix(0,length(treat)*nrow(input_guy11),5)
colnames(res_boot)=c("ID","treat","diff","pval","padj")
res_boot[,"treat"]=rep(treat,each=nrow(input_guy11))

index=gsub("_rep1|_rep2|_rep3","",colnames(input_guy11))

for (i in 1:length(treat)){
  temp_imp= cbind(input_guy11[,index%in%ctrl],input_guy11[,index%in%treat[i]])
  ###Limma
  group <- factor(c( "Control", "Control","Treatment", "Treatment"))
  design <- model.matrix(~ group)
  # Fit linear model (invertir control y treatment para uqe quede igual el signo de la comparacion)
  fit <- lmFit(temp_imp, design)
  # Apply empirical Bayes moderation
  fit <- eBayes(fit)
  # Extract top differentially expressed features
  results <- topTable(fit, coef = "groupTreatment", number = Inf)
  results=as.matrix(results[order(as.numeric(rownames(results))),])
  res_boot[res_boot[,"treat"]==treat[i],"ID"]=rownames(input_guy11)
  res_boot[res_boot[,"treat"]==treat[i],c("diff","pval","padj")]=results[,c("logFC","P.Value","adj.P.Val")]
}

res_boot_guy11=res_boot

####  DCDS1

ctrl="dcds1_moc_1h"
treat=colnames(input_dcds1)[grep(paste(ctrl,c("rep1","rep3"),sep="_",collapse="|"),colnames(input_dcds1),invert=T)]
treat=unique(gsub("_rep1|_rep2|_rep3","",treat))

res_boot=matrix(0,length(treat)*nrow(input_dcds1),5)
colnames(res_boot)=c("ID","treat","diff","pval","padj")
res_boot[,"treat"]=rep(treat,each=nrow(input_dcds1))

index=gsub("_rep1|_rep2|_rep3","",colnames(input_dcds1))

for (i in 1:length(treat)){
  temp_imp= cbind(input_dcds1[,index%in%ctrl],input_dcds1[,index%in%treat[i]])
  ###Limma
  group <- factor(c( "Control", "Control","Treatment", "Treatment"))
  design <- model.matrix(~ group)
  # Fit linear model (invertir control y treatment para uqe quede igual el signo de la comparacion)
  fit <- lmFit(temp_imp, design)
  # Apply empirical Bayes moderation
  fit <- eBayes(fit)
  # Extract top differentially expressed features
  results <- topTable(fit, coef = "groupTreatment", number = Inf)
  results=as.matrix(results[order(as.numeric(rownames(results))),])
  res_boot[res_boot[,"treat"]==treat[i],"ID"]=rownames(input_dcds1)
  res_boot[res_boot[,"treat"]==treat[i],c("diff","pval","padj")]=results[,c("logFC","P.Value","adj.P.Val")]
}

res_boot_dcds1=res_boot

###Rename rows
code_dcds1=paste(res_boot_dcds1[,"ID"],gsub("dcds1_","",res_boot_dcds1[,"treat"]),sep="-")
code_guy11=paste(res_boot_guy11[,"ID"],gsub("guy11_","",res_boot_guy11[,"treat"]),sep="-")
code_comb=unique(c(code_dcds1,code_guy11))

rownames(res_boot_dcds1)=code_dcds1
rownames(res_boot_guy11)=code_guy11

final=matrix("",length(code_comb),11)
colnames(final)=c("ID","treat","diff_guy11","pval_guy11","padj_guy11",
                  "diff_dcds1","pval_dcds1","padj_dcds1","filter_guy11","filter_dcds1","candidate_CDS1")
rownames(final)=code_comb
final[,"ID"]=unlist(strsplit(code_comb,"-"))[c(T,F)]
final[,"treat"]=unlist(strsplit(code_comb,"-"))[c(F,T)]

final[code_guy11,c("diff_guy11","pval_guy11","padj_guy11")]=res_boot_guy11[code_guy11,c("diff","pval","padj")]
final[code_dcds1,c("diff_dcds1","pval_dcds1","padj_dcds1")]=res_boot_dcds1[code_dcds1,c("diff","pval","padj")]

final[code_guy11,"filter_guy11"]="1"
final[code_dcds1,"filter_dcds1"]="1"

###Adding categorical column to look for candidate sites  

temp=final[final[,"filter_guy11"]=="1" & final[,"filter_dcds1"]=="1",]

sel1=as.numeric(temp[,"diff_guy11"])>1 & as.numeric(temp[,"padj_guy11"])<=0.05 & 
     as.numeric(temp[,"diff_dcds1"])<1 & as.numeric(temp[,"padj_dcds1"])>0.05

sel2=as.numeric(temp[,"diff_guy11"])<1 & as.numeric(temp[,"padj_guy11"])>0.05 &
     as.numeric(temp[,"diff_dcds1"])<=(-1) & as.numeric(temp[,"padj_dcds1"])<=0.05

final[rownames(temp)[sel1|sel2],"candidate_CDS1"]="1"

sel1=as.numeric(temp[,"diff_guy11"])>1 & as.numeric(temp[,"padj_guy11"])<=0.01 & 
  as.numeric(temp[,"diff_dcds1"])<1 & as.numeric(temp[,"padj_dcds1"])>0.05

sel2=as.numeric(temp[,"diff_guy11"])<1 & as.numeric(temp[,"padj_guy11"])>0.05 &
  as.numeric(temp[,"diff_dcds1"])<=(-1) & as.numeric(temp[,"padj_dcds1"])<=0.01

final[rownames(temp)[sel1|sel2],"candidate_CDS1"]="2"


######Saving all files as RDS (dependencies for the following scripts)

saveRDS(sites,"output/sites.rds")
saveRDS(input_dcds1,"output/input_dcds1.rds")
saveRDS(input_guy11,"output/input_guy11.rds")
saveRDS(tab_norm,"output/tab_norm.rds")
saveRDS(final,"output/limma_res.rds")





