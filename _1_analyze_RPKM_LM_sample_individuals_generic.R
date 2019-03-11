rm(list=ls())
fdir<-getwd();

#read in expression from TCGA: "normalized_count"(quantile normalized RSEM) value 
#from illumina hiseq/ga2 mRNAseq level_3 (v2) data released on 20160228.

RSEM_Gene_IlluminaHiSeq<-read.table("/CancerType.RSEM_genes_normalized.txt",header=TRUE,sep="\t");

AllData<-RSEM_Gene_IlluminaHiSeq[,2:ncol(RSEM_Gene_IlluminaHiSeq)];
names1<-RSEM_Gene_IlluminaHiSeq[2:20532,1];
names2=as.character(names1);
names=matrix(nrow=length(names2),ncol=1);

for (i in 1:length(names1))
{names3=strsplit(names2[i],"|",fixed=TRUE);
 names[i]=names3[[1]][1];
}


RPKM<-matrix(nrow=20531,ncol=ncol(RSEM_Gene_IlluminaHiSeq)-1);
sample_names<-matrix(nrow=ncol(RSEM_Gene_IlluminaHiSeq)-1,ncol=1);
for (i in 1:(ncol=ncol(RSEM_Gene_IlluminaHiSeq)-1))
{
RPKM[,i]<-as.numeric(as.matrix(AllData[2:20532,i])); 
sample_names[i]<-colnames(RSEM_Gene_IlluminaHiSeq)[i+1];
}
rownames(RPKM)<-unlist(names);
colnames(RPKM)<-sample_names;

new_sample_names=matrix(nrow=ncol(RSEM_Gene_IlluminaHiSeq)-1,ncol=1);
Case_Control=matrix(nrow=ncol(RSEM_Gene_IlluminaHiSeq)-1,ncol=1);
for (i in 1:ncol(RSEM_Gene_IlluminaHiSeq)-1)
{x=unlist(strsplit(sample_names[i], "[.]"));
 new_sample_names[i]=paste(x[1],x[2],x[3],sep=".");
 Case_Control[i]=x[4];
}
tumor_id=grep("01",Case_Control);
RPKM_log2=log2(RPKM+1);

#read in the immune cell signatures
TCell_Sign<-read.table("TCell_Sign.txt",header=FALSE,sep="\t");
BCell_Sign<-read.table("BCell_Sign.txt",header=FALSE,sep="\t");
DC_Cell_Sign<-read.table("DCCell_Sign.txt",header=FALSE,sep="\t");
Eos_Cell_Sign<-read.table("EosCell_Sign.txt",header=FALSE,sep="\t");
MPhage_Cell_Sign<-read.table("MPhageCell_Sign.txt",header=FALSE,sep="\t");
Neutrophil_Cell_Sign<-read.table("NeutrophilCell_Sign.txt",header=FALSE,sep="\t");
NK_Cell_Sign<-read.table("NKCell_Sign.txt",header=FALSE,sep="\t");
Mast_Cell_Sign<-read.table("MastCell_Sign.txt",header=FALSE,sep="\t");


#Calculate average gene expression among signature genes in tumor samples
TCell_Genes=intersect(unlist(TCell_Sign),names);
TCell_Genes_id=match(TCell_Genes,names);
keep=which(apply(RPKM_log2[TCell_Genes_id,tumor_id],1,mean)>1);
TCell_index=colMeans(RPKM_log2[TCell_Genes_id[keep],tumor_id]);
 
BCell_Genes=intersect(unlist(BCell_Sign),names);
BCell_Genes_id=match(BCell_Genes,names);
keep=which(apply(RPKM_log2[BCell_Genes_id,tumor_id],1,mean)>1);
BCell_index=colMeans(RPKM_log2[BCell_Genes_id[keep],tumor_id]);
 
DC_Genes=intersect(unlist(DC_Cell_Sign),names);
DC_Genes_id=match(DC_Genes,names);
keep=which(apply(RPKM_log2[DC_Genes_id,tumor_id],1,mean)>1);
DC_index=colMeans(RPKM_log2[DC_Genes_id[keep],tumor_id]);

Eos_Genes=intersect(unlist(Eos_Cell_Sign),names);
Eos_Genes_id=match(Eos_Genes,names);
keep=which(apply(RPKM_log2[Eos_Genes_id,tumor_id],1,mean)>1);
Eos_index=colMeans(RPKM_log2[Eos_Genes_id[keep],tumor_id]);
 
MPhage_Genes=intersect(unlist(MPhage_Cell_Sign),names);
MPhage_Genes_id=match(MPhage_Genes,names);
keep=which(apply(RPKM_log2[MPhage_Genes_id,tumor_id],1,mean)>1);
MPhage_index=colMeans(RPKM_log2[MPhage_Genes_id[keep],tumor_id]);
 
Neutrophil_Genes=intersect(unlist(Neutrophil_Cell_Sign),names);
Neutrophil_Genes_id=match(Neutrophil_Genes,names);
keep=which(apply(RPKM_log2[Neutrophil_Genes_id,tumor_id],1,mean)>1);
Neutrophil_index=colMeans(RPKM_log2[Neutrophil_Genes_id[keep],tumor_id]);
 
NK_Genes=intersect(unlist(NK_Cell_Sign),names);
NK_Genes_id=match(NK_Genes,names);
keep=which(apply(RPKM_log2[NK_Genes_id,tumor_id],1,mean)>1);
NK_index=colMeans(RPKM_log2[NK_Genes_id[keep],tumor_id]);


Mast_Genes=intersect(unlist(Mast_Cell_Sign),names);
Mast_Genes_id=match(Mast_Genes,names);
keep=which(apply(RPKM_log2[Mast_Genes_id,tumor_id],1,mean)>1);
Mast_index=colMeans(RPKM_log2[Mast_Genes_id[keep],tumor_id]);

# randomly select 100 samples 10 times and generate residual values using the full set of markers.

for (k in 1:10)
{
tumor_id1=sample(c(1:length(tumor_id)),100,replace=FALSE);
A1=rbind(RPKM_log2[,tumor_id],TCell_index,BCell_index,DC_index,Eos_index,MPhage_index,Neutrophil_index,NK_index,Mast_index);
A2=as.data.frame(t(A1[,tumor_id1]));

##########################
 
Residuals_Individual<-matrix(nrow=20502,ncol=length(tumor_id1));
rownames(Residuals_Individual)=rownames(RPKM)[30:20531];

X=gsub("-","",colnames(A2));
colnames(A2)=X;
 
 
New_RPKM1=RPKM[30:20531,];
New_RPKM=log2(New_RPKM1[,tumor_id[tumor_id1]]+1);

#claculate residuals of immune infiltrate regression
for (i in 30:20531)
 {
  form2 <- paste(colnames(A2)[i], "~","TCell_index+BCell_index+DC_index+Eos_index+MPhage_index+Neutrophil_index+NK_index+Mast_index",sep="");
  q2=lm(formula=form2,data=A2);
  Residuals_Individual[i-29,]=q2$residuals;
 }
 
 colnames(Residuals_Individual)=colnames(RPKM[,tumor_id[tumor_id1]]);

 save("New_RPKM","Residuals_Individual",file=paste("CancerType_Residuals_Log2_sample_v",k,".RData"));
}

