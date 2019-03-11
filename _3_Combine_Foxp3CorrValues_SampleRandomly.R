rm(list=ls())
fdir<-getwd();

# Load Residual Value of BRCA, COADREAD, LUAD and PAAD calculated using the full set of samples and markers.
load("BRCA_Residuals_Log2_Jan2018.RData");
BRCA_Before=New_RPKM;
BRCA_After=Residuals_Individual;
load("COADREAD_Residuals_Log2_Jan2018.RData");
COADREAD_Before=New_RPKM;
COADREAD_After=Residuals_Individual;
load("LUAD_Residuals_Log2_Jan2018.RData");
LUAD_Before=New_RPKM;
LUAD_After=Residuals_Individual;
load("PAAD_Residuals_Log2_Jan2018.RData");
PAAD_Before=New_RPKM;
PAAD_After=Residuals_Individual;
###########################################################
# Reshuffle data to calculate correlation as a control.

BRCA_After_Rnd=matrix(0,nrow=nrow(BRCA_After),ncol=ncol(BRCA_After));
COADREAD_After_Rnd=matrix(0,nrow=nrow(COADREAD_After),ncol=ncol(COADREAD_After));
LUAD_After_Rnd=matrix(0,nrow=nrow(LUAD_After),ncol=ncol(LUAD_After));
PAAD_After_Rnd=matrix(0,nrow=nrow(PAAD_After),ncol=ncol(PAAD_After));

for (i in 1:nrow(BRCA_After))
{
  BRCA_After_Rnd[i,]=sample(BRCA_After[i,]);
  COADREAD_After_Rnd[i,]=sample(COADREAD_After[i,]);
  LUAD_After_Rnd[i,]=sample(LUAD_After[i,]);
  PAAD_After_Rnd[i,]=sample(PAAD_After[i,]);
}

###########################################################
#Merging  correlation  values  across  datasets: In  order  to  have  a  comparable  set  of 
#correlation values for each of the 4 cancers we randomly selected a subset of equal number of 
#samples from each cancer dataset 10 times. "rnd" refers to the correlations calculated from reshuffled data-figure S4B.

a1=ncol(BRCA_Before);
a2=ncol(COADREAD_Before);
a3=ncol(LUAD_Before);
a4=ncol(PAAD_Before);

a=min(min(a1,a2),min(a3,a4));

id=which(rownames(Residuals_Individual)=="FOXP3");

BRCA_cor_before=matrix(nrow=20502,ncol=10);
BRCA_cor_after=matrix(nrow=20502,ncol=10);
BRCA_cor_after_rnd=matrix(nrow=20502,ncol=10);
for (i in 1:10)
{n=sample(a1,a,replace=FALSE);
BRCA_cor_before[,i]=cor(t(BRCA_Before[,n]),(BRCA_Before[id,n]),use="complete.obs");
BRCA_cor_after[,i]=cor(t(BRCA_After[,n]),BRCA_After[id,n],use="complete.obs");
BRCA_cor_after_rnd[,i]=cor(t(BRCA_After_Rnd[,n]),BRCA_After_Rnd[id,n],use="complete.obs");
}


COADREAD_cor_before=matrix(nrow=20502,ncol=10);
COADREAD_cor_after=matrix(nrow=20502,ncol=10);
COADREAD_cor_after_rnd=matrix(nrow=20502,ncol=10);

for (i in 1:10)
{n=sample(a2,a,replace=FALSE);
COADREAD_cor_before[,i]=cor(t(COADREAD_Before[,n]),(COADREAD_Before[id,n]),use="complete.obs");
COADREAD_cor_after[,i]=cor(t(COADREAD_After[,n]),COADREAD_After[id,n],use="complete.obs");
COADREAD_cor_after_rnd[,i]=cor(t(COADREAD_After_Rnd[,n]),COADREAD_After_Rnd[id,n],use="complete.obs");
}

LUAD_cor_before=matrix(nrow=20502,ncol=10);
LUAD_cor_after=matrix(nrow=20502,ncol=10);
LUAD_cor_after_rnd=matrix(nrow=20502,ncol=10);

for (i in 1:10)
{n=sample(a3,a,replace=FALSE);
LUAD_cor_before[,i]=cor(t(LUAD_Before[,n]),(LUAD_Before[id,n]),use="complete.obs");
LUAD_cor_after[,i]=cor(t(LUAD_After[,n]),LUAD_After[id,n],use="complete.obs");
LUAD_cor_after_rnd[,i]=cor(t(LUAD_After_Rnd[,n]),LUAD_After_Rnd[id,n],use="complete.obs");
}

PAAD_corB=cor(t(PAAD_Before),(PAAD_Before[id,]),use="complete.obs");
PAAD_corA=cor(t(PAAD_After),PAAD_After[id,],use="complete.obs");
PAAD_corA_rnd=cor(t(PAAD_After_Rnd),PAAD_After_Rnd[id,],use="complete.obs");

BRCA_corB=rowMeans(BRCA_cor_before);
BRCA_corA=rowMeans(BRCA_cor_after);
BRCA_corA_rnd=rowMeans(BRCA_cor_after_rnd);

COADREAD_corB=rowMeans(COADREAD_cor_before);
COADREAD_corA=rowMeans(COADREAD_cor_after);
COADREAD_corA_rnd=rowMeans(COADREAD_cor_after_rnd);

LUAD_corB=rowMeans(LUAD_cor_before);
LUAD_corA=rowMeans(LUAD_cor_after);
LUAD_corA_rnd=rowMeans(LUAD_cor_after_rnd);

Four_Total_Before=(BRCA_corB+COADREAD_corB+LUAD_corB+PAAD_corB)/4;
Four_Total_After=(BRCA_corA+COADREAD_corA+LUAD_corA+PAAD_corA)/4;
Four_Total_After_rnd=(BRCA_corA_rnd+COADREAD_corA_rnd+LUAD_corA_rnd+PAAD_corA_rnd)/4;



BRCA_cor=cbind(BRCA_corA,BRCA_corB);
rownames(BRCA_cor)=rownames(Four_Total_Before);
COADREAD_cor=cbind(COADREAD_corA,COADREAD_corB);
rownames(COADREAD_cor)=rownames(Four_Total_Before);
LUAD_cor=cbind(LUAD_corA,LUAD_corB);
rownames(LUAD_cor)=rownames(Four_Total_Before);
PAAD_cor=cbind(PAAD_corA,PAAD_corB);
rownames(PAAD_cor)=rownames(Four_Total_Before);
FourC_cor=cbind(Four_Total_After,Four_Total_Before);

