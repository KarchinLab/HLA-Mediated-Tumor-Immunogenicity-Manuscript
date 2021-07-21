library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(psych)
library(RColorBrewer)
library(gplots)
library(corrplot)
library(ggthemes)
library(rms)
library(gtools)
library(survival)
library(survminer)
library(rstatix)
library(finalfit)
library(readr)
library(pracma)
library(reshape2)
library(stats)
library(survcomp)
library(ggplot2)
library(viridis)
library(fdrtool)
library(plyr)
library(stringr)
library(gdata)
library(optparse)
library(openxlsx)
library(readr)
library(Rsubread)
library("org.Hs.eg.db")
library("AnnotationDbi")
library(stringr)
library(data.table)
library(CePa)
theme_set(theme_classic((base_size=20)))

# overview section (Mar 19 2021). Related to Figure 1 and S1.
# Figure 1A: see heatmap.R

# Figure 1B
overview_df = subset(df_immuneLandscape, df_immuneLandscape$Study!='DLBC' & df_immuneLandscape$Study!='MESO' & df_immuneLandscape$Study!='TGCT' & df_immuneLandscape$Study!='OV' & df_immuneLandscape$Study!='LAML' & df_immuneLandscape$Study!='UCS' & df_immuneLandscape$Study!='UVM' & df_immuneLandscape$Study!='ACC' & df_immuneLandscape$Study!='LGG')
overview_df$tumor = overview_df$geom_HLAI
overview_df$normal = overview_df$geom_HLAI_N
dat.m1 <- melt(overview_df,id.vars='Study', measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable
dat.m1$rankI=NA
dat.m1$rankI[which(dat.m1$Study=='KIRC')] = 1
dat.m1$rankI[which(dat.m1$Study=='CESC')] = 2
dat.m1$rankI[which(dat.m1$Study=='LUAD')] = 3
dat.m1$rankI[which(dat.m1$Study=='HNSC')] = 4
dat.m1$rankI[which(dat.m1$Study=='PAAD')] = 5
dat.m1$rankI[which(dat.m1$Study=='KIRP')] = 6
dat.m1$rankI[which(dat.m1$Study=='SKCM')] = 7
dat.m1$rankI[which(dat.m1$Study=='COAD')] = 8
dat.m1$rankI[which(dat.m1$Study=='READ')] = 9
dat.m1$rankI[which(dat.m1$Study=='THCA')] = 10
dat.m1$rankI[which(dat.m1$Study=='CHOL')] = 11
dat.m1$rankI[which(dat.m1$Study=='UCEC')] = 12
dat.m1$rankI[which(dat.m1$Study=='KICH')] = 13
dat.m1$rankI[which(dat.m1$Study=='BLCA')] = 14
dat.m1$rankI[which(dat.m1$Study=='LUSC')] = 15
dat.m1$rankI[which(dat.m1$Study=='SARC')] = 16
dat.m1$rankI[which(dat.m1$Study=='STAD')] = 17
dat.m1$rankI[which(dat.m1$Study=='PCPG')] = 18
dat.m1$rankI[which(dat.m1$Study=='THYM')] = 19
dat.m1$rankI[which(dat.m1$Study=='LIHC')] = 20
dat.m1$rankI[which(dat.m1$Study=='BRCA')] = 21
dat.m1$rankI[which(dat.m1$Study=='ESCA')] = 22
dat.m1$rankI[which(dat.m1$Study=='GBM')] = 23
dat.m1$rankI[which(dat.m1$Study=='PRAD')] = 24

# FDR: q values
pvalues=ddply(dat.m1, .(Study), function(x) { wilcox.test(value~Tissue, data=x, paired=FALSE)$p.value })
pvalues$rankI=NA
pvalues$rankI[which(pvalues$Study=='KIRC')] = 1
pvalues$rankI[which(pvalues$Study=='CESC')] = 2
pvalues$rankI[which(pvalues$Study=='LUAD')] = 3
pvalues$rankI[which(pvalues$Study=='HNSC')] = 4
pvalues$rankI[which(pvalues$Study=='PAAD')] = 5
pvalues$rankI[which(pvalues$Study=='KIRP')] = 6
pvalues$rankI[which(pvalues$Study=='SKCM')] = 7
pvalues$rankI[which(pvalues$Study=='COAD')] = 8
pvalues$rankI[which(pvalues$Study=='READ')] = 9
pvalues$rankI[which(pvalues$Study=='THCA')] = 10
pvalues$rankI[which(pvalues$Study=='CHOL')] = 11
pvalues$rankI[which(pvalues$Study=='UCEC')] = 12
pvalues$rankI[which(pvalues$Study=='KICH')] = 13
pvalues$rankI[which(pvalues$Study=='BLCA')] = 14
pvalues$rankI[which(pvalues$Study=='LUSC')] = 15
pvalues$rankI[which(pvalues$Study=='SARC')] = 16
pvalues$rankI[which(pvalues$Study=='STAD')] = 17
pvalues$rankI[which(pvalues$Study=='PCPG')] = 18
pvalues$rankI[which(pvalues$Study=='THYM')] = 19
pvalues$rankI[which(pvalues$Study=='LIHC')] = 20
pvalues$rankI[which(pvalues$Study=='BRCA')] = 21
pvalues$rankI[which(pvalues$Study=='ESCA')] = 22
pvalues$rankI[which(pvalues$Study=='GBM')] = 23
pvalues$rankI[which(pvalues$Study=='PRAD')] = 24
pvalues2=pvalues %>% arrange(rankI)
pVal_vec=pvalues2$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=rev(stars.pval(qVal))
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
calibrate=0.1
y<-rep(c(14.7-calibrate),each=24)
y[which(str_length(stars_vec) == 1)]=14.8-calibrate
y[which(str_length(stars_vec) == 3)]=14.6-calibrate
x<-rep(1:24)+0.15

y2<-rep(c(14.8),each=24)
x2<-rep(1:24)+0.01
format_pVals=function(x) {
  if (x < 0.01) {
    formatC(x, format = "e", digits = 2)
  } else {
    round(x,2)
  }
}
pVal_format=unlist(sapply(pVal_vec, format_pVals))
pVal_format=rev(pVal_format)
pVal_format=paste('p = ',pVal_format,sep='')

# Image: 650*700
ggplot(dat.m1, aes(x=reorder(Study, (rankI)*(-1)), y=value, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
  xlab('') + ylab('HLA Class I Expression') + ylim(c(6.5,16))
  coord_flip() + theme_bw() + 
  annotate(geom = "text", y = c(y,y2), x = c(x,x2), label = c(stars_vec,pVal_format), size = 3.5, hjust = 0)


# Figure 1C
overview_df = subset(df_immuneLandscape, df_immuneLandscape$Study!='DLBC' & df_immuneLandscape$Study!='MESO' & df_immuneLandscape$Study!='TGCT' & df_immuneLandscape$Study!='OV' & df_immuneLandscape$Study!='LAML' & df_immuneLandscape$Study!='UCS' & df_immuneLandscape$Study!='UVM' & df_immuneLandscape$Study!='ACC' & df_immuneLandscape$Study!='LGG')
overview_df$tumor = overview_df$geom_HLAII
overview_df$normal = overview_df$geom_HLAII_N
dat.m1 <- melt(overview_df,id.vars='Study', measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable
dat.m1$rankII=NA
dat.m1$rankII[which(dat.m1$Study=='KIRC')] = 1
dat.m1$rankII[which(dat.m1$Study=='LUAD')] = 2
dat.m1$rankII[which(dat.m1$Study=='THYM')] = 3
dat.m1$rankII[which(dat.m1$Study=='PAAD')] = 4
dat.m1$rankII[which(dat.m1$Study=='BRCA')] = 5
dat.m1$rankII[which(dat.m1$Study=='SKCM')] = 6
dat.m1$rankII[which(dat.m1$Study=='KIRP')] = 7
dat.m1$rankII[which(dat.m1$Study=='STAD')] = 8
dat.m1$rankII[which(dat.m1$Study=='CESC')] = 9
dat.m1$rankII[which(dat.m1$Study=='THCA')] = 10
dat.m1$rankII[which(dat.m1$Study=='LUSC')] = 11
dat.m1$rankII[which(dat.m1$Study=='GBM')] = 12
dat.m1$rankII[which(dat.m1$Study=='SARC')] = 13
dat.m1$rankII[which(dat.m1$Study=='HNSC')] = 14
dat.m1$rankII[which(dat.m1$Study=='UCEC')] = 15
dat.m1$rankII[which(dat.m1$Study=='CHOL')] = 16
dat.m1$rankII[which(dat.m1$Study=='BLCA')] = 17
dat.m1$rankII[which(dat.m1$Study=='KICH')] = 18
dat.m1$rankII[which(dat.m1$Study=='COAD')] = 19
dat.m1$rankII[which(dat.m1$Study=='ESCA')] = 20
dat.m1$rankII[which(dat.m1$Study=='PRAD')] = 21
dat.m1$rankII[which(dat.m1$Study=='READ')] = 22
dat.m1$rankII[which(dat.m1$Study=='PCPG')] = 23
dat.m1$rankII[which(dat.m1$Study=='LIHC')] = 24
dat.m1=dat.m1 %>% arrange(rankII)

# FDR: q values
pvalues=ddply(dat.m1, .(Study), function(x) { wilcox.test(value~Tissue, data=x, paired=FALSE)$p.value })
pvalues$rankII=NA
pvalues$rankII[which(pvalues$Study=='KIRC')] = 1
pvalues$rankII[which(pvalues$Study=='LUAD')] = 2
pvalues$rankII[which(pvalues$Study=='THYM')] = 3
pvalues$rankII[which(pvalues$Study=='PAAD')] = 4
pvalues$rankII[which(pvalues$Study=='BRCA')] = 5
pvalues$rankII[which(pvalues$Study=='SKCM')] = 6
pvalues$rankII[which(pvalues$Study=='KIRP')] = 7
pvalues$rankII[which(pvalues$Study=='STAD')] = 8
pvalues$rankII[which(pvalues$Study=='CESC')] = 9
pvalues$rankII[which(pvalues$Study=='THCA')] = 10
pvalues$rankII[which(pvalues$Study=='LUSC')] = 11
pvalues$rankII[which(pvalues$Study=='GBM')] = 12
pvalues$rankII[which(pvalues$Study=='SARC')] = 13
pvalues$rankII[which(pvalues$Study=='HNSC')] = 14
pvalues$rankII[which(pvalues$Study=='UCEC')] = 15
pvalues$rankII[which(pvalues$Study=='CHOL')] = 16
pvalues$rankII[which(pvalues$Study=='BLCA')] = 17
pvalues$rankII[which(pvalues$Study=='KICH')] = 18
pvalues$rankII[which(pvalues$Study=='COAD')] = 19
pvalues$rankII[which(pvalues$Study=='ESCA')] = 20
pvalues$rankII[which(pvalues$Study=='PRAD')] = 21
pvalues$rankII[which(pvalues$Study=='READ')] = 22
pvalues$rankII[which(pvalues$Study=='PCPG')] = 23
pvalues$rankII[which(pvalues$Study=='LIHC')] = 24
pvalues2=pvalues %>% arrange(rankII)
pVal_vec=pvalues2$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=rev(stars.pval(qVal))
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
x=0.4
y<-rep(c(14.7-x),each=24)
y[which(str_length(stars_vec) == 1)]=14.85-x
y[which(str_length(stars_vec) == 3)]=14.55-x
x<-rep(1:24)+0.15

y2<-rep(c(14.55),each=24)
x2<-rep(1:24)+0.01
format_pVals=function(x) {
  if (x < 0.01) {
    formatC(x, format = "e", digits = 2)
  } else {
    round(x,2)
  }
}
pVal_format=unlist(sapply(pVal_vec, format_pVals))
pVal_format=rev(pVal_format)
pVal_format=paste('p = ',pVal_format,sep='')

# Image: 650*700
ggplot(dat.m1, aes(x=reorder(Study, (rankII)*(-1)), y=value, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
  xlab('') + ylab('HLA Class II Expression') + ylim(c(3,16.2)) +
  coord_flip() + theme_bw() + 
  annotate(geom = "text", y = c(y,y2), x = c(x,x2), label = c(stars_vec,pVal_format), size = 3.5, hjust = 0)

# Figure S1A-D: correlation between HLA expression and TMB/neoantigen loads
# 650*430
ggplot(df_immuneLandscape, aes(x=geom_HLAI, y=mutationrate_nonsilent_per_Mb)) + 
  geom_point() + geom_smooth(method=loess) +
  scale_color_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_cor(method = "spearman", label.x = c(10, 0.4), label.y = c(400, 1.8)) +
  xlab('HLA-I Expression') + ylab('TMB') + theme_bw()

ggplot(df_immuneLandscape, aes(x=geom_HLAII, y=mutationrate_nonsilent_per_Mb)) + 
  geom_point() + geom_smooth(method=loess) +
  scale_color_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_cor(method = "spearman", label.x = c(6.5, 0.4), label.y = c(400, 1.8)) +
  xlab('HLA-II Expression') + ylab('TMB') + theme_bw()

ggplot(df_immuneLandscape, aes(x=geom_HLAI, y=numberOfBindingExpressedPMHC)) + 
  geom_point() + geom_smooth(method=loess) +
  scale_color_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_cor(method = "spearman", label.x = c(9, 0.4), label.y = c(25000, 1.8)) +
  xlab('HLA-I Expression') + ylab('Neoantigen Load') + theme_bw()

ggplot(df_immuneLandscape, aes(x=geom_HLAII, y=numberOfBindingExpressedPMHC)) + 
  geom_point() + geom_smooth(method=loess) +
  scale_color_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_cor(method = "spearman", label.x = c(6, 0.4), label.y = c(25000, 1.8)) +
  xlab('HLA-II Expression') + ylab('Neoantigen Load') + theme_bw()

ggplot(df_immuneLandscape, aes(x=mutationrate_nonsilent_per_Mb, y=numberOfBindingExpressedPMHC)) + 
  geom_point() + geom_smooth(method=loess) +
  scale_color_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_cor(method = "spearman", label.x = c(250, 0.4), label.y = c(25000, 1.8)) +
  xlab('TMB') + ylab('Neoantigen Load') + theme_bw()


# Figure S1E-F: see GTEx-figures.R
# immunogeneity section (Mar 19 2021). Related to Figure 2 and S2.
# Figure 2A: see heatmap.R

# Figure S2A: CYT tumor vs. normal
overview_df = subset(df_immuneLandscape, df_immuneLandscape$Study!='DLBC' & df_immuneLandscape$Study!='MESO' & df_immuneLandscape$Study!='TGCT' & df_immuneLandscape$Study!='OV' & df_immuneLandscape$Study!='LAML' & df_immuneLandscape$Study!='UCS' & df_immuneLandscape$Study!='UVM' & df_immuneLandscape$Study!='ACC' & df_immuneLandscape$Study!='LGG')
overview_df$CYT=log2(overview_df$tumor_GZMA + 1) + log2(overview_df$tumor_PRF1 + 1) / 2
overview_df$CYT_N=log2(overview_df$normal_GZMA + 1) + log2(overview_df$normal_PRF1 + 1) / 2
overview_df$tumor = overview_df$CYT
overview_df$normal = overview_df$CYT_N
dat.m1 <- melt(overview_df,id.vars='Study', measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable
dat.m1$rankI=NA
dat.m1$rankI[which(dat.m1$Study=='KIRC')] = 1
dat.m1$rankI[which(dat.m1$Study=='CESC')] = 2
dat.m1$rankI[which(dat.m1$Study=='LUAD')] = 3
dat.m1$rankI[which(dat.m1$Study=='THYM')] = 4
dat.m1$rankI[which(dat.m1$Study=='LUSC')] = 5
dat.m1$rankI[which(dat.m1$Study=='HNSC')] = 6
dat.m1$rankI[which(dat.m1$Study=='STAD')] = 7
dat.m1$rankI[which(dat.m1$Study=='PAAD')] = 8
dat.m1$rankI[which(dat.m1$Study=='CHOL')] = 9
dat.m1$rankI[which(dat.m1$Study=='UCEC')] = 10
dat.m1$rankI[which(dat.m1$Study=='SKCM')] = 11
dat.m1$rankI[which(dat.m1$Study=='SARC')] = 12
dat.m1$rankI[which(dat.m1$Study=='COAD')] = 13
dat.m1$rankI[which(dat.m1$Study=='BRCA')] = 14
dat.m1$rankI[which(dat.m1$Study=='READ')] = 15
dat.m1$rankI[which(dat.m1$Study=='KIRP')] = 16
dat.m1$rankI[which(dat.m1$Study=='BLCA')] = 17
dat.m1$rankI[which(dat.m1$Study=='ESCA')] = 18
dat.m1$rankI[which(dat.m1$Study=='LIHC')] = 19
dat.m1$rankI[which(dat.m1$Study=='KICH')] = 20
dat.m1$rankI[which(dat.m1$Study=='THCA')] = 21
dat.m1$rankI[which(dat.m1$Study=='PRAD')] = 22
dat.m1$rankI[which(dat.m1$Study=='GBM')] = 23
dat.m1$rankI[which(dat.m1$Study=='PCPG')] = 24

# FDR: q values
pvalues=ddply(dat.m1, .(Study), function(x) { wilcox.test(value~Tissue, data=x, paired=FALSE)$p.value })
pvalues$rankI=NA
pvalues$rankI[which(pvalues$Study=='KIRC')] = 1
pvalues$rankI[which(pvalues$Study=='CESC')] = 2
pvalues$rankI[which(pvalues$Study=='LUAD')] = 3
pvalues$rankI[which(pvalues$Study=='THYM')] = 4
pvalues$rankI[which(pvalues$Study=='LUSC')] = 5
pvalues$rankI[which(pvalues$Study=='HNSC')] = 6
pvalues$rankI[which(pvalues$Study=='STAD')] = 7
pvalues$rankI[which(pvalues$Study=='PAAD')] = 8
pvalues$rankI[which(pvalues$Study=='CHOL')] = 9
pvalues$rankI[which(pvalues$Study=='UCEC')] = 10
pvalues$rankI[which(pvalues$Study=='SKCM')] = 11
pvalues$rankI[which(pvalues$Study=='SARC')] = 12
pvalues$rankI[which(pvalues$Study=='COAD')] = 13
pvalues$rankI[which(pvalues$Study=='BRCA')] = 14
pvalues$rankI[which(pvalues$Study=='READ')] = 15
pvalues$rankI[which(pvalues$Study=='KIRP')] = 16
pvalues$rankI[which(pvalues$Study=='BLCA')] = 17
pvalues$rankI[which(pvalues$Study=='ESCA')] = 18
pvalues$rankI[which(pvalues$Study=='LIHC')] = 19
pvalues$rankI[which(pvalues$Study=='KICH')] = 20
pvalues$rankI[which(pvalues$Study=='THCA')] = 21
pvalues$rankI[which(pvalues$Study=='PRAD')] = 22
pvalues$rankI[which(pvalues$Study=='GBM')] = 23
pvalues$rankI[which(pvalues$Study=='PCPG')] = 24
pvalues2=pvalues %>% arrange(rankI)
pVal_vec=pvalues2$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=rev(stars.pval(qVal))
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
calibrate=0.3
y<-rep(c(14.7-calibrate),each=24)
y[which(str_length(stars_vec) == 1)]=14.9-calibrate
y[which(str_length(stars_vec) == 3)]=14.5-calibrate
x<-rep(1:24)+0.15

y2<-rep(c(14.8),each=24)
x2<-rep(1:24)+0.01
format_pVals=function(x) {
  if (x < 0.01) {
    formatC(x, format = "e", digits = 2)
  } else {
    round(x,2)
  }
}
pVal_format=unlist(sapply(pVal_vec, format_pVals))
pVal_format=rev(pVal_format)
pVal_format=paste('p = ',pVal_format,sep='')

# Image: 650*700
ggplot(dat.m1, aes(x=reorder(Study, (rankI)*(-1)), y=value, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
  xlab('') + ylab('Cytolytic Activity: GZMA, PRF1') +
  coord_flip() + theme_bw() + ylim(c(0,17)) +
  annotate(geom = "text", y = c(y,y2), x = c(x,x2), label = c(stars_vec,pVal_format), size = 3.5, hjust = 0)


# Figure S2B: see heatmap.R for correlation heatmap of comparison between TMB, neoantigens and HLA genes

# Figure 2B: scatterplot between supertypes and CYT
immune_df=df_immuneLandscape
immune_df$B2M=immune_df$B2M_log
immune_df$`HLA-A`=immune_df$HLA_A_log
immune_df$`HLA-B`=immune_df$HLA_B_log
immune_df$`HLA-C`=immune_df$HLA_C_log
immune_df$`HLA-E`=immune_df$HLA_E_log
immune_df$`HLA-G`=immune_df$HLA_G_log
immune_df$`HLA-DP`=immune_df$HLA_DP_log
immune_df$`HLA-DQ`=immune_df$HLA_DQ_log
immune_df$`HLA-DR`=immune_df$HLA_DR_log
immune_df$`HLA-DM`=immune_df$HLA_DM_log
immune_df$`HLA-DO`=immune_df$HLA_DO_log
dat.m1 <- melt(immune_df,id.vars='PRF1_GZMA', measure.vars=c('B2M','HLA-A','HLA-B','HLA-C','HLA-E','HLA-G','HLA-DP','HLA-DQ','HLA-DR','HLA-DM','HLA-DO'))
dat.m1$Gene=dat.m1$variable
# Image: 890 * 580
# Image: 840 * 520
ggplot(dat.m1, aes(x=value, y=PRF1_GZMA, color=Gene)) + 
  geom_point(size=0.00005) + theme_bw() +
  geom_smooth(aes(color = Gene, fill = Gene), method = "loess", se=F) + 
  scale_color_viridis(discrete = TRUE, option = 'D')+
  scale_fill_viridis(discrete = TRUE) +
  stat_cor(method = "spearman", label.x = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 3.7, 3.7, 3.7, 3.7, 3.7), label.y = c(14.5,14,13.5, 13,12.5,12,14.5,14,13.5, 13,12.5, 12)) +
  xlab('Expression') + ylab('CYT')

# Figure 2C: PCA, see PCA.R

# Figure 2SC-D: see 0711-UPDATE-highCYT-classifiers.ipynb



# HED section (Mar 18 2021). Related to Figure 3 and S3
# Figure S3A
# count number of heterozygous alleles
df_immuneLandscape$numHeterozygous=0
for (i in (1:nrow(df_immuneLandscape))) {
  if (is.na(df_immuneLandscape$HED[i])) {
    df_immuneLandscape$numHeterozygous[i]=NA
    next
  }
  if (df_immuneLandscape$grantham_A[i] != 0) {
    df_immuneLandscape$numHeterozygous[i] = df_immuneLandscape$numHeterozygous[i] + 1
  }
  if (df_immuneLandscape$grantham_B[i] != 0) {
    df_immuneLandscape$numHeterozygous[i] = df_immuneLandscape$numHeterozygous[i] + 1
  }
  if (df_immuneLandscape$grantham_C[i] != 0) {
    df_immuneLandscape$numHeterozygous[i] = df_immuneLandscape$numHeterozygous[i] + 1
  }
  if (df_immuneLandscape$HED[i] != 0 & df_immuneLandscape$numHeterozygous[i] == 0) {
    message(df_immuneLandscape$ParticipantBarcode[i])
  }
}
df_immuneLandscape$numHeterozygous=factor(df_immuneLandscape$numHeterozygous, levels(factor(df_immuneLandscape$numHeterozygous)))
df_immuneLandscape$`Heterozygous Loci` = df_immuneLandscape$numHeterozygous
# Density plot: heterozygous loci
# 810 * 520
ggplot(df_immuneLandscape %>% filter(!is.na(HED))) + 
  geom_density(aes(x = HED, fill = `Heterozygous Loci`), alpha = 0.7) + 
  scale_fill_viridis(discrete = T, option = 'A')

# Figure 3A: histogram of all tumor samples
median_HED=round(median(df_immuneLandscape$HED,na.rm=TRUE),2)
# Image: 690 * 450
ggplot(df_immuneLandscape) +
  geom_histogram(aes(x = HED, y = ..density.., fill = ..x..), binwidth = 0.1) + 
  scale_fill_gradient2(low='darkgreen', mid='goldenrod', high='red', midpoint=6.65, name= 'HED') +
  xlab('HED score') + theme_bw() + theme(legend.position = 'none') +
  annotate(geom = "text", y=0.25, x = median_HED-1, label = paste("median HED = ",median_HED,sep=''), size = 4, hjust = 0)

# Figure S3B: HED correlation matrix with immune signatures 
HED_immuneLandscape=df_immuneLandscape
df=subset(HED_immuneLandscape,select=(c("HED","T.cells.CD8","T.cells.CD4.memory.activated","Th1.cells","Th2.cells","Th17.cells","T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
colnames(df)=c("HED","CD8 Tc1","CD4 Tmem","Th1","Th2","Th17","Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
matrix=corr.test(df, method = "spearman", adjust = 'fdr')
# Image: 790 * 640
corrplot(matrix$r, tl.col = "black", tl.srt = 40,method="color",p.mat = matrix$p, sig.level = 0.01,
         col=colorRampPalette(c("deepskyblue4","white","firebrick3"))(200))

# Figure 3B: quadrants
HED_immuneLandscape$highCYT=HED_immuneLandscape$PRF1_GZMA > unname(quantile(HED_immuneLandscape$PRF1_GZMA,na.rm=TRUE)["75%"])
HED_immuneLandscape$highCYT=factor(HED_immuneLandscape$highCYT, rev(levels(factor(HED_immuneLandscape$highCYT))))
HED_CYT_corr=cor(HED_immuneLandscape$HED, HED_immuneLandscape$PRF1_GZMA,  method = "spearman", use = "complete.obs")
HLAI_CYT_corr=cor(HED_immuneLandscape$geom_HLAI, HED_immuneLandscape$PRF1_GZMA,  method = "spearman", use = "complete.obs")

# Image: 660 * 520
ggplot(HED_immuneLandscape %>% filter(!is.na(highCYT)), aes(x=HED, y=geom_HLAI,color=highCYT)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c("gold2", "forestgreen"), 
                    name="CYT",
                    labels=c("high", "low")) +
  geom_vline(xintercept = unname(quantile(HED_immuneLandscape$HED,na.rm=TRUE)["50%"])) +
  geom_hline(yintercept = unname(quantile(HED_immuneLandscape$geom_HLAI,na.rm=TRUE)["75%"])) +
  xlab('HED') + ylab('HLA-I expression') + theme_bw() +
  annotate(geom = "text", y = c(13,13,4.5,4.5), x = c(11,0,0,11), label = c("I","II","III","IV"), size = 7, hjust = 0, parse = TRUE)

# Figure 3C: all tumor samples; highHED: top 50% (ns); highHLAI: top 25%.
HED_immuneLandscape$highHED=HED_immuneLandscape$HED > unname(quantile(HED_immuneLandscape$HED,na.rm=TRUE)["50%"])
HED_immuneLandscape$highHLAI75=HED_immuneLandscape$geom_HLAI > unname(quantile(HED_immuneLandscape$geom_HLAI,na.rm=TRUE, c(0.5,0.7,0.75,0.8,0.85,0.9))["75%"])

HED_immuneLandscape$highHED_HLAI=HED_immuneLandscape$highHED==T & HED_immuneLandscape$highHLAI75==T
fit <- survfit(Surv(PFI_time_1,PFI_1) ~highHED_HLAI,data = HED_immuneLandscape)
ggsurv=ggsurvplot(fit, data = HED_immuneLandscape, palette = c('darkgreen','darkgoldenrod2'),  
           conf.int = F, 
           pval = F,  
           risk.table = TRUE, 
           risk.table.col = "strata",
           legend.labs = c('low HLAI or HED','high HLAI and HED'),
           legend.title = '',
           xlab = '', ylab = 'Survival',
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank()))
ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                            x = c(600, 800, 1700), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                            label = c("p<0.001", "HR = 0.80", "95% CI = 0.72-0.90"), size = 5)
# Image: 710 * 500
ggsurv
coxph(Surv(PFI_time_1, PFI_1) ~ highHED_HLAI, data = HED_immuneLandscape) %>% gtsummary::tbl_regression(exp = TRUE)

# Figure 3D:
highHED_highHLAI=subset(HED_immuneLandscape,HED_immuneLandscape$highHED==T & HED_immuneLandscape$highHLAI75==T)
highHED_highHLAI$Groups='high HED high HLAI\n(Quandrant I)'
highHED_lowHLAI=subset(HED_immuneLandscape,HED_immuneLandscape$highHED==T & HED_immuneLandscape$highHLAI75==F)
highHED_lowHLAI$Groups='high HED low HLAI\n(Quandrant II)'
lowHED_highHLAI=subset(HED_immuneLandscape,HED_immuneLandscape$highHED==F & HED_immuneLandscape$highHLAI75==T)
lowHED_highHLAI$Groups='low HED high HLAI\n(Quandrant IV)'
lowHED_lowHLAI=subset(HED_immuneLandscape,HED_immuneLandscape$highHED==F & HED_immuneLandscape$highHLAI75==F)
lowHED_lowHLAI$Groups='low HED low HLAI\n(Quandrant III)'
df=rbind(highHED_highHLAI, highHED_lowHLAI, lowHED_highHLAI, lowHED_lowHLAI)
df$Groups = factor(df$Groups, rev(levels(factor(df$Groups))))
model <- coxph( Surv(PFI_time_1, PFI_1) ~ Groups, data = df)
# Image: 870 * 500
ggforest(model, main = "",
         cpositions = c(0.02, 0.18, 0.4),
         fontsize = 0.5,
         refLabel = "reference",
         noDigits = 2)



# HLA expression fold change section (Mar 19 2021). Related to Figure 4 and S4
# Figure 4A-B: HLA expression in tumor vs. normal across immune subtypes
HLAupregulation_df=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based!='')
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C1']='C1\nwound-healing'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C2']='C2\nIFN-γ dominant'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C3']='C3\ninflammatory'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C4']='C4\nlymphocyte-depleted'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C5']='C5\nimmunologically quiet'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C6']='C6\nTGF-β dominant'

HLAupregulation_df$tumor = HLAupregulation_df$geom_HLAI
HLAupregulation_df$normal = HLAupregulation_df$geom_HLAI_N
dat.m1 <- melt(HLAupregulation_df,id.vars='immuneSubtype', measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable
# FDR: q values
pvalues=ddply(dat.m1, .(immuneSubtype), function(x) { wilcox.test(value~Tissue, data=x, paired=FALSE)$p.value })
pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=stars.pval(qVal)
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
y<-rep(c(15.15),each=6)
x<-c(0.63, 1.64, 2.63, 3.7, 4.71, 5.67)

# Image: 900 * 500
ggplot(dat.m1, aes(x=immuneSubtype, y=value, fill=Tissue)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
  stat_compare_means(aes(group = Tissue, label = paste0("p = ", ..p.format..)), label.x = 1, label.y =15, size=3) +
  ylab('HLA Class I Expression') + xlab('') +
  theme_bw() + 
  annotate(geom = "text", y = c(y), x = c(x), label = c(stars_vec), size = 3.5, hjust = 0)


HLAupregulation_df$tumor = HLAupregulation_df$geom_HLAII
HLAupregulation_df$normal = HLAupregulation_df$geom_HLAII_N
dat.m1 <- melt(HLAupregulation_df,id.vars='immuneSubtype', measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable
# FDR: q values
pvalues=ddply(dat.m1, .(immuneSubtype), function(x) { wilcox.test(value~Tissue, data=x, paired=FALSE)$p.value })
pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=stars.pval(qVal)
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
y<-rep(c(15.2),each=6)
x<-c(0.64, 2.0, 2.63, 3.7, 4.73, 5.67)

# Image: 770 * 450
ggplot(dat.m1, aes(x=immuneSubtype, y=value, fill=Tissue)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
  stat_compare_means(aes(group = Tissue, label = paste0("p = ", ..p.format..)), label.x = 1, label.y = 15, size=3) +
  xlab('') + ylab('HLA Class II Expression') +
  theme_bw() +
  annotate(geom = "text", y = c(y), x = c(x), label = c(stars_vec), size = 3.5, hjust = 0)

# Figure S4C
HLAupregulation_df=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based!='')
HLAupregulation_df$geom_high_HLAIupreg=HLAupregulation_df$geom_HLAI_upregulation > unname(quantile(HLAupregulation_df$geom_HLAI_upregulation,na.rm=TRUE)["50%"])
HLAupregulation_df$geom_high_HLAIIupreg=HLAupregulation_df$geom_HLAII_upregulation > unname(quantile(HLAupregulation_df$geom_HLAII_upregulation,na.rm=TRUE)["50%"])
HLAupregulation_df$geom_high_HLAI=HLAupregulation_df$geom_HLAI > unname(quantile(HLAupregulation_df$geom_HLAI,na.rm=TRUE)["50%"])
HLAupregulation_df$geom_high_HLAII=HLAupregulation_df$geom_HLAII > unname(quantile(HLAupregulation_df$geom_HLAII,na.rm=TRUE)["50%"])
# pairwise boxplots: exclude tumor samples with undefined immune subtypes
high_HLAI_high_HLAII=subset(HLAupregulation_df,HLAupregulation_df$geom_high_HLAIupreg==T & HLAupregulation_df$geom_high_HLAIIupreg==T)
label1='high HLAI high HLAII'
high_HLAI_high_HLAII$HLAI_HLAII=label1
high_HLAI_high_HLAII$Groups='high HLAI high HLAII'

high_HLAI_low_HLAII=subset(HLAupregulation_df,HLAupregulation_df$geom_high_HLAIupreg==T & HLAupregulation_df$geom_high_HLAIIupreg==F)
label2='high HLAI low HLAII'
high_HLAI_low_HLAII$HLAI_HLAII=label2
high_HLAI_low_HLAII$Groups='high HLAI low HLAII'

low_HLAI_low_HLAII=subset(HLAupregulation_df,HLAupregulation_df$geom_high_HLAIupreg==F & HLAupregulation_df$geom_high_HLAIIupreg==F)
label3='low HLAI low HLAII'
low_HLAI_low_HLAII$HLAI_HLAII=label3
low_HLAI_low_HLAII$Groups='low HLAI low HLAII'

low_HLAI_high_HLAII=subset(HLAupregulation_df,HLAupregulation_df$geom_high_HLAIupreg==F & HLAupregulation_df$geom_high_HLAIIupreg==T)
label4='low HLAI high HLAII'
low_HLAI_high_HLAII$HLAI_HLAII=label4
low_HLAI_high_HLAII$Groups='low HLAI high HLAII'

df=rbind(low_HLAI_low_HLAII, high_HLAI_low_HLAII, low_HLAI_high_HLAII, high_HLAI_high_HLAII)
# strata survival: HR
# 614*463
df$Groups = factor(df$Groups, rev(levels(factor(df$Groups))))
model <- coxph( Surv(PFI_time_1, PFI_1) ~ Groups, data = df)
ggforest(model, main = "",
         cpositions = c(0.02, 0.18, 0.4),
         fontsize = 0.5,
         refLabel = "reference",
         noDigits = 2)

# Figure 4C: fold change heatmap across immune subtypes, see heatmap.R

# Figure 4D, S4A-B
# pairwise boxplots
df$HLAI_HLAII = factor(df$HLAI_HLAII, rev(levels(factor(df$HLAI_HLAII))))
df$Groups = factor(df$Groups, rev(levels(factor(df$Groups))))
theme_set(theme_classic((base_size=15)))
my_comparisons <- list( c(label1, label4), c(label1, label2), c(label2, label3), c(label4, label3) )
# 870 * 520
ggplot(df, aes(x= HLAI_HLAII, y=Th1.cells, fill=Groups)) + 
  geom_boxplot() + scale_fill_brewer(palette="PuBuGn") + 
  stat_compare_means(label.x = 1, label.y = 1200) + 
  stat_compare_means(comparisons = my_comparisons) + 
  xlab('HLA-I vs. HLA-II fold change') + ylab('Th1 cells') +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16))
ggplot(df, aes(x= HLAI_HLAII, y=T.cells.CD8, fill=Groups)) + 
  geom_boxplot() + scale_fill_brewer(palette="PuBuGn") + 
  stat_compare_means(label.x = 1, label.y = 0.8) + 
  stat_compare_means(comparisons = my_comparisons) + 
  xlab('HLA-I vs. HLA-II  fold change') + ylab('CD8 Cytotoxic T cells') +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16))
ggplot(df, aes(x= HLAI_HLAII, y=PRF1_GZMA, fill=Groups)) + 
  geom_boxplot() + scale_fill_brewer(palette="PuBuGn") + 
  stat_compare_means(label.x = 3, label.y = 14) + 
  stat_compare_means(comparisons = my_comparisons) + 
  xlab('HLA-I vs. HLA-II  fold change') + ylab('CYT') +
  theme(axis.text=element_text(size=10),axis.title=element_text(size=16))

# Figure 4E: neural network, see 0711-UPDATE-multiclass-neural-network.ipynb

# Figure 4F, S4D: hazard ratio, see HR_plots.R


# methylation section (Mar 18 2021). Related to Figure 5 and S5
# exclude undefined immune subtypes
methylation_df = subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based != '')
methylation_df$geom_high_HLAIupreg=methylation_df$geom_HLAI_upregulation > unname(quantile(methylation_df$geom_HLAI_upregulation,na.rm=TRUE)["50%"])
methylation_df$geom_high_HLAIIupreg=methylation_df$geom_HLAII_upregulation > unname(quantile(methylation_df$geom_HLAII_upregulation,na.rm=TRUE)["50%"])

# Figure 5A
hypermethylation=1.25
hypomethylation=0.75
methylation_df$highCYT75=methylation_df$PRF1_GZMA > unname(quantile(methylation_df$PRF1_GZMA,na.rm=TRUE)["75%"])
highCYT_df=subset(methylation_df,methylation_df$highCYT75==T)
highCYT_df$CYT='high'
lowCYT_df=subset(methylation_df,methylation_df$highCYT75==F)
lowCYT_df$CYT='low'
df=rbind(highCYT_df,lowCYT_df)

# 720 * 500
ggplot(df, aes(x=normalized_met_I, y=geom_HLAI)) + 
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = CYT)) +
  geom_vline(xintercept = hypomethylation,linetype=2) +
  geom_vline(xintercept = hypermethylation,linetype=2) + 
  theme_bw() +
  scale_fill_manual(values = c("firebrick2", "deepskyblue2")) +
  stat_cor(method = "spearman", label.x = c(1.1), label.y = c(12)) +
  guides(alpha = F) +
  xlab('HLA I Methylation Level') + ylab('HLA I Expression')
  

# Figure 5B: methylation across six immune subtypes
my_comparisons <- list( c('C2','C3'), c('C2', 'C1'), c('C2', 'C4') )
methylation_df$`immune subtype`=methylation_df$Subtype_Immune_Model_Based
# Image: 670 * 450
ggplot(methylation_df %>% filter((methylation_df$Subtype_Immune_Model_Based != '') & !is.na(methylation_df$normalized_met_I)), aes(x= Subtype_Immune_Model_Based, y=normalized_met_I, fill=`immune subtype`)) + 
  geom_violin() + 
  stat_compare_means(label.x = 4, label.y = 2.7) +
  stat_compare_means(comparisons = my_comparisons) + 
  geom_boxplot(width=0.1) + 
  theme_bw() +
  guides(fill = F) +
  scale_fill_brewer(palette="Blues") + 
  xlab('Immune Subtype') + ylab('HLA I Methylation Level')

# Figure S5A: HLA-I and HLA-II fold change in C2 vs. C3
immunogenic=subset(methylation_df,methylation_df$Subtype_Immune_Model_Based == 'C2' | methylation_df$Subtype_Immune_Model_Based == 'C3')
immunogenic$`immune subtype`=immunogenic$Subtype_Immune_Model_Based
immunogenic$`HLA I`=immunogenic$geom_HLAI_upregulation
immunogenic$`HLA II`=immunogenic$geom_HLAII_upregulation
dat.m1 <- melt(immunogenic,id.vars='immune subtype', measure.vars=c('HLA I','HLA II'))
dat.m1$expression=dat.m1$value
# 644 * 545
ggplot(dat.m1, aes(x='immuhe subtype', y=expression, fill=`immune subtype`)) +
  geom_boxplot() + facet_grid(. ~variable ) + 
  scale_fill_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_compare_means(label.x = 1, label.y = 3.5, size=5) +
  xlab('') + ylab('Expression Fold Change')

# Figure 5C: CYT vs. methylation among immunogenic tumors (C2 + C3)
immunogenic=subset(methylation_df,methylation_df$Subtype_Immune_Model_Based == 'C2' | methylation_df$Subtype_Immune_Model_Based == 'C3')
immunogenic$`Immune Subtype`=immunogenic$Subtype_Immune_Model_Based
# Image: 690 * 500
ggplot(immunogenic, aes(x=PRF1_GZMA, y=normalized_met_I, color=`Immune Subtype`)) + 
  geom_point() + geom_smooth(method=loess) +
  scale_color_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_cor(method = "spearman", label.x = c(0.4, 0.4), label.y = c(1.9, 1.8)) +
  xlab('CYT') + ylab('HLA I Methylation Level') +
  guides(color = F)

immunogenic$highCYT=immunogenic$PRF1_GZMA > unname(quantile(immunogenic$PRF1_GZMA,na.rm=TRUE)["50%"])
highCYT_immunogenic=subset(immunogenic,highCYT == T)
lowCYT_immunogenic=subset(immunogenic,highCYT == F)
# Image: 540 * 530
ggplot(highCYT_immunogenic, aes(x=Subtype_Immune_Model_Based, y=normalized_met_I, fill=`Immune Subtype`)) + 
  geom_violin() + 
  stat_compare_means(label.x = 1.7, label.y = 1.6)+
  geom_boxplot(width=0.1) + 
  scale_fill_manual(values=c("goldenrod2", "deepskyblue3")) +
  xlab('Immune Subtype') + ylab('HLA I Methylation Level') +
  guides(fill = F)

# no significance in lowCYT_immunogenic
ggplot(lowCYT_immunogenic, aes(x=Subtype_Immune_Model_Based, y=normalized_met_I, fill=`Immune Subtype`)) + 
  geom_violin() + 
  stat_compare_means(label.x = 1, label.y = 1.5)+
  geom_boxplot(width=0.1) + 
  scale_fill_manual(values=c("goldenrod2", "deepskyblue3")) +
  geom_hline(aes(yintercept=1),linetype=2) +
  xlab('immune subtype') + ylab('HLA-I methylation')

# Figure 5D: fishers: HLA-II fold change vs. HLA I methylation
df=methylation_df
df$hypermethylation=df$normalized_met_I > 1.25
df$hypomethylation=df$normalized_met_I < 0.75

highHLAIIupreg=subset(df,df$geom_high_HLAIIupreg==T)
lowHLAIIupreg=subset(df,df$geom_high_HLAIIupreg==F)

piechart_df1=data.frame(c("HLA I\nHypermethylation", "HLA I\nHypomethylation"),c(nrow(subset(highHLAIIupreg,highHLAIIupreg$hypermethylation == T)),nrow(subset(highHLAIIupreg,highHLAIIupreg$hypomethylation == T))))
colnames(piechart_df1)=c("Methyl_state","num")
piechart_df1 <- piechart_df1 %>% 
  arrange(desc(Methyl_state)) %>%
  mutate(prop = num / sum(piechart_df1$num) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

piechart_df1$label=paste(piechart_df1$Methyl_state,'\n',piechart_df1$num,sep='')
# Image: 650 * 470
ggplot(piechart_df1, aes(x="", y=prop, fill=Methyl_state)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none", plot.title = element_text(size = 15)) +
  geom_text(aes(y = ypos, label = label), color = "white", size=5) +
  scale_fill_brewer(palette="Set1") + ggtitle("high HLA-II fold change (top 50%)") +
  theme(plot.title = element_text(hjust = 0.5))

piechart_df2=data.frame(c("HLA I\nHypermethylation", "HLA I\nHypomethylation"),c(nrow(subset(lowHLAIIupreg,lowHLAIIupreg$hypermethylation == T)),nrow(subset(lowHLAIIupreg,lowHLAIIupreg$hypomethylation == T))))
colnames(piechart_df2)=c("Methyl_state","num")
piechart_df2 <- piechart_df2 %>% 
  arrange(desc(Methyl_state)) %>%
  mutate(prop = num / sum(piechart_df2$num) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
piechart_df2$label=paste(piechart_df2$Methyl_state,'\n',piechart_df2$num,sep='')
# Image: 650 * 470
ggplot(piechart_df2, aes(x="", y=prop, fill=Methyl_state)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none", plot.title = element_text(size = 15)) +
  geom_text(aes(y = ypos, label = label), color = "white", size=5) +
  scale_fill_brewer(palette="Set1") + ggtitle("low HLA-II fold change (bottom 50%)") +
  theme(plot.title = element_text(hjust = 0.5))

df_highUpregII=data.frame('demethylation'= piechart_df1$num[1],'methylation'= piechart_df1$num[2])
df_lowUpregII=data.frame('demethylation'= piechart_df2$num[1],'methylation'= piechart_df2$num[2])
df_combined=rbind(df_highUpregII,df_lowUpregII)
rownames(df_combined)=c('upreg_HLAII','no.upreg_HLAII')
fisher.test(df_combined,alternative = "two.sided")


# HLA LOH section. Related to Figure 6 and S6. See HLALOH_plots.R


