# Figure 3A-B
HLAupregulation_df=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based!='')
HLAupregulation_df=HLAupregulation_df[which(!is.na(HLAupregulation_df$norm_T_HLAI) & !is.na(HLAupregulation_df$norm_N_HLAI)),]
immuneSubtype=unique(HLAupregulation_df$Subtype_Immune_Model_Based)
keep_rows=c()
for (type in immuneSubtype) {
  df=subset(HLAupregulation_df,HLAupregulation_df$Subtype_Immune_Model_Based == type)
  if (nrow(df) >= 5) {
    keep_rows=append(keep_rows,rownames(df))
  } else {
    message(type)
  }
}
HLAupregulation_df=HLAupregulation_df[keep_rows,]

HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C1']='C1 (n=171)\nwound-healing'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C2']='C2 (n=144)\nIFN-γ dominant'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C3']='C3 (n=256)\ninflammatory'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C4']='C4 (n=67)\nlymphocyte-depleted'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C5']='C5 (n=0)\nimmunologically quiet'
HLAupregulation_df$immuneSubtype[HLAupregulation_df$Subtype_Immune_Model_Based=='C6']='C6 (n=18)\nTGF-β dominant'

HLAupregulation_df$tumor = HLAupregulation_df$norm_T_HLAI
HLAupregulation_df$normal = HLAupregulation_df$norm_N_HLAI
dat.m1 <- melt(HLAupregulation_df,id.vars='immuneSubtype', measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable
# FDR: q values
pvalues=ddply(dat.m1, .(immuneSubtype), function(x) { wilcox.test(value~Tissue, data=x, paired=T)$p.value })
pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH") 
stars_vec=stars.pval(qVal)
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
y<-rep(c(15.15),each=5)
x<-c(0.62, 1.67, 2.71, 3.73, 4.75)-0.05

# Image: 900 * 500
ggplot(dat.m1, aes(x=immuneSubtype, y=value, fill=Tissue)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
  stat_compare_means(aes(group = Tissue, label = ifelse( p < 2e-16, "p < 2e-16", 
                                                         paste0("p = ", ..p.format..))), paired=T, label.x = 1, label.y =15, size=3) +
  ylab('HLA Class I Expression') + xlab('') +
  theme_bw() +
  annotate(geom = "text", y = c(y), x = c(x), label = c(stars_vec), size = 3.5, hjust = 0)
ggsave("immuneSubtype_HLAI.png",height = 4.5,width =8,dpi = 300)


HLAupregulation_df$tumor = HLAupregulation_df$norm_HLAII
HLAupregulation_df$normal = HLAupregulation_df$norm_HLAII_N
dat.m1 <- melt(HLAupregulation_df,id.vars='immuneSubtype', measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable
# FDR: q values
pvalues=ddply(dat.m1, .(immuneSubtype), function(x) { wilcox.test(value~Tissue, data=x, paired=T)$p.value })
pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=stars.pval(qVal)
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
y<-rep(c(15.2),each=5)
x<-c(0.69, 2.0, 2.70, 3.76, 4.74)-0.07

# Image: 770 * 450
ggplot(dat.m1, aes(x=immuneSubtype, y=value, fill=Tissue)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
  stat_compare_means(aes(group = Tissue, label = paste0("p = ", ..p.format..)), paired=T, label.x = 1, label.y = 15, size=3) +
  xlab('') + ylab('HLA Class II Expression') +
  theme_bw() +
  annotate(geom = "text", y = c(y), x = c(x), label = c(stars_vec), size = 3.5, hjust = 0)
ggsave("immuneSubtype_HLAII.png",height = 4.5,width =8,dpi = 300)


# Figure 3C: HR plot
HLAupregulation_df=df_immuneLandscape
HLAupregulation_df=subset(HLAupregulation_df,!is.na(HLAupregulation_df$DE_HLAI))

HLAupregulation_df$highCYT=HLAupregulation_df$CYT > unname(quantile(HLAupregulation_df$CYT,na.rm=TRUE)["50%"])
HLAupregulation_df$high_HLAI_fold=HLAupregulation_df$DE_HLAI > unname(quantile(HLAupregulation_df$DE_HLAI,na.rm=TRUE)["50%"])
HLAupregulation_df$high_HLAII_fold=HLAupregulation_df$DE_HLAII > unname(quantile(HLAupregulation_df$DE_HLAII,na.rm=TRUE)["50%"])

explanatory = c("highCYT","high_HLAI_fold","high_HLAII_fold")
dependent = "Surv(PFI_time_1, PFI_1)"

png(filename = "HR.png",width = 2300, height = 1000,res = 300)
HLAupregulation_df %>%
  hr_plot(dependent, explanatory, remove_ref = T,
          column_space = c(-0.4, -0.2, 0.8, 1), dependent_label = "", prefix = "TCGA", suffix = "", table_text_size = 5, title_text_size = 18)
dev.off()


# Figure 3D
dat.m1 <- melt(df_immuneLandscape,id.vars=c('norm_T_HLAI','norm_HLAII','DE_HLAI','DE_HLAII'), 
               measure.vars=c('pathway_Cell.Cycle','pathway_HIPPO','pathway_MYC','pathway_NOTCH',
                              'pathway_NRF2','pathway_PI3K','pathway_RTK.RAS','pathway_TGF.Beta',
                              'pathway_TP53','pathway_WNT'))
dat.m1=subset(dat.m1,!is.na(dat.m1$value))
dat.m1$value=as.character(dat.m1$value)
dat.m1$Pathway=NA
dat.m1$Pathway[which(dat.m1$value == '1')]='Altered'
dat.m1$Pathway[which(dat.m1$value == '0')]='Not Altered'
dat.m1$pathway_name=substring(dat.m1$variable,regexpr('_',dat.m1$variable)+1)
dat.m1$pathway_name[which(dat.m1$pathway_name == 'TGF.Beta')]='TGFβ'
dat.m1$pathway_name[which(dat.m1$pathway_name == 'RTK.RAS')]='RTK/RAS'
dat.m1$pathway_name[which(dat.m1$pathway_name == 'Cell.Cycle')]='Cell Cycle'
dat.m1$pathway_name=factor(dat.m1$pathway_name,rev(levels(factor(dat.m1$pathway_name))))

pvalues=ddply(dat.m1, .(variable), function(x) { wilcox.test(DE_HLAII~value, data=x, paired=FALSE)$p.value })
pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=rev(stars.pval(qVal))
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
caliberate=0.05
y<-rep(c(2.185+5.5-caliberate),each=length(qVal))
y[which(str_length(stars_vec) == 1)]=2.21+5.7-caliberate
y[which(str_length(stars_vec) == 3)]=2.155+5.4-caliberate
x<-rep(1:length(qVal))+0.1

y2<-rep(c(2.2+6),each=length(qVal))
x2<-rep(1:length(qVal))+0.01
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

ggplot(dat.m1, aes(x=pathway_name, y=dat.m1$DE_HLAII, fill=Pathway)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("steelblue3", "thistle3")) +
  xlab('Oncogenic Pathway') + ylab('HLA II Log2FC') + ylim(c(-5,12)) +
  coord_flip() + theme_bw() + 
  annotate(geom = "text", y = c(y,y2), x = c(x,x2), label = c(stars_vec,pVal_format), size = 3.5, hjust = 0)
ggsave("pathway_HLAIIDE.png",height = 5,width =6,dpi = 300)




# Figure S3A-B
overview_df=df_immuneLandscape[which(!is.na(df_immuneLandscape$norm_T_HLAI) & !is.na(df_immuneLandscape$norm_N_HLAI)),]
tumortypes=unique(overview_df$Study)
keep_rows=c()
overview_df$rankI=NA
for (type in tumortypes) {
  df=subset(overview_df,overview_df$Study == type)
  if (nrow(df) >= 5) {
    keep_rows=append(keep_rows,rownames(df))
    overview_df[rownames(df),]$rankI=median(df$norm_T_HLAI)
  } else {
    message(type)
  }
}
overview_df=overview_df[keep_rows,]

# create new labels
overview_df$label=NA
tumortypes=unique(overview_df$Study)
for (type in tumortypes) {
  overview_df$label[which(overview_df$Study == type)]=paste0(type,' ','(n=',sum(overview_df$Study == type),')')
}

overview_df$tumor = overview_df$norm_T_HLAI
overview_df$normal = overview_df$norm_N_HLAI
dat.m1 <- melt(overview_df,id.vars=c('label','rankI'), measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable

pvalues=ddply(dat.m1, .(label), function(x) { wilcox.test(value~Tissue, data=x, paired=T)$p.value })
pvalues$rankI=NA
for (i in (1:nrow(pvalues))) {
  cancerType=pvalues$label[i]
  index=which(dat.m1$label == cancerType)[1]
  pvalues$rankI[i]=dat.m1$rankI[index]
}
pvalues2=pvalues %>% arrange(-rankI)
pVal_vec=pvalues2$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=rev(stars.pval(qVal))
stars_vec=replace(stars_vec, stars_vec=='.', ' ')

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

ggplot(dat.m1, aes(x=reorder(label, rankI), y=value, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("steelblue3", "thistle3")) +
  xlab('') + ylab('HLA Class I Expression') + ylim(c(5.7,16)) +
  coord_flip() + theme_bw()
ggsave("HLAI-overall.png",height = 8.5,width =7.5,dpi = 300)


# Figure 1C
overview_df=df_immuneLandscape[which(!is.na(df_immuneLandscape$norm_HLAII) & !is.na(df_immuneLandscape$norm_HLAII_N)),]
tumortypes=unique(overview_df$Study)
keep_rows=c()
overview_df$rankI=NA
for (type in tumortypes) {
  df=subset(overview_df,overview_df$Study == type)
  if (nrow(df) >= 5) {
    keep_rows=append(keep_rows,rownames(df))
    overview_df[rownames(df),]$rankI=median(df$norm_HLAII)
  } else {
    message(type)
  }
}
overview_df=overview_df[keep_rows,]

# create new labels
overview_df$label=NA
tumortypes=unique(overview_df$Study)
for (type in tumortypes) {
  overview_df$label[which(overview_df$Study == type)]=paste0(type,' ','(n=',sum(overview_df$Study == type),')')
}

overview_df$tumor = overview_df$norm_HLAII
overview_df$normal = overview_df$norm_HLAII_N
dat.m1 <- melt(overview_df,id.vars=c('label','rankI'), measure.vars=c('tumor','normal'))
dat.m1$Tissue = dat.m1$variable

pvalues=ddply(dat.m1, .(label), function(x) { wilcox.test(value~Tissue, data=x, paired=T)$p.value })
pvalues$rankI=NA
pvalues$rankI=NA
for (i in (1:nrow(pvalues))) {
  cancerType=pvalues$label[i]
  index=which(dat.m1$label == cancerType)[1]
  pvalues$rankI[i]=dat.m1$rankI[index]
}
pvalues2=pvalues %>% arrange(-rankI)
pVal_vec=pvalues2$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=rev(stars.pval(qVal))
stars_vec=replace(stars_vec, stars_vec=='.', ' ')

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

ggplot(dat.m1, aes(x=reorder(label, rankI), y=value, fill=Tissue)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("steelblue3", "thistle3")) +
  xlab('') + ylab('HLA Class II Expression') + ylim(c(2,16)) +
  coord_flip() + theme_bw()
ggsave("HLAII-overall.png",height = 8.5,width =7.5,dpi = 300)



# Figure S3C-D
HLAupregulation_df=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based!='')
HLAupregulation_df=HLAupregulation_df[which(!is.na(HLAupregulation_df$norm_T_HLAI) & !is.na(HLAupregulation_df$norm_N_HLAI)),]

keep_rows=c()
tumorTypes=unique(HLAupregulation_df$Study)
for (tumorTypes in tumorTypes) {
  temp=subset(HLAupregulation_df,HLAupregulation_df$Study == tumorTypes)
  immuneSubtype=unique(temp$Subtype_Immune_Model_Based)
  
  for (type in immuneSubtype) {
    df=subset(temp,temp$Subtype_Immune_Model_Based == type)
    if (nrow(df) >= 5) {
      keep_rows=append(keep_rows,rownames(df))
    }
  }
}
HLAupregulation_facet_df=HLAupregulation_df[keep_rows,]


# HLA I 
for (TumorType in sort(tumorTypes)) {
  tumorType_df=subset(HLAupregulation_facet_df,HLAupregulation_facet_df$Study == TumorType)
  
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C1']=paste0('C1 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C1'),')\nwound-healing')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C2']=paste0('C2 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C2'),')\nIFN-γ dominant')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C3']=paste0('C3 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C3'),')\ninflammatory')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C4']=paste0('C4 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C4'),')\nlymphocyte-depleted')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C5']=paste0('C5 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C5'),')\nimmunologically quiet')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C6']=paste0('C6 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C6'),')\nTGF-β dominant')
  
  tumorType_df$tumor = tumorType_df$norm_T_HLAI
  tumorType_df$normal = tumorType_df$norm_N_HLAI
  dat.m1 <- melt(tumorType_df,id.vars=c('immuneSubtype','Study'), measure.vars=c('tumor','normal'))
  dat.m1$Tissue = dat.m1$variable
  
  pVal_format=subset(pVal_df,pVal_df$tumorType == TumorType)$pVal_format
  stars_vec=subset(pVal_df,pVal_df$tumorType == TumorType)$stars
  
  image_height=3.5
  image_width=7

  # Image: 900 * 500
  ggplot(dat.m1, aes(x=immuneSubtype, y=value, fill=Tissue)) +
    geom_boxplot() + 
    scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
    ylab('') + xlab('') +
    theme_bw() +
    facet_grid(cols = vars(Study),scales = "free") + theme(legend.position = "none")
  ggsave(paste0(TumorType,"immuneSubtype_HLAI.png"),height = image_height,width =image_width,dpi = 300)
}



# HLA II
for (TumorType in sort(tumorTypes)) {
  tumorType_df=subset(HLAupregulation_facet_df,HLAupregulation_facet_df$Study == TumorType)
  
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C1']=paste0('C1 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C1'),')\nwound-healing')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C2']=paste0('C2 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C2'),')\nIFN-γ dominant')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C3']=paste0('C3 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C3'),')\ninflammatory')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C4']=paste0('C4 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C4'),')\nlymphocyte-depleted')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C5']=paste0('C5 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C5'),')\nimmunologically quiet')
  tumorType_df$immuneSubtype[tumorType_df$Subtype_Immune_Model_Based=='C6']=paste0('C6 (n=',sum(tumorType_df$Subtype_Immune_Model_Based=='C6'),')\nTGF-β dominant')
  
  tumorType_df$tumor = tumorType_df$norm_HLAII
  tumorType_df$normal = tumorType_df$norm_HLAII_N
  dat.m1 <- melt(tumorType_df,id.vars=c('immuneSubtype','Study'), measure.vars=c('tumor','normal'))
  dat.m1$Tissue = dat.m1$variable
  
  pVal_format=subset(pVal_df,pVal_df$tumorType == TumorType)$pVal_format
  stars_vec=subset(pVal_df,pVal_df$tumorType == TumorType)$stars
  
  image_height=3.5
  image_width=7
  
  # Image: 900 * 500
  ggplot(dat.m1, aes(x=immuneSubtype, y=value, fill=Tissue)) +
    geom_boxplot() + 
    scale_fill_manual(values = c("firebrick", "deepskyblue3")) +
    ylab('') + xlab('') +
    theme_bw() +
    annotate(geom = "text", x = c(x,x2), y = c(y,y2), label = c(pVal_format,stars_vec), size = 4, hjust = 0) +
    facet_grid(cols = vars(Study),scales = "free") + theme(legend.position = "none")
  ggsave(paste0(TumorType,"immuneSubtype_HLAII.png"),height = image_height,width =image_width,dpi = 300)
}


# Figure S3E-F
HLAupregulation_df$high_HLAI_fold=HLAupregulation_df$DE_HLAI > unname(quantile(HLAupregulation_df$DE_HLAI,na.rm=TRUE)["50%"])
HLAupregulation_df$high_HLAII_fold=HLAupregulation_df$DE_HLAII > unname(quantile(HLAupregulation_df$DE_HLAII,na.rm=TRUE)["50%"])

HLAupregulation_df$Groups=NA
HLAupregulation_df$Groups[HLAupregulation_df$high_HLAI_fold == T & HLAupregulation_df$high_HLAII_fold ==T]='high HLAI high HLAII'
HLAupregulation_df$Groups[HLAupregulation_df$high_HLAI_fold == T & HLAupregulation_df$high_HLAII_fold ==F]='high HLAI low HLAII'
HLAupregulation_df$Groups[HLAupregulation_df$high_HLAI_fold == F & HLAupregulation_df$high_HLAII_fold ==T]='low HLAI high HLAII'
HLAupregulation_df$Groups[HLAupregulation_df$high_HLAI_fold == F & HLAupregulation_df$high_HLAII_fold ==F]='low HLAI low HLAII'
HLAupregulation_df=subset(HLAupregulation_df,!is.na(HLAupregulation_df$Groups))
HLAupregulation_df$Groups = factor(HLAupregulation_df$Groups, rev(levels(factor(HLAupregulation_df$Groups))))
my_comparisons <- list( c('high HLAI high HLAII', 'low HLAI low HLAII'), c('high HLAI high HLAII', 'high HLAI low HLAII'), c('high HLAI low HLAII', 'low HLAI high HLAII'), c('low HLAI low HLAII', 'low HLAI high HLAII') )
ggplot(HLAupregulation_df, aes(x= Groups, y=Th1.cells, fill=Groups)) + 
  geom_boxplot() + scale_fill_brewer(palette="PuBuGn") + 
  stat_compare_means(label.x = 1, label.y = 1200) + 
  stat_compare_means(comparisons = my_comparisons) + 
  xlab('HLA-I vs. HLA-II fold change') + ylab('Th1 signature score') +
  theme_bw() + theme(legend.position = "none")
ggsave("Th1_HLAfolds.png",height = 4,width =6,dpi = 300)


ggplot(HLAupregulation_df, aes(x= Groups, y=T.cells.CD8, fill=Groups)) + 
  geom_boxplot() + scale_fill_brewer(palette="PuBuGn") + 
  stat_compare_means(label.x = 1, label.y = 0.8) +
  stat_compare_means(comparisons = my_comparisons) + 
  xlab('HLA-I vs. HLA-II fold change') + ylab('CD8 T cells') +
  theme_bw() + theme(legend.position = "none")
ggsave("CD8T_HLAfolds.png",height = 4,width =6,dpi = 300)


# Figure S3G
dat.m1 <- melt(df_immuneLandscape,id.vars=c('norm_T_HLAI','norm_HLAII','DE_HLAI','DE_HLAII'), 
               measure.vars=c('pathway_Cell.Cycle','pathway_HIPPO','pathway_MYC','pathway_NOTCH',
                              'pathway_NRF2','pathway_PI3K','pathway_RTK.RAS','pathway_TGF.Beta',
                              'pathway_TP53','pathway_WNT'))
dat.m1=subset(dat.m1,!is.na(dat.m1$value))
dat.m1$value=as.character(dat.m1$value)
dat.m1$Pathway=NA
dat.m1$Pathway[which(dat.m1$value == '1')]='Altered'
dat.m1$Pathway[which(dat.m1$value == '0')]='Not Altered'
dat.m1$pathway_name=substring(dat.m1$variable,regexpr('_',dat.m1$variable)+1)
dat.m1$pathway_name[which(dat.m1$pathway_name == 'TGF.Beta')]='TGFβ'
dat.m1$pathway_name[which(dat.m1$pathway_name == 'RTK.RAS')]='RTK/RAS'
dat.m1$pathway_name[which(dat.m1$pathway_name == 'Cell.Cycle')]='Cell Cycle'
dat.m1$pathway_name=factor(dat.m1$pathway_name,rev(levels(factor(dat.m1$pathway_name))))

pvalues=ddply(dat.m1, .(variable), function(x) { wilcox.test(DE_HLAI~value, data=x, paired=FALSE)$p.value })
pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH")
stars_vec=rev(stars.pval(qVal))
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
caliberate=0.05
y<-rep(c(2.185+5.5-caliberate),each=length(qVal))
y[which(str_length(stars_vec) == 1)]=2.21+5.7-caliberate
y[which(str_length(stars_vec) == 3)]=2.155+5.4-caliberate
x<-rep(1:length(qVal))+0.1

y2<-rep(c(2.2+6),each=length(qVal))
x2<-rep(1:length(qVal))+0.01
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

ggplot(dat.m1, aes(x=pathway_name, y=dat.m1$DE_HLAI, fill=Pathway)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("steelblue3", "thistle3")) +
  xlab('Oncogenic Pathway') + ylab('HLA I Log2FC') + ylim(c(-5,12)) +
  coord_flip() + theme_bw() + 
  annotate(geom = "text", y = c(y,y2), x = c(x,x2), label = c(stars_vec,pVal_format), size = 3.5, hjust = 0)
ggsave("pathway_HLAIDE.png",height = 5,width =6,dpi = 300)

