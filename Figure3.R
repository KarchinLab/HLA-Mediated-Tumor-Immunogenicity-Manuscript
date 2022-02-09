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


# Figure 3E
fisher.stats=function(df) {
  df1=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C1' & df$pathway == 1)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C1' & df$pathway == 1)))       
  df2=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C1'  & df$pathway == 0)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C1' & df$pathway == 0)))       
  df_combined=rbind(df1,df2)
  stats1=fisher.test(df_combined,alternative = "two.sided")
  
  df1=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C2' & df$pathway == 1)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C2' & df$pathway == 1)))       
  df2=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C2'  & df$pathway == 0)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C2' & df$pathway == 0)))       
  df_combined=rbind(df1,df2)
  stats2=fisher.test(df_combined,alternative = "two.sided")
  
  df1=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C3' & df$pathway == 1)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C3' & df$pathway == 1)))       
  df2=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C3'  & df$pathway == 0)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C3' & df$pathway == 0)))       
  df_combined=rbind(df1,df2)
  stats3=fisher.test(df_combined,alternative = "two.sided")
  
  df1=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C4' & df$pathway == 1)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C4' & df$pathway == 1)))       
  df2=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C4'  & df$pathway == 0)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C4' & df$pathway == 0)))       
  df_combined=rbind(df1,df2)
  stats4=fisher.test(df_combined,alternative = "two.sided")
  
  df1=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C5' & df$pathway == 1)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C5' & df$pathway == 1)))       
  df2=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C5'  & df$pathway == 0)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C5' & df$pathway == 0)))       
  df_combined=rbind(df1,df2)
  stats5=fisher.test(df_combined,alternative = "two.sided")
  
  df1=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C6' & df$pathway == 1)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C6' & df$pathway == 1)))       
  df2=data.frame('subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based == 'C6'  & df$pathway == 0)),
                 'non.subtype'= nrow(subset(df, df$Subtype_Immune_Model_Based != 'C6' & df$pathway == 0)))       
  df_combined=rbind(df1,df2)
  stats6=fisher.test(df_combined,alternative = "two.sided")
  
  return(list(c(stats1$estimate,stats2$estimate,stats3$estimate,stats4$estimate,stats5$estimate,stats6$estimate),
              c(stats1$conf.int[1],stats2$conf.int[1],stats3$conf.int[1],stats4$conf.int[1],stats5$conf.int[1],stats6$conf.int[1]),
              c(stats1$conf.int[2],stats2$conf.int[2],stats3$conf.int[2],stats4$conf.int[2],stats5$conf.int[2],stats6$conf.int[2]),
              c(stats1$p.value,stats2$p.value,stats3$p.value,stats4$p.value,stats5$p.value,stats6$p.value)))
}
df=subset(df_immuneLandscape,!is.na(df_immuneLandscape$Subtype_Immune_Model_Based) & df_immuneLandscape$Subtype_Immune_Model_Based != '')

immuneSubtypes=c('C1','C2','C3','C4','C5','C6')
df$pathway=df$pathway_Cell.Cycle
stats.cell.cycle=fisher.stats(df)
stats.cell.cycle <- data.frame(
  Pathway='Cell Cycle',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.cell.cycle[[1]],
  boxCILow = stats.cell.cycle[[2]],
  boxCIHigh = stats.cell.cycle[[3]],
  pVal = stats.cell.cycle[[4]]
)

df$pathway=df$pathway_HIPPO
stats.hippo=fisher.stats(df)
stats.hippo <- data.frame(
  Pathway='Hippo',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.hippo[[1]],
  boxCILow = stats.hippo[[2]],
  boxCIHigh = stats.hippo[[3]],
  pVal = stats.hippo[[4]]
)

df$pathway=df$pathway_MYC
stats.myc=fisher.stats(df)
stats.myc <- data.frame(
  Pathway='Myc',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.myc[[1]],
  boxCILow = stats.myc[[2]],
  boxCIHigh = stats.myc[[3]],
  pVal = stats.myc[[4]]
)

df$pathway=df$pathway_NOTCH
stats.notch=fisher.stats(df)
stats.notch <- data.frame(
  Pathway='Notch',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.notch[[1]],
  boxCILow = stats.notch[[2]],
  boxCIHigh = stats.notch[[3]],
  pVal = stats.notch[[4]]
)

df$pathway=df$pathway_NRF2
stats.nrf2=fisher.stats(df)
stats.nrf2 <- data.frame(
  Pathway='Nrf2',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.nrf2[[1]],
  boxCILow = stats.nrf2[[2]],
  boxCIHigh = stats.nrf2[[3]],
  pVal = stats.nrf2[[4]]
)

df$pathway=df$pathway_PI3K
stats.pi3k=fisher.stats(df)
stats.pi3k <- data.frame(
  Pathway='PI3K',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.pi3k[[1]],
  boxCILow = stats.pi3k[[2]],
  boxCIHigh = stats.pi3k[[3]],
  pVal = stats.pi3k[[4]]
)

df$pathway=df$pathway_RTK.RAS
stats.ras=fisher.stats(df)
stats.ras <- data.frame(
  Pathway='RTK/RAS',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.ras[[1]],
  boxCILow = stats.ras[[2]],
  boxCIHigh = stats.ras[[3]],
  pVal = stats.ras[[4]]
)

df$pathway=df$pathway_TGF.Beta
stats.tgfb=fisher.stats(df)
stats.tgfb <- data.frame(
  Pathway='TGFβ',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.tgfb[[1]],
  boxCILow = stats.tgfb[[2]],
  boxCIHigh = stats.tgfb[[3]],
  pVal = stats.tgfb[[4]]
)

df$pathway=df$pathway_TP53
stats.tp53=fisher.stats(df)
stats.tp53 <- data.frame(
  Pathway='TP53',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.tp53[[1]],
  boxCILow = stats.tp53[[2]],
  boxCIHigh = stats.tp53[[3]],
  pVal = stats.tp53[[4]]
)

df$pathway=df$pathway_WNT
stats.wnt=fisher.stats(df)
stats.wnt <- data.frame(
  Pathway='Wnt',
  'Immune Subtype'=immuneSubtypes,
  boxOdds = stats.wnt[[1]],
  boxCILow = stats.wnt[[2]],
  boxCIHigh = stats.wnt[[3]],
  pVal = stats.wnt[[4]]
)

stats.combined=rbind(stats.cell.cycle,stats.hippo,stats.myc,stats.notch,stats.nrf2,stats.pi3k,stats.ras,stats.tgfb,stats.tp53,stats.wnt)
stats.combined$logOR=log10(stats.combined$boxOdds)
plot.data <- stats.combined[,c('Pathway','Immune.Subtype','logOR')]
plot.pVal <- stats.combined[,c('Pathway','Immune.Subtype','pVal')]
qVal=p.adjust(plot.pVal$pVal,method="BH") 
stars_vec=stars.pval(qVal)
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
plot.data$stars=stars_vec
plot.data$Pathway=factor(plot.data$Pathway,rev(levels(factor(plot.data$Pathway))))

p <- ggplot(aes(x=Pathway, y=Immune.Subtype, fill=logOR), data=plot.data)
p + geom_tile() + scale_fill_gradient2(low='darkorchid4',mid='white',high='forestgreen') + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="log10 (OR)") + theme_bw() + coord_flip() + 
  xlab('Oncogenic Pathway') + ylab('Immune Subtype') +
  theme(axis.text.x=element_text(angle = -55, hjust = 0),legend.title=element_text(size=8)) +
  guides(fill = guide_colourbar(barwidth = 1,
                                barheight = 10))
ggsave("DE_HLAII_pathways.png",height = 5,width = 5.5,dpi = 300)



# Figure S3A-B
# pairwise boxplots
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


# Figure S3C
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

