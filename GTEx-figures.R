df_gtex=read.xlsx("/Users/katana/PycharmProjects/Data-revision/Supplementary-Tables/Supplementary Table 2.xlsx","HLA/B2M Expression (GTEx)")
# starts from here
### TCGA tumor vs. GTEx normal
gtex_kidney=subset(df_gtex, df_gtex$Primary.Site == 'Kidney')
kidney_df1=data.frame('Study'=rep(c('KIRC'),each=nrow(gtex_kidney)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_kidney)),'value'=gtex_kidney$geom_HLAI)
kidney_df2=data.frame('Study'=rep(c('KIRP'),each=nrow(gtex_kidney)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_kidney)),'value'=gtex_kidney$geom_HLAI)
kidney_df3=data.frame('Study'=rep(c('KICH'),each=nrow(gtex_kidney)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_kidney)),'value'=gtex_kidney$geom_HLAI)

gtex_cervix=subset(df_gtex, df_gtex$Primary.Site == 'Cervix Uteri')
cervix_df=data.frame('Study'=rep(c('CESC'),each=nrow(gtex_cervix)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_cervix)),'value'=gtex_cervix$geom_HLAI)

gtex_lung=subset(df_gtex, df_gtex$Primary.Site == 'Lung')
lung_df1=data.frame('Study'=rep(c('LUAD'),each=nrow(gtex_lung)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_lung)),'value'=gtex_lung$geom_HLAI)
lung_df2=data.frame('Study'=rep(c('LUSC'),each=nrow(gtex_lung)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_lung)),'value'=gtex_lung$geom_HLAI)

gtex_skin=subset(df_gtex, df_gtex$Primary.Site == 'Skin')
skin_df=data.frame('Study'=rep(c('SKCM'),each=nrow(gtex_skin)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_skin)),'value'=gtex_skin$geom_HLAI)

gtex_pancreas=subset(df_gtex, df_gtex$Primary.Site == 'Pancreas')
pancreas_df=data.frame('Study'=rep(c('PAAD'),each=nrow(gtex_pancreas)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_pancreas)),'value'=gtex_pancreas$geom_HLAI)

gtex_colon=subset(df_gtex, df_gtex$Primary.Site == 'Colon')
colon_df=data.frame('Study'=rep(c('COAD'),each=nrow(gtex_colon)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_colon)),'value'=gtex_colon$geom_HLAI)

gtex_thyroid=subset(df_gtex, df_gtex$Primary.Site == 'Thyroid')
thyroid_df=data.frame('Study'=rep(c('THCA'),each=nrow(gtex_thyroid)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_thyroid)),'value'=gtex_thyroid$geom_HLAI)

gtex_uterus=subset(df_gtex, df_gtex$Primary.Site == 'Uterus')
uterus_df=data.frame('Study'=rep(c('UCEC'),each=nrow(gtex_uterus)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_uterus)),'value'=gtex_uterus$geom_HLAI)

gtex_bladder=subset(df_gtex, df_gtex$Primary.Site == 'Bladder')
bladder_df=data.frame('Study'=rep(c('BLCA'),each=nrow(gtex_bladder)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_bladder)),'value'=gtex_bladder$geom_HLAI)

gtex_stomach=subset(df_gtex, df_gtex$Primary.Site == 'Stomach')
stomach_df=data.frame('Study'=rep(c('STAD'),each=nrow(gtex_stomach)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_stomach)),'value'=gtex_stomach$geom_HLAI)

gtex_nerve=subset(df_gtex, df_gtex$Primary.Site == 'Nerve')
nerve_df=data.frame('Study'=rep(c('PCPG'),each=nrow(gtex_nerve)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_nerve)),'value'=gtex_nerve$geom_HLAI)

gtex_liver=subset(df_gtex, df_gtex$Primary.Site == 'Liver')
liver_df=data.frame('Study'=rep(c('LIHC'),each=nrow(gtex_liver)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_liver)),'value'=gtex_liver$geom_HLAI)

gtex_breast=subset(df_gtex, df_gtex$Primary.Site == 'Breast')
breast_df=data.frame('Study'=rep(c('BRCA'),each=nrow(gtex_breast)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_breast)),'value'=gtex_breast$geom_HLAI)

gtex_esophagus=subset(df_gtex, df_gtex$Primary.Site == 'Esophagus')
esophagus_df=data.frame('Study'=rep(c('ESCA'),each=nrow(gtex_esophagus)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_esophagus)),'value'=gtex_esophagus$geom_HLAI)

gtex_brain=subset(df_gtex, df_gtex$Primary.Site == 'Brain')
brain_df=data.frame('Study'=rep(c('GBM'),each=nrow(gtex_brain)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_brain)),'value'=gtex_brain$geom_HLAI)

gtex_prostate=subset(df_gtex, df_gtex$Primary.Site == 'Prostate')
prostate_df=data.frame('Study'=rep(c('PRAD'),each=nrow(gtex_prostate)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_prostate)),'value'=gtex_prostate$geom_HLAI)

overview_df = subset(df_immuneLandscape, df_immuneLandscape$Study!='DLBC' & df_immuneLandscape$Study!='MESO' & df_immuneLandscape$Study!='TGCT' & df_immuneLandscape$Study!='OV' & df_immuneLandscape$Study!='LAML' & df_immuneLandscape$Study!='UCS' & df_immuneLandscape$Study!='UVM' & df_immuneLandscape$Study!='ACC' & df_immuneLandscape$Study!='LGG')
overview_df$`TCGA-tumor` = overview_df$geom_HLAI
overview_df_filtered=subset(overview_df,overview_df$Study != 'THYM' & overview_df$Study != 'SARC' & overview_df$Study != 'READ' & overview_df$Study != 'HNSC' & overview_df$Study != 'CHOL')
dat.m1 <- melt(overview_df_filtered,id.vars='Study', measure.vars=c('TCGA-tumor'))
dat.m1=rbind(dat.m1,kidney_df1,kidney_df2,kidney_df3,cervix_df,lung_df1,lung_df2,skin_df,pancreas_df,
             colon_df,thyroid_df,uterus_df,bladder_df,stomach_df,nerve_df,liver_df,breast_df,esophagus_df,
             brain_df,prostate_df)
dat.m1$Source = dat.m1$variable
dat.m1$Source=factor(dat.m1$Source)
dat.m1=dat.m1 %>% filter(!is.na(dat.m1$Source))
# FDR: q values
pvalues=ddply(dat.m1, .(Study), function(x) { wilcox.test(value~Source, data=x, paired=FALSE)$p.value })

pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH")
# fdr_df=fdrtool(pVal_vec,statistic = 'pvalue')
# qVal=fdr_df$qval
# stars_vec=rev(stars.pval(qVal))
stars_vec=stars.pval(qVal)
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
calibrate=0.1
y<-rep(c(14.7-calibrate),each=length(stars_vec))
y[which(str_length(stars_vec) == 1)]=14.9-calibrate
y[which(str_length(stars_vec) == 3)]=14.5-calibrate
x<-rep(1:length(stars_vec))+0.15
ggplot(dat.m1, aes(x=Study, y=value, fill=Source)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("firebrick", "slategray3")) + 
  stat_compare_means(aes(group = Source, label = paste0("p = ", ..p.format..)), label.y=16) +
  xlab('') + ylab('HL                                                                                                                                                                                                                    A Class I Expression') +
  coord_flip(ylim=c(5,17)) + theme_bw() +
  annotate(geom = "text", y = c(y), x = c(x), label = c(stars_vec), size = 3.5, hjust = 0)







gtex_kidney=subset(df_gtex, df_gtex$Primary.Site == 'Kidney')
kidney_df1=data.frame('Study'=rep(c('KIRC'),each=nrow(gtex_kidney)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_kidney)),'value'=gtex_kidney$geom_HLAII)
kidney_df2=data.frame('Study'=rep(c('KIRP'),each=nrow(gtex_kidney)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_kidney)),'value'=gtex_kidney$geom_HLAII)
kidney_df3=data.frame('Study'=rep(c('KICH'),each=nrow(gtex_kidney)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_kidney)),'value'=gtex_kidney$geom_HLAII)

gtex_cervix=subset(df_gtex, df_gtex$Primary.Site == 'Cervix Uteri')
cervix_df=data.frame('Study'=rep(c('CESC'),each=nrow(gtex_cervix)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_cervix)),'value'=gtex_cervix$geom_HLAII)

gtex_lung=subset(df_gtex, df_gtex$Primary.Site == 'Lung')
lung_df1=data.frame('Study'=rep(c('LUAD'),each=nrow(gtex_lung)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_lung)),'value'=gtex_lung$geom_HLAII)
lung_df2=data.frame('Study'=rep(c('LUSC'),each=nrow(gtex_lung)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_lung)),'value'=gtex_lung$geom_HLAII)

gtex_skin=subset(df_gtex, df_gtex$Primary.Site == 'Skin')
skin_df=data.frame('Study'=rep(c('SKCM'),each=nrow(gtex_skin)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_skin)),'value'=gtex_skin$geom_HLAII)

gtex_pancreas=subset(df_gtex, df_gtex$Primary.Site == 'Pancreas')
pancreas_df=data.frame('Study'=rep(c('PAAD'),each=nrow(gtex_pancreas)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_pancreas)),'value'=gtex_pancreas$geom_HLAII)

gtex_colon=subset(df_gtex, df_gtex$Primary.Site == 'Colon')
colon_df=data.frame('Study'=rep(c('COAD'),each=nrow(gtex_colon)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_colon)),'value'=gtex_colon$geom_HLAII)

gtex_thyroid=subset(df_gtex, df_gtex$Primary.Site == 'Thyroid')
thyroid_df=data.frame('Study'=rep(c('THCA'),each=nrow(gtex_thyroid)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_thyroid)),'value'=gtex_thyroid$geom_HLAII)

gtex_uterus=subset(df_gtex, df_gtex$Primary.Site == 'Uterus')
uterus_df=data.frame('Study'=rep(c('UCEC'),each=nrow(gtex_uterus)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_uterus)),'value'=gtex_uterus$geom_HLAII)

gtex_bladder=subset(df_gtex, df_gtex$Primary.Site == 'Bladder')
bladder_df=data.frame('Study'=rep(c('BLCA'),each=nrow(gtex_bladder)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_bladder)),'value'=gtex_bladder$geom_HLAII)

gtex_stomach=subset(df_gtex, df_gtex$Primary.Site == 'Stomach')
stomach_df=data.frame('Study'=rep(c('STAD'),each=nrow(gtex_stomach)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_stomach)),'value'=gtex_stomach$geom_HLAII)

gtex_nerve=subset(df_gtex, df_gtex$Primary.Site == 'Nerve')
nerve_df=data.frame('Study'=rep(c('PCPG'),each=nrow(gtex_nerve)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_nerve)),'value'=gtex_nerve$geom_HLAII)

gtex_liver=subset(df_gtex, df_gtex$Primary.Site == 'Liver')
liver_df=data.frame('Study'=rep(c('LIHC'),each=nrow(gtex_liver)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_liver)),'value'=gtex_liver$geom_HLAII)

gtex_breast=subset(df_gtex, df_gtex$Primary.Site == 'Breast')
breast_df=data.frame('Study'=rep(c('BRCA'),each=nrow(gtex_breast)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_breast)),'value'=gtex_breast$geom_HLAII)

gtex_esophagus=subset(df_gtex, df_gtex$Primary.Site == 'Esophagus')
esophagus_df=data.frame('Study'=rep(c('ESCA'),each=nrow(gtex_esophagus)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_esophagus)),'value'=gtex_esophagus$geom_HLAII)

gtex_brain=subset(df_gtex, df_gtex$Primary.Site == 'Brain')
brain_df=data.frame('Study'=rep(c('GBM'),each=nrow(gtex_brain)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_brain)),'value'=gtex_brain$geom_HLAII)

gtex_prostate=subset(df_gtex, df_gtex$Primary.Site == 'Prostate')
prostate_df=data.frame('Study'=rep(c('PRAD'),each=nrow(gtex_prostate)),'variable'=rep(c('GTEx-normal'),each=nrow(gtex_prostate)),'value'=gtex_prostate$geom_HLAII)

overview_df$`TCGA-tumor` = overview_df$geom_HLAII
overview_df_filtered=subset(overview_df,overview_df$Study != 'THYM' & overview_df$Study != 'SARC' & overview_df$Study != 'READ' & overview_df$Study != 'HNSC' & overview_df$Study != 'CHOL')
dat.m1 <- melt(overview_df_filtered,id.vars='Study', measure.vars=c('TCGA-tumor'))
dat.m1=rbind(dat.m1,kidney_df1,kidney_df2,kidney_df3,cervix_df,lung_df1,lung_df2,skin_df,pancreas_df,
             colon_df,thyroid_df,uterus_df,bladder_df,stomach_df,nerve_df,liver_df,breast_df,esophagus_df,
             brain_df,prostate_df)
dat.m1$Source = dat.m1$variable
# FDR: q values
pvalues=ddply(dat.m1, .(Study), function(x) { wilcox.test(value~Source, data=x, paired=FALSE)$p.value })

pVal_vec=pvalues$V1
qVal=p.adjust(pVal_vec,method="BH")
# fdr_df=fdrtool(pVal_vec,statistic = 'pvalue')
# qVal=fdr_df$qval
# stars_vec=rev(stars.pval(qVal))
stars_vec=stars.pval(qVal)
stars_vec=replace(stars_vec, stars_vec=='.', ' ')
calibrate=0.8
y<-rep(c(14.8-calibrate),each=length(stars_vec))
y[which(str_length(stars_vec) == 1)]=15.1-calibrate
y[which(str_length(stars_vec) == 3)]=14.5-calibrate
x<-rep(1:length(stars_vec))+0.15
ggplot(dat.m1, aes(x=Study, y=value, fill=Source)) +
  geom_boxplot(outlier.shape=NA) + 
  scale_fill_manual(values = c("firebrick", "slategray3")) + 
  stat_compare_means(aes(group = Source, label = paste0("p = ", ..p.format..)), label.y=16) +
  xlab('') + ylab('HLA Class II Expression') +
  coord_flip(ylim=c(0,17)) + theme_bw() +
  annotate(geom = "text", y = c(y), x = c(x), label = c(stars_vec), size = 3.5, hjust = 0)
















