# Figure 6C
fit <- survfit(Surv(progression_free,progression) ~ Prediction,data = ici_df)
ggsurv=ggsurvplot(fit, data = ici_df,
                  conf.int = F, 
                  pval = T,  
                  risk.table = TRUE,
                  legend = 'none',
                  xlab = '', ylab = 'Progression Free Survival', 
                  ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank()))
png(filename = "predicted_surv3.png",width = 1400, height = 1700,res = 300)
ggsurv
dev.off()
ici_df=ici_df
write.csv(ici_df,'/Users/katana/PycharmProjects/MCB-revision/panModel4_cohort122_HED_predicted_immuneSubytpes.csv')



# Figure 6B
ici_df$BOR3=ici_df$BOR
ici_df$BOR3[which(ici_df$BOR == 'CR' | ici_df$BOR == 'PR')]='CRPR'
features=c('norm_B2M_log','norm_HLA_A_log','norm_HLA_B_log','norm_HLA_C_log','norm_HLA_E_log','norm_HLA_G_log','HLA_DRA_log', 'HLA_DRB1_log', 'HLA_DQA1_log', 'HLA_DQB1_log', 'HLA_DPA1_log', 'HLA_DPB1_log', 'HLA_DMA_log', 'HLA_DMB_log', 'HLA_DOA_log', 'HLA_DOB_log')
df=subset(ici_df,select=(c('patientId','BOR3','Prediction',features)))
C1=subset(df,df$Prediction == 'C1')
C2=subset(df,df$Prediction == 'C2')
C3=subset(df,df$Prediction == 'C3')
C4=subset(df,df$Prediction == 'C4')
df=rbind(C1,C2,C3,C4)
patientList=df$patientId
rownames(df)=patientList
BORList=df$BOR3
immuneSubtypeList=df$Prediction

df=df[4:ncol(df)]
colnames(df)=c('B2M','HLA-A','HLA-B','HLA-C','HLA-E','HLA-G',"HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA",'HLA-DMB',"HLA-DOA",'HLA-DOB')
data=as.matrix(df)


annotation_col <- data.frame(row.names = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-G",'B2M',"HLA-DRA", "HLA-DRB1","HLA-DQA1","HLA-DQB1","HLA-DPA1","HLA-DPB1","HLA-DMA",'HLA-DMB',"HLA-DOA",'HLA-DOB'), 
                             HLA = c(rep("Class I", 6), rep("Class II", 10)))
annotation_row <- data.frame(row.names = patientList, 
                             'BOR' = BORList,
                             'Prediction'=immuneSubtypeList)

annotation_row$BOR=factor(annotation_row$BOR, c('CRPR','SD','PD','MR'))


BOR_colors=rev(magma(4))
names(BOR_colors)=c('CRPR','SD','PD','MR')
prediction_colors=wes_palette("Moonrise3", n = 4)
names(prediction_colors)=c('C1','C2','C3','C4')
HLA_colors=wes_palette("Darjeeling2", n = 2)
names(HLA_colors)=c('Class I','Class II')
annotation_colors = list(
  'Prediction' = prediction_colors,
  'BOR' = BOR_colors,
  'HLA' = HLA_colors
)
my_palette <- colorRampPalette(c("deepskyblue4","white","goldenrod2"))(200)

png(filename = "HLA_heatmap_ici.png",width = 1300, height = 1700,res = 300)
pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=F, clustering_distance_rows = "euclidean",
                   clustering_method = 'complete', 
                   annotation_row=annotation_row, annotation_col=annotation_col,annotation_colors=annotation_colors,
                   annotation_names_row = F,annotation_names_col = F, show_rownames = F,
                   legend=T, fontsize = 8, angle_col = "315",
                   main='Liu Melanoma anti-PD1')
dev.off()


# Figure 6D
ici_df$immuneSubtype[ici_df$Prediction=='C1']='C1\nwound-healing'
ici_df$immuneSubtype[ici_df$Prediction=='C2']='C2\nIFN-γ dominant'
ici_df$immuneSubtype[ici_df$Prediction=='C3']='C3\ninflammatory'
ici_df$immuneSubtype[ici_df$Prediction=='C4']='C4\nlymphocyte-depleted'
ici_df$immuneSubtype[ici_df$Prediction=='C5']='C5\nimmunologically quiet'
ici_df$immuneSubtype[ici_df$Prediction=='C6']='C6\nTGF-β dominant'
ici_df$exhaustion=(log2(ici_df$HAVCR2 + 1) + log2(ici_df$TIGIT + 1) + log2(ici_df$LAG3 + 1) + log2(ici_df$IDO1 + 1))/4
ggplot(ici_df, aes(x=immuneSubtype, y=exhaustion, fill=immuneSubtype)) +
  geom_violin() +
  stat_compare_means(label.x = 3, label.y=7, size=6) + 
  geom_boxplot(width=0.1, fill="white") + 
  guides(fill=FALSE) + 
  theme_bw() + xlab('Predicted Immune Subtypes') + ylab('Immune Exhaustion Score')
ggsave("predicted_immuneSubtypes_exhaustion.png",height = 5,width =8,dpi = 300)



# Figure S6A
ici_df$CYT=(log2(ici_df$PRF1 + 1) + log2(ici_df$GZMA + 1))/2
ici_df$immuneSubtype[ici_df$Prediction=='C1']='C1\nwound-healing'
ici_df$immuneSubtype[ici_df$Prediction=='C2']='C2\nIFN-γ dominant'
ici_df$immuneSubtype[ici_df$Prediction=='C3']='C3\ninflammatory'
ici_df$immuneSubtype[ici_df$Prediction=='C4']='C4\nlymphocyte-depleted'
ici_df$immuneSubtype[ici_df$Prediction=='C5']='C5\nimmunologically quiet'
ici_df$immuneSubtype[ici_df$Prediction=='C6']='C6\nTGF-β dominant'

ici_df$norm_HLAI = log2(nthroot(ici_df$norm_HLA_A * ici_df$norm_HLA_B * ici_df$norm_HLA_C * ici_df$norm_B2M,4))
ici_df$norm_HLAII = log2(nthroot(ici_df$HLA_DRA * ici_df$HLA_DRB1 * ici_df$HLA_DQA1 * ici_df$HLA_DQB1 * ici_df$HLA_DPA1 * ici_df$HLA_DPB1, 6))
ici_df$norm_HLAII[which(is.infinite(ici_df$norm_HLAII))]=0
ggplot(ici_df, aes(x=immuneSubtype, y=CYT, fill=immuneSubtype)) +
  geom_violin() +
  stat_compare_means(label.x = 3, label.y=7, size=6) + 
  geom_boxplot(width=0.1, fill="white") + 
  guides(fill=FALSE) + 
  theme_bw() + xlab('Predicted Immune Subtypes') + ylab('Cytolytic Score')
ggsave("predicted_immuneSubtypes_CYT.png",height = 5,width =8,dpi = 300)





