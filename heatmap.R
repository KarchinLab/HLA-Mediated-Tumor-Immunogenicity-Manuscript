library(gtools)
library(pheatmap)

ACC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'ACC')
BLCA_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'BLCA')
BRCA_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'BRCA')
CESC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'CESC')
CHOL_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'CHOL')
COAD_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'COAD')
DLBC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'DLBC')
ESCA_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'ESCA')
GBM_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'GBM')
HNSC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'HNSC')
KICH_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'KICH')
KIRC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'KIRC')
KIRP_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'KIRP')
LAML_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'LAML')
LGG_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'LGG')
LIHC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'LIHC')
LUAD_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'LUAD')
LUSC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'LUSC')
MESO_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'MESO')
OV_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'OV')
PAAD_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'PAAD')
PCPG_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'PCPG')
PRAD_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'PRAD')
READ_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'READ')
SARC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'SARC')
SKCM_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'SKCM')
STAD_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'STAD')
TGCT_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'TGCT')
THCA_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'THCA')
THYM_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'THYM')
UCEC_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'UCEC')
UCS_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'UCS')
UVM_immuneLandscape=subset(df_immuneLandscape,df_immuneLandscape$Study == 'UVM')


# Figure 1A
# HLA supertypes heatmap (included nonclassical ones)
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, DLBC_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LAML_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, THYM_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
HLA_supertypes=lapply(dfList, function(x) {
  #x = x %>% filter(tumor_HLA_G != 0) %>% filter(tumor_HLA_DOB != 0)
  #supertypes = data.frame(mean(log2(x$tumor_HLA_A + 1),na.rm=T), mean(log2(x$tumor_HLA_B + 1),na.rm=T), mean(log2(x$tumor_HLA_C + 1),na.rm=T), mean(log2(x$tumor_HLA_E + 1),na.rm=T), mean(log2(x$tumor_HLA_G + 1),na.rm=T), mean(log2(x$tumor_HLA_DR + 1),na.rm=T), mean(log2(x$tumor_HLA_DP + 1),na.rm=T), mean(log2(x$tumor_HLA_DQ + 1),na.rm=T), mean(log2(x$tumor_HLA_DM + 1),na.rm=T), mean(log2(x$tumor_HLA_DO + 1),na.rm=T))
  supertypes = data.frame(mean(x$HLA_A_log,na.rm=T), mean(x$HLA_B_log,na.rm=T), mean(x$HLA_C_log,na.rm=T), mean(x$HLA_E_log,na.rm=T), mean(x$HLA_G_log,na.rm=T), mean(x$B2M_log,na.rm=T), mean(x$HLA_DRA_log,na.rm=T), mean(x$HLA_DRB1_log,na.rm=T), mean(x$HLA_DQA1_log,na.rm=T), mean(x$HLA_DQB1_log,na.rm=T), mean(x$HLA_DPA1_log,na.rm=T), mean(x$HLA_DPB1_log,na.rm=T), mean(x$HLA_DMA_log,na.rm=T), mean(x$HLA_DMB_log,na.rm=T), mean(x$HLA_DOA_log,na.rm=T), mean(x$HLA_DOB_log,na.rm=T))
  
  colnames(supertypes)=c("A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB')
  return(supertypes)
} )
combined_matrix <- data.frame(matrix(unlist(HLA_supertypes), nrow=length(HLA_supertypes), byrow=TRUE))
colnames(combined_matrix)=c("A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB')
rownames(combined_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP",'LAML',"LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
data=as.matrix(combined_matrix)

annotation_col <- data.frame(row.names = c("A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB'), 
                             HLA = c(rep("Class I", 6), rep("Class II", 10)))


annotation_colors = list(
  HLA = c("Class I" = "thistle2","Class II" = 'forestgreen')
)
# Image: 590 * 570
my_palette <- colorRampPalette(c("deepskyblue4","white","goldenrod2"))(200)
# heatmap.2(data,col= my_palette, scale = 'none', trace = "none", density.info = "none",key.xlab="Expression",keysize = 1)
# 660 * 570
pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=T, annotation_colors = annotation_colors,
                   annotation_col=annotation_col, annotation_legend = T, 
                   annotation_names_row = F,annotation_names_col = F, show_rownames = T,
                   border_color=NA, legend=T, angle_col = "315")




# Figure 2A
# correlation heatmap for HLA I, II, TMB:
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, DLBC_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, THYM_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
row=lapply(dfList, function(x) {
  df=subset(x,select=(c("mutationrate_nonsilent_per_Mb","T.cells.CD8","T.cells.CD4.memory.activated",'Th1.cells','Th2.cells','Th17.cells',"T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
  # missings=sapply(df, function(x)all(is.na(x)))
  # if (any(missings)) {
  #   for (index in which(missings)) {
  #     df[index] = 1000  # replace NA column with nonsense values
  #   }
  # }
  matrix_df=corr.test(df, method = "spearman",use = "complete", adjust = 'fdr')
  matrix=matrix_df$r
  matrix_p=matrix_df$p
  mat=data.frame(matrix)
  colnames(mat)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
  row=mat[1,]
  return(row)
} )
TMB_matrix <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
colnames(TMB_matrix)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
rownames(TMB_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
TMB_matrix=TMB_matrix[2:ncol(TMB_matrix)]

# correlation heatmap:
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, DLBC_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, THYM_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
row=lapply(dfList, function(x) {
  df=subset(x,select=(c("geom_HLAI","T.cells.CD8","T.cells.CD4.memory.activated",'Th1.cells','Th2.cells','Th17.cells',"T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
  # missings=sapply(df, function(x)all(is.na(x)))
  # if (any(missings)) {
  #   for (index in which(missings)) {
  #     df[index] = 1000  # replace NA column with nonsense values
  #   }
  # }
  matrix_df=corr.test(df, method = "spearman",use = "complete", adjust = 'fdr')
  matrix=matrix_df$r
  matrix_p=matrix_df$p
  mat=data.frame(matrix)
  colnames(mat)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
  row=mat[1,]
  return(row)
} )
HLAI_matrix <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
colnames(HLAI_matrix)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
rownames(HLAI_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HLAI_matrix=HLAI_matrix[2:ncol(HLAI_matrix)]


# correlation heatmpa:
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, DLBC_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, THYM_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
row=lapply(dfList, function(x) {
  df=subset(x,select=(c("geom_HLAII","T.cells.CD8","T.cells.CD4.memory.activated",'Th1.cells','Th2.cells','Th17.cells',"T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
  # missings=sapply(df, function(x)all(is.na(x)))
  # if (any(missings)) {
  #   for (index in which(missings)) {
  #     df[index] = 1000  # replace NA column with nonsense values
  #   }
  # }
  matrix_df=corr.test(df, method = "spearman",use = "complete", adjust = 'fdr')
  matrix=matrix_df$r
  matrix_p=matrix_df$p
  mat=data.frame(matrix)
  colnames(mat)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
  row=mat[1,]
  return(row)
} )
HLAII_matrix <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
colnames(HLAII_matrix)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4", "LAG3","TIGIT",'HAVCR2')
rownames(HLAII_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HLAII_matrix=HLAII_matrix[2:ncol(HLAII_matrix)]

# coul <- colorRampPalette(brewer.pal(8, "RdBu"))(25)
my_palette <- colorRampPalette(c("deepskyblue4","white","firebrick3"))(200)
# heatmap.2(data, scale = "none", col=my_palette, margins=c(8,8), trace = "none", density.info = "none",key.xlab="TMB\nSpearman Corr.",keysize = 1) 

annotation_col <- data.frame(row.names = c("CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4", "LAG3","TIGIT",'HAVCR2'), 
                             'I' = c(rep("Immune Infiltration", 9), rep("Proinflammatory Genes", 5), rep("Immune Checkpoints", 6)))
annotation_row <- data.frame(row.names = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM","ACC1","BLCA1","BRCA1","CESC1","CHOL1","COAD1","DLBC1","ESCA1","GBM1","HNSC1","KICH1","KIRC1","KIRP1","LGG1","LIHC1","LUAD1","LUSC1","MESO1","OV1","PAAD1","PCPG1","PRAD1","READ1","SARC1","SKCM1","STAD1","TGCT1","THCA1","THYM1","UCEC1","UCS1","UVM1","ACC2","BLCA2","BRCA2","CESC2","CHOL2","COAD2","DLBC2","ESCA2","GBM2","HNSC2","KICH2","KIRC2","KIRP2","LGG2","LIHC2","LUAD2","LUSC2","MESO2","OV2","PAAD2","PCPG2","PRAD2","READ2","SARC2","SKCM2","STAD2","TGCT2","THCA2","THYM2","UCEC2","UCS2","UVM2"), 
                            'II' = c(rep("HLA I", 32), rep("HLA II", 32), rep("TMB", 32)))


annotation_colors = list(
  'I' = c("Immune Infiltration" = "darkorange1","Proinflammatory Genes" = 'deepskyblue3',"Immune Checkpoints" = 'olivedrab'),
  'II' = c("HLA I" = "dodgerblue3","HLA II" = 'goldenrod2',"TMB" = 'gray35')
)
mat=rbind(HLAI_matrix, HLAII_matrix, TMB_matrix)
data2=as.matrix(mat)
# Image: 590 * 630
# Image: 730 * 900
out=pheatmap::pheatmap(data2, color = my_palette, cluster_cols=T,cluster_rows=F, annotation_colors = annotation_colors,
                   annotation_col=annotation_col, annotation_row=annotation_row, annotation_legend = T, gaps_col = c(9, 14),gaps_row = c(32, 64),
                   annotation_names_row = F,annotation_names_col = F, show_rownames = F, angle_col = "315")

# dendrogram: 700 * 430
plot(out$tree_col)




# Figure S2B: spearman correlation between CYT and HLAI, HLAII, pMHC, TMB
# DLBC, LAML, THYM excluded
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
row=lapply(dfList, function(x) {
  HLAI_corr=cor.test(x$geom_HLAI,x$PRF1_GZMA,method='spearman')
  HLAII_corr=cor.test(x$geom_HLAII,x$PRF1_GZMA,method='spearman')
  tmb_corr=cor.test(x$mutationrate_nonsilent_per_Mb,x$PRF1_GZMA,method='spearman')
  pMHC_corr=cor.test(x$numberOfBindingExpressedPMHC,x$PRF1_GZMA,method='spearman')
  corr_vec=c(HLAI_corr$estimate,HLAI_corr$p.value,HLAII_corr$estimate,HLAII_corr$p.value,tmb_corr$estimate,tmb_corr$p.value,pMHC_corr$estimate,pMHC_corr$p.value)
  corr_df=data.frame(corr_vec)
  
  return(corr_df)
} )
corr_matrix <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
colnames(corr_matrix)=c("HLA I","HLAI p-value","HLA II","HLAII p-value",'TMB','TMB p-value',"neoantigen\nload","neoantigen load p-value")
rownames(corr_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS","UVM")
corr_matrix$histology=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS","UVM")

plot.data <- melt(corr_matrix,id.vars='histology', measure.vars=c("HLA I","HLA II",'TMB',"neoantigen\nload"))
plot.pVal <- melt(corr_matrix,id.vars='histology', measure.vars=c("HLAI p-value","HLAII p-value",'TMB p-value',"neoantigen load p-value"))
plot.data$stars <- cut(plot.pVal$value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels

# 550 * 620
p <- ggplot(aes(x=histology, y=variable, fill=value), data=plot.data)
p + geom_tile() + scale_fill_gradient2() + 
  geom_text(aes(label=stars), color="black", size=5) + 
  labs(y=NULL, x=NULL, fill="Spearman Rho") + coord_flip() + theme_bw() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0))




# Figure 4C
C1=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based=='C1')
C2=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based=='C2')
C3=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based=='C3')
C4=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based=='C4')
C5=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based=='C5')
C6=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based=='C6')

dfList=list(C1,C2,C3,C4,C5,C6)
dfList=list(C1,C2,C3,C4,C5,C6)
HLA_supertypes=lapply(dfList, function(x) {
  supertypes = data.frame(mean(x$HLA_A_log/x$HLA_A_log_N,na.rm=T), mean(x$HLA_B_log/x$HLA_B_log_N,na.rm=T), mean(x$HLA_C_log/x$HLA_C_log_N,na.rm=T), mean(x$HLA_E_log/x$HLA_E_log_N,na.rm=T), mean(x$HLA_G_log/x$HLA_G_log_N,na.rm=T), mean(x$B2M_log/x$B2M_log_N,na.rm=T), mean(x$HLA_DRA_log/x$HLA_DRA_log_N,na.rm=T), mean(x$HLA_DRB1_log/x$HLA_DRB1_log_N,na.rm=T), mean(x$HLA_DQA1_log/x$HLA_DQA1_log_N,na.rm=T), mean(x$HLA_DQB1_log/x$HLA_DQB1_log_N,na.rm=T), mean(x$HLA_DPA1_log/x$HLA_DPA1_log_N,na.rm=T), mean(x$HLA_DPB1_log/x$HLA_DPB1_log_N,na.rm=T), mean(x$HLA_DMA_log/x$HLA_DMA_log_N,na.rm=T), mean(x$HLA_DMB_log/x$HLA_DMB_log_N,na.rm=T), mean(x$HLA_DOA_log/x$HLA_DOA_log_N,na.rm=T), mean(x$HLA_DOB_log/x$HLA_DOB_log_N,na.rm=T))
  colnames(supertypes)=c("A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB')
  return(supertypes)
} )
combined_matrix <- data.frame(matrix(unlist(HLA_supertypes), nrow=length(HLA_supertypes), byrow=TRUE))
colnames(combined_matrix)=c("A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB')
rownames(combined_matrix)=c('C1','C2','C3','C4','C5','C6')
data=as.matrix(combined_matrix)


annotation_col <- data.frame(row.names = c("A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB'), 
                             HLA = c(rep("Class I", 6), rep("Class II", 10)))

annotation_colors = list(
  HLA = c("Class I" = "thistle2","Class II" = 'forestgreen')
)
# Image: 590 * 570
my_palette <- colorRampPalette(c("deepskyblue4","white","firebrick3"))(200)
# heatmap.2(data,col= my_palette, scale = 'none', trace = "none", density.info = "none",key.xlab="Expression",keysize = 1)
# 550 * 530
pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=F, annotation_colors = annotation_colors,
                   annotation_col=annotation_col, annotation_legend = T, 
                   annotation_names_row = F,annotation_names_col = F, show_rownames = T,
                   border_color=NA, legend=T, 
                   scale='column', angle_col = "315")




