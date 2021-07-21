library(FactoMineR)
library(factoextra)

# Figure 2C
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, DLBC_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LAML_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, THYM_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
row=lapply(dfList, function(x) {
  x$tumor_B2M = log2(x$tumor_B2M + 1)
  x$tumor_HLA_A = log2(x$tumor_HLA_A + 1)
  x$tumor_HLA_B = log2(x$tumor_HLA_B + 1)
  x$tumor_HLA_C = log2(x$tumor_HLA_C + 1)
  x$tumor_HLA_E = log2(x$tumor_HLA_E + 1)
  x$tumor_HLA_G = log2(x$tumor_HLA_G + 1)
  x$tumor_HLA_DRA = log2(x$tumor_HLA_DRA + 1)
  x$tumor_HLA_DRB1 = log2(x$tumor_HLA_DRB1 + 1)
  x$tumor_HLA_DQA1 = log2(x$tumor_HLA_DQA1 + 1)
  x$tumor_HLA_DQB1 = log2(x$tumor_HLA_DQB1 + 1)
  x$tumor_HLA_DPA1 = log2(x$tumor_HLA_DPA1 + 1)
  x$tumor_HLA_DPB1 = log2(x$tumor_HLA_DPB1 + 1)
  x$tumor_HLA_DMA = log2(x$tumor_HLA_DMA + 1)
  x$tumor_HLA_DMB = log2(x$tumor_HLA_DMB + 1)
  x$tumor_HLA_DOA = log2(x$tumor_HLA_DOA + 1)
  x$tumor_HLA_DOB = log2(x$tumor_HLA_DOB + 1)
  x=do.call(data.frame,lapply(x, function(x) replace(x, is.infinite(x),NA)))
  hlas_df=x[c('tumor_B2M','tumor_HLA_A','tumor_HLA_B','tumor_HLA_C','tumor_HLA_E','tumor_HLA_G','tumor_HLA_DRA','tumor_HLA_DRB1','tumor_HLA_DQA1','tumor_HLA_DQB1','tumor_HLA_DPA1','tumor_HLA_DPB1','tumor_HLA_DMA','tumor_HLA_DMB','tumor_HLA_DOA','tumor_HLA_DOB','PRF1_GZMA')]
  hlas_mean=colMeans(hlas_df, na.rm = T, dims = 1)
  return(hlas_mean)
} )
combined_matrix2 <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
rownames(combined_matrix2)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP",'LAML',"LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
# center the data
colMeans=colMeans(combined_matrix2)
combined_matrix2$X1=combined_matrix2$X1 - colMeans[1]
combined_matrix2$X2=combined_matrix2$X2 - colMeans[2]
combined_matrix2$X3=combined_matrix2$X3 - colMeans[3]
combined_matrix2$X4=combined_matrix2$X4 - colMeans[4]
combined_matrix2$X5=combined_matrix2$X5 - colMeans[5]
combined_matrix2$X6=combined_matrix2$X6 - colMeans[6]
combined_matrix2$X7=combined_matrix2$X7 - colMeans[7]
combined_matrix2$X8=combined_matrix2$X8 - colMeans[8]
combined_matrix2$X9=combined_matrix2$X9 - colMeans[9]
combined_matrix2$X10=combined_matrix2$X10 - colMeans[10]
combined_matrix2$X11=combined_matrix2$X11 - colMeans[11]
combined_matrix2$X12=combined_matrix2$X12 - colMeans[12]
combined_matrix2$X13=combined_matrix2$X13 - colMeans[13]
combined_matrix2$X14=combined_matrix2$X14 - colMeans[14]
combined_matrix2$X15=combined_matrix2$X15 - colMeans[15]
combined_matrix2$X16=combined_matrix2$X16 - colMeans[16]
combined_matrix2$X17=combined_matrix2$X17 - colMeans[17]

combined_matrix = combined_matrix2
combined_matrix$highCYT=combined_matrix$X17 > unname(quantile(combined_matrix$X17,na.rm=TRUE)["75%"])
highCYT=subset(combined_matrix,combined_matrix$highCYT==T)
highCYT$CYT='high'
lowCYT=subset(combined_matrix,combined_matrix$highCYT==F)
lowCYT$CYT='low'
combined_matrix=rbind(highCYT,lowCYT)
combined_matrix=combined_matrix[c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10','X11','X12','X13','X14','X15','X16','X17',"CYT")]
res.pca <- PCA(combined_matrix,quanti.sup = 17,quali.sup=18)
fviz_eig(res.pca)  # scree plot
# Image: 680 * 582
fviz_pca_ind(res.pca, habillage=18, repel = TRUE, title='', addEllipses = F, axes=c(1,2)) + 
  scale_color_manual(values=c("firebrick3", "deepskyblue3"))


