# df_immuneLandscape: PanCancer matrix

# Figure 1A: heatmap
heatmap_df=subset(df_immuneLandscape,!is.na(df_immuneLandscape$molecule.subtype) & !is.na(df_immuneLandscape$tumor_HLA_A))
heatmap_df$Study_molecularSubtype=paste0(heatmap_df$Study,' ',heatmap_df$molecule.subtype)
heatmap_df$Study_molecularSubtype[which(heatmap_df$molecule.subtype == 'Not_Applicable')]=heatmap_df$Study[which(heatmap_df$molecule.subtype == 'Not_Applicable')]
Study_molecularSubtype_vec=sort(unique(heatmap_df$Study_molecularSubtype))

ACC=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'ACC')
BLCA=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'BLCA')
BRCA.Basal=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'BRCA Basal')
BRCA.Her2=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'BRCA Her2')
BRCA.LumA=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'BRCA LumA')

BRCA.LumB=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'BRCA LumB')
BRCA.Normal=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'BRCA Normal')
CESC.AdenoCarcinoma=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'CESC AdenoCarcinoma')
CESC.SquamousCarcinoma=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'CESC SquamousCarcinoma')
CHOL=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'CHOL')
COAD.CIN=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'COAD CIN')
COAD.GS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'COAD GS')
COAD.MSI=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'COAD MSI')
COAD.POLE=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'COAD POLE')

ESCA.CIN=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'ESCA CIN')
ESCA.ESCC=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'ESCA ESCC')
ESCA.MSI=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'ESCA MSI')
ESCA.POLE=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'ESCA POLE')

GBM.IDHmut.non.codel=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'GBM IDHmut-non-codel')
GBM.IDHwt=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'GBM IDHwt')
HNSC.HPV.NEG=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'HNSC HPV-')

HNSC.HPV.POS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'HNSC HPV+')
KICH=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'KICH')
KIRC=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'KIRC')
KIRP=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'KIRP')
LGG.IDHmut.codel=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'LGG IDHmut-codel')
LGG.IDHmut.non.codel=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'LGG IDHmut-non-codel')

LGG.IDHwt=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'LGG IDHwt')
LIHC=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'LIHC')
LUAD=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'LUAD')
LUSC=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'LUSC')
MESO=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'MESO')
OV=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'OV')
PAAD=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'PAAD')
PCPG=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'PCPG')
PRAD=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'PRAD')
READ.CIN=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'READ CIN')
READ.GS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'READ GS')
READ.MSI=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'READ MSI')

READ.POLE=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'READ POLE')
SARC.DDLPS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'SARC DDLPS')
SARC.LMS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'SARC LMS')
SARC.MFS.UPS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'SARC MFS/UPS')
SARC.Other=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'SARC Other')
SKCM=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'SKCM')
STAD.CIN=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'STAD CIN')
STAD.EBV=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'STAD EBV')
STAD.GS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'STAD GS')
STAD.MSI=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'STAD MSI')
STAD.POLE=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'STAD POLE')
TGCT.non.seminoma=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'TGCT non-seminoma')

TGCT.seminoma=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'TGCT seminoma')
THCA=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'THCA')
THYM=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'THYM')
UCEC.CN_HIGH=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'UCEC CN_HIGH')
UCEC.CN_LOW=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'UCEC CN_LOW')
UCEC.MSI=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'UCEC MSI')
UCEC.POLE=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'UCEC POLE')
UCS=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'UCS')
UVM=subset(heatmap_df,heatmap_df$Study_molecularSubtype == 'UVM')

# Figure 1A
# HLA supertypes heatmap (included nonclassical ones)
dfList=list(GBM.IDHwt,GBM.IDHmut.non.codel,LGG.IDHwt,LGG.IDHmut.codel,LGG.IDHmut.non.codel,
            UVM,
            HNSC.HPV.POS,HNSC.HPV.NEG,
            THCA,ACC,PCPG,
            THYM,
            LUAD,LUSC,MESO,
            BRCA.Basal,BRCA.Her2,BRCA.LumA,BRCA.LumB,BRCA.Normal,
            ESCA.CIN,ESCA.ESCC,ESCA.MSI,ESCA.POLE,STAD.CIN,STAD.EBV,STAD.GS,STAD.MSI,STAD.POLE,COAD.CIN,COAD.GS,COAD.MSI,COAD.POLE,READ.CIN,READ.GS,READ.MSI,READ.POLE,
            LIHC,CHOL,PAAD,
            KICH,KIRC,KIRP,BLCA,PRAD,TGCT.non.seminoma,TGCT.seminoma,
            OV,UCEC.CN_HIGH,UCEC.CN_LOW,UCEC.MSI,UCEC.POLE,UCS,CESC.AdenoCarcinoma,CESC.SquamousCarcinoma,
            SKCM,
            SARC.DDLPS,SARC.LMS,SARC.MFS.UPS,SARC.Other)
molecular.subtypes.names=c('GBM.IDHwt (n=133)','GBM.IDHmut.non.codel (n=7)','LGG.IDHwt (n=92)','LGG.IDHmut.codel (n=166)','LGG.IDHmut.non.codel (n=244)',
                           'UVM (n=80)',
                           'HNSC.HPV.POS (n=69)','HNSC.HPV.NEG (412)',
                           'THCA (n=477)','ACC (n=76)','PCPG (n=161)',
                           'THYM (n=118)',
                           'LUAD (n=499)','LUSC (n=464)','MESO (n=81)',
                           'BRCA.Basal (n=171)','BRCA.Her2 (n=78)','BRCA.LumA (n=497)','BRCA.LumB (n=196)','BRCA.Normal (n=36)',
                           'ESCA.CIN (n=67)','ESCA.ESCC (n=78)','ESCA.MSI (n=2)','ESCA.POLE (n=2)','STAD.CIN (n=207)','STAD.EBV (n=27)','STAD.GS (n=45)','STAD.MSI (n=61)','STAD.POLE (n=7)','COAD.CIN (n=225)','COAD.GS (n=49)','COAD.MSI (n=60)','COAD.POLE (n=6)','READ.CIN (n=102)','READ.GS (n=9)','READ.MSI (n=3)','READ.POLE (n=4)',
                           'LIHC (n=348)','CHOL (n=36)','PAAD (n=151)',
                           'KICH (n=64)','KIRC (n=350)','KIRP (n=269)','BLCA (n=399)','PRAD (n=477)','TGCT.non.seminoma (n=82)','TGCT.seminoma (n=62)',
                           'OV (n=156)','UCEC.CN_HIGH (n=162)','UCEC.CN_LOW (n=146)','UCEC.MSI (n=148)','UCEC.POLE (n=49)','UCS (n=55)','CESC.AdenoCarcinoma (n=43)','CESC.SquamousCarcinoma (n=229)',
                           'SKCM (n=362)',
                           'SARC.DDLPS (n=46)','SARC.LMS (n=83)','SARC.MFS.UPS (n=80)','SARC.Other (n=20)')
count=0
HLA_supertypes=lapply(dfList, function(x) {
  supertypes = data.frame(mean(x$norm_T_HLAI,na.rm=T),mean(x$norm_HLAII,na.rm=T),mean(x$norm_T_HLA_A_log,na.rm=T), mean(x$norm_T_HLA_B_log,na.rm=T), mean(x$norm_T_HLA_C_log,na.rm=T), mean(x$norm_T_HLA_E_log,na.rm=T), mean(x$norm_T_HLA_G_log,na.rm=T), mean(x$norm_T_B2M_log,na.rm=T), mean(x$HLA_DRA_log,na.rm=T), mean(x$HLA_DRB1_log,na.rm=T), mean(x$HLA_DQA1_log,na.rm=T), mean(x$HLA_DQB1_log,na.rm=T), mean(x$HLA_DPA1_log,na.rm=T), mean(x$HLA_DPB1_log,na.rm=T), mean(x$HLA_DMA_log,na.rm=T), mean(x$HLA_DMB_log,na.rm=T), mean(x$HLA_DOA_log,na.rm=T), mean(x$HLA_DOB_log,na.rm=T))
  colnames(supertypes)=c('HLAI','HLAII',"A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB')
  count=count+1
  return(supertypes)
} )
combined_matrix <- data.frame(matrix(unlist(HLA_supertypes), nrow=length(HLA_supertypes), byrow=TRUE))
colnames(combined_matrix)=c('HLAI','HLAII',"A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB')
rownames(combined_matrix)=molecular.subtypes.names
data=as.matrix(combined_matrix[3:ncol(combined_matrix)])

combined_matrix$molecular.subtype=rownames(combined_matrix)
HLAI_sort=combined_matrix %>% arrange(combined_matrix$HLAI)
rownames(HLAI_sort)=HLAI_sort$molecular.subtype
HLAI_sort$HLAI_rank=NA
HLAI_sort$HLAI_rank[45:60]='75-100th'
HLAI_sort$HLAI_rank[30:44]='50-75th'
HLAI_sort$HLAI_rank[15:29]='25-50th'
HLAI_sort$HLAI_rank[1:14]='0-25th'
HLAI_rank=HLAI_sort[molecular.subtypes.names,'HLAI_rank']

HLAII_sort=combined_matrix %>% arrange(combined_matrix$HLAII)
rownames(HLAII_sort)=HLAII_sort$molecular.subtype
HLAII_sort$HLAII_rank=NA
HLAII_sort$HLAII_rank[45:60]='75-100th'
HLAII_sort$HLAII_rank[30:44]='50-75th'
HLAII_sort$HLAII_rank[15:29]='25-50th'
HLAII_sort$HLAII_rank[1:14]='0-25th'
HLAII_rank=HLAII_sort[molecular.subtypes.names,'HLAII_rank']

tissue_vec=c(rep("CNS", 5),rep("EYE", 1),rep("Head and Neck", 2),rep("Endocrine", 3),rep("Thymus", 1),
             rep("Thoracic", 3),rep("Breast", 5),rep("Core Gastrointestinal", 17),rep("Developmental GI Tract", 3),
             rep("Genitourinary", 7),rep("Gynecologic", 8),rep("Skin", 1),rep("Soft Tissue", 4))
annotation_col <- data.frame(row.names = c("A","B","C","E","G",'B2M',"DRA", "DRB1","DQA1","DQB1","DPA1","DPB1","DMA",'DMB',"DOA",'DOB'), 
                             Class = c(rep("HLA-I", 6), rep("HLA-II", 10)))
annotation_row <- data.frame(row.names = molecular.subtypes.names,
                             Tissue = tissue_vec,
                             HLA.II.Expression.Quantile = HLAII_rank,
                             HLA.I.Expression.Quantile = HLAI_rank)
my_palette <- colorRampPalette(c("deepskyblue4","white","goldenrod2"))(200)

paletteFunc <- colorRampPalette(c('thistle2', 'forestgreen'))
palette_HLAI     <- paletteFunc(4)
names(palette_HLAI)=c('0-25th','25-50th','50-75th','75-100th')
paletteFunc <- colorRampPalette(c('lightsteelblue2', 'mediumpurple3'))
palette_HLAII     <- paletteFunc(4)
names(palette_HLAII)=c('0-25th','25-50th','50-75th','75-100th')
palette_subtypes=c(wes_palette('Moonrise3',5),wes_palette('FantasticFox1',5),wes_palette('Darjeeling2',3))
names(palette_subtypes)=unique(tissue_vec)

annotation_colors = list(
  Class = c("HLA-I" = "goldenrod2","HLA-II" = 'deepskyblue4'),
  HLA.I.Expression.Quantile = palette_HLAI,
  HLA.II.Expression.Quantile = palette_HLAII,
  Tissue = palette_subtypes
)
png(filename = "HLA_heatmap.png",width = 2700, height = 2600,res = 300)
pheatmap::pheatmap(data, color = my_palette, cluster_cols=T,cluster_rows=F, 
                   annotation_row=annotation_row, annotation_col=annotation_col,annotation_colors=annotation_colors,
                   annotation_legend = F, 
                   annotation_names_row = F,annotation_names_col = F, show_rownames = T,
                   border_color=NA, legend=T, angle_col = "315",fontsize = 7)
dev.off()



# Figure 1B: heatmap
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, THYM_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
row=lapply(dfList, function(x) {
  df=subset(x,select=(c("mutationrate_nonsilent_per_Mb","T.cells.CD8","T.cells.CD4.memory.activated",'Th1.cells','Th2.cells','Th17.cells',"T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
  matrix_df=corr.test(df, method = "spearman",use = "complete", adjust = 'fdr')
  matrix=matrix_df$r
  matrix_p=matrix_df$p
  mat=data.frame(matrix)
  colnames(mat)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
  row=mat[1,]
  return(row)
} )

TMB_matrix <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
colnames(TMB_matrix)=c("TMB","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
rownames(TMB_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
TMB_matrix=TMB_matrix[2:ncol(TMB_matrix)]

row=lapply(dfList, function(x) {
  df=subset(x,select=(c("numberOfBindingExpressedPMHC","T.cells.CD8","T.cells.CD4.memory.activated",'Th1.cells','Th2.cells','Th17.cells',"T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
  matrix_df=corr.test(df, method = "spearman",use = "complete", adjust = 'fdr')
  matrix=matrix_df$r
  matrix_p=matrix_df$p
  mat=data.frame(matrix)
  colnames(mat)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
  row=mat[1,]
  return(row)
} )
pMHC_matrix <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
colnames(pMHC_matrix)=c("pMHC","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
rownames(pMHC_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
pMHC_matrix=pMHC_matrix[2:ncol(pMHC_matrix)]

row=lapply(dfList, function(x) {
  df=subset(x,select=(c("norm_T_HLAI","T.cells.CD8","T.cells.CD4.memory.activated",'Th1.cells','Th2.cells','Th17.cells',"T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
  matrix_df=corr.test(df, method = "spearman",use = "complete", adjust = 'fdr')
  matrix=matrix_df$r
  matrix_p=matrix_df$p
  mat=data.frame(matrix)
  colnames(mat)=c("HLA-II","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
  row=mat[1,]
  return(row)
} )
HLAI_matrix <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
colnames(HLAI_matrix)=c("HLA-I","CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4","LAG3","TIGIT",'HAVCR2')
rownames(HLAI_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HLAI_matrix=HLAI_matrix[2:ncol(HLAI_matrix)]


row=lapply(dfList, function(x) {
  df=subset(x,select=(c("norm_HLAII","T.cells.CD8","T.cells.CD4.memory.activated",'Th1.cells','Th2.cells','Th17.cells',"T.cells.regulatory.Tregs","NK.cells.activated","Macrophages.M1","Macrophages.M2","IDO1_log",'CXCL9_log','CXCL10_log','STAT1_log','IFNG_log',"PDCD1_log","CD274_log","CTLA4_log","LAG3_log","TIGIT_log","HAVCR2_log")))
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
rownames(HLAII_matrix)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")
HLAII_matrix=HLAII_matrix[2:ncol(HLAII_matrix)]

my_palette <- colorRampPalette(c("deepskyblue4","white","firebrick3"))(200)

annotation_col <- data.frame(row.names = c("CD8 Tc1","CD4 Tmem",'Th1','Th2','Th17',"Tregs","NK","M1","M2","IDO1",'CXCL9','CXCL10','STAT1','IFNG',"PDCD1","CD274","CTLA4", "LAG3","TIGIT",'HAVCR2'), 
                             'Immune Characteristics' = c(rep("Immune Infiltration", 9), rep("Proinflammatory Genes", 5), rep("Immune Checkpoints", 6)))
annotation_row <- data.frame(row.names = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM","ACC1","BLCA1","BRCA1","CESC1","CHOL1","COAD1","ESCA1","GBM1","HNSC1","KICH1","KIRC1","KIRP1","LGG1","LIHC1","LUAD1","LUSC1","MESO1","OV1","PAAD1","PCPG1","PRAD1","READ1","SARC1","SKCM1","STAD1","TGCT1","THCA1","THYM1","UCEC1","UCS1","UVM1","ACC2","BLCA2","BRCA2","CESC2","CHOL2","COAD2","ESCA2","GBM2","HNSC2","KICH2","KIRC2","KIRP2","LGG2","LIHC2","LUAD2","LUSC2","MESO2","OV2","PAAD2","PCPG2","PRAD2","READ2","SARC2","SKCM2","STAD2","TGCT2","THCA2","THYM2","UCEC2","UCS2","UVM2","ACC3","BLCA3","BRCA3","CESC3","CHOL3","COAD3","ESCA3","GBM3","HNSC3","KICH3","KIRC3","KIRP3","LGG3","LIHC3","LUAD3","LUSC3","MESO3","OV3","PAAD3","PCPG3","PRAD3","READ3","SARC3","SKCM3","STAD3","TGCT3","THCA3","THYM3","UCEC3","UCS3","UVM3"), 
                             'Variables' = c(rep("HLA I", 31), rep("HLA II", 31), rep("TMB", 31), rep('Neoantigen Load', 31)))

annotation_colors = list(
  'Immune.Characteristics' = c("Immune Infiltration" = "#85D4E3","Proinflammatory Genes" = '#F4B5BD',"Immune Checkpoints" = '#9C964A'),
  'Variables' = c("HLA I" = "#0073C2FF","HLA II" = '#EFC000FF',"TMB" = '#868686FF', 'Neoantigen Load' = '#CD534CFF')
)
mat=rbind(HLAI_matrix, HLAII_matrix, TMB_matrix, pMHC_matrix)
data2=as.matrix(mat)

png(filename = "immunogenicity-heatmap.png",width = 2000, height = 2300,res = 300)
out=pheatmap::pheatmap(data2, color = my_palette, cluster_cols=T,cluster_rows=F, annotation_colors = annotation_colors,
                       annotation_col=annotation_col, annotation_row=annotation_row, annotation_legend = T, gaps_col = c(9, 14),gaps_row = c(31, 62, 95),
                       annotation_names_row = F,annotation_names_col = F, show_rownames = F, angle_col = "315", border_color = 'black')
dev.off()


# Figure 1C: ici cohort HR plot
df=cohort122_HED
df$highHLAI=df$norm_HLAI > unname(quantile(df$norm_HLAI,na.rm=TRUE)["50%"])
df$highHLAII=df$norm_HLAII > unname(quantile(df$norm_HLAII,na.rm=TRUE)["50%"])
df$CYT=(log2(df$PRF1+1)+log2(df$GZMA+1))/2
df$highCYT=df$CYT > unname(quantile(df$CYT,na.rm=TRUE)["50%"])
df$PDCD1_log=log2(df$PDCD1+1)
df$highPD1=df$PDCD1_log > unname(quantile(df$PDCD1_log,na.rm=TRUE)["50%"])
df$HLAI = df$highHLAI
df$HLAII = df$highHLAII
df$TMB = df$nonsyn_muts > unname(quantile(df$nonsyn_muts,na.rm=TRUE)["50%"])

explanatory = c('highCYT','highPD1','TMB',"HLAI", "HLAII")
dependent = "Surv(progression_free, progression)"

png(filename = "HR_ici.png",width = 2300, height = 1200,res = 300)
df %>%
  hr_plot(dependent, explanatory, remove_ref = T,
          column_space = c(-0.4, -0.2, 0.8, 1), dependent_label = "", prefix = "Liu Anti-PD1 Melanoma (n=122)", suffix = "", table_text_size = 5, title_text_size = 10)
dev.off()


# Figure 1D-E
# HR high HLAI vs. low HLAI
df=subset(df_immuneLandscape,!is.na(df_immuneLandscape$MSI_status) & (df_immuneLandscape$Study == 'COAD' | df_immuneLandscape$Study == 'STAD' | df_immuneLandscape$Study == 'UCEC'))
df$MSI.high=F
df$MSI.high[df$MSI_status == 'MSI-H']=T
df$highHLAI=df$norm_T_HLAI > unname(quantile(df$norm_T_HLAI,na.rm=TRUE)["50%"])
df$highHLAII=df$norm_HLAII > unname(quantile(df$norm_HLAII,na.rm=TRUE)["50%"])
highHLAI_df=subset(df, df$highHLAI == T)
lowHLAI_df=subset(df, df$highHLAI == F)

dfList=list(highHLAI_df, lowHLAI_df)
row=lapply(dfList, function(x) {
  res.cox=coxph(Surv(PFI_time_1, PFI_1) ~ MSI.high, data = x)
  summary_cox = summary(res.cox)
  list=list(summary_cox$conf.int[1],summary_cox$conf.int[3],summary_cox$conf.int[4],summary(res.cox)$sctest[3])
  return(list)
} )
row=unlist(row)
coeff=c()
lower_CI=c()
upper_CI=c()
pVal=c()
for (i in (1:(length(row) / 4))) {
  coeff=append(coeff,row[(i-1)*4+1])
  lower_CI=append(lower_CI,row[(i-1)*4+2])
  upper_CI=append(upper_CI,row[(i-1)*4+3])
  pVal=append(pVal,row[(i-1)*4+4])
}
df_OR <- data.frame(
  HLA=c('MSI-H\n(high HLAI group)','MSI-H\n(low HLAI group)'),
  boxOdds = coeff,
  boxCILow = lower_CI,
  boxCIHigh = upper_CI,
  pVal = pVal
)
df_OR$HLA = factor(df_OR$HLA, c('MSI-H\n(high HLAI group)','MSI-H\n(low HLAI group)'))
stars=stars.pval(pVal)

ggplot(df_OR, aes(x = boxOdds, y = HLA)) +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") + 
  geom_point(size = 3.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0,7,1) ) +
  ylab("") +
  xlab("Odds ratio (log scale)") + ggtitle('TCGA COAD, STAD, UCEC') +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + coord_cartesian(xlim = c(0, 2.5)) +
  annotate(geom = "text",label = c('HR=0.72\nns','HR=0.45\n**p = 0.007'),y =c(1,2), x = c(2,2),size=4)
ggsave("HR_MSIh_I.png",height = 3,width =6,dpi = 300)

# HR high HLAII vs. low HLAII
df=subset(df_immuneLandscape,!is.na(df_immuneLandscape$MSI_status) & (df_immuneLandscape$Study == 'COAD' | df_immuneLandscape$Study == 'STAD' | df_immuneLandscape$Study == 'UCEC'))
df$MSI.high=F
df$MSI.high[df$MSI_status == 'MSI-H']=T
df$highHLAI=df$norm_T_HLAI > unname(quantile(df$norm_T_HLAI,na.rm=TRUE)["50%"])
df$highHLAII=df$norm_HLAII > unname(quantile(df$norm_HLAII,na.rm=TRUE)["50%"])
highHLAII_df=subset(df, df$highHLAII == T)
lowHLAII_df=subset(df, df$highHLAII == F)

dfList=list(highHLAII_df, lowHLAII_df)
row=lapply(dfList, function(x) {
  res.cox=coxph(Surv(PFI_time_1, PFI_1) ~ MSI.high, data = x)
  summary_cox = summary(res.cox)
  list=list(summary_cox$conf.int[1],summary_cox$conf.int[3],summary_cox$conf.int[4],summary(res.cox)$sctest[3])
  return(list)
} )
row=unlist(row)
coeff=c()
lower_CI=c()
upper_CI=c()
pVal=c()
for (i in (1:(length(row) / 4))) {
  coeff=append(coeff,row[(i-1)*4+1])
  lower_CI=append(lower_CI,row[(i-1)*4+2])
  upper_CI=append(upper_CI,row[(i-1)*4+3])
  pVal=append(pVal,row[(i-1)*4+4])
}
df_OR <- data.frame(
  HLA=c('MSI-H\n(high HLAII group)','MSI-H\n(low HLAII group)'),
  boxOdds = coeff,
  boxCILow = lower_CI,
  boxCIHigh = upper_CI,
  pVal = pVal
)
df_OR$HLA = factor(df_OR$HLA, c('MSI-H\n(high HLAII group)','MSI-H\n(low HLAII group)'))
stars=stars.pval(pVal)

ggplot(df_OR, aes(x = boxOdds, y = HLA)) +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") + 
  geom_point(size = 3.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0,7,1) ) +
  ylab("") +
  xlab("Odds ratio (log scale)") + ggtitle('TCGA COAD, STAD, UCEC') +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + coord_cartesian(xlim = c(0, 2.5)) +
  annotate(geom = "text",label = c('HR=0.45\n**p = 0.005','HR=0.76\nns'),y =c(1,2), x = c(2,2),size=4)
ggsave("HR_MSIh_II.png",height = 3,width =6,dpi = 300)




# Figure S1A: scatterplot between supertypes and CYT
immune_df=df_immuneLandscape %>% filter(df_immuneLandscape$Study != 'DLBC' & df_immuneLandscape$Study != 'LAML')
immune_df$B2M=immune_df$norm_T_B2M_log
immune_df$`HLA-A`=immune_df$norm_T_HLA_A_log
immune_df$`HLA-B`=immune_df$norm_T_HLA_B_log
immune_df$`HLA-C`=immune_df$norm_T_HLA_C_log
immune_df$`HLA-E`=immune_df$norm_T_HLA_E_log
immune_df$`HLA-G`=immune_df$norm_T_HLA_G_log
immune_df$`HLA-DP`=immune_df$HLA_DP_log
immune_df$`HLA-DQ`=immune_df$HLA_DQ_log
immune_df$`HLA-DR`=immune_df$HLA_DR_log
immune_df$`HLA-DM`=immune_df$HLA_DM_log
immune_df$`HLA-DO`=immune_df$HLA_DO_log
dat.m1 <- melt(immune_df,id.vars='CYT', measure.vars=c('B2M','HLA-A','HLA-B','HLA-C','HLA-E','HLA-G','HLA-DP','HLA-DQ','HLA-DR','HLA-DM','HLA-DO'))
dat.m1$Gene=dat.m1$variable

ggplot(dat.m1, aes(x=value, y=CYT, color=Gene)) + 
  geom_point(size=0.00005) + theme_bw() +
  geom_smooth(aes(color = Gene, fill = Gene), method = "loess", se=F) + 
  scale_color_viridis(discrete = TRUE, option = 'D')+
  scale_fill_viridis(discrete = TRUE) +
  stat_cor(method = "spearman", label.x = c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 4.1, 4.1, 4.1, 4.1, 4.1), label.y = c(rev(12:17),rev(13:17))) +
  xlab('Expression') + ylab('CYT')
ggsave("HLA_CYT_scatter.png",height = 4,width =7,dpi = 300)


# Figure S1B: PCA
dfList=list(ACC_immuneLandscape, BLCA_immuneLandscape, BRCA_immuneLandscape, CESC_immuneLandscape, CHOL_immuneLandscape, COAD_immuneLandscape, ESCA_immuneLandscape, GBM_immuneLandscape, HNSC_immuneLandscape, KICH_immuneLandscape, KIRC_immuneLandscape, KIRP_immuneLandscape, LGG_immuneLandscape, LIHC_immuneLandscape, LUAD_immuneLandscape, LUSC_immuneLandscape, MESO_immuneLandscape, OV_immuneLandscape, PAAD_immuneLandscape, PCPG_immuneLandscape, PRAD_immuneLandscape, READ_immuneLandscape, SARC_immuneLandscape, SKCM_immuneLandscape, STAD_immuneLandscape, TGCT_immuneLandscape, THCA_immuneLandscape, THYM_immuneLandscape, UCEC_immuneLandscape, UCS_immuneLandscape, UVM_immuneLandscape)
row=lapply(dfList, function(x) {
  x=do.call(data.frame,lapply(x, function(x) replace(x, is.infinite(x),NA)))
  hlas_df=x[c('norm_T_B2M_log','norm_T_HLA_A_log','norm_T_HLA_B_log','norm_T_HLA_C_log','norm_T_HLA_E_log','norm_T_HLA_G_log','HLA_DRA_log','HLA_DRB1_log','HLA_DQA1_log','HLA_DQB1_log','HLA_DPA1_log','HLA_DPB1_log','HLA_DMA_log','HLA_DMB_log','HLA_DOA_log','HLA_DOB_log','CYT')]
  hlas_mean=colMeans(hlas_df, na.rm = T, dims = 1)
  return(hlas_mean)
} )
combined_matrix2 <- data.frame(matrix(unlist(row), nrow=length(row), byrow=TRUE))
rownames(combined_matrix2)=c("ACC","BLCA","BRCA","CESC","CHOL","COAD","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","TGCT","THCA","THYM","UCEC","UCS","UVM")

combined_matrix = combined_matrix2
combined_matrix$CYT='low'
combined_matrix$CYT[combined_matrix$X17 > unname(quantile(combined_matrix$X17,na.rm=TRUE)["50%"])]='high'

res.pca <- PCA(combined_matrix,quanti.sup = 17,quali.sup=18)
png(filename = "PCA.png",width = 250, height = 200,res = 300)
fviz_pca_ind(res.pca, habillage=18, repel = TRUE, title='', addEllipses = T, axes=c(1,2), ellipse.type = "confidence", ellipse.level = 0.95) + 
  scale_color_manual(values=c("firebrick3", "deepskyblue3")) + theme_bw()
ggsave("pca.png",height = 5,width =6,dpi = 300)

# one-sided permutation test
pc.df=data.frame(res.pca[5])
pc.df$labels=as.numeric(as.factor(combined_matrix$CYT))
cancor=cancor(pc.df[,c('x.PC1','x.PC2')],pc.df[,'labels'])
obs.cor = cancor$cor

B=10000
new.idx = t(sapply(1:B,function(b){
  set.seed(b)
  vec = sample(1:nrow(pc.df),replace = F)
  if(identical(vec,1:nrow(pc.df))) vec = NULL
  
  return(vec)
  
})) %>% unique
print(paste0('#unique permuted rows: ',nrow(new.idx)))

new.cor = sapply(1:nrow(new.idx),function(i){
  
  newdat = pc.df
  newdat[,'labels'] = pc.df[new.idx[i,],'labels']
  cancor(newdat[,c('x.PC1','x.PC2')],newdat[,'labels'])$cor
  
})
p.value = sum(new.cor>=obs.cor)/nrow(new.idx)
print(paste0('=== Observed cancor: ',obs.cor))
print(paste0('=== P-value: ',p.value))


# Figure S1C-D
df=subset(df_immuneLandscape,!is.na(df_immuneLandscape$MSI_status) & (df_immuneLandscape$Study == 'COAD' | df_immuneLandscape$Study == 'STAD' | df_immuneLandscape$Study == 'UCEC'))
df$MSI_status=factor(df$MSI_status,rev(levels(factor(df$MSI_status))))
ggplot(df, aes(x=MSI_status, y=df$norm_T_HLAI, fill=MSI_status)) +
  geom_violin() +
  stat_compare_means(label.y=12) + 
  geom_boxplot(width=0.1, fill="white") + 
  guides(fill=FALSE) +
  facet_grid(cols = vars(Study),scales = "free") +
  scale_fill_brewer(palette="Blues") +
  theme_bw() + xlab('MSI Category') + ylab('HLA I Expression')
ggsave("MSI_HLAI.png",height = 4,width =7,dpi = 300)

ggplot(df, aes(x=MSI_status, y=df$norm_HLAII, fill=MSI_status)) +
  geom_violin() +
  stat_compare_means(label.y=12) + 
  geom_boxplot(width=0.1, fill="white") + 
  guides(fill=FALSE) +
  facet_grid(cols = vars(Study),scales = "free") +
  scale_fill_brewer(palette="Blues") +
  theme_bw() + xlab('MSI Category') + ylab('HLA II Expression')
ggsave("MSI_HLAII.png",height = 4,width =7,dpi = 300)










