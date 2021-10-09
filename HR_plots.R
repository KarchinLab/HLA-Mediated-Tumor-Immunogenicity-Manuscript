library(gtools)

# Figure 4F: overall HLAI, HLAII fold change, CYT
HLAupregulation_df=subset(df_immuneLandscape,df_immuneLandscape$Subtype_Immune_Model_Based!='')
HLAupregulation_df$geom_high_HLAIupreg=HLAupregulation_df$geom_HLAI_upregulation > unname(quantile(HLAupregulation_df$geom_HLAI_upregulation,na.rm=TRUE)["50%"])
HLAupregulation_df$geom_high_HLAIIupreg=HLAupregulation_df$geom_HLAII_upregulation > unname(quantile(HLAupregulation_df$geom_HLAII_upregulation,na.rm=TRUE)["50%"])
HLAupregulation_df$CYT=HLAupregulation_df$PRF1_GZMA > unname(quantile(HLAupregulation_df$PRF1_GZMA,na.rm=TRUE)["75%"])
HLAupregulation_df$HLAI = HLAupregulation_df$geom_high_HLAIupreg
HLAupregulation_df$HLAII = HLAupregulation_df$geom_high_HLAIIupreg
HLAupregulation_df$TMB = HLAupregulation_df$mutationrate_nonsilent_per_Mb > unname(quantile(HLAupregulation_df$mutationrate_nonsilent_per_Mb,na.rm=TRUE)["75%"])
HLAupregulation_df$Neoantigens = HLAupregulation_df$numberOfBindingExpressedPMHC > unname(quantile(HLAupregulation_df$numberOfBindingExpressedPMHC,na.rm=TRUE)["75%"])

df=HLAupregulation_df
explanatory = c('CYT','TMB','Neoantigens',"HLAI", "HLAII")
dependent = "Surv(PFI_time_1, PFI_1)"
df %>%
  hr_plot(dependent, explanatory, remove_ref = T,
          column_space = c(-0.4, -0.2, 0.8, 1), dependent_label = "", prefix = "", suffix = "", table_text_size = 3, title_text_size = 18)

# overall HLAI, HLAII expression, CYT
HLAupregulation_df$geom_high_HLAI=HLAupregulation_df$geom_HLAI > unname(quantile(HLAupregulation_df$geom_HLAI,na.rm=TRUE)["50%"])
HLAupregulation_df$geom_high_HLAII=HLAupregulation_df$geom_HLAII > unname(quantile(HLAupregulation_df$geom_HLAII,na.rm=TRUE)["50%"])
HLAupregulation_df$HLAI = HLAupregulation_df$geom_high_HLAI
HLAupregulation_df$HLAII = HLAupregulation_df$geom_high_HLAII
explanatory = c('CYT',"HLAI", "HLAII")
dependent = "Surv(PFI_time_1, PFI_1)"
HLAupregulation_df %>%
  hr_plot(dependent, explanatory, remove_ref = T,
          column_space = c(-0.4, -0.2, 0.8, 1), dependent_label = "", prefix = "", suffix = "", table_text_size = 3, title_text_size = 18)


# Figure S4C: HLAII upregulation in different tumors
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

# Note: high HLAII upregulation is computed within each tumor type, not in all tumors
dfList=list(LUAD_immuneLandscape, BRCA_immuneLandscape, COAD_immuneLandscape, THCA_immuneLandscape, ESCA_immuneLandscape, UCEC_immuneLandscape, HNSC_immuneLandscape, LIHC_immuneLandscape, BLCA_immuneLandscape, KIRP_immuneLandscape, KIRC_immuneLandscape, CHOL_immuneLandscape, LUSC_immuneLandscape)
row=lapply(dfList, function(x) {
  x$high_HLAIIupreg=x$geom_HLAII_upregulation > unname(quantile(x$geom_HLAII_upregulation,na.rm=TRUE)["50%"])
  res.cox=coxph(Surv(PFI_time_1, PFI_1) ~ high_HLAIIupreg, data = x)
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
df <- data.frame(
  HLA=c('LUAD','BRCA',"COAD","THCA",'ESCA','UCEC','HNSC','LIHC',"BLCA","KIRP",'KIRC','CHOL','LUSC'),
  boxOdds = coeff,
  boxCILow = lower_CI,
  boxCIHigh = upper_CI,
  pVal = pVal
)
df$HLA = factor(df$HLA, c('LUAD','BRCA','KIRC',"THCA",'COAD','UCEC','ESCA',"HNSC","LIHC",'BLCA','LUSC',"KIRP",'CHOL'))
stars=stars.pval(pVal)
# Plot
# 575 * 570
ggplot(df, aes(x = boxOdds, y = HLA)) +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") + 
  geom_point(size = 3.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0,7,1) ) +
  ylab("") +
  xlab("Hazard ratio (log scale)") +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + coord_cartesian(xlim = c(-1, 3)) +
  annotate(geom = "text",label = c('*p = 0.011','*p = 0.048'),y =c(1,2), x = c(-0.6,-0.6),size=4)

# hazard ratio: all tumor types 
model <- coxph( Surv(PFI_time_1, PFI_1) ~ Study,
                data = df_immuneLandscape )
ggforest(model)




