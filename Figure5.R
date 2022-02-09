# Figure 5A: fishers CYT ~ HLA LOH
df_fullyHeterozygous$highCYT=df_fullyHeterozygous$CYT > unname(quantile(df_fullyHeterozygous$CYT,na.rm=TRUE)["50%"])
df_highCYT=subset(df_fullyHeterozygous,df_fullyHeterozygous$highCYT==T)
piechart_df1=data.frame(c("LOH", "no LOH"),c(nrow(subset(df_highCYT,df_highCYT$`HLA LOH` == T)),nrow(subset(df_highCYT,df_highCYT$`HLA LOH` == F))))
colnames(piechart_df1)=c("HLALOH","ifLOH")
piechart_df1 <- piechart_df1 %>% 
  arrange(rev(HLALOH)) %>%
  mutate(prop = ifLOH / sum(piechart_df1$ifLOH) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
piechart_df1$label=paste(piechart_df1$HLALOH,'\n',piechart_df1$ifLOH,sep='')
ggplot(piechart_df1, aes(x="", y=prop, fill=HLALOH)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = label), color = "white", size=6) +
  scale_fill_brewer(palette="Set1") + ggtitle("high CYT (top 50%)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("HLALOH_pie.png",height = 3.5,width =3.5,dpi = 300)

df_lowCYT=subset(df_fullyHeterozygous,df_fullyHeterozygous$highCYT==F)
piechart_df2=data.frame(c("LOH", "no LOH"),c(nrow(subset(df_lowCYT,df_lowCYT$`HLA LOH` == T)),nrow(subset(df_lowCYT,df_lowCYT$`HLA LOH` == F))))
colnames(piechart_df2)=c("HLALOH","ifLOH")
piechart_df2 <- piechart_df2 %>% 
  arrange(rev(HLALOH)) %>%
  mutate(prop = ifLOH / sum(piechart_df2$ifLOH) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
piechart_df2$label=paste(piechart_df2$HLALOH,'\n',piechart_df2$ifLOH,sep='')
ggplot(piechart_df2, aes(x="", y=prop, fill=HLALOH)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = label), color = "white", size=6) +
  scale_fill_brewer(palette="Set1") + ggtitle("low CYT (bottom 50%)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("HLALOH_pie2.png",height = 3.5,width =3.5,dpi = 300)

df_highUpregII=data.frame('no.LOH'= piechart_df1$ifLOH[1],'LOH'= piechart_df1$ifLOH[2])
df_lowUpregII=data.frame('no.LOH'= piechart_df2$ifLOH[1],'LOH'= piechart_df2$ifLOH[2])
df_combined=rbind(df_highUpregII,df_lowUpregII)
rownames(df_combined)=c('upreg_HLAII','no.upreg_HLAII')
fisher.test(df_combined,alternative = "two.sided")


# Figure 5B: facet BRCA LUAD SKCM bar graph
df=df_HLALOH
df$immuneSubtype=df$Subtype_Immune_Model_Based
facet_df=df %>% filter(immuneSubtype != '') %>% filter(!is.na(`HLA LOH`))
brca_c1=subset(facet_df,facet_df$Study=='BRCA' & facet_df$Subtype_Immune_Model_Based=='C1')
brca_c1_LOHpercent=paste(round(length(which(brca_c1$HLA.LOH == T)) / nrow(brca_c1) * 100,0), '%',sep='')
brca_c2=subset(facet_df,facet_df$Study=='BRCA' & facet_df$Subtype_Immune_Model_Based=='C2')
brca_c2_LOHpercent=paste(round(length(which(brca_c2$HLA.LOH == T)) / nrow(brca_c2) * 100,0), '%',sep='')
brca_c3=subset(facet_df,facet_df$Study=='BRCA' & facet_df$Subtype_Immune_Model_Based=='C3')
brca_c3_LOHpercent=paste(round(length(which(brca_c3$HLA.LOH == T)) / nrow(brca_c3) * 100,0), '%',sep='')
brca_c4=subset(facet_df,facet_df$Study=='BRCA' & facet_df$Subtype_Immune_Model_Based=='C4')
brca_c4_LOHpercent=paste(round(length(which(brca_c4$HLA.LOH == T)) / nrow(brca_c4) * 100,0), '%',sep='')
brca_c6=subset(facet_df,facet_df$Study=='BRCA' & facet_df$Subtype_Immune_Model_Based=='C6')
brca_c6_LOHpercent=paste(round(length(which(brca_c6$HLA.LOH == T)) / nrow(brca_c6) * 100,0), '%',sep='')

luad_c1=subset(facet_df,facet_df$Study=='LUAD' & facet_df$Subtype_Immune_Model_Based=='C1')
luad_c1_LOHpercent=paste(round(length(which(luad_c1$HLA.LOH == T)) / nrow(luad_c1) * 100,0), '%',sep='')
luad_c2=subset(facet_df,facet_df$Study=='LUAD' & facet_df$Subtype_Immune_Model_Based=='C2')
luad_c2_LOHpercent=paste(round(length(which(luad_c2$HLA.LOH == T)) / nrow(luad_c2) * 100,0), '%',sep='')
luad_c3=subset(facet_df,facet_df$Study=='LUAD' & facet_df$Subtype_Immune_Model_Based=='C3')
luad_c3_LOHpercent=paste(round(length(which(luad_c3$HLA.LOH == T)) / nrow(luad_c3) * 100,0), '%',sep='')
luad_c4=subset(facet_df,facet_df$Study=='LUAD' & facet_df$Subtype_Immune_Model_Based=='C4')
luad_c4_LOHpercent=paste(round(length(which(luad_c4$HLA.LOH == T)) / nrow(luad_c4) * 100,0), '%',sep='')
luad_c6=subset(facet_df,facet_df$Study=='LUAD' & facet_df$Subtype_Immune_Model_Based=='C6')
luad_c6_LOHpercent=paste(round(length(which(luad_c6$HLA.LOH == T)) / nrow(luad_c6) * 100,0), '%',sep='')

skcm_c1=subset(facet_df,facet_df$Study=='SKCM' & facet_df$Subtype_Immune_Model_Based=='C1')
skcm_c1_LOHpercent=paste(round(length(which(skcm_c1$HLA.LOH == T)) / nrow(skcm_c1) * 100,0), '%',sep='')
skcm_c2=subset(facet_df,facet_df$Study=='SKCM' & facet_df$Subtype_Immune_Model_Based=='C2')
skcm_c2_LOHpercent=paste(round(length(which(skcm_c2$HLA.LOH == T)) / nrow(skcm_c2) * 100,0), '%',sep='')
skcm_c3=subset(facet_df,facet_df$Study=='SKCM' & facet_df$Subtype_Immune_Model_Based=='C3')
skcm_c3_LOHpercent=paste(round(length(which(skcm_c3$HLA.LOH == T)) / nrow(skcm_c3) * 100,0), '%',sep='')
skcm_c4=subset(facet_df,facet_df$Study=='SKCM' & facet_df$Subtype_Immune_Model_Based=='C4')
skcm_c4_LOHpercent=paste(round(length(which(skcm_c4$HLA.LOH == T)) / nrow(skcm_c4) * 100,0), '%',sep='')
skcm_c6=subset(facet_df,facet_df$Study=='SKCM' & facet_df$Subtype_Immune_Model_Based=='C6')
skcm_c6_LOHpercent=paste(round(length(which(skcm_c6$HLA.LOH == T)) / nrow(skcm_c6) * 100,0), '%',sep='')

anno <- data.frame(xstar = c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5), 
                   ystar = c(nrow(brca_c1),nrow(luad_c1),nrow(skcm_c1),nrow(brca_c2),nrow(luad_c2),nrow(skcm_c2),nrow(brca_c3),nrow(luad_c3),nrow(skcm_c3),nrow(brca_c4),nrow(luad_c4),nrow(skcm_c4),nrow(brca_c6),nrow(luad_c6),nrow(skcm_c6)) + 13,
                   lab = c(brca_c1_LOHpercent,luad_c1_LOHpercent,skcm_c1_LOHpercent,brca_c2_LOHpercent,luad_c2_LOHpercent,skcm_c2_LOHpercent,brca_c3_LOHpercent,luad_c3_LOHpercent,skcm_c3_LOHpercent,brca_c4_LOHpercent,luad_c4_LOHpercent,skcm_c4_LOHpercent,brca_c6_LOHpercent,luad_c6_LOHpercent,skcm_c6_LOHpercent),
                   Study = c("BRCA", "LUAD", 'SKCM'))
ggplot(facet_df, aes(immuneSubtype)) + 
  geom_bar(aes(fill=`HLA LOH`)) +
  scale_fill_viridis(discrete = T, option = "E") +
  facet_wrap(~Study) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab)) +
  xlab("") + ylab('Count') +
  ggtitle('Percent HLA LOH among Immune Subtypes')
ggsave("bar_hlaloh.png",height = 4,width =7,dpi = 300)


# Figure 5C
df_fullyHeterozygous$`Allelic Loss`=df_fullyHeterozygous$totalLoss
ggplot(df_fullyHeterozygous %>% filter(!is.na(totalLoss)), aes(x=meanCN, y=norm_T_HLAI, color=`Allelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(6), label.y = c(12)) +
  theme_bw() +
  xlab('HLA-I Mean Copy Number') + ylab('HLA-I Expression')
ggsave("CN_expression.png",height = 3.5,width =6.5,dpi = 300)


# Figure 5D
LUAD_HLALOH=subset(df_HLALOH,df_HLALOH$Study == 'LUAD')
LUAD_HLALOH$highHLAI=LUAD_HLALOH$norm_T_HLAI > unname(quantile(LUAD_HLALOH$norm_T_HLAI,na.rm=TRUE)["50%"])

highHLAI_df=subset(LUAD_HLALOH,LUAD_HLALOH$highHLAI == T)
lowHLAI_df=subset(LUAD_HLALOH,LUAD_HLALOH$highHLAI == F)

dfList=list(highHLAI_df, lowHLAI_df)
row=lapply(dfList, function(x) {
  res.cox=coxph(Surv(PFI_time_1, PFI_1) ~ HLA.LOH, data = x)
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
  HLA=c('HLA LOH\n(high HLAI group)','HLA LOH\n(low HLAI group)'),
  boxOdds = coeff,
  boxCILow = lower_CI,
  boxCIHigh = upper_CI,
  pVal = pVal
)
df$HLA = factor(df$HLA, c('HLA LOH\n(low HLAI group)','HLA LOH\n(high HLAI group)'))
stars=stars.pval(pVal)

ggplot(df, aes(x = boxOdds, y = HLA)) +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") + 
  geom_point(size = 3.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0,7,1) ) +
  ylab("") +
  xlab("Odds ratio (log scale)") + ggtitle('TCGA LUAD') +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + coord_cartesian(xlim = c(0, 2.5)) +
  annotate(geom = "text",label = c('HR=1.82\n**p = 0.004','HR=1.06\nns'),y =c(1,2), x = c(0.3,0.3),size=4)
ggsave("HR_LUAD_LOH.png",height = 3,width =6,dpi = 300)


# Figure S6A-C
fit <- survfit(Surv(PFI_time_1,PFI_1) ~ HLA.LOH,data = SKCM_HLALOH)
ggsurv=ggsurvplot(fit, data = SKCM_HLALOH, palette = c("deepskyblue3", "firebrick3"),
                  conf.int = FALSE, 
                  pval = F,  
                  risk.table = TRUE, 
                  title = 'TCGA SKCM',
                  risk.table.col = "strata",
                  ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank())) + ylab("PFI")
coxph(Surv(PFI_time_1, PFI_1) ~ HLA.LOH, data = SKCM_HLALOH) %>% gtsummary::tbl_regression(exp = TRUE)

ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                            x = c(150, 180, 335), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                            label = c("p = 0.20", "HR = 0.36", "95% CI = 0.09-1.50"), size = 4)
png(filename = "skcm_hlaloh_surv.png",width = 1800, height = 1500,res = 300)
ggsurv
dev.off()



fit <- survfit(Surv(PFI_time_1,PFI_1) ~ HLA.LOH,data = BRCA_HLALOH)
ggsurv=ggsurvplot(fit, data = BRCA_HLALOH, palette = c("deepskyblue3", "firebrick3"),
                  conf.int = FALSE, 
                  pval = F,  
                  risk.table = TRUE, 
                  title = 'TCGA BRCA',
                  risk.table.col = "strata",
                  ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank())) + ylab("PFI")

ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                            x = c(690, 800, 1440), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                            label = c("p = 0.10", "HR = 0.66", "95% CI = 0.40-1.09"), size = 4)
png(filename = "brca_hlaloh_surv.png",width = 1800, height = 1500,res = 300)
ggsurv
dev.off()
coxph(Surv(PFI_time_1, PFI_1) ~ HLA.LOH, data = BRCA_HLALOH) %>% gtsummary::tbl_regression(exp = TRUE)

fit <- survfit(Surv(PFI_time_1,PFI_1) ~ HLA.LOH,data = LUAD_HLALOH)
ggsurv=ggsurvplot(fit, data = LUAD_HLALOH, palette = c("deepskyblue3", "firebrick3"),
                  conf.int = FALSE, 
                  pval = F,  
                  risk.table = TRUE, 
                  risk.table.col = "strata",
                  title = 'TCGA LUAD',
                  ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank())) + ylab("PFI")

ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                            x = c(640, 680, 1280), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                            label = c("p = 0.020", "HR = 1.42", "95% CI = 1.06-1.91"), size = 4)
png(filename = "luad_hlaloh_surv.png",width = 1800, height = 1500,res = 300)
ggsurv
dev.off()
coxph(Surv(PFI_time_1, PFI_1) ~ HLA.LOH, data = LUAD_HLALOH) %>% gtsummary::tbl_regression(exp = TRUE)




# Figure S5D-F
df_fullyHeterozygous$`HLA-A\nAllelic Loss`=df_fullyHeterozygous$lostA
df_fullyHeterozygous$meanCN_A=(df_fullyHeterozygous$HLA_A1 + df_fullyHeterozygous$HLA_A2)/2

ggplot(df_fullyHeterozygous %>% filter(!is.na(lostA)), aes(x=meanCN_A, y=HLA_A_log, color=`HLA-A\nAllelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(7), label.y = c(16)) +
  theme_bw() +
  xlab('HLA-A Mean Copy Number') + ylab('HLA-A Expression')
ggsave("HLAA_CN.png",height = 3,width =5,dpi = 300)

df_fullyHeterozygous$`HLA-B\nAllelic Loss`=df_fullyHeterozygous$lostB
df_fullyHeterozygous$meanCN_B=(df_fullyHeterozygous$HLA_B1 + df_fullyHeterozygous$HLA_B2)/2

ggplot(df_fullyHeterozygous %>% filter(!is.na(lostB)), aes(x=meanCN_B, y=HLA_B_log, color=`HLA-B\nAllelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(5), label.y = c(16)) +
  theme_bw() +
  xlab('HLA-B Mean Copy Number') + ylab('HLA-B Expression')
ggsave("HLAB_CN.png",height = 3,width =5,dpi = 300)

df_fullyHeterozygous$`HLA-C\nAllelic Loss`=df_fullyHeterozygous$lostC
df_fullyHeterozygous$meanCN_C=(df_fullyHeterozygous$HLA_C1 + df_fullyHeterozygous$HLA_C2)/2

ggplot(df_fullyHeterozygous %>% filter(!is.na(lostC)), aes(x=meanCN_C, y=HLA_C_log, color=`HLA-C\nAllelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(4), label.y = c(16)) +
  theme_bw() +
  xlab('HLA-C Mean Copy Number') + ylab('HLA-C Expression')
ggsave("HLAC_CN.png",height = 3,width =5,dpi = 300)


# Figure S5G-H
ggplot(df_fullyHeterozygous %>% filter(!is.na(`HLA LOH`)), aes(x= `HLA LOH`, y=geom_HLAI)) + 
  geom_violin(aes(fill=`HLA LOH`)) + 
  scale_fill_manual(values=c("deepskyblue4", "firebrick")) + 
  stat_compare_means(label.x = 1, label.y = 11.5) +
  geom_boxplot(width=0.1) + 
  facet_wrap(~Study) +
  theme_bw() +
  ylab('HLA-I Expression')
ggsave("HLAI_LOH.png",height = 5,width =7,dpi = 300)

# SKCM has no available normal sampls
ggplot(df_fullyHeterozygous %>% filter(!is.na(`HLA LOH`) & df_fullyHeterozygous$Study != 'SKCM'), aes(x= `HLA LOH`, y=geom_HLAI_upregulation)) + 
  geom_violin(aes(fill=`HLA LOH`)) + 
  scale_fill_manual(values=c("deepskyblue4", "firebrick")) + 
  stat_compare_means(label.x = 1, label.y = 1.4) +
  geom_boxplot(width=0.1) + 
  facet_wrap(~Study) +
  theme_bw() +
  ylab('HLA-I Fold Change')
ggsave("HLAI_DE_LOH.png",height = 3,width =5,dpi = 300)





















