BRCA_HLALOH=read.csv("/Users/katana/PycharmProjects/LOHHLA_analysis/BRCA/BRCA_LOHHLA_immunelandscape.csv")
LUAD_HLALOH=read.csv("/Users/katana/PycharmProjects/LOHHLA_analysis/LUAD/LUAD_LOHHLA_immunelandscape.csv")
SKCM_HLALOH=read.csv("/Users/katana/PycharmProjects/LOHHLA_analysis/SKCM/SKCM_LOHHLA_immunelandscape.csv")

BRCA_HLALOH$high_HLAI=BRCA_HLALOH$geom_HLAI > unname(quantile(BRCA_HLALOH$geom_HLAI,na.rm=TRUE)["50%"])
LUAD_HLALOH$high_HLAI=LUAD_HLALOH$geom_HLAI > unname(quantile(LUAD_HLALOH$geom_HLAI,na.rm=TRUE)["50%"])
SKCM_HLALOH$high_HLAI=SKCM_HLALOH$geom_HLAI > unname(quantile(SKCM_HLALOH$geom_HLAI,na.rm=TRUE)["50%"])

BRCA_HLALOH$`HLA LOH` = BRCA_HLALOH$HLA.LOH
LUAD_HLALOH$`HLA LOH` = LUAD_HLALOH$HLA.LOH
SKCM_HLALOH$`HLA LOH` = SKCM_HLALOH$HLA.LOH

BRCA_HLALOH$meanCN=(BRCA_HLALOH$HLA_A1 + BRCA_HLALOH$HLA_A2 + BRCA_HLALOH$HLA_B1 + BRCA_HLALOH$HLA_B2 + BRCA_HLALOH$HLA_C1 + BRCA_HLALOH$HLA_C2) / 6
LUAD_HLALOH$meanCN=(LUAD_HLALOH$HLA_A1 + LUAD_HLALOH$HLA_A2 + LUAD_HLALOH$HLA_B1 + LUAD_HLALOH$HLA_B2 + LUAD_HLALOH$HLA_C1 + LUAD_HLALOH$HLA_C2) / 6
SKCM_HLALOH$meanCN=(SKCM_HLALOH$HLA_A1 + SKCM_HLALOH$HLA_A2 + SKCM_HLALOH$HLA_B1 + SKCM_HLALOH$HLA_B2 + SKCM_HLALOH$HLA_C1 + SKCM_HLALOH$HLA_C2) / 6

df_HLALOH=rbind(BRCA_HLALOH,LUAD_HLALOH,SKCM_HLALOH)
# fully heterozygous
df_fullyHeterozygous=subset(df_HLALOH,!is.na(df_HLALOH$lostA) & !is.na(df_HLALOH$lostB) & !is.na(df_HLALOH$lostC) )
df_fullyHeterozygous=subset(df_HLALOH,df_HLALOH$fully_heterozygous == T)
df_fullyHeterozygous$immuneSubtype=df_fullyHeterozygous$Subtype_Immune_Model_Based
df_fullyHeterozygous$`HLA LOH` = df_fullyHeterozygous$HLA.LOH
df_fullyHeterozygous$totalLoss=df_fullyHeterozygous$lostA + df_fullyHeterozygous$lostB + df_fullyHeterozygous$lostC

# Figure 6A: fishers CYT ~ HLA LOH
df_fullyHeterozygous$highCYT=df_fullyHeterozygous$PRF1_GZMA > unname(quantile(df_fullyHeterozygous$PRF1_GZMA,na.rm=TRUE)["50%"])
df_highCYT=subset(df_fullyHeterozygous,df_fullyHeterozygous$highCYT==T)
piechart_df1=data.frame(c("LOH", "no LOH"),c(nrow(subset(df_highCYT,df_highCYT$`HLA LOH` == T)),nrow(subset(df_highCYT,df_highCYT$`HLA LOH` == F))))
colnames(piechart_df1)=c("HLALOH","ifLOH")
piechart_df1 <- piechart_df1 %>% 
  arrange(desc(HLALOH)) %>%
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

df_lowCYT=subset(df_fullyHeterozygous,df_fullyHeterozygous$highCYT==F)
piechart_df2=data.frame(c("LOH", "no LOH"),c(nrow(subset(df_lowCYT,df_lowCYT$`HLA LOH` == T)),nrow(subset(df_lowCYT,df_lowCYT$`HLA LOH` == F))))
colnames(piechart_df2)=c("HLALOH","ifLOH")
piechart_df2 <- piechart_df2 %>% 
  arrange(desc(HLALOH)) %>%
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

df_highUpregII=data.frame('no.LOH'= piechart_df1$ifLOH[1],'LOH'= piechart_df1$ifLOH[2])
df_lowUpregII=data.frame('no.LOH'= piechart_df2$ifLOH[1],'LOH'= piechart_df2$ifLOH[2])
df_combined=rbind(df_highUpregII,df_lowUpregII)
rownames(df_combined)=c('upreg_HLAII','no.upreg_HLAII')
fisher.test(df_combined,alternative = "two.sided")


# Figure 6B
# scatterplot: CN vs. RNA
# Image: 770 * 470
df_fullyHeterozygous$`Allelic Loss`=df_fullyHeterozygous$totalLoss
ggplot(df_fullyHeterozygous %>% filter(!is.na(totalLoss)), aes(x=meanCN, y=geom_HLAI, color=`Allelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(7), label.y = c(12)) +
  theme_bw() +
  xlab('HLA-I Mean Copy Number') + ylab('HLA-I Expression')

# Figure S6A-D
fit <- survfit(Surv(PFI_time_1,PFI_1) ~ HLA.LOH,data = SKCM_HLALOH)
ggsurv=ggsurvplot(fit, data = SKCM_HLALOH, palette = c("deepskyblue3", "firebrick3"),
                  conf.int = FALSE, 
                  pval = F,  
                  risk.table = TRUE, 
                  title = 'TCGA SKCM',
                  risk.table.col = "strata",
                  ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                             panel.background = element_blank())) + ylab("PFI")
# Image: 760 * 550
ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                            x = c(130, 180, 315), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                            label = c("p = 0.2", "HR = 0.36", "95% CI = 0.09-1.50"), size = 5)
ggsurv
coxph(Surv(PFI_time_1, PFI_1) ~ HLA.LOH, data = SKCM_HLALOH) %>% gtsummary::tbl_regression(exp = TRUE)


fit <- survfit(Surv(PFI_time_1,PFI_1) ~ HLA.LOH,data = BRCA_HLALOH)
ggsurv=ggsurvplot(fit, data = BRCA_HLALOH, palette = c("deepskyblue3", "firebrick3"),
           conf.int = FALSE, 
           pval = F,  
           risk.table = TRUE, 
           title = 'TCGA BRCA',
           risk.table.col = "strata",
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank())) + ylab("PFI")
# Image: 760 * 550
ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                x = c(690, 800, 1360), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                label = c("p = 0.11", "HR = 0.66", "95% CI = 0.40-1.10"), size = 5)
ggsurv
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
# Image: 760 * 550
ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                            x = c(640, 680, 1240), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                            label = c("p = 0.020", "HR = 1.42", "95% CI = 1.06-1.91"), size = 5)
ggsurv
coxph(Surv(PFI_time_1, PFI_1) ~ HLA.LOH, data = LUAD_HLALOH) %>% gtsummary::tbl_regression(exp = TRUE)

fit <- survfit(Surv(PFI_time_1,PFI_1) ~ HLA.LOH,data = df_HLALOH)
ggsurv=ggsurvplot(fit, data = df_HLALOH, palette = c("deepskyblue3", "firebrick3"),
           conf.int = FALSE, 
           pval = F,  
           risk.table = TRUE, 
           title = 'TCGA BRCA, LUAD, SKCM Combined',
           risk.table.col = "strata",
           ggtheme = theme_bw()+theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank())) + ylab("PFI")
# Image: 760 * 550
ggsurv$plot=ggsurv$plot + ggplot2::annotate("text", 
                                            x = c(640, 680, 1240), y = c(0.3, 0.2, 0.1), # x and y coordinates of the text
                                            label = c("p = 0.043", "HR = 1.28", "95% CI = 1.01-1.62"), size = 5)
ggsurv
coxph(Surv(PFI_time_1, PFI_1) ~ HLA.LOH, data = df_HLALOH) %>% gtsummary::tbl_regression(exp = TRUE)


# Figure 6C: facet BRCA LUAD bar graph
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
                   ystar = c(nrow(brca_c1),nrow(luad_c1),nrow(skcm_c1),nrow(brca_c2),nrow(luad_c2),nrow(skcm_c2),nrow(brca_c3),nrow(luad_c3),nrow(skcm_c3),nrow(brca_c4),nrow(luad_c4),nrow(skcm_c4),nrow(brca_c6),nrow(luad_c6),nrow(skcm_c6)) + 10,
                   lab = c(brca_c1_LOHpercent,luad_c1_LOHpercent,skcm_c1_LOHpercent,brca_c2_LOHpercent,luad_c2_LOHpercent,skcm_c2_LOHpercent,brca_c3_LOHpercent,luad_c3_LOHpercent,skcm_c3_LOHpercent,brca_c4_LOHpercent,luad_c4_LOHpercent,skcm_c4_LOHpercent,brca_c6_LOHpercent,luad_c6_LOHpercent,skcm_c6_LOHpercent),
                   Study = c("BRCA", "LUAD", 'SKCM'))
# Image: 870 * 570
ggplot(facet_df, aes(immuneSubtype)) + 
  geom_bar(aes(fill=`HLA LOH`)) +
  scale_fill_viridis(discrete = T, option = "E") +
  facet_wrap(~Study) +
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab)) +
  xlab("") + ylab('Count') +
  ggtitle('Percent HLA LOH among Immune Subtypes')

# Fisher's: C2 vs. C3 HLA LOH
df=df_HLALOH
BRCA_C2=subset(df, df$Subtype_Immune_Model_Based == 'C2' & df$Study == 'BRCA')
BRCA_C3=subset(df, df$Subtype_Immune_Model_Based == 'C3' & df$Study == 'BRCA')
LUAD_C2=subset(df, df$Subtype_Immune_Model_Based == 'C2' & df$Study == 'LUAD')
LUAD_C3=subset(df, df$Subtype_Immune_Model_Based == 'C3' & df$Study == 'LUAD')
SKCM_C2=subset(df, df$Subtype_Immune_Model_Based == 'C2' & df$Study == 'SKCM')
SKCM_C3=subset(df, df$Subtype_Immune_Model_Based == 'C3' & df$Study == 'SKCM')

BRCA_C2_temp=data.frame('no.LOH'= length(which(BRCA_C2$HLA.LOH == F)),'LOH'= length(which(BRCA_C2$HLA.LOH == T)))
BRCA_C3_temp=data.frame('no.LOH'= length(which(BRCA_C3$HLA.LOH == F)),'LOH'= length(which(BRCA_C3$HLA.LOH == T)))
BRCA_df_combined=rbind(BRCA_C2_temp,BRCA_C3_temp)
rownames(BRCA_df_combined)=c('C2','C3')
fisher.test(BRCA_df_combined,alternative = "two.sided")

LUAD_C2_temp=data.frame('no.LOH'= length(which(LUAD_C2$HLA.LOH == F)),'LOH'= length(which(LUAD_C2$HLA.LOH == T)))
LUAD_C3_temp=data.frame('no.LOH'= length(which(LUAD_C3$HLA.LOH == F)),'LOH'= length(which(LUAD_C3$HLA.LOH == T)))
LUAD_df_combined=rbind(LUAD_C2_temp,LUAD_C3_temp)
rownames(LUAD_df_combined)=c('C2','C3')
fisher.test(LUAD_df_combined,alternative = "two.sided")

SKCM_C2_temp=data.frame('no.LOH'= length(which(SKCM_C2$HLA.LOH == F)),'LOH'= length(which(SKCM_C2$HLA.LOH == T)))
SKCM_C3_temp=data.frame('no.LOH'= length(which(SKCM_C3$HLA.LOH == F)),'LOH'= length(which(SKCM_C3$HLA.LOH == T)))
SKCM_df_combined=rbind(SKCM_C2_temp,SKCM_C3_temp)
rownames(SKCM_df_combined)=c('C2','C3')
fisher.test(SKCM_df_combined,alternative = "two.sided")


# Figure 6D: HR plot, high HLAI vs. low HLA I expression ~ LOH
highHLAI_df=subset(LUAD_HLALOH,LUAD_HLALOH$high_HLAI == T)
lowHLAI_df=subset(LUAD_HLALOH,LUAD_HLALOH$high_HLAI == F)

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
  HLA=c('high HLAI','low HLAI'),
  boxOdds = coeff,
  boxCILow = lower_CI,
  boxCIHigh = upper_CI,
  pVal = pVal
)
df$HLA = factor(df$HLA, c('low HLAI','high HLAI'))
stars=stars.pval(pVal)
# Plot
# 610*330
ggplot(df, aes(x = boxOdds, y = HLA)) +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") + 
  geom_point(size = 3.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(0,7,1) ) +
  ylab("") +
  xlab("Odds ratio (log scale)") + ggtitle('TCGA LUAD') +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + coord_cartesian(xlim = c(0, 2.5)) +
  annotate(geom = "text",label = c('HR=1.99\n***p < 0.001','HR=0.90\nns'),y =c(1,2), x = c(0.3,0.3),size=4)


# Figure S6E-G: HLA gene cn vs. expression scatterplot
df_fullyHeterozygous$`HLA-A\nAllelic Loss`=df_fullyHeterozygous$lostA
df_fullyHeterozygous$meanCN_A=(df_fullyHeterozygous$HLA_A1 + df_fullyHeterozygous$HLA_A2)/2
# 696 * 406
ggplot(df_fullyHeterozygous %>% filter(!is.na(lostA)), aes(x=meanCN_A, y=HLA_A_log, color=`HLA-A\nAllelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(7), label.y = c(16)) +
  theme_bw() +
  xlab('HLA-A Mean Copy Number') + ylab('HLA-A Expression')

df_fullyHeterozygous$`HLA-B\nAllelic Loss`=df_fullyHeterozygous$lostB
df_fullyHeterozygous$meanCN_B=(df_fullyHeterozygous$HLA_B1 + df_fullyHeterozygous$HLA_B2)/2
# 696 * 406
ggplot(df_fullyHeterozygous %>% filter(!is.na(lostB)), aes(x=meanCN_B, y=HLA_B_log, color=`HLA-B\nAllelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(7), label.y = c(16)) +
  theme_bw() +
  xlab('HLA-B Mean Copy Number') + ylab('HLA-B Expression')

df_fullyHeterozygous$`HLA-C\nAllelic Loss`=df_fullyHeterozygous$lostC
df_fullyHeterozygous$meanCN_C=(df_fullyHeterozygous$HLA_C1 + df_fullyHeterozygous$HLA_C2)/2
# 696 * 406
ggplot(df_fullyHeterozygous %>% filter(!is.na(lostC)), aes(x=meanCN_C, y=HLA_C_log, color=`HLA-C\nAllelic Loss`)) + 
  geom_point() + geom_smooth(method=loess,se=F) + 
  scale_color_viridis(discrete = F, option = 'B') +
  stat_cor(method = "spearman", label.x = c(7), label.y = c(16)) +
  theme_bw() +
  xlab('HLA-C Mean Copy Number') + ylab('HLA-C Expression')


# Figure S6H-I
# facet LOH ~ expression
# Image: 750 * 480
ggplot(df_fullyHeterozygous %>% filter(!is.na(`HLA LOH`)), aes(x= `HLA LOH`, y=geom_HLAI)) + 
  geom_violin(aes(fill=`HLA LOH`)) + 
  scale_fill_manual(values=c("deepskyblue4", "firebrick")) + 
  stat_compare_means(label.x = 1, label.y = 12) +
  geom_boxplot(width=0.1) + 
  facet_wrap(~Study) +
  theme_bw() +
  ylab('HLA-I Expression')

# SKCM has no available normal sampls
# 500 * 480
ggplot(df_fullyHeterozygous %>% filter(!is.na(`HLA LOH`) & df_fullyHeterozygous$Study != 'SKCM'), aes(x= `HLA LOH`, y=geom_HLAI_upregulation)) + 
  geom_violin(aes(fill=`HLA LOH`)) + 
  scale_fill_manual(values=c("deepskyblue4", "firebrick")) + 
  stat_compare_means(label.x = 1, label.y = 1.6) +
  geom_boxplot(width=0.1) + 
  facet_wrap(~Study) +
  theme_bw() +
  ylab('HLA-I Fold Change')




