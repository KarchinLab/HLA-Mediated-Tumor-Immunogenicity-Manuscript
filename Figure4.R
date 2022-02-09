methylation_df = subset(df_immuneLandscape,!is.na(df_immuneLandscape$met_diff_HLAI400))

# Figure 4A: promoter methylation
hypermethylation=unname(quantile(methylation_df$met_diff_HLAI400,na.rm=TRUE)["75%"])
hypomethylation=unname(quantile(methylation_df$met_diff_HLAI400,na.rm=TRUE)["25%"])
methylation_df$highCYT=methylation_df$CYT > unname(quantile(methylation_df$CYT,na.rm=TRUE)["50%"])
highCYT_df=subset(methylation_df,methylation_df$highCYT==T)
highCYT_df$CYT='high'
lowCYT_df=subset(methylation_df,methylation_df$highCYT==F)
lowCYT_df$CYT='low'
df=rbind(highCYT_df,lowCYT_df)

ggplot(df, aes(x=met_diff_HLAI400, y=norm_T_HLAI)) + 
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = CYT)) +
  geom_vline(xintercept = hypomethylation,linetype=2) +
  geom_vline(xintercept = hypermethylation,linetype=2) + 
  theme_bw() +
  scale_fill_manual(values = c("firebrick2", "deepskyblue2")) +
  stat_cor(method = "spearman", label.x = c(-0.23), label.y = c(12)) +
  guides(alpha = F) +
  xlab('HLA I Promoter Methylation') + ylab('HLA I Expression')
ggsave("met_expression.png",height = 4,width =6,dpi = 300)


# Figure 4B: methylation across six immune subtypes
# promoter
my_comparisons <- list( c('C2','C3'), c('C2', 'C1'), c('C2', 'C4') )
methylation_df$`immune subtype`=methylation_df$Subtype_Immune_Model_Based

ggplot(methylation_df %>% filter((methylation_df$Subtype_Immune_Model_Based != '') & !is.na(methylation_df$met_diff_HLAI400)), aes(x= Subtype_Immune_Model_Based, y=met_diff_HLAI400, fill=`immune subtype`)) + 
  geom_violin() + 
  stat_compare_means(label.x = 4, label.y = 0.5) +
  stat_compare_means(comparisons = my_comparisons) + 
  geom_boxplot(width=0.1) + 
  theme_bw() +
  guides(fill = F) +
  scale_fill_brewer(palette="Blues") + 
  xlab('Immune Subtype') + ylab('HLA I Promoter Methylation')
ggsave("met_immuneSubtype.png",height = 4,width =6,dpi = 300)


# Figure 4C: CYT vs. methylation among immunogenic tumors (C2 + C3)
immunogenic=subset(methylation_df,(methylation_df$Subtype_Immune_Model_Based == 'C2' | methylation_df$Subtype_Immune_Model_Based == 'C3') & !is.na(methylation_df$met_diff_HLAI400))
immunogenic$`Immune Subtype`=immunogenic$Subtype_Immune_Model_Based
# promoter
ggplot(immunogenic, aes(x=CYT, y=met_diff_HLAI400, color=`Immune Subtype`)) + 
  geom_point() + geom_smooth(method=loess) +
  scale_color_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_cor(method = "spearman", label.x = c(0.4, 0.4), label.y = c(0.33, 0.3)) +
  xlab('CYT') + ylab('HLA I Promoter Methylation') +
  guides(color = F)
ggsave("scatter_C2_C3.png",height = 4,width =5,dpi = 300)

immunogenic$highCYT=immunogenic$CYT > unname(quantile(immunogenic$CYT,na.rm=TRUE)["50%"])
highCYT_immunogenic=subset(immunogenic,highCYT == T)
lowCYT_immunogenic=subset(immunogenic,highCYT == F)
ggplot(highCYT_immunogenic, aes(x=Subtype_Immune_Model_Based, y=met_diff_HLAI400, fill=`Immune Subtype`)) + 
  geom_violin() + 
  stat_compare_means(label.x = 1, label.y = 0.18)+
  geom_boxplot(width=0.1) + 
  theme_bw() +
  scale_fill_manual(values=c("goldenrod2", "deepskyblue3")) +
  xlab('Immune Subtype') + ylab('HLA I Promoter Methylation') +
  guides(fill = F)
ggsave("met_C2_C3.png",height = 3,width =3,dpi = 300)


# no significance in lowCYT_immunogenic
ggplot(lowCYT_immunogenic, aes(x=Subtype_Immune_Model_Based, y=met_diff_HLAI400, fill=`Immune Subtype`)) + 
  geom_violin() + 
  stat_compare_means(label.x = 1, label.y = 0.3)+
  geom_boxplot(width=0.1) + 
  theme_bw() +
  scale_fill_manual(values=c("goldenrod2", "deepskyblue3")) +
  xlab('Immune Subtype') + ylab('HLA I Promoter Methylation')
ggsave("met_lowCYT_C2_C3.png",height = 5.8,width =5.8,dpi = 300)

# Figure 4D: fishers: HLA-II fold change vs. HLA I methylation
methylation_df = subset(df_immuneLandscape,!is.na(df_immuneLandscape$met_diff_HLAI400))
df=methylation_df
hypermethylation=unname(quantile(df$met_diff_HLAI_body500,na.rm=TRUE)["75%"])
hypomethylation=unname(quantile(df$met_diff_HLAI_body500,na.rm=TRUE)["25%"])
df$hypermethylation=df$met_diff_HLAI400 > hypermethylation
df$hypomethylation=df$met_diff_HLAI400 < hypomethylation

df$highDE_HLAI=df$DE_HLAI > unname(quantile(df$DE_HLAI,na.rm=TRUE)["50%"])
df$highDE_HLAII=df$DE_HLAII > unname(quantile(df$DE_HLAII,na.rm=TRUE)["50%"])

highHLAIIupreg=subset(df,df$highDE_HLAII==T)
lowHLAIIupreg=subset(df,df$highDE_HLAII==F)

piechart_df1=data.frame(c("HLA I\nHypermethylation", "HLA I\nHypomethylation"),c(nrow(subset(highHLAIIupreg,highHLAIIupreg$hypermethylation == T)),nrow(subset(highHLAIIupreg,highHLAIIupreg$hypomethylation == T))))
colnames(piechart_df1)=c("Methyl_state","num")
piechart_df1 <- piechart_df1 %>% 
  arrange(rev(Methyl_state)) %>%
  mutate(prop = num / sum(piechart_df1$num) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
piechart_df1$label=paste(piechart_df1$Methyl_state,'\n',piechart_df1$num,sep='')

ggplot(piechart_df1, aes(x="", y=prop, fill=Methyl_state)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none", plot.title = element_text(size = 15)) +
  geom_text(aes(y = ypos, label = label), color = "white", size=5) +
  scale_fill_brewer(palette="Set1") + ggtitle("high HLA-II fold change (top 50%)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("met_pie.png",height = 4.5,width =4.5,dpi = 300)

piechart_df2=data.frame(c("HLA I\nHypermethylation", "HLA I\nHypomethylation"),c(nrow(subset(lowHLAIIupreg,lowHLAIIupreg$hypermethylation == T)),nrow(subset(lowHLAIIupreg,lowHLAIIupreg$hypomethylation == T))))
colnames(piechart_df2)=c("Methyl_state","num")
piechart_df2 <- piechart_df2 %>% 
  arrange(rev(Methyl_state)) %>%
  mutate(prop = num / sum(piechart_df2$num) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
piechart_df2$label=paste(piechart_df2$Methyl_state,'\n',piechart_df2$num,sep='')

ggplot(piechart_df2, aes(x="", y=prop, fill=Methyl_state)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none", plot.title = element_text(size = 15)) +
  geom_text(aes(y = ypos, label = label), color = "white", size=5) +
  scale_fill_brewer(palette="Set1") + ggtitle("low HLA-II fold change (bottom 50%)") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("met_pie2.png",height = 4.5,width =4.5,dpi = 300)

df_highUpregII=data.frame('demethylation'= piechart_df1$num[1],'methylation'= piechart_df1$num[2])
df_lowUpregII=data.frame('demethylation'= piechart_df2$num[1],'methylation'= piechart_df2$num[2])
df_combined=rbind(df_highUpregII,df_lowUpregII)
rownames(df_combined)=c('upreg_HLAII','no.upreg_HLAII')
fisher.test(df_combined,alternative = "two.sided")






# Figure S4A: gene body methylation
hypermethylation=unname(quantile(methylation_df$met_diff_HLAI_body500,na.rm=TRUE)["75%"])
hypomethylation=unname(quantile(methylation_df$met_diff_HLAI_body500,na.rm=TRUE)["25%"])
methylation_df$highCYT=methylation_df$CYT > unname(quantile(methylation_df$CYT,na.rm=TRUE)["50%"])
highCYT_df=subset(methylation_df,methylation_df$highCYT==T)
highCYT_df$CYT='high'
lowCYT_df=subset(methylation_df,methylation_df$highCYT==F)
lowCYT_df$CYT='low'
df=rbind(highCYT_df,lowCYT_df)

ggplot(df, aes(x=met_diff_HLAI_body500, y=norm_T_HLAI)) + 
  stat_density_2d(geom = "polygon", aes(alpha = ..level.., fill = CYT)) +
  geom_vline(xintercept = hypomethylation,linetype=2) +
  geom_vline(xintercept = hypermethylation,linetype=2) + 
  theme_bw() +
  scale_fill_manual(values = c("firebrick2", "deepskyblue2")) +
  stat_cor(method = "spearman", label.x = c(-0.08), label.y = c(12)) +
  guides(alpha = F) +
  xlab('HLA I Gene body Methylation') + ylab('HLA I Expression')
ggsave("met_expression_body.png",height = 4,width =6,dpi = 300)



# Figure S4B: HLA-I and HLA-II fold change in C2 vs. C3
methylation_df = df_immuneLandscape
immunogenic=subset(methylation_df,methylation_df$Subtype_Immune_Model_Based == 'C2' | methylation_df$Subtype_Immune_Model_Based == 'C3')
immunogenic$`immune subtype`=immunogenic$Subtype_Immune_Model_Based
immunogenic$`HLA I`=immunogenic$DE_HLAI
immunogenic$`HLA II`=immunogenic$DE_HLAII
dat.m1 <- melt(immunogenic,id.vars='immune subtype', measure.vars=c('HLA I','HLA II'))
dat.m1$expression=dat.m1$value

ggplot(dat.m1, aes(x='immuhe subtype', y=expression, fill=`immune subtype`)) +
  geom_boxplot() + facet_grid(. ~variable ) + 
  scale_fill_manual(values = c("goldenrod2", "deepskyblue3")) + 
  stat_compare_means(label.x = 1, label.y = 5) +
  xlab('') + ylab('Log2FC')
ggsave("HLAII_C2_C3.png",height = 5,width =5,dpi = 300)




