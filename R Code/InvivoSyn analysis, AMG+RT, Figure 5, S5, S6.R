#Load required packages
library(invivoSyn)
library(tidyverse)

#Plots and synergy analyses for TGI data

##Read in csv file with tumor growth data 
tv <- read_tv("JHUEM2, AMG + RT, D17-26.csv")

##Boxplot with only Day 7 and Day 9 (Figure 5B)
tv_79 <- tv[which(tv$Day== 7 | tv$Day== 9),]
cmps=as.list(as.data.frame(combn(as.character(levels(tv_79$Group)),2)))
n_day=length(unique(tv_79$Day))
col_n=ceiling(sqrt(n_day))
row_n=ceiling(n_day/col_n)
p=ggpubr::ggboxplot(tv_79,x='Group',y='TV',add=c('jitter'),color='Treatment',add.params = list(shape=17),xlab='Day',repel = TRUE)+facet_wrap(~Day,scales = 'free_y')+
  ggpubr::stat_compare_means(method='t.test',comparisons = cmps,aes(label = paste0("p=", ..p.format..)),label.x.npc = 'left',label.y.npc = 'top')+
  theme_Publication()+ggsci::scale_color_npg()

p <- p + labs(y= expression(bold(paste("Tumor Volume (", mm^3, ")"))))

png("Boxplots, D7&9.png", width= 1700, height= 900, units= "px", res= 150)
print(p)
dev.off()

##Generate tumor growth plots (Supplementary Figure 5)
p1=tv %>% ggplot()+aes(Day,TV,group=Mouse,color=Group) + labs(y= expression(bold(paste("Tumor Volume (", mm^3, ")"))))+geom_point()+geom_line()+
  theme_Publication()+facet_wrap(~Treatment)+
  ggsci::scale_color_jco() + theme(axis.title.x=element_text(face = "bold",size = rel(1.6)),axis.text=element_text(size=rel(1.8)),axis.ticks=element_line(size=rel(5)))+
  theme(legend.position = 'none')
ggsave("Tumor growth curve_AMGRT.png",width=10,height=8,dpi=300)


##Generate density plot for HSA Synergy Scores with Day 9 TGI (Figure 5C)
TGI_lst <- getTGI(tv, sel_day = 9)
method <- "HSA"
bsTGI_all=TGI_lst$bsTGI_r$t
colnames(bsTGI_all)=c('D1','D2','Combo')
bsTGI_all = bsTGI_all %>% as.data.frame() %>% mutate('expected_TGI'=switch(method,"Bliss"=D1+D2-D1*D2,"HSA"=pmax(D1,D2),"RA"=D1+D2))
bsTGI_all = bsTGI_all %>% mutate(Synergy_score=Combo-expected_TGI)
bsTGI_all = bsTGI_all*100 #in percentage
pval=round(1-sum(bsTGI_all$Combo>bsTGI_all$expected_TGI)/nrow(bsTGI_all),4)

p2= plot_density(bsTGI_all, sel_var = "Synergy_score", pval= pval)
p2 <- p2 + xlab("HSA Synergy Score") + theme_Publication()

png("Day 9 HSA TGI Density.png", width= 1200, height= 720, units= "px", res= 150)
print(p2)
dev.off()

##Generate dot plot for bootstrap analysis on Day 9 TGI using HSA model (Supplementary Figure 6)
HSA_Synergy <- TGI_synergy(TGI_lst, method = "HSA")


##Generate density plot for Bliss Synergy Scores with Day 9 TGI (Figure 5C)
TGI_lst <- getTGI(tv, sel_day = 9)
method <- "Bliss"
bsTGI_all=TGI_lst$bsTGI_r$t
colnames(bsTGI_all)=c('D1','D2','Combo')
bsTGI_all = bsTGI_all %>% as.data.frame() %>% mutate('expected_TGI'=switch(method,"Bliss"=D1+D2-D1*D2,"HSA"=pmax(D1,D2),"RA"=D1+D2))
bsTGI_all = bsTGI_all %>% mutate(Synergy_score=Combo-expected_TGI)
bsTGI_all = bsTGI_all*100 #in percentage
pval=round(1-sum(bsTGI_all$Combo>bsTGI_all$expected_TGI)/nrow(bsTGI_all),4)

p3= plot_density(bsTGI_all, sel_var = "Synergy_score", pval= pval)
p3 <- p3 + xlab("Bliss Synergy Score") + theme_Publication()

png("Day 9 Bliss TGI Density.png", width= 1200, height= 720, units= "px", res= 150)
print(p3)
dev.off()

##Generate dot plot for bootstrap analysis on Day 9 TGI using Bliss model (Supplementary Figure 6)
Bliss_Synergy <- TGI_synergy(TGI_lst, method = "Bliss")


#HSA AUC Synergy density plots (Figure 5D)
tv_AUC <- read_tv("JHUEM2, AMG+RT, D14-26.csv")
AUC_lst <- get_mAUCr(tv_AUC)
AUC_Synergy <- AUC_synergy(AUC_lst, t= 12, method= "HSA")

  
  
  
  