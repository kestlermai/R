library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggthemes)
library(ggprism)
library(ggpattern)
library(agricolae)
library(dplyr)
#install.packages("ggbreak")
#install_github("JanCoUnchained/ggunchained")
#BiocManager::install("vioplot")
setwd("C:/Users/maihuanzhuo/Desktop/R/alpha_diversity")
alpha <- read.delim("alpha_diversity.txt",sep="\t",header=T)
head(alpha)
alpha$group <- factor(alpha$group,levels=c("H1","H2","H3"))

#多组间做wilcox.test或者t.test肯定是不合理的，要经过bonferroni校正
data_sig <- aov(Observed ~ group,data = alpha)
summary(data_sig)
sig_lsd <- LSD.test(data_sig,"group",p.adj = "bonferroni")
sig_lsd[[5]] %>%
  as_tibble(rownames = "x") %>%
  left_join(alpha %>%
              group_by(group) %>%
              summarize(Observed_adj = max(Observed)),
            by = c("x" = "group")) -> sig_lsd_adj

#箱线图
p <- ggplot(alpha,aes(x=group,y=Observed),color = "black")+
  stat_boxplot(geom = "errorbar",position = position_dodge(width = 1),width = 0.2,
               color = "black",size = 1) +
  geom_boxplot_pattern(aes(fill=group), size=1, 
                       key_glyph=draw_key_rect, # 图例的形状
                       pattern='stripe', # 线条的样式
                       pattern_spacing = 0.01,# 线条之间的间距
                       pattern_density=0.01, # 线条的密度
                       pattern_angle= 45, 
                       position=position_dodge(width=1), # 对齐方式
                       outlier.shape = NA)+ # 异常值的形状（NA表示不显示）
  geom_text(data = sig_lsd_adj, aes(x = x, y = Observed_adj + 5 ,label = groups),fontface = "bold")+
  scale_fill_nejm() + # 设置颜色映射 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100),
                     breaks = seq(0, 100, 20),guide = "prism_minor") +
  labs(x = "", y = "") +  # 设置x轴和y轴标签
  theme(panel.grid.major.y = element_line(linetype = 3, color = "black"),  # 设置y轴主网格线的样式
        panel.grid.major.x = element_blank(),  # 去除x轴主网格线
        panel.grid.minor = element_blank(),   # 去除次要网格线
        panel.background = element_rect(color = "grey"), 
        axis.text = element_text(color = "black", face = "bold", size = 10), 
        axis.line = element_line(color = "black", size = 0.8), 
        axis.ticks = element_line(color = "black", size = 0.8),  
        axis.ticks.x = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(color = "black", face = "bold", size = 8),
        legend.key.height = unit(0.2, 'in')) +
  guides(fill = guide_legend(override.aes = list(alpha = 1),
                             nrow = 4, byrow = TRUE),  # 设置图例分为4行，并显示所有图例项
         keywidth = unit(4, units = "mm"))  # 设置图例项宽度
p
ggsave("Observed+boxplot.pdf", p, units = "cm", width = 15, height = 10,dpi = 300)

####后面两两比较没有校正，p值假阳性
#Observed
p<-ggplot(alpha,aes(x=group,y=Observed))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Observed.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Chao1
p<-ggplot(alpha,aes(x=group,y=Chao1))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Chao1.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#ACE
p<-ggplot(alpha,aes(x=group,y=ACE))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("ACE.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Shannon
p<-ggplot(alpha,aes(x=group,y=Shannon))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Shannon.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Simpson
p<-ggplot(alpha,aes(x=group,y=Simpson))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Simpson.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#InvSimpson
p<-ggplot(alpha,aes(x=group,y=InvSimpson))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("InvSimpson.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Fisher
p<-ggplot(alpha,aes(x=group,y=Fisher))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Fisher.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#PD
p<-ggplot(alpha,aes(x=group,y=PD))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  #scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  scale_fill_manual(values = c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#5ac7a2","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","Obesity"),
                                      c("Healthy Control","NAFL"),
                                      c("Healthy Control","NASH"),
                                      c("Obesity","NAFL"),
                                      c("Obesity","NASH"),
                                      c("NAFL","NASH")))+
  
  labs(title = 'Kordy K 2021')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("PD.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

