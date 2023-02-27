library(ggplot2)
library(ggpubr)
#install.packages("ggbreak")
#install_github("JanCoUnchained/ggunchained")
#BiocManager::install("vioplot")
setwd("C:/Users/11146/Desktop/NAFLD2.0/11 Sonja L 2020/alpha_diversity2")
alpha <- read.delim("alpha_diversity.txt",sep="\t",header=T)
head(alpha)
#Observed
p<-ggplot(alpha,aes(x=group,y=Observed))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
                                      #c("Healthy Control","NAFLD-6months later"),
                                      #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Observed.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Chao1
p<-ggplot(alpha,aes(x=group,y=Chao1))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
  #c("Healthy Control","NAFLD-6months later"),
  #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Chao1.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#ACE
p<-ggplot(alpha,aes(x=group,y=ACE))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
  #c("Healthy Control","NAFLD-6months later"),
  #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("ACE.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Shannon
p<-ggplot(alpha,aes(x=group,y=Shannon))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
  #c("Healthy Control","NAFLD-6months later"),
  #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Shannon.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Simpson
p<-ggplot(alpha,aes(x=group,y=Simpson))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
  #c("Healthy Control","NAFLD-6months later"),
  #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Simpson.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#InvSimpson
p<-ggplot(alpha,aes(x=group,y=InvSimpson))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
  #c("Healthy Control","NAFLD-6months later"),
  #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("InvSimpson.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#Fisher
p<-ggplot(alpha,aes(x=group,y=Fisher))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
  #c("Healthy Control","NAFLD-6months later"),
  #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("Fisher.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

#PD
p<-ggplot(alpha,aes(x=group,y=PD))+
  geom_violin(aes(fill=group),show.legend = T,scale = "width",color=NA,alpha=0.5)+
  scale_fill_manual(values = c("#7697CB","#F47E62"))+#"BD8EC0"))+
  #scale_fill_manual(values = c("#7697CB","#F47E62","#BD8EC0"))+
  geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62"))+#,"BD8EC0"))+
  #geom_boxplot(width=0.3,outlier.shape = NA,position = position_dodge(0.9),fill=c("#7697CB","#F47E62","#BD8EC0"))+
  stat_compare_means(method = "wilcox.test",
                     #method ="kruskal.test",
                     label="p.signif",
                     comparisons=list(c("Healthy Control","NAFLD")))+
  #c("Healthy Control","NAFLD-6months later"),
  #c("NAFLD","NAFLD-6months later")))+
  labs(title = 'Sonja L 2020')+
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) 
p
ggsave("PD.PDF", p, units = "cm", width = 15, height = 12,dpi = 1200)

