library(openxlsx)
library(ggplot2)
library(ggsci)
setwd("C:/Users/11146/Desktop")
data <- read.xlsx("data.xlsx", check.names = FALSE,sep = "\t")
#facet自定义排序
data$Taxonomy <- factor(data$Taxonomy,levels = c("Phylum","Class","Order","Family","Genus","Species"))
#柱状图
p1 <- ggplot(data,aes(x=mean_relative_abundance,y=ID,fill=Taxonomy))+
  geom_bar(stat="identity", width = 0.8, position = "dodge")+
  labs(x="mean relative abundance(%)",y="")+
  scale_fill_discrete(name="Taxonomy",limits=c("Phylum","Class","Order","Family","Genus","Species"))+
  theme(axis.text.y=element_text(face="italic"),
        panel.background = element_rect(fill="white",color="black"),
        panel.grid.major = element_line(color = "gray90",size=0.1))
#theme(axis.title.y=element_text(face="italic"))
#strip.position = "top",
p1
ggsave("bar-有益菌.pdf", p1, units = "cm", width = 20, height = 35,dpi = 600)

#柱状图facet
p1 <- ggplot(data,aes(x=mean_relative_abundance,y=ID,fill=Taxonomy))+
  geom_bar(stat="identity", width = 0.8, position = "dodge")+
  labs(x="mean relative abundance(%)",y="")+
  scale_fill_discrete(name="Taxonomy",limits=c("Family","Genus"))+
  facet_grid(Taxonomy ~ ., scales = "free_y")+
  theme(axis.text.y=element_text(face="italic"),
        panel.background = element_rect(fill="white",color="black"),
        strip.background = element_rect(fill=c("#99D8CF")),#BD8EC0
        strip.text = element_text(size=9 ,colour = "white"),
        panel.grid.major = element_line(color = "gray90",size=0.1))
  #theme(axis.title.y=element_text(face="italic"))
  #strip.position = "top",
p1
ggsave("bar-facet-科属有益菌.pdf", p1, units = "cm", width = 20, height = 30,dpi = 600)

#LDA气泡图
p2 <- ggplot(data, aes(x=LDA,y=ID)) +
  geom_point(aes(color = FDR, size = Significance))+
  scale_color_gradient(high = 'blue', low = 'red') +
  labs(x="LDA",y="")+
  theme(axis.text.y=element_text(face="italic"),
        panel.background = element_rect(fill="white",color="black"),
        panel.grid.major = element_line(color = "gray90",size=0.1))
p2
ggsave("bubble-有害菌.pdf", p2, units = "cm", width = 20, height = 25,dpi = 600)

#facet气泡图
p2 <- ggplot(data, aes(x=LDA,y=ID)) +
  geom_point(aes(color = FDR, size = Significance))+
  scale_color_gradient(high = 'blue', low = 'red') +
  labs(x="LDA",y="")+
  facet_grid(Taxonomy ~ ., scales = "free_y")+#x ~ .纵向排布  . ~ X横向排布 
  #scales控制坐标轴，fixed固定xy坐标轴均一致，free_y即y轴跟随子图
  #facet_wrap跟facet_grid一样，但wrap可以指定行数跟列数ncol=2,nrow=2
  theme(axis.text.y=element_text(face="italic"),
        panel.background = element_rect(fill="white",color="black"),
        strip.background = element_rect(fill=c("#99D8CF")),#子图标题背景颜色
        strip.text = element_text(size=9 ,colour = "white", face="bold"),#子图标题文字颜色
        panel.grid.major = element_line(color = "gray90",size=0.1))
p2
ggsave("bubble-facet-科属有害菌.pdf", p2, units = "cm", width = 20, height = 30,dpi = 600)

#丰度气泡图
data <- read.xlsx("相对丰度.xlsx", check.names = FALSE,sep = "\t")
p3 <- ggplot(data, aes(x=group,y=ID, z=mean_relative_abundance,fill=group)) +
  geom_point(aes(size= size))+
  scale_color_gradient(high = 'blue', low = 'red') +
  labs(x="",y="")+
  theme(axis.text.y=element_text(face="italic"),
        panel.background = element_rect(fill="white",color="black"),
        panel.grid.major = element_line(color = "gray90",size=0.1))
p3
ggsave("bubble2.pdf", p3, units = "cm", width = 20, height = 25,dpi = 600)


#KEGG
library(ggplot2)
setwd("C:/Users/11146/Desktop/NAFLD2.0/17 Vincent W 2013")
KEGG <- read.xlsx("kegg_level3-pre-nafld.xlsx",sep = "\t")

#Level1太长，实现自动换行
library(reshape2)
swr = function(string, nwrap = 12){
  paste(strwrap(string,width = nwrap),collapse = "\n")
}
swr = Vectorize(swr)
KEGG$group <- swr(KEGG$group)

p <- ggplot(KEGG,aes(mean,level3)) +
  geom_bar(aes(fill = group),stat = "identity",width = 0.6) +
  xlab("Mean Relative Abundance (%)") +
  ylab("KEGG Pathway level3") +
  theme(panel.background = element_rect(fill = "white",colour='black'),
        panel.grid.major = element_line(color = "grey",linetype = "dotted",size = 0.3),
        panel.grid.minor = element_line(color = "grey",linetype = "dotted",size = 0.3),
        axis.ticks.length = unit(0.4,"lines"),
        #axis.ticks = element_line(color='black'),
        #axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size = 8,face = "bold"),
        axis.title.y=element_text(colour='black', size = 10),
        axis.text.x=element_text(colour='black',size = 6 ),
        axis.text.y = element_text(color = "black",size = 5),
        legend.position = "none",
        strip.text.y = element_text(angle = 0,size = 6.5,face = "bold")) +#分面字体
        facet_grid(group~.,scales = "free_y")
p
pal_npg(palette = c("nrc"),alpha = 0.6)(10)
p1 <- p + scale_color_npg()
p2 <- p1 + scale_fill_npg()
p2
ggsave("KEGG_level3-pre-nafld条形图.pdf", p, width = 18, height = 30, units = "cm")

#lefse-LDA
data <- read.xlsx("data.xlsx", check.names = FALSE,sep = "\t")

p <- ggplot(data,aes(x=reorder(ID,LDA),y=LDA,fill=Group)) + 
  geom_bar(position = "dodge",stat = "identity",color = "black",width = 0.8,show.legend=T)+
  labs(y="LDA score(log10)",x="")+#标签名
  #scale_fill_discrete(name="Group",limits=c("Healthy Control","NAFLD"))+
  theme(axis.title.x=element_text(size=10))+
  theme(axis.text.y=element_text(size=10))+
  scale_fill_nejm()+ 
  theme(axis.text.y=element_text(face="italic"),
        panel.background = element_rect(fill="white",color="black"),
        panel.grid.major = element_line(color = "gray90",size=0.1))+
  labs(title = 'Wang J 2017')+
  coord_flip()
p
ggsave("LDA.pdf", p, units = "cm", width = 20, height = 10, dpi = 600)

