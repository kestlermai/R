setwd("C:/Users/11146/Desktop")
#导入
abundance <- read.delim('PICRUSt2相对丰度表.txt', sep="\t", header=T,row.names = 1)
head(abundance)
# 将绝对丰度转化为百分比形式的相对丰度
abundance <-t(abundance)
abundance <- abundance /rowSums(abundance)
head(abundance)
#导出
write.table(abundance, file ="picrust2相对丰度表.txt", sep ="\t", row.names = T)

# KEGG分面条形图绘制
library(ggplot2)
library(ggsci)
library(ggpubr)
library(scales)
library(ggbreak)
library(openxlsx)

KEGG <- read.xlsx('PICRUSt2.xlsx', sep="\t")

library(reshape2)
swr = function(string, nwrap = 12){
  paste(strwrap(string,width = nwrap),collapse = "\n")
}
swr = Vectorize(swr)
KEGG$level1 <- swr(KEGG$level1)

p <- ggplot(data=KEGG,aes(x=mean,y=level2,fill=group))+
  geom_bar(position=position_dodge(),
           width = 0.8, # position设置柱子位置，width设置柱子宽度
           stat = "identity")+
  geom_errorbar(aes(xmin =mean-std,xmax =mean+std),
                width = 0.3,
                color="black",size=0.75,
                position = position_dodge(0.8))+
  facet_grid(level1 ~ ., scales = 'free_y',space="free_y")+ # scale和space设置分面的轴范围和空间大小，默认每个分面的轴范围和空间大小都一样。
  #scale_fill_manual(values = mypal)+
  #geom_text(aes(x=mean+mean*0.5,label = sig,group=group),
            #vjust = 0.8,size = 5,fontface = "bold")+
  theme_grey()+
  labs(#title="KEGG Chart", # 设置图形标题
    x = "Mean Relative Abundance%",y= "KEGG Level II")+
  theme(legend.position="top",
        legend.title = element_text(face = "bold", 
                                    size =12,color="black"),
        legend.text = element_text(face = "bold", 
                                   size =12,color="black"))+
  theme(axis.title = element_text(face = "bold", 
                                  size = 14,colour = "black"),
        axis.ticks.y = element_blank(),
        axis.line.x =element_blank(),
        axis.text = element_text(face = "bold", 
                                 size = 12,color="black"))+
  theme(#plot.title = element_text(hjust =0.5),# 设置图形标题顶部居中对齐
    panel.grid=element_blank(),
  )+
  theme(strip.text.y = element_text(size = 14,face = "bold",angle=360)) #strip.*参数用于调节分面图形细节，这里设置图分面标题水平放置。
p
p1 <- p+scale_x_break(breaks = c(5, 35), scales = 0.5, ticklabels = c(38))
p1
p2 <- p1 + scale_color_nejm()+scale_fill_nejm()
p2
ggsave("KEGG-yx2.pdf", p2, units = "cm", width = 28, height = 25, dpi = 600)
