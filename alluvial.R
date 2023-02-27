rm(list=ls())
setwd("C:/Users/11146/Desktop")
Class <- read.table('纲水平丰度表.txt', sep="\t", header=T, row.names=1)
#Unclassifieds先在excel相加合并相对丰度，再导入R
head(Class)
#microeco包是算好了相对丰度
Class.ave <- apply(Class, 1, FUN=mean)
Class.sort <- cbind(Class, Class.ave)[order(-Class.ave),] #排序
# 选择丰度最高的10个纲 剩下的放入others里
Class.sort <- subset(Class.sort[1:10,], select=-Class.ave)
# 统计others丰度
Class.sort <- rbind(Class.sort, Others=apply(Class.sort, 2, function(x){1-sum(x)}))
# 加一列行名 便于后续的长宽转换
Class.sort <- cbind(ClassID=row.names(Class.sort), Class.sort)
library(reshape2)
# 分组排序，长宽数据转换
Class.sort$ClassID <- factor(Class.sort$ClassID, levels = rev(Class.sort$ClassID))
Class.group <- melt(Class.sort, id.vars="ClassID", variable.name="SampleID", value.name="Abundance")
head(Class.group)
# 导入分组
group <- read.table("group.txt", sep="\t", header=T)
#匹配分组
merged_Class <- merge(group,Class.group,by="SampleID")

#单个样本堆积图
library(wesanderson)
library(colortools)
library(ggplot2)
library(ggpubr)
p <- ggbarplot(merged_Class, x = "SampleID", y="Abundance", color="black", fill="ClassID",
               legend="right", 
               legend.title="Top10 Class", #main="Relative counts per Class",
               font.main = c(14, "bold", "black"), font.x = c(12, "bold"), 
               font.y=c(12,"bold")) + 
  theme_bw() +
  rotate_x_text() + 
  scale_fill_manual(values=c("gray",rev(wheel("#5A8BB4FF")[1:10]))) + # 颜色设置"skyblue3"
  facet_grid(~ group, scales = "free_x", space='free') + 
  labs(x = "Sample", y = "Relative Abundance") + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(face = "bold"), 
        plot.title = element_text(face = "bold"), 
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "italic")) 
p
ggsave("Class.pdf", p, units = "cm", width = 20, height = 12, dpi = 600)

#计算每组的平均丰度
library(dplyr)
merged_Class_ave <- merged_Class %>%
  group_by(ClassID,group) %>%#对数据进行分组，分组变量分别为ClassID和group
  summarize(mean_abundance = mean(Abundance))#用summarize函数计算均值
head(merged_Class_ave)
#变回矩形，这里与下面无关
#library(reshape2)
#merged_Class_ave <- dcast(merged_Class_ave, ClassID ~ group,value.var = "mean_abundance")
#write.table(merged_Class_ave,"merged_Class_ave.csv",row.names = FALSE, sep = ',', quote = FALSE)
#dcast函数将长数据转成宽数据，ClassID作为新框的行名，group为列命名
#dcast(data, formula, fun.aggregate, value.var)
#formula是公式，用来描述新数据框中行和列的组合，fun.aggregate是对值进行聚合的函数，value.var是需要聚合的列。
#使用~分隔的公式左右两部分，左侧为第一列，右侧作为第二列
#melt函数是将宽数据转换成长数据，melt(data, id.vars, measure.vars, variable.name, value.name)
#id.vars是需要保留的列，measure.vars是需要转换的列，variable.name是新列的名称，value.name是新值的名称。

#整合样本冲击堆积图
library(wesanderson)
library(colortools)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(forcats)
library(shiny)
library(ggThemeAssist)
# 使用 reorder() 函数对 ClassID 变量重新排序，按照 mean_abundance 变量的大小从大到小排序
merged_Class_ave$ClassID <- reorder(merged_Class_ave$ClassID, -merged_Class_ave$mean_abundance)

p <- ggplot(merged_Class_ave,aes(group,mean_abundance,stratum = ClassID, alluvium = ClassID))+
  geom_alluvium(aes(fill = ClassID),color="black",alpha = 0.7,width = 0.5)+ 
  geom_stratum(aes(fill = ClassID),width = 0.5)+
  scale_fill_manual(values=c("skyblue3",rev(wheel("skyblue3")[1:11]))) +# 颜色设置"skyblue3""#5A8BB4FF" 
  labs(x = "group", y = "Relative Abundance", fill="Top10 Class") + 
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.ticks.length=unit(-0.25, "cm"), 
        axis.text.x = element_text(face = "bold",margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(face = "bold",margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+
  theme_bw()+
  theme(legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(face = "italic",size = 6),
        axis.text = element_text(face = "bold"),
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA),
        legend.key.size = unit(0.5, "cm"))+
  guides(fill = guide_legend(reverse = T))#倒序显示
p
ggsave("Class.pdf", p, units = "cm", width = 12, height = 10, dpi = 600)

#简化成函数
library(dplyr)
library(reshape2)
library(wesanderson)
library(colortools)
library(ggplot2)
library(ggpubr)
library(ggalluvial)
library(forcats)
setwd("C:/Users/11146/Desktop")
#导入相对丰度表
#没有比对的菌群用Unclassifieds合并
Class <- read.table('Class_abundance.txt', sep="\t", header=T, row.names=1)
# 导入分组
group <- read.table("group.txt", sep="\t", header=T)

draw_alluvial <- function(Class, group) {
  Class.ave <- apply(Class, 1, FUN=mean)
  Class.sort <- cbind(Class, Class.ave)[order(-Class.ave),] #排序
  # 选择丰度最高的10个纲 剩下的放入others里
  Class.sort <- subset(Class.sort[1:10,], select=-Class.ave)
  # 统计others丰度
  Class.sort <- rbind(Class.sort, Others=apply(Class.sort, 2, function(x){1-sum(x)}))
  # 加一列行名 便于后续的长宽转换
  Class.sort <- cbind(ClassID=row.names(Class.sort), Class.sort)
  
  # 分组排序，长宽数据转换
  Class.sort$ClassID <- factor(Class.sort$ClassID, levels = rev(Class.sort$ClassID))
  Class.group <- melt(Class.sort, id.vars="ClassID", variable.name="SampleID", value.name="Abundance")
  
  #匹配分组
  merged_Class <- merge(group,Class.group,by="SampleID")
  merged_Class_ave <- merged_Class %>%
    group_by(ClassID,group) %>%
    summarize(mean_abundance = mean(Abundance))
  
  merged_Class_ave$ClassID <- reorder(merged_Class_ave$ClassID, -merged_Class_ave$mean_abundance)
  
  p <- ggplot(merged_Class_ave,aes(group,mean_abundance,stratum = ClassID, alluvium = ClassID))+
    geom_alluvium(aes(fill = ClassID),color="black",alpha = 0.7,width = 0.5)+ 
    geom_stratum(aes(fill = ClassID),width = 0.5)+
    scale_fill_manual(values=c("skyblue3",rev(wheel("#5A8BB4FF" )[1:11]))) +
    labs(x = "group", y = "Relative Abundance", fill="Top10 Class") + 
    theme(panel.background = element_blank(),
          axis.line = element_line(colour = "black", size = 0.5),
          axis.ticks.length=unit(-0.25, "cm"), 
          axis.text.x = element_text(face = "bold",margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(face = "bold",margin=unit(c(0.5,0.5,0.5,0.5), "cm")))+
    theme_bw()+
    theme(legend.title = element_text(face = "bold", size = 8),
          legend.text = element_text(face = "italic",size = 6),
          axis.text = element_text(face = "bold"),
          panel.background = element_rect(fill = NA),
          plot.background = element_rect(colour = NA),
          legend.key.size = unit(0.5, "cm"))+
    guides(fill = guide_legend(reverse = T))
  ggsave("Class.pdf", p, units = "cm", width = 12, height = 10, dpi = 600)
  return(p)
}

draw_alluvial(Class,group)
