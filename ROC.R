library(dplyr)
library(pROC)
library(ggplot2)
library(verification)
library(RColorBrewer)
#install.packages("verification")
setwd("C:/Users/11146/Desktop")

df <- read.table("test.txt",row.names = 1, header = T, sep = "\t")
# 提取除第一列外的所有列名
colnames <- names(df)[-1]
# 将分组列转换为因子类型
df$group <- as.factor(df$group)

# 创建一个空数据框，用于存储结果
result <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(result) <- c("Genus", "AUC", "P.value")

# 循环计算每列对应的AUC值和P值
for(col in colnames) {
  roc <- roc(group ~ df[[col]], df, levels = c("Healthy Control", "NAFLD"))
  auc <- roc$auc
  row <- data.frame(Variable = col, AUC = auc)
  if(length(levels(df$group)) == 2){
    df$result <- ifelse(df$group == "Healthy Control", 0, 1)
    pval <- roc.area(df$result, roc$predictor)$p.value
    row$`P.value` <- pval
  }
  result <- rbind(result, row)
}

# 过滤
filtered_result <- result %>% 
  filter(AUC > 0.7 & P.value < 0.05)
write.csv(filtered_result, "result.csv", row.names = FALSE)  


roc_list <- list()
auc_list <- list()
for (col in colnames) {
  roc <- roc(group ~ df[[col]], df, levels = c("Healthy Control", "NAFLD"))
  auc <- roc$auc
  if (auc > 0.7) {
    df$result <- ifelse(df$group == "Healthy Control", 0, 1)
    pval <- roc.area(df$result, roc$predictor)$p.value
    if (pval < 0.05) {
      roc_list[[col]] <- roc
      auc_list[[col]] <- auc
    }
  }
}

# 去除前缀 "g__"，创建标签
labels <- paste0(sub("^g__", "", names(roc_list)),", AUC = ", round(sapply(auc_list, function(x) round(x, 3)), 3))
head(labels)

p <- ggroc(roc_list, legacy.axes = TRUE,linewidth = 1)+ 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), #添加对角线
               color="red",linewidth = 1,linetype="dashed")+
  ggtitle("ROC Curves of Genera") + xlab("1 - Specificity") + ylab("Sensitivity") +
  theme_minimal()+
  theme_bw()+
  scale_color_brewer(labels = labels, type = "qual", palette = "Set2")+#RColorBrewer::display.brewer.all() 
  theme(legend.position = c(0.75,0.25),legend.key.width = unit(0.5,"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(color = guide_legend(ncol = 1, keyheight = unit(0.8, "lines"), 
                              label.position = "right",label.hjust = 0))+
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 14),
        axis.text = element_text(size = 12), legend.title = element_blank(),
        legend.text = element_text(size = 8))
p
ggsave("roc.png", p, width = 16, height = 12, dpi = 600, units = "cm")




#单个菌ROC
df <- read.table("test.txt",row.names = 1, header = T, sep = "\t")

df$group <- as.factor(df$group)
roc <- roc(group ~ g__.Clostridium._methylpentosum_group,df,levels=c("Healthy Control","NAFLD"))
df$result[df$group =="Healthy Control"] <- 0
df$result[df$group =="NAFLD"] <- 1
df$result
roc.area(df$result,roc$predictor)


plot(roc,print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2),
     grid.col=c("green","red"), max.auc.polygon=TRUE,
     auc.polygon.col="lightblue",print.thres=TRUE)


ggroc(roc, legacy.axes = TRUE,linewidth = 1,color = "blue") +
  #scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0, 0))+#(0,0)相交
  xlab("1-Specificity") + ylab("Sensitivity") +
  theme_minimal() + 
  theme_bw()+
  theme(legend.position = c(0.8,0.3))+
  ggtitle("My ROC curve") + 
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), #添加对角线
               color="red",linewidth = 1)#linetype="dashed"
