setwd("C:/Users/11146/Desktop/NAFLD2.0/9 Tien S.D 2020")
otutab<- read.delim('otu_table_rare.txt', header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F)
metadata <- read.delim('9group1.txt', header=T, row.names= 1, sep="\t",comment.char = "")
taxonomy <- read.delim('taxonomy1.txt', header=T, row.names= 1, sep="\t",comment.char = "", stringsAsFactors = F)

format2lefse <- function(otutab, taxonomy, metadata, thre = 0.01, groupID = "Group", output = "LEfSe.txt"){

  # 数据交叉筛选
  idx = rownames(otutab) %in% rownames(taxonomy)
  otutab = otutab[idx,]
  tax = taxonomy[rownames(otutab),]
  
  # 依赖关系检测与安装
  p_list = c("dplyr")
  for(p in p_list){
    if (!requireNamespace(p)){
      install.packages(p)}
    library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}
  
  # 标准化，并筛选高丰度菌均值最小百万分之一0.0001%
  norm = t(t(otutab)/colSums(otutab,na=T))*100
  # dim(norm)
  # colSums(norm)
  idx = rowMeans(norm) > 0.01
  HA = norm[idx,]
  # dim(HA)
  # colSums(HA)
  # 数据筛选并排序，要求每个OTU必须的注释，可以为空
  tax = tax[rownames(HA),]
  
  # 转换为等级|连接格式
  tax$Phylum=paste(tax$Kingdom,tax$Phylum,sep = "|")
  tax$Class=paste(tax$Phylum,tax$Class,sep = "|")
  tax$Order=paste(tax$Class,tax$Order,sep = "|")
  tax$Family=paste(tax$Order,tax$Family,sep = "|")
  tax$Genus=paste(tax$Family,tax$Genus,sep = "|")
  tax$Species=paste(tax$Genus,tax$Species,sep = "|")
  # head(tax)
  
  # 按Kingdom合并
  grp <- tax[rownames(tax), "Kingdom", drop=F]
  merge=cbind(HA, grp)
  HA_Kingdom = merge %>% group_by(Kingdom) %>% summarise_all(sum)
  colnames(HA_Kingdom)[1]="Class"
  
  # 按Phylum合并
  grp <- tax[rownames(tax), "Phylum", drop=F]
  merge=cbind(HA, grp)
  HA_Phylum = merge %>% group_by(Phylum) %>% summarise_all(sum)
  colnames(HA_Phylum)[1]="Class"
  
  # 按Class合并
  grp <- tax[rownames(tax), "Class", drop=F]
  merge=cbind(HA, grp)
  HA_Class = merge %>% group_by(Class) %>% summarise_all(sum)
  colnames(HA_Class)[1]="Class"
  
  # 按Order合并
  grp <- tax[rownames(tax), "Order", drop=F]
  merge=cbind(HA, grp)
  HA_Order = merge %>% group_by(Order) %>% summarise_all(sum)
  colnames(HA_Order)[1]="Class"
  
  # 按Family合并
  grp <- tax[rownames(tax), "Family", drop=F]
  merge=cbind(HA, grp)
  HA_Family = merge %>% group_by(Family) %>% summarise_all(sum)
  colnames(HA_Family)[1]="Class"
  
  # 按Genus合并
  grp <- tax[rownames(tax), "Genus", drop=F]
  merge=cbind(HA, grp)
  HA_Genus = merge %>% group_by(Genus) %>% summarise_all(sum)
  colnames(HA_Genus)[1]="Class"
  
  # 合并6个分类级
  all = rbind(HA_Kingdom, HA_Phylum, HA_Class, HA_Order, HA_Family, HA_Genus)
  
  # 将选定的分组列统一命名为group
  metadata$group = metadata[, groupID]

  # 修改样品名为组名，用实验设计中的group列替换样品名
  colnames(all)[2:dim(all)[2]] = as.character(metadata[colnames(all)[2:dim(all)[2]],]$group)
  
  # 保存结果为lefse
  write.table(all, file=paste(output, sep = ""), append = FALSE, sep="\t", quote=F, row.names=F, col.names=T)
}

format2lefse(otutab=otutab,taxonomy=taxonomy,metadata=metadata)
