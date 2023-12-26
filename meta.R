#install.packages("meta")
library(meta)

setwd("C:/Users/maihuanzhuo/Desktop/NAFLD2.0/metaforest")

meta <- read.delim("meta-NAFLD-control.txt",check.names = FALSE)
meta <- read.delim('C:/Users/maihuanzhuo/Desktop/adverse events.txt',check.names = F)
#两分类变量metabin；连续资料metacont
metaresult <- metabin(event.e, n.e,event.c,n.c,
                      data = meta,
                      byvar = subgroup,
                      print.byvar = F,
                      studlab = paste(Study),
                      #sm = "OR",
                      method="MH",#Mantel-Haenszel、方差倒数Inverse、Peto
                      comb.random = FALSE)#不用随机效应模型进行效应合并，而是固定效应模型

metaresult <- metacont(n.e,mean.e,sd.e,n.c,mean.c,sd.c,
                       data = meta,
                       sm = "SMD",
                       byvar = subDiversity,
                       print.byvar = F,
                       studlab = paste(Study),
                       comb.fixed = F,
                       comb.random = T)
summary(metaresult)

#影响性分析，即每个研究对总估计效应的影响大小
metainf(metaresult, pooled = "random")#fixed

#森林图
settings.meta('revman5')#则可绘制出RevMan 5风格的图表
#settings.meta('JAMA')
pdf("meta-NAFLD-HC-test.pdf", width = 10,height = 30)
pdf("test.pdf", width = 15,height = 10)
p <- forest(metaresult,
            type.study = "square",#diamond,每个研究的形状正方形或者是菱形
            type.random = "diamond",#汇总的形状
            type.subgroup = "diamond",#每个亚组汇总的形状
            col.square = "gray",#col.square表示森林图中方框的颜色；
            col.square.lines = "gray",#col.square表示森林图中外框的颜色；
            col.diamond.random  = "grey",#col.diamond表示森林图中菱形的颜色;
            col.diamond.lines.random  = "black",#col.diamond.lines表示菱形外框的颜色；
            #leftcols = "Study",
            lwd = 2.5,#设置线粗
            xlim = c(0.5, 10),#设置x轴范围
            colgap.forest.left = "3cm",#调整森林图与左边的间隙
            squaresize = 1,#控制方块大小
            digits.pval = 3,#控制P值的有效数
            digits.pval.Q = 3, #控制Q值的有效数
            digits.sd = 2, #控制sd值的有效数
            col.random = "red",#虚线颜色
            col.predict.lines = "red",
            label.e = "NAFLD",#experimental表头更改
            label.c = "Healthy Control",#control表头更改
            label.left = "",#添加左坐标轴标签
            label.right = "",#添加右坐标轴标签
            hetstat = TRUE)#hetatat=TRUE表示汇报异质性。
dev.off()

#发表偏倚的检验由于纳入的研究个数小于10个，在cochrane手册中不建议做test，只要求做一个漏斗图。
pdf("funnel-NAFLD-control.pdf", width = 5, height = 5)
p1 <- funnel(metaresult)
dev.off()
#当研究个数大于等于10个的时候，对于这种四格表资料不建议使用egger检验或begg检验，可以使用Peters检验
metabias(metaresult,method.bias = "Egger")#有subgroup的话会报错

#使用trim and filled或者copas模型等进行校正
tf1 <- trimfill(metaresult, comb.fixed = TRUE)#没意义就用固定效应模型
summary(tf1)
funnel(tf1)

#敏感性分析
metainf(metaresult, pooled = "fixed")#pooled="random"改用随机效应模型
forest(metainf(metaresult), comb.fixed = TRUE)

#meta回归
#回归需要取对数进行
reg_age <- metareg(metaresult, ~ age)#~age为协变量，多个协变量可以用+连在一起，~age+sex+country+.....
bubble(reg_age)


