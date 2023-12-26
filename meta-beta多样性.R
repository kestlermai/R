#meta-beta多样性
setwd("C:/Users/11146/Desktop")
library(ggplot2)
library(ggsci)
library(ggprism)

df <- read.table('beta多样性统计.txt', header = T, sep = "\t")
df$study <- factor(df$study,levels = rev(c("Monga Kravetz, et al. 2020","Jiang, et al. 2015","Lang, et al. 2021",
                                       "Zhang, et al. 2019","Ahmed, et al. 2021","Baumann, et al. 2021","Wang, et al. 2017",
                                       "Liang, et al. 2022","Pan, et al. 2021","Dong, et al. 2020","Caussy, et al. 2019",
                                       "Wong, et al. 2013","Kordy, et al. 2021")))
df$distance <- factor(df$distance,levels = c("Bray Curtis","Weighted UniFrac","Unweighted UniFrac","Jaccard"))

p <- ggplot(df,aes(x=R2,y=study,fill=distance))+
  geom_point(aes(colour=distance), size = 5, shape = 16,alpha = 0.8, position = position_jitter(height = 0))+
  #scale_fill_locuszoom(palette = c("default"), alpha = 0.7) +
  labs(x = "R^2", y = "Study",title = "PERMANOVA")+
  theme_prism()+
  theme(legend.text = element_text(face = "bold"),legend.spacing = unit(1, "cm"))
p
ggsave("C:/Users/11146/Desktop/beta多样性-PERMANOVA.pdf", p, width = 20, height = 15, dpi = 600, units = "cm")
