#BiocManager::install("dada2")
#https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html
library(dada2);packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(gridExtra); packageVersion("gridExtra")
library(phyloseq); packageVersion("phyloseq")

path <- "C:/Users/maihuanzhuo/Desktop/XIE/data" 
fns <- list.files(path, pattern = "fastq.gz", full.names = TRUE)
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
path.out <- "C:/Users/maihuanzhuo/Desktop/XIE/data/Figures/"
path.rds <- "C:/Users/maihuanzhuo/Desktop/XIE/data/RDS/"
theme_set(theme_bw())

#Remove primers and orient reads:
nops <- file.path(path, "noprimers", basename(fns))
prim <- removePrimers(fns, nops, primer.fwd = F27, primer.rev = dada2:::rc(R1492), orient = TRUE)

#Higher loss here than in the mock communities, but still very reasonable. 
#Probably not surprising its a bit higher given this is a real sample rather than mock DNA, 
#and the DNA extraction step was done by CoreBiome without special care to ensuring higher molecular weight DNA.
#这里的损失比模拟社区更高，但仍然非常合理。这可能并不奇怪，因为这是一个真实的样本而不是模拟DNA，
#而且DNA提取步骤是由CoreBiome在没有特别注意确保更高分子量DNA的情况下完成的。

#Inspect length distribution.检查长度分布。
lens.fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

#Filter:
filts <- file.path(path, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ = 3, #默认为0。截断后，包含质量分数小于minQ的读取将被丢弃。
                       minLen = 995, 
                       maxLen = 1800, 
                       maxN = 0, #截断后，超过maxN个Ns的序列将被丢弃。注意，dada不允许Ns。
                       rm.phix = FALSE, #如果为TRUE，则丢弃与isPhiX确定的phiX基因组匹配的读数。
                       maxEE = 2#截断后，“预期错误”高于maxEE的读取将被丢弃。预期误差根据质量分数的标称定义计算：EE=总和（10^（-Q/10））
                       )
track
dim(track)
## Run DADA2
#DADA2算法是一种分裂式分割算法。首先，将每个reads全部看作单独的单元，sequence相同的reads被纳入一个sequence，
#reads个数即成为该sequence的丰度（abundance）；其次，计算每个sequence丰度的p-value，当最小的p-value低于设定的阈值时，
#将产生一个新的partition。每一个sequence将会被归入最可能生成该sequence的partition；最后，依次类推，完成分割归并。
#Dereplicate:
drp <- derepFastq(filts, verbose = TRUE)

#通过在样本推断和误差率估计之间交替来学习误差率，直到收敛。样本推断由dada函数执行。
#错误率估计由errorEstimationFunction执行。此函数的输出用作作为err参数的dada函数调用的输入。
#Learn errors:
err <- learnErrors(drp, errorEstimationFunction = PacBioErrfun, 
                   BAND_SIZE = 32, #BAND_SIZE参数指定了在去噪算法中使用的最大嵌合物大小（以bp为单位）
                   multithread = TRUE)
saveRDS(err, file.path(path.rds, "err.rds"))

#Inspect errors:
plotErrors(err)
plotErrors(err, nominalQ = TRUE)
#图形分布对称，则表明降噪质量好

#Denoise:
dd <- dada(drp, err = err, BAND_SIZE = 32, multithread = TRUE)
saveRDS(dd, file.path(path.rds, "dd.rds"))

#生成一个数据框，其中包含了CCS序列、引物序列、经过质量过滤的序列数量以及每个样本的去噪后的序列数量。
#评估dada2流程中每个步骤的效果和样本的质量.
#Read tracking:
after_dada2 <- cbind(ccs = prim[,1], primers = prim[,2], filtered = track[,2], 
                     denoised = sapply(dd, function(x) sum(x$denoised)))
write.csv(after_dada2, file = "C:/Users/maihuanzhuo/Desktop/XIE/data/after_dada2.csv")

#Sequence table:构建ASV表
st <- makeSequenceTable(dd)
dim(st)

## Taxonomy and Chimeras嵌合体
#http://benjjneb.github.io/dada2/training.html，下载参考数据库
#Assign taxonomy:
tax <- assignTaxonomy(st, "C:/Users/maihuanzhuo/Desktop/R包/silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE) # Slowest part
# 完成分类后，用参考序列数据包，对应填充数据信息
tax <- addSpecies(tax, "C:/Users/maihuanzhuo/Desktop/R包/silva_species_assignment_v138.1.fa.gz")
# 输出taxa文件
write.csv(tax,file="C:/Users/ASUS/Desktop/R/rawdata/filtered/taxa.CSV",append = FALSE, quote = FALSE , 
          sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE, 
          qmethod = c("escape", "double"),fileEncoding = "")
# 另存物种注释变量，去除序列名，只显示物种信息
# Removing sequence rownames for display only
tax.print <- tax
rownames(tax.print) <- NULL
head(tax.print)
# 输出tax.print文件
write.csv(taxa.print,file="C:/Users/ASUS/Desktop/XIE/data/filtered/taxa.print.CSV",append = FALSE, 
          quote = FALSE , sep = " ",eol = "\n", na = "NA", dec = ".", row.names = TRUE,
          col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")
#将所有包含"Escherichia/Shigella"的字符串替换为"Escherichia"
#tax[,"Genus"] <- gsub("Escherichia/Shigella", "Escherichia", tax[,"Genus"]) # Reformat to be compatible with other data sources
head(unname(tax))

#Check chimeras:
bim <- isBimeraDenovo(st, 
                      minFoldParentOverAbundance = 3.5,# Only sequences greater than this-fold more abundant than a sequence can be its "parents".
                      multithread = 1)#Windows不支持'mc.cores' > 1
table(bim)
#bim
#FALSE  TRUE 
#5928   940 
sum(st[,bim])/sum(st)#0.04368147
#去除嵌合体
st.nochim <- removeBimeraDenovo(st, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(st.nochim)
sum(st.nochim)/sum(st)
# 导出生成的ASV表（seqtab.nochim）
write.csv(seqtab.nochim,file="C:/Users/maihuanzhuo/Desktop/XIE/data/filtered/ASV.CSV",
          append = FALSE, quote = FALSE , sep = " ",
          eol = "\n", na = "NA", dec = ".", 
          row.names = TRUE,col.names = TRUE, qmethod = c("escape", "double"),fileEncoding = "")

#Extract sample names from the filenames:
sample.names <- sapply(strsplit(fns, "_"), function(x) paste(x[3], x[4], sep="_"))
sample.names <- gsub(".ccs.fastq.gz", "", sample.names)
rownames(st) <- sample.names
sample.names

#Save processed objects for future analysis.
saveRDS(st, file.path(path.rds, "st.rds"))
saveRDS(bim, file.path(path.rds, "bim.rds"))
saveRDS(tax, file.path(path.rds, "tax_Silva138.1.rds"))

#Reload processed data objects (can run code below from here in absence of input sequence data):
st <- readRDS(file.path(path.rds, "st.rds"))
tax <- readRDS(file.path(path.rds, "tax_Silva138.1.rds"))

## Sample Metadta

#Import the metadata for these samples, which is just the subject and 
#the time ordering of the sample from that subject (only relevant for the two subjects with multiple samples).
ft <- sweep(st, 1, rowSums(st), "/")
df <- read.table("Docs/Fecal_Metadata.csv", header = T, sep = "\t", stringsAsFactors = FALSE)
df$SampleID <- gsub("_", ".", df$X)
df$SampleID <- gsub("^D", "D_", df$SampleID)
df$SampleID <- gsub("^R", "R_", df$SampleID)
rownames(df) <- df$SampleID
df <- df[sample.names2,-1]
head(df)

