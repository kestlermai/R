setwd("D:/Rdocument")
library(Tax4Fun2)
#物种注释
#指定 OTU 代表序列、Tax4Fun2 库的位置、参考数据库版本、序列比对（blastn）线程数等
runRefBlast(path_to_otus = 'rep-seqs.fasta', 
            path_to_reference_data = './Tax4Fun2_ReferenceData_v2', 
            path_to_temp_folder = 'otu_Ref99NR', database_mode = 'Ref99NR', 
            use_force = TRUE, num_threads = 4)
#预测群落功能
#指定 OTU 丰度表、Tax4Fun2 库的位置、参考数据库版本、上步的物种注释结果路径等
makeFunctionalPrediction(path_to_otu_table = 'otu_table.txt', #导入otu原始表，后面丰度标化true
                         path_to_reference_data = './Tax4Fun2_ReferenceData_v2', 
                         path_to_temp_folder = 'otu_Ref99NR', 
                         database_mode = 'Ref99NR', 
                         normalize_by_copy_number = TRUE, 
                         min_identity_to_reference = 0.97, normalize_pathways = FALSE)

#或者
makeFunctionalPrediction(path_to_otu_table = 'KELP_otu_table.txt', 
                         path_to_reference_data = './Tax4Fun2_ReferenceData_v2', 
                         path_to_temp_folder = 'Kelp_Ref99NR', 
                         database_mode = 'Ref99NR', 
                         normalize_by_copy_number = TRUE, 
                         min_identity_to_reference = 0.97, normalize_pathways = TRUE)