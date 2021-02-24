#install.packages("dplyr", lib="/scratch/bell/dewoody/rlibs/", repos='http://cran.us.r-project.org')
.libPaths("/scratch/bell/dewoody/rlibs/") 
library(dplyr)

pwd = commandArgs(trailingOnly=TRUE)[2]
g_s = commandArgs(trailingOnly=TRUE)[3]
acc = commandArgs(trailingOnly=TRUE)[4]

setwd(pwd)

#update scaffold names in annotation file to match  assembly names

# read data and name files
df <- read.table(paste("/scratch/bell/dewoody/theta/", g_s, "/", acc, "_gtf/gtf.gtf", sep=""), header = F , comment.char = '#', sep="\t")
colnames(df) = c('ID')

df.name <- read.table(paste("/scratch/bell/dewoody/theta/", g_s, "/", acc, "_gtf/ID.txt", sep=""), header =F)
colnames(df.name) = c('ID','scaffold')

# find replace ID with scaffold name
df1 <- merge(df.name, df, by = "ID") 
df1$ID = NULL
write.table(df1, file = paste("/scratch/bell/dewoody/theta/", g_s, "/", acc, "_gtf/gtf_edited.body", sep=""), sep = "\t",row.names = F, col.names = F)
