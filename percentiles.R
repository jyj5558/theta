# percentile sites to exclude
#in R

# read data and name files
df <- read.table("strict.depthGlobal", header = F)

# coverage
names(df) <- seq(0, 600, by=1)

# 1% and 99% percentile
df1 <- quantile(df, c(.01, .99))
df2 <- round(df1, digits = 0)

write.table(df2, file = "strict.percentile", sep = "\t",row.names = F, col.names = F)

# END

