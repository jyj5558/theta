# estimate number of sites with mappability == 1
#in R

a <- read.table("mappability.bed", header = F)
# the three columns represent sites with mappability == 1

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
okmap <- sum(a$diff)

write.table(okmap, file = "okmap.txt", row.names = F, col.names = F)

# END
