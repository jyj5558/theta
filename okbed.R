# estimate number of sites with ok mappability and no repeats
#in R

a <- read.table("ok.bed", header = F)
# the three columns represent sites with mappability == 1

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
ok <- sum(a$diff)

write.table(ok, file = "okbed.txt", row.names = F, col.names = F)

# END
