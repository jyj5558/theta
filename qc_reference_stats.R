# get stats for the reference qc
#in R

#################
# estimate number of sites with no repeats
#################

a <- read.table("nonrepeat.bed", header = F)

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
ok <- sum(a$diff)

write.table(ok, file = "norepeat.txt", row.names = F, col.names = F)

#################
# estimate number of sites with mappability == 1
#################

a <- read.table("mappability.bed", header = F)

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
okmap <- sum(a$diff)

write.table(okmap, file = "okmap.txt", row.names = F, col.names = F)

################
# estimate number of sites with ok mappability and no repeats
#################

a <- read.table("ok.bed", header = F)

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
ok <- sum(a$diff)

write.table(ok, file = "okbed.txt", row.names = F, col.names = F)

# END
