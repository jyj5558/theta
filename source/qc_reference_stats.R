pwd = commandArgs(trailingOnly=TRUE)[2]
g_s = commandArgs(trailingOnly=TRUE)[3]
acc = commandArgs(trailingOnly=TRUE)[4]

setwd(pwd)

#################
# estimate number of sites with no repeats
#################



a <- read.table(paste("/scratch/bell/dewoody/theta/", g_s, "/", acc, "_rm/nonrepeat.bed", sep=""), header = F)

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
ok <- sum(a$diff)

write.table(ok, file = paste("/scratch/bell/dewoody/theta/", g_s, "/", "norepeat.txt", sep=""), row.names = F, col.names = F)

#################
# estimate number of sites with mappability == 1
#################

a <- read.table(paste("/scratch/bell/dewoody/theta/", g_s, "/mappability/mappability.bed", sep=""), header = F)

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
okmap <- sum(a$diff)

write.table(okmap, file = paste("/scratch/bell/dewoody/theta/", g_s, "/", "okmap.txt", sep=""), row.names = F, col.names = F)

################
# estimate number of sites with ok mappability and no repeats
#################

a <- read.table(paste("/scratch/bell/dewoody/theta/", g_s, "/", "ok.bed", sep=""), header = F)

# make header
names(a) <- c("chr","start","stop")
 
# subtract $3 - $2 = sites to get the number of OK sites
a$diff <- (a$stop - a$start)

# sum colum 4
ok <- sum(a$diff)

write.table(ok, file = paste("/scratch/bell/dewoody/theta/", g_s, "/", "okbed.txt", sep=""), row.names = F, col.names = F)

# END
