# find thresholds for global depth
       
depth.global <- as.matrix(read.table("strict.depthGlobal", col.names =F))[,1]
depth.global.cumsum <- cumsum(as.numeric(depth.global))

# Find thresholds.
q99 <- depth.global.cumsum[length(depth.global.cumsum)]*0.99
q01 <- depth.global.cumsum[length(depth.global.cumsum)]*0.01
q99.threshold <- min(which(depth.global.cumsum > q99))
q01.threshold <- min(which(depth.global.cumsum > q01))

# Plot.
# Full range.
png('fullrange.png', width = 600, height = 300)
plot(depth.global, type='l',
     xlab = 'Global depth (count)', ylab = 'Number of sites')
abline(v = q99.threshold, col='red', lty=2); abline(v = q01.threshold, col='red', lty=2)
legend('topright', legend = c(paste('1% thres:', q01.threshold, sep=''),
      paste('99% thres: ', q99.threshold, sep='')),  bty = "n")
dev.off()

png('zoomin.png', width = 600, height = 300)
plot(depth.global, type='l',
     xlab = 'Global depth (count)', ylab = 'Number of sites', xlim = c(0, 30))
abline(v = q99.threshold, col='red', lty=2); abline(v = q01.threshold, col='red', lty=2)
legend('topright', legend = c(paste('1% thres:', q01.threshold, sep=''),
      paste('99% thres: ', q99.threshold, sep='')),  bty = "n")
dev.off()

write.table(q99.threshold, file = "q99.threshold.txt", row.names = F, col.names = F)
write.table(q01.threshold, file = "q01.threshold.txt", row.names = F, col.names = F)

# END