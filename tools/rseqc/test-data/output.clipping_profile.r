pdf("output.clipping_profile.pdf")
read_pos=c(0,1,2,3,4,5,6,7,8,9,44,45,46,47,48,49,50)
count=c(16,12,11,8,6,5,1,1,1,1,1,2,2,2,3,4,4)
plot(read_pos,1-(count/40),col="blue",main="clipping profile",xlab="Position of reads",ylab="Mappability",type="b")
dev.off()
