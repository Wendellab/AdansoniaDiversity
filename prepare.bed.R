df <- read.table("bombax.bed", sep="\t")
df$V5 <- paste(df$V4, df$V1, sep="-")
bed <- (merge((aggregate(V2~V5, df, min)),(aggregate(V3~V5, df, max)), by = "V5"))
bed$diff <- bed$V3-bed$V2
bed<-bed[bed$diff<10001,]
bed$diff <- NULL
bed$name <- gsub("-.*","",bed$V5)
write.table(bed, file="bombax.genes.bed", sep="\t", quote=F)
