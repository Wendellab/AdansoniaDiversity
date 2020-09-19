args = commandArgs(trailingOnly=TRUE)

snp <- read.table(args[1], col.names=c("allele","position","GT","ref.reads","alt1.reads","alt2.reads"),fill=TRUE)
length <- read.table(args[2], col.names=c("allele","length"))
allele <- merge(snp,length,by="allele")

numsnp <- as.data.frame(table(allele$allele))
names(numsnp) <- c("allele","var_sites")

variants <- merge(allele, numsnp, by="allele")
variants$perc_ref <- 100*round(variants$ref.reads/rowSums(variants[,c("ref.reads","alt1.reads","alt2.reads")],na.rm=T))
variants$perc_alt1 <- 100*round(variants$alt1.reads/rowSums(variants[,c("ref.reads","alt1.reads","alt2.reads")],na.rm=T))
variants$perc_alt2 <- 100*round(variants$alt1.reads/rowSums(variants[,c("ref.reads","alt1.reads","alt2.reads")],na.rm=T))

variants$perc_variant <- 100*variants$var_sites/variants$length


write.table(variants, file=args[3], quote=F, sep="\t", row.names=F)