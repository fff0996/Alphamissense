#만명 gene별로 pathogenic missense variant 추출하기 



library(data.table)
hg38gene <- read.csv("/BiO/Hyein/Alphafold/pathogenicSCORE/hg38/hg38.ensGene.protein.coding.txt",sep="\t")
pathogenicbim <- fread("../bedbimfam/pathogenic.ukb.wes.merge.bim")
for ( i in 1:nrow(hg38gene) ){
geneid <- hg38gene[i,]$ensembl_gene_id
chr <- hg38gene[i,]$chromosome_name
start.position <- hg38gene[i,]$start_position - 1
end.position <- hg38gene[i,]$end_position + 1
geneA <- pathogenicbim[ (pathogenicbim$V1 == chr & pathogenicbim$V4 >= start.position & pathogenicbim$V4 <= end.position),]
if( nrow(geneA) == 0 ){
cat(paste("Gene - ",geneid," don't have pathogenic variants","\n",sep=""),file="error.txt",append=TRUE)
#next
}else{
geneA <- geneA[,c(2,5)]
geneA <- data.frame(geneA)
geneA["SCORE"] <- 1
names(geneA) <- c("SNP","A1","SCORE")

fwrite(geneA,paste("/BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/gene.",geneid,".score",sep=""),quote=FALSE,row.names=FALSE,sep="\t")

scorecmd <- paste("plink --bfile /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/bedbimfam/pathogenic.ukb.wes.merge --extract /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/bedbimfam/white.freq.extract.1E-05.v2 --score /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/gene.",geneid,".score sum --keep /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/tmp_MC4R/ZZZZZZZ.fam --out /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/gene.",geneid,sep="")
try(system(scorecmd),silent=TRUE)
}
}


#gene sampling bedbimfam만들기
for ( i in 1:500 ){
geneid <- hg38gene[i,]$ensembl_gene_id
chr <- hg38gene[i,]$chromosome_name
start.position <- hg38gene[i,]$start_position - 1
end.position <- hg38gene[i,]$end_position + 1
if( file.exists(paste("/BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/gene.",geneid,".profile",sep="")) ){
geneSCORESUM <- read.table(paste("/BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/gene.",geneid,".profile",sep=""),header=T)
genelist <- read.csv(paste("/BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/gene.",geneid,".score",sep=""),sep="\t")
genelist <- data.frame(genelist$SNP)
names(genelist)[1] <- c("SNP")
write.table(genelist,paste("/BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/genevariantlist/",geneid,".variantlist.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
fix <- geneSCORESUM[geneSCORESUM$SCORESUM >=1,]
unfix <- geneSCORESUM[geneSCORESUM$SCORESUM == 0,]
sample_rows <- sample(1:nrow(unfix),10000,replace=FALSE)
tmp <- unfix[sample_rows,]
tmp2 <- rbind(fix,tmp)
tmp2 <- tmp2[,c("FID","IID")]
write.table(tmp2,paste("/BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/genescoresumsampling10000+/",geneid,".sampling10000.txt",sep=""),sep="\t",row.names=F,col.names=F,quote=F)
bedbimfamcmd <- paste("plink --bfile /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/bedbimfam/pathogenic.ukb.wes.merge  --keep /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/genescoresumsampling10000+/",
geneid,".sampling10000.txt"," --extract /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/genevariantlist/",geneid,".variantlist.txt --make-bed --out /BiO/Hyein/Alphafold/pathogenicSCORE/hg38/gene_h2/geneSCORESUM/genescoresumsampling10000+/bedbimfam/",geneid,".sampling",sep="")
system(bedbimfamcmd)
}else{
#next
}
}
