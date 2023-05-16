library(coloc)
source('Coloc_locus_zoom2.R')
#if(!exists('l2')){
#    l2<-read.table('LocusZooms/Example.ld',as.is=T,header=T)
#}
#if(!exists('l3')){
#    l3<-read.table('LocusZooms/UCSC_GRCh37_Genes_UniqueList2021.txt',as.is=T,header=T,sep='\t')
#}

#if(!exists('genetic.map')){
#    genetic.map<-read.table('/Users/shu/Dropbox/MD/coloc/data/genetic_map_chr16_GRCh37.txt',as.is=T,header=T)
#}
#leadSNP = "rs11646512"
#secondary = c("rs2540781", "rs8053279")
#if(!exists('data.list')){
#    data.list<-readRDS('test_suise.rds')
#}

#if(!exists('ld.mat')){
#    ld.mat<-readRDS('ld_matrix_generate.rds')
#}

generate.coloc.list<-function(data.list=data.list,lp=ld.mat){
    lp.od<-lp$SNP_B
    nmt<-length(lp.od)
    mat<-diag(1,nmt)
    mat[1,2:nmt]<-lp$R2[-1]
    for(i in 2:(nmt-1)){
        for(j in (i-1):1){
            mat[i,j]<-mat[j,i]
        }
        for(j in (i+1):nmt){
            mat[i,j]<-mat[i-1,i]*mat[i-1,j]
        }
    }
    rownames(mat)<-colnames(mat)<-lp.od
    set.seed(2021)
    MAF<-c(0.12*c(1,0.9,0.84,0.63),runif(nmt-4,0.05,0.5))
    names(MAF)<-colnames(mat)
    d1<-data.list[[1]]
    d2<-data.list[[2]]
    d1<-d1[match(lp.od,d2$SNP),]
    d2<-d2[match(lp.od,d2$SNP),]
    dv1<-list(beta=d1$BETA,varbeta=d1$sd^2,N=5000,sdY=1,type='quant',MAF=MAF,LD=mat,snp=lp.od,position=d1$BP)
    dv2<-list(beta=d2$BETA,varbeta=d2$sd^2,N=7000,sdY=1,type='quant',MAF=MAF,LD=mat,snp=lp.od,position=d2$BP)
    dall<-list(D1=dv1,D2=dv2)
}

make.abf.sensitivity<-function(coloc.res,filename){
    png(filename)

    dev.off()
}

plot.lz<-function(des=data.list,ld.dt,genetic.map,GeneReg,SNP1,secondary=NA,filename='Test_suise.jpg'){
    locus.zoom(data = des, snp = SNP1, ld.file = ld.dt, offset = 200000, genes.data = GeneReg, plot.title = c("BMI1","BMI2","BMI3"), nominal = 6, significant = 7.3, file.name = filename, secondary.snp = secondary, population = "EUR", sig.type = "P",nplots=T,geneticMAP=genetic.map)
}
