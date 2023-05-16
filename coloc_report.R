library(knitr)
library(markdown)
library(imager)
args<-commandArgs(T)
source('run_suise_coloc.R')
output.dir<-args[1]
author=args[2]
study.name=args[3]
template.Rmd.file=args[4]
dt.file=args[5]
LDmat.file=args[6]
plotLD.file=args[7]
genereg.file=args[8]
genetic.mapfile=args[9]
leadSNP=args[10]

print(paste0("output.dir=",output.dir))
print(paste0("author=",author))
print(paste0("study.name=",study.name))
print(paste0("template.Rmd.file=",template.Rmd.file))
print(paste0("dt.file=",dt.file))
print(paste0("LDmat.file=",LDmat.file))

if(!exists('l2')){
    print(paste0("plotLD.file=",plotLD.file))
    l2<-read.table(plotLD.file,as.is=T,header=T)
}

if(!exists('l3')){
    print(paste0("genereg.file=",genereg.file))
    l3<-read.table(genereg.file,as.is=T,header=T,sep='\t')
}

if(!exists('genetic.map')){
    genetic.map<-read.table(genetic.mapfile,as.is=T,header=T)
}

outputfile=paste0(output.dir,'/coloc_report.pdf')
out.tmp.dir<-paste0(output.dir,'/tmp/')
if(!dir.exists(out.tmp.dir)){
    dir.create(out.tmp.dir)
}
if(!exists('coloc.dt.list')){
    coloc.dt.list<-readRDS(dt.file)
}
if(!exists('ld.dt')){
    ld.dt<-readRDS(LDmat.file)
}
print("running coloc(suise) analysis")
# Knit report using working environment
# Warning: It is quite likely that this will be called within an Rmd file
# implying a recursive call to \code{knit()}. This will generate "duplicate label"
# errors for unlabelled chunks. To avoid this, all code chunks
# in your Rmd file should be named.
# Supposedly this error can also be avoided by setting the following option:
# \code{options(knitr.duplicate.label = 'allow')}.
# input_filename Rmd file.
# output_filename Markdown or HTML output file.  An HTML file
# is specified using the .htm, .html, .HTM or .HTML file extension.
# When html is specified, a similarly named markdown file is also
# generated.
# All output files including cache and figures will appear in the
# same folder as \code{output_filename}.
# ... Arguments to be passed to \code{\link[knitr:knit]{knitr::knit}}
# @return NULL
# @keywords internal
knit_report <- function(input_filename, output_filename, ...)
{
	requireNamespace("knitr", quietly = TRUE)
    requireNamespace("markdown", quietly = TRUE)
    output_filename <- normalizePath(output_filename)

    output_dir <- dirname(output_filename)
    if (!file.exists(output_dir))
        dir.create(output_dir)

    current_dir <- getwd()
    on.exit(setwd(current_dir))
    setwd(output_dir)

    name <- gsub("\\.[^.]+$", "", basename(output_filename))
    suffix <- gsub(".*\\.([^.]+)$", "\\1", output_filename)

    is.html <- tolower(suffix) %in% c("htm","html")
    is.pdf <- tolower(suffix) == "pdf"
    is.docx <- tolower(suffix) %in% c("doc", "docx", "word")
    is.md <- tolower(suffix) %in% c("md", "markdown")

    if (is.html)
        return(knitr::knit2html(input_filename, output=paste0(name, ".html"), envir=parent.frame(), ...))
    else if (is.md)
        return(knitr::knit(input_filename, output=paste0(name, ".md"), envir=parent.frame(), ...))
    else if (is.pdf)
    {        
        requireNamespace("rmarkdown", quietly = TRUE)
        return(rmarkdown::render(input_filename, rmarkdown::pdf_document(), intermediates_dir=getwd(), output_dir=getwd(), output_file=paste0(name, ".pdf"), clean = TRUE, envir=parent.frame(), ...))
    }
    else if (is.docx)
    {        
        requireNamespace("rmarkdown", quietly = TRUE)
        return(rmarkdown::render(input_filename, rmarkdown::word_document(), intermediates_dir=getwd(), output_dir=getwd(), output_file=paste0(name, ".docx"), clean = TRUE, envir=parent.frame(), ...))
    }
    else
        stop("Please choose a filename with pdf, html, docx or md suffix")
}


# Generate coloc report
#
# Using the output from the \code{coloc} function this report will generate a report containing tables and graphs summarising the results.
# A separate report is produced for each exposure - outcome pair that was analysed.
# dat Output from \code{\link{harmonise_data}}
# output_path Directory in which reports should be saved.
# output_type Choose `"html"` or `"md"`. Default is `"html"`.
# All output files including cache and figures will appear in the
# folder specified in \code{output_path}.
# author Author name.
# study Study title.
# path The filepath to the report template.
# ... Extra options to be passed to \code{\link[knitr:knit]{knitr::knit}}.
# @return NULL
coloc_report <- function(output_path = ".", output_type = "pdf", author = "Analyst", study = "Coloc Demo", template.Rmd='./coloc_report.Rmd', ...)
{
    message("Writing report as ", output_type, " file to ", output_path)

    message("Performing analysis")
    gcl<-generate.coloc.list(coloc.dt.list,ld.dt)
    S1<-runsusie(gcl[[1]])
    S2<-runsusie(gcl[[2]])
    c.abf<-coloc.abf(gcl[[1]],gcl[[2]])
    #c.sui<-coloc.susie(gcl[[1]],gcl[[2]])
    c.sui<-coloc.susie(S1,S2)
    ppz=c.abf$summary[["PP.H4.abf"]]
    pps.c1<-c.sui$summary
    pps.c1<-pps.c1[order(pps.c1$PP.H4.abf,decreasing=T),]
    #pps=0.957
    if(all(pps.c1$PP.H4.abf<=0.8)){
        pps=pps.c1$PP.H4.abf[1]
        pps.snp1<-pps.c1[1,'hit1']
        pps.snp2<-pps.c1[1,'hit2']
        pps.snp11.pv<-coloc.dt.list[[1]]$P[coloc.dt.list[[1]]$SNP==pps.snp1]
        pps.snp12.pv<-coloc.dt.list[[2]]$P[coloc.dt.list[[1]]$SNP==pps.snp1]
        pps.snp21.pv<-coloc.dt.list[[1]]$P[coloc.dt.list[[1]]$SNP==pps.snp2]
        pps.snp22.pv<-coloc.dt.list[[2]]$P[coloc.dt.list[[1]]$SNP==pps.snp2]
        vv<-matrix(c(pps.snp11.pv,pps.snp12.pv,pps.snp21.pv,pps.snp22.pv),nrow=1,byrow=T)
        colnames(vv)<-c('snp1.trait1.pv','snp1.trait2.pv','snp2.trait1.pv','snp2.trait2.pv')
        vv<-cbind(data.frame(snp1=pps.snp1,snp2=pps.snp2,stringsAsFactors=F),vv)
    }
    else{
        l.pps<-pps.c1$PP.H4.abf>0.8
        pps.snp1<-pps.c1$hit1[l.pps]
        pps.snp2<-pps.c1$hit2[l.pps]
        print(c(pps.snp1,pps.snp2))
        pps.snp11.pv<-coloc.dt.list[[1]]$P[match(pps.snp1,coloc.dt.list[[1]]$SNP)]
        print(pps.snp11.pv)
        pps.snp12.pv<-coloc.dt.list[[2]]$P[match(pps.snp1,coloc.dt.list[[2]]$SNP)] 
        pps.snp21.pv<-coloc.dt.list[[1]]$P[match(pps.snp2,coloc.dt.list[[1]]$SNP)] 
        pps.snp22.pv<-coloc.dt.list[[2]]$P[match(pps.snp2,coloc.dt.list[[2]]$SNP)] 
        vv<-data.frame(snp1=pps.snp1,snp2=pps.snp2,snp1.trait1.pv=pps.snp11.pv,snp1.trait2.pv=pps.snp12.pv,snp2.trait1.pv=pps.snp21.pv,snp2.trait2.pv=pps.snp22.pv)
        pps=pps.c1$PP.H4.abf[l.pps]

    }
    l.share<-apply(as.matrix(vv[,-(1:2)]),1,function(x)all(x<0.05))
    vv<-vv[l.share,]
    pps=pps[l.share]
    vz<-sapply(1:length(pps),function(i)which(c.sui$summary$hit1==vv$snp1[i] &c.sui$summary$hit2==vv$snp2[i]))
    print(vz) 
    #abf.sense=sensitivity(c.abf,rule='H4>0.8')
    #suise.sense=sensitivity(c.sui,rule='H4>0.8')
    c.abf.tb<-do.call(rbind,list(c.abf$summary))
    c.sui.tb<-c.sui$summary
    m <- list(
        colocAbf = c.abf.tb,
        colocSuise = c.sui.tb
    )

    message("Generating graphs")
    plot.lz(des=coloc.dt.list,ld.dt=l2,GeneReg=l3,genetic.map=genetic.map,SNP1=leadSNP,filename=paste0(out.tmp.dir,'/test_coloc_report.jpg'))
    coloc.p<-load.image(paste0(out.tmp.dir,'/test_coloc_report.jpg'))
    #output_file <- array("", nrow(combinations))
    title<-'Coloc(Suise)_demo'
    output_file <- file.path(output_path, paste(study, sanitise_string(title), output_type, sep="."))
    #output_file[i] <- knit_report(template.Rmd, output_file[i], ...)
    #output_file<-knit_report(template.Rmd,file.path('coloc_report.pdf'), ...)
    output_file<-knit_report(template.Rmd,outputfile, ...)
    return(output_file)
}

sanitise_string <- function(x)
{
    gsub(" ", "_", gsub("[^[:alnum:] ]", "", x))
}

my.report<-coloc_report(output_path=output.dir,author=author,study=study.name,template.Rmd=template.Rmd.file)
q(save='no')
