######################################################################################
# CRG 
# Hana SUSAK
# date: 23/06/2016
#------------------------------------------------------------------------------------
# Checking calles from inhouse data
#------------------------------------------------------------------------------------
# 1st Argument:   MAF file   
# 2nd Argument:   Prior probability of having cancer 
# 3th Argument:   Gold set of genes
# 4th Argument:   If plots should be created 
# 5th Argument:   Outpit folder 
######################################################################################
rm(list=ls())

######################################################################################
# handling arguments
if(! 'rjson' %in% installed.packages()) {
   install.packages('rjson', verbose=F, repos='http://cran.us.r-project.org')
}

if(! 'argparse' %in% installed.packages()) {
   install.packages('argparse', verbose=F, repos='http://cran.us.r-project.org', type = "source")
}
suppressPackageStartupMessages(library(argparse))

parser <-  ArgumentParser(formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument("-m", "--maf_file", type="character", help="Input MAF-like file, with additional mandatory VAF column", metavar="file_path", nargs=1, required=TRUE)
parser$add_argument("-p", "--prior_probability", type="double", help="prior probability of cancer (cancer incidence) \\n\ Default value is 5e-05", metavar="numeric_value[0,1]", nargs=1, default=1/20000)
parser$add_argument("-c", "--gold_standard_genes_file", type="character", help="Cancer Gold Standard genes, as one column without header. \\n\ HUGO gene names should be used, same as in MAF file ", metavar="file_path", nargs=1)
parser$add_argument("-f", "--figures", type="logical", help="If figures should be created. Default is set to T (True).", metavar="T/F", nargs=1, default=T)
parser$add_argument("-o", "--output_folder", type="character", help="Path to output folder. If not provided current directory is used.", metavar="folder_path", nargs=1)
parser$add_argument("--genes_length", type="character", help="Path to file with 2 columns (no header) \\n\ First one is Hugo Symbols of gsnes; \\n\ Second one is length of genes (sum of  exons lengts in gene) ", metavar="file_path", nargs=1)
#parser$add_argument("--sequneced_genes_list", type="character", help="Path to file with (first) column of all sequenced genes. No header.", metavar="file_path", nargs=1)


args <- parser$parse_args()

file.maf <- args$maf_file
cancer.prior <- as.numeric(args$prior_probability)
gold.set.file <- args$gold_standard_genes_file
to.plot <- args$figures
out.folder <- args$output_folder
#sequneced.genes.list.file <- args$sequneced_genes_list
genes.length.file <- args$genes_length

# file.maf <- 'CLL_maf_like_file.txt'
# cancer.prior <- 5/100000
# gold.set.file <- 'CLL_Gold_Standard_Genes.txt'
# sequneced.genes.list.file <- NULL #option not specified
# genes.length.file <- NULL #option not specified
# sequneced.genes.list.file <- 'CLL_Genes_length.txt'
# genes.length.file <- 'CLL_Genes_length.txt'
# out.folder <- '~/test_cDriver_tool/' 
# to.plot <- T
######################################################################################
if (is.null(out.folder)){
   out.folder <-  getwd()
} 
setwd(out.folder)

logFile <- paste0(getwd(),paste0('/cDriver_log_file_',format(Sys.time(), format = "%Y-%j-%H%M%S") , '.txt'))
cat(timestamp(), file=logFile, append=FALSE, sep = "\n")

if(! 'ggplot2' %in% installed.packages()) {
   install.packages('ggplot2', verbose=F,  repos='http://cran.us.r-project.org')
}
suppressPackageStartupMessages(library(ggplot2))

if(!'cDriver' %in% installed.packages()) {
   if(!'devtools' %in% installed.packages()){
      install.packages("devtools")
   }      
   suppressPackageStartupMessages(library(devtools))   
   print(paste0("Installing cDriver package"))
   cat(paste0("Installing cDriver package ..."), file=logFile, append=TRUE, sep = "\n")
   install_github("hanasusak/cDriver",quiet=T)   
}
suppressPackageStartupMessages(library(cDriver))

if(! 'RCurl' %in% installed.packages()) {
   install.packages('RCurl', verbose=F,  repos='http://cran.us.r-project.org')
}
suppressPackageStartupMessages(library(RCurl))

link <- getURL("https://raw.githubusercontent.com/hanasusak/cDriver/master/DESCRIPTION")
lines <- unlist(strsplit(link, '\n'))
vers.line <- lines[grepl('^Version:', lines)]
vers.line <- gsub('Version:', '', vers.line)
vers.line <- trimws(vers.line)

if (packageVersion('cDriver') != vers.line) {
   if(!'devtools' %in% installed.packages()){
      install.packages("devtools")
   }      
   suppressPackageStartupMessages(library(devtools))   
   print(paste0("Installing new version of cDriver package"))
   cat(paste0("Installing new version of cDriver package ..."), file=logFile, append=TRUE, sep = "\n")
   install_github("hanasusak/cDriver",quiet=T)     
}


print(paste('Input Maf File: ',file.maf))
print(paste('Prior risk of having cancer is: ',cancer.prior))
print(paste('Plots: ',to.plot))
cat(paste('Input Maf File full path: ',file.maf), file=logFile, append=TRUE, sep = "\n")
cat(paste('Prior risk of having cancer is: ',cancer.prior), file=logFile, append=TRUE, sep = "\n")
cat(paste('Plots: ',to.plot), file=logFile, append=TRUE, sep = "\n")

if(!is.null(gold.set.file)){
   print(paste('Gold set file is: ',gold.set.file))
   cat(paste('Gold set file is: ',gold.set.file), file=logFile, append=TRUE, sep = "\n") 
   gold.set <- read.table(gold.set.file, header=F, sep='\t', quote="")
   gold.set <- as.character(gold.set[,1])   
   print(paste('There are',length(unique(gold.set)), 'genes as gold standard genes.'))
   cat(paste('There are',length(unique(gold.set)), 'genes as gold standard genes.'),
             file=logFile, append=TRUE, sep = "\n")   
} else {
   gold.set <- driver.genes.concensus
   txt <- paste('You did not provide set of genes to be used as Canceer Gold Standard. \n We will use set of ',length(unique(gold.set)), 'genes previously associated with cancer (in character vector driver.genes.concensus which comes as part of cDriver package).')
   print(txt)
   cat(txt, file=logFile, append=TRUE, sep = "\n")     
}


if(!is.null(out.folder)){
   print(paste('Out folder where all results will be written: ',out.folder))
   cat(paste('Out folder where all results will be written: ',out.folder), file=logFile, append=TRUE, sep = "\n")
} else {
   out.folder <- getwd()
   txt <- paste('You did not specify output folder. Results will be written in current directory: ',out.folder)
   print(txt)
   cat(txt, file=logFile, append=TRUE, sep = "\n")    
}


if(!is.null(genes.length.file)){
   sequneced.genes.list <- read.table(genes.length.file, header=F, sep='\t', quote="")
   sequneced.genes.list <- as.character(sequneced.genes.list[,1])   
   genes.length <- read.table(genes.length.file, header=F, sep='\t', quote="", row.names=NULL)
   genes.length <- (genes.length[,1:2])   
   rownames(genes.length) <- genes.length[,1]
   colnames(genes.length)  <- c('Hugo_Symbol','Length')
   txt <- paste('File with genes and it\'s lengths: ',genes.length.file, '. There are ',nrow(genes.length),' genes in the table.', sep="")
   print(txt)
   cat(txt, file=logFile, append=TRUE, sep = "\n")  
} else {
   sequneced.genes.list <- NULL
   txt <- paste('You did not provide list of genes which were sequenced. 
   Union of all protein coding genes in table all.genes.lengths (comes as part of cDriver package) and your genes is taken')
   #print(txt)
   cat(txt, file=logFile, append=TRUE, sep = "\n")  
   
   genes.length <- NULL
   txt <- paste('You did not provide lenght of sequenced genes, so sum of lenghts of all exons for each gene is taken from table all.genes.lengths (comes as part of cDriver package)')
   #print(txt)
   cat(txt, file=logFile, append=TRUE, sep = "\n")  
}
######################################################################################

cat(paste('Loading MAF file ...'), file=logFile, append=TRUE, sep = "\n")   
print(paste('Loading MAF file ...'))
cancer.maf <- read.table(file.maf, header=T, quote="", sep='\t', stringsAsFactors=F)


cat(paste('Calculating CCF ...'), file=logFile, append=TRUE, sep = "\n")   
print(paste('Calculating CCF ...'))
cancer.maf <- CCF(cancer.maf)


cat(paste('Calculating background ...'), file=logFile, append=TRUE, sep = "\n")   
print(paste('Calculating background ...'))
cancer.bcgr <- bcgr.combine(cancer.maf, genes=sequneced.genes.list, lengthGenes=genes.length[,2])


cat(paste('Running Bayesian Risk (Hazard) Inference model  ...'), file=logFile, append=TRUE, sep = "\n")   
print(paste('Running Bayesian Risk (Hazard) Inference model ...'))
res1 <- bayes.risk(cancer.maf, cancer.bcgr, prior.sick=cancer.prior, genes=sequneced.genes.list)

cat(paste('Running Bayesian Driver Inference model  ...'), file=logFile, append=TRUE, sep = "\n")   
print(paste('Running Bayesian Driver Inference model  ...'))
res2 <- bayes.driver(cancer.maf, cancer.bcgr, driver.genes= gold.set, genes=sequneced.genes.list)

cat(paste('Combinig two results ...'), file=logFile, append=TRUE, sep = "\n")   
print(paste('Combinig two results ...'))
res.final <- combine.ranking(list(res1, res2), genes=sequneced.genes.list, min.mut = 2 )
#head(res.final)

out.file <- unlist(strsplit(file.maf, '/'))
out.file <- out.file[length(out.file)]
out.file <- sub("([^.]+)(\\.[[:alnum:]]+$)", "\\1", out.file)
if(!is.null(gold.set.file)){
   out.file2 <- unlist(strsplit(gold.set.file, '/'))
   out.file2 <- out.file2[length(out.file2)]
   out.file2 <- sub("([^.]+)(\\.[[:alnum:]]+$)", "\\1", out.file2)
} else {
   out.file2 <-paste0(length(driver.genes.concensus),'_genes')
}
out.file <- paste(out.file,'_GS_', out.file2,collapse="", sep="")
out.file <- paste(c(getwd(),out.file), sep='/', collapse="/")

cat(paste('Writing results to a file: ', paste0(out.file,'_cDrive.tsv')), file=logFile, append=TRUE, sep = "\n")   
print(paste('Writing results to a file: ', paste0(out.file,'_cDrive.tsv')))
write.table(res.final, file=paste0(out.file,'_cDrive.tsv'), sep='\t', row.names=F, quote=F)   

if(to.plot){
   # Multiple plot function
   #
   # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
   # - cols:   Number of columns in layout
   # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
   #
   # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
   # then plot 1 will go in the upper left, 2 will go in the upper right, and
   # 3 will go all the way across the bottom.
   #
   multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      library(grid)
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
         # Make the panel
         # ncol: Number of columns of plots
         # nrow: Number of rows needed, calculated from # of cols
         layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                          ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
         print(plots[[1]])
         
      } else {
         # Set up the page
         grid.newpage()
         pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
         
         # Make each plot, in the correct location
         for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
         }
      }
   }
   p1 <- plotSamplesMut(sample.mutations=cancer.maf)
   p1 <- p1 + ggtitle('Number of mutations per sample') + theme(plot.title = element_text(size = 20, face = "bold"))
   p2 <- plotSamplesMut(sample.mutations=cancer.maf, fill=T)
   p2 <- p2 + ggtitle('Type of mutation proportion') + theme(plot.title = element_text(size = 20, face = "bold"))
   p3 <- plotMutChange(sample.mutations=cancer.maf)
   p3 <- p3 + ggtitle('Type of mutation change proportion') + theme(plot.title = element_text(size = 20, face = "bold"))
   #p4 <- boxplotCCF.mutations(sample.mutations=cancer.maf, result.df=res.final, color="Variant_Classification")
   p5 <- boxplotCCF.mutations(sample.mutations=cancer.maf, result.df=res.final, color="Damage_score", shape="Variant_Type")
   p5 <- p5 + ggtitle('Top 20 cDriver ranked genes \n Mutations CCF boxplots ') + theme(plot.title = element_text(size = 20, face = "bold"))
   p6 <- boxplotCCF.patients(sample.mutations=cancer.maf, result.df=res.final)
   p6 <- p6 + ggtitle('Top 20 cDriver ranked genes \n Patient-Gene max CCF boxplots ') + theme(plot.title = element_text(size = 20, face = "bold"))
   p7 <- plotStaircase(sample.mutations=cancer.maf, result.df=res.final, allSamples=T)
   p7 <- p7 + ggtitle('Staircase plot for top 40 cDriver ranked genes ') + theme(plot.title = element_text(size = 20, face = "bold"))
   
   cat(paste('Saving plots to PDF: ', paste0(out.file,'_report_figures.pdf')), file=logFile, append=TRUE, sep = "\n")   
   print(paste('Saving plots to PDF: ', paste0(out.file,'_report_figures.pdf')))  
   
   pdf(file=paste0(out.file,'_report_figures.pdf'),width=20, height=15)
   multiplot(p1, p2, p3, p5, p6, p7, cols=2)
   dev.off()
}

cat(paste('Everything done!'), file=logFile, append=TRUE, sep = "\n")   
print(paste('Everything done!'))
cat(timestamp(), file=logFile, append=TRUE, sep = "\n")
