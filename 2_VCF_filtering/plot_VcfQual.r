#!/usr/bin/Rscript

#########################################
# Usage: Rscript plot_VcfQualiControl.r -i <file prefix, Obtained from vcftools calculation programme "cmd.VcfQualiContr.sh">
#
# The R script makes vcf quality control related plots and basic statistics from the file lists of vcftools calculation results.
# like samples.frq; samples.het; samples.idepth; samples.imiss; samples.ldepth.mean; samples.lmiss
#
# Authors Tianpeng Wang, July 2020, wangtianpeng19@hotmail.com
########################################

# load the relevant packages
library(optparse)
library(readr)
library(dplyr)
library(ggplot2)

# Read in the arguments
option_list <-  list(
  make_option(c("-i","--input"),type = "character",default = NULL,help = "prefix name of the vcf quality control files",metavar = "character"),
  make_option(c("-o","--output"),type = "character",default = "default",help = "prefix name of outputfiles",metavar = "character")
)
opt_parse <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parse)

# check the required arguments
if (is.null(opt$input)) {
  print_help(opt_parse)
  stop("please provide the vcf quality control input name prefix",call. = F)
}

# set the output prefix as the input if it is not given.
if(opt$output=="default"){opt$output=opt$input}
prefix <- opt$input


## 1. variant quality plot (phred encoded)
var_qual <- read_delim(paste0(prefix,".lqual"),delim = "\t",col_names = c("chr","pos","qual"),skip = 1)
p_qual <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+theme_light()
ggsave(p_qual,filename = paste0(prefix,".lqual",".png"))


## 2. mean depth of each variants
var_depth <- read_delim(paste0(prefix,".ldepth.mean"),delim = "\t",col_names = c("chr","pos","mean_depth","var_depth"),skip = 1)
p_depth <- ggplot(var_depth,aes(mean_depth))+geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+theme_light()+xlim(0,50)
ggsave(p_depth,filename = paste0(prefix,".ldepth.mean",".png"))
### rule of thumb:  5 and 95% quantiles, paralogous or repetitive regions

## 3. Variant missingness
var_miss <- read_delim(paste0(prefix,".lmiss"),delim = "\t",col_names = c("chr","pos","nchr","nfiltered","nmiss","fmiss"),skip = 1)
p_miss <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+theme_light()
ggsave(p_miss,filename = paste0(prefix,".lmiss",".png"))
### rule of thumb: 75-95%

## 4. Minor allele frequency distribution plot
var_freq <- read_delim(paste0(prefix,".frq"),delim="\t",col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))
p_freq <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)+theme_light()
ggsave(p_freq,filename = paste0(prefix,".frq",".png"))

# individual based statistics
## 5. Mean depth per individual
ind_depth <- read_delim(paste0(prefix,".idepth"),delim = "\t",col_names = c("ind", "nsites", "depth"),skip = 1)
p_idepth <- ggplot(ind_depth, aes(depth))+geom_histogram(fill="dodgerblue1",colour="black",alpha=0.3)+theme_light()
ggsave(p_idepth,filename = paste0(prefix,"idepth",".png"))


## 6. missing data per individual
ind_miss  <- read_delim(paste0(prefix,".imiss"),delim = "\t",col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"),skip = 1)
p_imiss <- ggplot(ind_miss,aes(fmiss))+geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)+theme_light()
ggsave(p_imiss,filename = paste0(prefix,"imiss",".png"))


## 7. heterozygosity and inbreeding coeffient per individual
ind_het <- read_delim(paste0(prefix,".het"),delim = "\t",col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
p_ihet <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)+theme_light()
ggsave(p_ihet,filename = paste0(prefix,"ihet.png"))


# calculate the basic statistics.
statis <- function(x){
  ### quantiles calculation
  probs <- c(0.001,0.01,0.05,seq(0.1,0.9,by=0.1),0.95, 0.99,0.999)
  out_stat <- vector(mode = "list",length = length(probs))
  for (i in seq_along(probs)){
    out_stat[[i]]=quantile(x,probs = probs[[i]])
  }
  out_stat <- bind_cols(out_stat)
  probname <- paste("quantile",probs,sep = "_")
  colnames(out_stat) <- probname
  ### others calculation
  tab1 <- tibble(
         min=min(x),
         max=max(x),
         mean=mean(x),
         median=median(x)
         )
  bind_cols(out_stat,tab1)
}

## for all seven quality parameters

quali <- c("var_qual","var_depth","var_miss","var_freq","ind_depth","ind_miss","ind_het")

out_quali <- vector("list",length = length(quali))
out_quali[[1]] <- statis(var_qual$qual)
out_quali[[2]] <- statis(var_depth$mean_depth)
out_quali[[3]] <- statis(var_miss$fmiss)
out_quali[[4]] <- statis(var_freq$maf)
out_quali[[5]] <- statis(ind_depth$depth)
out_quali[[6]] <- statis(ind_miss$fmiss)
out_quali[[7]] <- statis(ind_het$f)

out_quali <- bind_rows(out_quali)
out_quali <- mutate(out_quali,name=quali) %>% select(name,everything())

write_tsv(out_quali,path = paste0(prefix,"statistics.tsv"),col_names = T)




