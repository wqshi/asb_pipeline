#Given a TF read depth file, call the ASB events.

source('s_function.R')
source('test/t_get_epi_from_database3.r')
source('s_library.R')
source('s_file_function.R')
source('s_mlr_learn_fun.r')
library('futile.logger')
library(logging)
library(Biostrings)

library("optparse")
option_list = list(
    make_option(c("--base_dir"),      type="character", default=NULL, help="the dir name of the data", metavar="character"),
    make_option(c("--data_version"),      type="character", default='', help="the sub folder name under the data_dir", metavar="character"),    
    make_option(c("--cell"),      type="character", default=NULL, help="Name of the cell", metavar="character"),
    make_option(c("--tf"),      type="character", default=NULL, help="Name of tf", metavar="character"),
    make_option(c("--dp_threshold"),      type="integer", default=10, help="Minimum read depth", metavar="character"),
    make_option(c("--fdr"),      type="double", default=0.05, help="fdr for the ASB call", metavar="character")
    
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
str(opt)

if(TRUE){
    #Test case1
    base_dir = "/homed/home/shi/projects/asb_pipeline//data/raw_data/fastq/"
    dp_threshold = 10
    fdr = 0.05
    db_version = 'asb_q30'
    cell_loc = 'gm12878'
    tf = 'ctcf'


    
}else if(is.null(opt$base_dir)){
    base_dir = "./data/test/"
    dp_threshold = 10
    fdr = 0.05
    db_version = 'asb_q30'
    cell_loc = 'gm12878'
    tf = 'batf'

}else{
    base_dir=opt$base_dir
    dp_threshold = opt$dp_threshold
    db_version= opt$db_version
    fdr = opt$fdr
    cell_loc = opt$cell
    tf = opt$tf
}
lab_list = c('uw','sydh','uta','haib','embl')
diffpeak_method = 'PePr' #Not used any more.
loc_rep = '.'


guest_data = f_get_epi_from_paird_database(cell_loc, cell_loc, tf, base_dir, dp_threshold = opt$dp_threshold, opt$fdr, strong_allele_exchange = FALSE,
                                               cell_filter = cell_loc, het_filter=c("het"), add_het = TRUE,
                                               db_data = NULL, rep=loc_rep, labs= lab_list, target_lab = NA, diffpeak_method = diffpeak_method,
                                           debug=F)


write.table(x=guest_data, file = f_p('%s/%s-%s.asb.txt',base_dir, cell_loc, tf ), quote=F, row.names = F, sep ='\t' )
##write.table(x = subset(combined_data, ! tf %in% c('runx3')), file = './data/oriol/combined.20160715.txt', quote = F, sep = '\t', row.names = F)
























