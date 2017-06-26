t_asb_calling <- function(){
      #Des: The test file for f_asb_calling
#By Wenqiang, 2016-04-27

	    source('s_asb_meta_data.R')
      #write the test data
      #The data are from t_get_epi_from_database3.
      #write.table( data_replicates, file='./data/test/t_asb_calling.txt')  
      test_data=read.table(file='./data/test/t_asb_calling.txt')
      target_lab = 'uw'
      #lab_cols =
      input_data = test_data
      loc_tf = 'ctcf'
      labs = c('sydh', 'uw', 'haib', 'uta')
      rep = '.'
      lab_pattern= f_p( '.*%s.*_(%s)[1-9]*', loc_tf ,paste(labs,collapse = '|')) 
      rep_pattern = f_p('(%s).*(%s%s)$', loc_tf, target_lab, rep )
      
      f_check_data_features(test_data, necessary = c(paste0('lab_',target_lab), lab_pattern, rep_pattern, 'peak_.*_dis','start', 'peak.*_pwm_ref_score', 'pwm_ref_score'))
      test_data$ref_simulate_dp_broad1
      
      #Filter only one lab data according to the lab name and 
      #necessary_list = c(paste0('lab_',target_lab), )
      #f_check_data_features(db_data, necessary_list)
      
      
      #ASB calling method, update ref_tf_dp, alt_tf_dp, p_value, FDR. Also might subset the db_data.
      lab_cols=grep(lab_pattern, colnames(test_data), value = TRUE)
      keep_cols = grep(rep_pattern, lab_cols, value = T)
      
      
      sum(test_data$ref_simulate_dp_broad1 > 72)/nrow(test_data)
      
      
      #debug onece
      debug(f_asb_calling)
      #undebug(f_asb_calling)
      #mtrace(f_asb_calling)
      
      binomial_data = f_asb_calling(test_data, 'ctcf', target_lab, keep_cols, calling_method = 'binomial', fdr = 0.05, dp_threshold = 10)
      table(binomial_data$ASB)
	    sum( binomial_data$ASB == 'ASB' &  binomial_data$ref_tf_dp/(binomial_data$ref_tf_dp + binomial_data$alt_tf_dp) > 0.6)
	    
	    likelihood_data = f_asb_calling(test_data, 'ctcf', target_lab, keep_cols, calling_method = 'likelihood', fdr = 0.05, dp_threshold = 10)
	    table(likelihood_data$ASB)
      sum( likelihood_data$ASB == 'ASB' &  likelihood_data$ref_tf_dp/(likelihood_data$ref_tf_dp + likelihood_data$alt_tf_dp) > 0.6)
	    
      
	    test_data$ref_simulate_dp_broad1 = NULL
      library(edgeR)
      edgeR_data = f_asb_calling(test_data, 'ctcf', target_lab, keep_cols, calling_method = 'edgeR', fdr = 0.05, dp_threshold = 10)
      
      
	    sum(test_data[rownames(binomial_data), 'ref_simulate_dp_broad1'] > 72)/nrow(binomial_data)
      ASB_data = (subset(binomial_data, ASB == 'ASB'))
      
	    sum(test_data[rownames(ASB_data), 'ref_simulate_dp_broad1'] > 72)/nrow(ASB_data)
      
      checkTrue(sum(binomial_data$ASB == 'ASB') == 683, 'Wrong binomial test number')
      checkTrue(sum(edgeR_data$ASB == 'ASB') == 170, 'Wrong edgeR test number')
      
}



f_edgeR_call <- function(x, group){
  y <- DGEList(counts=x,group=group)
  y <- calcNormFactors(y)
  y <- estimateCommonDisp(y)
  y <- estimateTagwiseDisp(y)
  et <- exactTest(y)
  edgeR_predictions=topTags(et,adjust.method = 'BH',n = nrow(x),)$table
  edgeR_predictions$ref_counts = rowSums(y$pseudo.counts[, grep('ref', colnames(y$pseudo.counts), value = TRUE)])
  edgeR_predictions$alt_counts = rowSums(y$pseudo.counts[, grep('alt', colnames(y$pseudo.counts), value = TRUE)])
  
  return (edgeR_predictions)
  
}


f_asb_calling <- function(db_data, loc_tf , target_lab, lab_cols, calling_method = 'edgeR', fdr, dp_threshold){
  #Call ASB events by multiple methods
  #edgeR, binomial
  
  #Filter the biased sites.
  known_bias = NULL
  if('ref_simulate_dp_broad1' %in% colnames(db_data)){
    original_count = nrow(db_data)
    db_data$known_bias = (db_data[['ref_simulate_dp_broad1']] + 1)/(db_data[['ref_simulate_dp_broad1']] + db_data[['alt_simulate_dp_broad1']] + 2)
    db_data = subset(db_data, known_bias <= 0.6 &  known_bias >= 0.4)
    cat('Reference bias filtering:', nrow(db_data), 'of ', original_count )
  } 
  
  
  #####Filter the sites with less than 10 reads
  if( length(grep('ref', lab_cols, value = T)) ==1)
  {
    #Rep 1 or 2, this is only one available
    #flog.info('Select the replicate %s', rep)
    db_data[["ref_tf_dp"]]=db_data[,grep('ref', lab_cols, value = T)]
    db_data[["alt_tf_dp"]]=db_data[[grep('alt', lab_cols, value = T)]]
  }else
  { #Sum the rep
    flog.info('Sum all the replicates')
    db_data[["ref_tf_dp"]]=rowSums(db_data[,grep('ref', lab_cols, value = T)])
    db_data[["alt_tf_dp"]]=rowSums(db_data[,grep('alt', lab_cols, value = T)])
  }
  loc_data_sub=subset(db_data, ref_tf_dp + alt_tf_dp >= dp_threshold)
  cat("Depth filter (>=", dp_threshold ,"bp) keep ", nrow(loc_data_sub), "of ", nrow(db_data),"points\n")
  db_data = loc_data_sub
  f_assert(nrow(db_data) > 0, message = 'After the depth filter')
  
  
  db_data$ASB=0
  db_data$p_value=0
  
  if(calling_method == 'binomial'){
    
    if (calling_method == 'binomial'){
      if(!is.null(db_data$known_bias)){
        db_data[,c("p_value","FDR")]=f_binomial_test_adjust_ref_bias(db_data[,c("ref_tf_dp","alt_tf_dp")], fdr, normalize_flag = F, db_data$known_bias ) 
      }else{
        db_data[,c("p_value","FDR")]=f_binomial_test(db_data[,c("ref_tf_dp","alt_tf_dp")],fdr, q_value = T)
      }
      
      
      #db_data = subset(db_data, FDR < fdr | FDR > 0.25)
      db_data$ASB = db_data$FDR <= fdr
      #db_data$FDR = NULL
      
      db_data$fold_change = ((db_data$ref_tf_dp + 0.01)/(db_data$alt_tf_dp + 0.01) <= 0.66 | (db_data$ref_tf_dp + 0.01)/(db_data$alt_tf_dp + 0.01) >= 1.5)
      db_data$ASB[db_data$ASB]="ASB"
      db_data$ASB[db_data$ASB==FALSE | db_data$fold_change == FALSE ]="nonASB"
      #db_data$ASB[db_data$ASB==FALSE]="nonASB"
      db_data$fold_change = NULL      
    }
        
  }else if(calling_method == 'edgeR'){
    #edgeR
    
    best_lab = target_lab
    f_select_highest_target_lab_and_rep(db_data = db_data,tf = loc_tf,lab_list = c('haib','sydh','uw'))
    ref_group = sort(grep(f_p('ref_%s_dp_%s[0-9]$', loc_tf, best_lab), colnames(db_data), value = TRUE))
    alt_group = sort(grep(f_p('alt_%s_dp_%s[0-9]$', loc_tf, best_lab), colnames(db_data), value = TRUE))
    
    pasillaCountTable = db_data[, c(ref_group, alt_group)]
    group = rep(1:2, each = length(ref_group))  
    
    edgeR_predictions = f_edgeR_call(pasillaCountTable, group)
    edgeR_predictions = edgeR_predictions[rownames(db_data),]

    db_data[,c("p_value","FDR")] = edgeR_predictions[,c('PValue', "FDR")]
    db_data[,c("ref_tf_dp","alt_tf_dp")] = edgeR_predictions[,c('ref_counts', 'alt_counts')]
    db_data$ASB = db_data$FDR < fdr
    db_data$ASB[db_data$ASB]="ASB"
    db_data$ASB[db_data$ASB==FALSE]="nonASB"
  }else if(calling_method == 'likelihood'){
    
    best_lab = target_lab
    f_select_highest_target_lab_and_rep(db_data = db_data,tf = loc_tf,lab_list = c('haib','sydh','uw'))
    ref_group = sort(grep(f_p('ref_%s_dp_%s[0-9]$', loc_tf, best_lab), colnames(db_data), value = TRUE))
    alt_group = sort(grep(f_p('alt_%s_dp_%s[0-9]$', loc_tf, best_lab), colnames(db_data), value = TRUE))
    
    db_data[,"p_value"] = f_likelihood_call_ASB(db_data[, ref_group[1]], db_data[, ref_group[2]], rowSums(db_data[, c(ref_group[1], alt_group[1])]), rowSums(db_data[, c(ref_group[2], alt_group[2])]) )
    db_data$FDR = p.adjust(db_data[,'p_value'], method = 'BH')

    db_data$ASB = db_data$FDR < fdr
    db_data$ASB[db_data$ASB]="ASB"
    db_data$ASB[db_data$ASB==FALSE]="nonASB"
    
  }else{
    
    flog.error('Unkown ASB calling method %s', calling_method)
    
  }
  
  return (db_data)
  
}


options(run.main=FALSE)
if (getOption('run.main', default=TRUE)) {
rm(list = ls())
runTestFile('./test/t_asb_calling.r',testFuncRegexp = '^t_.*')
}