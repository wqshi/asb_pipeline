f_get_epi_from_paird_database_rep <- function(cell, guest_cell, tf, output_dir, dp_threshold, fdr, cell_filter, het_filter, add_het=FALSE, 
                                              cord_name=TRUE, strong_allele_exchange=TRUE, raw_data=FALSE, db_data=NULL, quiet=TRUE,
                                              labs = c('uw'), target_lab = 'uw', rep ='1', file_suffix = 'database$', diffpeak_mehtod='pepr', host_alt_exchange = F,
                                              diffpeak_method = 'PePr', return_opt = '', pwm_fdr = 0.001)
{ #cord_name: name each row as chr-start, default as chr-start-cell
  #strong_allele_exchange: For the ASB events in which strong allele is ALT, exchange the value of ALT with Ref. 
  #As a result, ref is strong alllele consistantly.
  
  if(is.null(db_data))
  {
    db_file=f_create_paired_database_name(cell, guest_cell, tf, file_suffix)
    
    db_data_raw=f_read_bed_table_remove_dup_by_keys(f_find_file(output_dir,db_file),key_cols = c("chr","start","cell"), quiet)
    
    dim(db_data_raw)
    
    #filter by diploid data
    if ('diploid' %in% colnames(db_data_raw) & cell == cell_filter)
    {
      db_data_raw = subset(db_data_raw, diploid != '.')
    }
    
    #low_quality_rows = rownames(subset(subset(db_data_raw, wgs_alt_dp <=5 | wgs_ref_dp <=5), het_type == 'het'))
    #db_data_raw = db_data_raw[!rownames(db_data_raw) %in% low_quality_rows,]
    
    #change the data format to numeric
    db_data=subset(db_data_raw, cell == cell_filter & het_type %in% het_filter)
    dim(db_data)
    
    if(nrow(db_data) == 0)
    {
      cat('Empty dataframe after filtering cell type\n')
      return (0)
    }
    
    colnames(db_data)
    
    numeric_cols=grep("strand",
                      grep("score|dp|cons|peak_.*dis|pwm|hct|value|count|depth|start|log_|minp_", colnames(db_data), value = TRUE), 
                      value=TRUE, invert=TRUE
    )
    table(db_data$het_type)
    
    for (col_name in numeric_cols)
    {
      if (any(db_data[[col_name]]=="."))
      {
        db_data[db_data[[col_name]]==".",col_name]="0"
      }
      
      db_data[[col_name]] = as.numeric(db_data[[col_name]])
    }
    
    target_tf = tf
    
    
    cofactor_list = str_replace_all(grep( 'peak|gene', grep('[a-z].*_dis', colnames(db_data),value = TRUE), value = TRUE, invert = TRUE), pattern = '_dis', replacement = '')  
    for (cofactor_tf in cofactor_list){
      db_data[[paste0(tolower(cofactor_tf),'_dis')]] = as.numeric(as.numeric(db_data[[paste0(cofactor_tf,'_dis')]]) > 0)
      db_data[[paste0(toupper(cofactor_tf),'_dis')]] = NULL
    }
    
    
    #snp_cols=grep(".*strand$|.+start$|p_value|tffm.*|neighbour|distance|type|cell",x = colnames(db_data),value = TRUE, invert = TRUE)
    
    if (cord_name==TRUE)
    {
      db_data=f_rename_by_genome_cordinate(db_data)
    }
    
    
    
    if (raw_data == TRUE)
    {
      return (db_data)
    }
  }
  
  colnames(db_data)
  
  #Change the pwm_ref_start
  pwm_information_content = read.table('./data/server/snv/pwm_information_content.list',header = T)
  rownames(pwm_information_content) = tolower(pwm_information_content$tf)
  pwm_len = pwm_information_content[tf, 'length']
  if(!is.na(pwm_len)){
    forward_strand = db_data$pwm_ref_strand == "+"
    db_data$pwm_position = db_data$pwm_ref_start
    db_data$pwm_position[forward_strand]=pwm_len - db_data$pwm_ref_start[forward_strand] -1 #1 based
  }
  
  
  
  if (is.na(target_lab)){
    
    result_list=f_select_highest_target_lab_and_rep(db_data, tf, labs)
    
    if (rep != '.'){
      target_lab = result_list$best[2] #best replicates
      rep = result_list$best[3]      
    }else{
      target_lab = result_list$best[4] #best lab
    }
    
    print(result_list$stats)
  }
  
  
  lab_pattern= f_p( '.*%s.*_(%s)[1-9]*', tf ,paste(labs,collapse = '|')) 
  rep_pattern = f_p('(%s).*(%s%s)', tf, target_lab, rep )
  
  f_check_data_features(db_data, necessary = c(paste0('lab_',target_lab), lab_pattern, rep_pattern, 'peak_.*_dis','start', 'peak.*_pwm_ref_score', 'pwm_ref_score'))
  
  #Filter only one lab data according to the lab name and 
  #necessary_list = c(paste0('lab_',target_lab), )
  #f_check_data_features(db_data, necessary_list)
  
  db_data$wgs_ref_dp=NULL
  db_data$wgs_alt_dp=NULL
  
  colnames(db_data)
  db_data = db_data[db_data[[paste0('lab_',target_lab)]] != '.',] 
  f_assert(nrow(db_data) > 0, message = 'Selecting after targeted lab')
  if(nrow(db_data)==0){
    return (NULL)
  }
  
  
  
  lab_cols=grep(lab_pattern, colnames(db_data), value = TRUE)
  
  
  keep_cols = grep(rep_pattern, lab_cols, value = T)
  
  #Change the TF's read depth into ref_tf_dp and alt_tf_dp
  if ( length(grep('ref', keep_cols, value = T)) ==1)
  {
    #Rep 1 or 2
    flog.info('Select the replicate %s', rep)
    db_data[["ref_tf_dp"]]=db_data[,grep('ref', keep_cols, value = T)]
    db_data[["alt_tf_dp"]]=db_data[[grep('alt', keep_cols, value = T)]]
  }else
  { #Sum the rep
    flog.info('Sum all the replicates')
    db_data[["ref_tf_dp"]]=rowSums(db_data[,grep('ref', keep_cols, value = T)])
    db_data[["alt_tf_dp"]]=rowSums(db_data[,grep('alt', keep_cols, value = T)])
  }
  
  f_assert(all(db_data[["ref_tf_dp"]] <= 100), 'Total read depth is above 100bp, so it is not rmdup data')
  
  
  #Change the epigenetic features
  #feature_list = c("dnase", 'h2az', 'h3k27ac', 'h3k27me3', 'h3k36me3', 'h3k4me1', "h3k4me2","h3k4me3", 'h3k79me2', "h3k9ac","h3k9me3","h4k20me1","input")
  dp_cols = grep('.*_dp_.*', colnames(db_data), value = TRUE)
  feature_cols = grep(tf, dp_cols, value = TRUE, invert = TRUE)
  
  if(length(feature_cols) > 0){
    feature_list = unique(matrix(unlist(str_split(string = feature_cols,pattern = '_')), byrow = TRUE, ncol = 4)[,2]) 
  }else{
    feature_list = list()
  }
  
  
  
  for (feature in feature_list){
    result_list=f_select_highest_target_lab_and_rep(db_data, feature, c('broad', labs))
    
    if (is.null(result_list)) next
    
    db_data[[f_p("ref_%s_dp", feature)]]=db_data[,f_p("ref_%s_dp_%s%s", feature, result_list$best[2], result_list$best[3])]
    db_data[[f_p("alt_%s_dp", feature)]]=db_data[,f_p("alt_%s_dp_%s%s", feature, result_list$best[2], result_list$best[3])]
    
  }
  
  
  
  
  
  
  
  #Peak_lab related features change back to peak_, choosing by lab name.
  peak_cols = grep(f_p('peak.*%s', target_lab), colnames(db_data), value = TRUE)
  renamed_peak_cols = str_replace_all(peak_cols, pattern = paste0('_', target_lab), replacement = '')
  db_data[,renamed_peak_cols] = db_data[,peak_cols] 
  
  
  #HCT colnames
  hct_cols = grep(f_p('hct.*%s', target_lab), colnames(db_data), value = TRUE)
  renamed_hct_cols = str_replace_all(hct_cols, pattern = paste0('_', target_lab), replacement = '')
  db_data[,renamed_hct_cols] = db_data[,hct_cols]
  
  #Remove the renamed colnunms: TF ChIP and TF's peak colunms.
  TF_pattern = f_p('(%s)_dp_(%s)', tf, paste(labs,collapse = '|'))
  tf_cols_remove = grep(TF_pattern, colnames(db_data), value = TRUE)
  
  peak_pattern = f_p('peak_(%s)', paste(labs,collapse = '|'))
  peak_cols_remove = grep(peak_pattern, colnames(db_data), value = TRUE)
  
  hct_pattern = f_p('hct.*(%s)', paste(labs,collapse = '|'))
  hct_cols_remove = grep(hct_pattern, colnames(db_data), value = TRUE)
  
  if(cell == 'helas3'){
    db_data$filter = as.numeric(db_data$filter)
    flog.info(f_p('Hela-S3 filter: %s ', nrow(subset(db_data, filter < 30))))
    db_data = subset(db_data, filter > 30)
  }
  
  
  db_data = db_data[, !(names(db_data) %in% c(feature_cols, tf_cols_remove, peak_cols_remove, hct_cols_remove, 'filter'))]
  colnames(db_data)
  
  
  
  
  #Yes, this is necessary
  db_data$peak_dis=db_data$peak_dis - db_data$start
  #str(db_data)
  
  db_data$ref_simulate_dp
  #Depth filter
  db_data$best_match="nonBest"
  db_data$best_match[abs(db_data$peak_pwm_ref_score - db_data$pwm_ref_score) < 0.0000001] = "Best"
  db_data$best_match = as.factor(db_data$best_match)
  #snp_cols=grep(".*strand$|.+start$|p_value|tffm.*|neighbour|distance|type|cell",x = colnames(db_data),value = TRUE, invert = TRUE)
  
  
  
  loc_data_sub=subset(db_data, ref_tf_dp + alt_tf_dp >= dp_threshold)
  cat("Depth filter (>=", dp_threshold ,"bp) keep ", nrow(loc_data_sub), "of ", nrow(db_data),"points\n")
  db_data = loc_data_sub
  f_assert(nrow(db_data) > 0, message = 'After the depth filter')
  
  db_data$ASB=0
  db_data$p_value=0
  
  if('ref_simulate_dp' %in% colnames(db_data)){
    
    known_bias = (db_data[['ref_simulate_dp']] + 1)/(db_data[['ref_simulate_dp']] + db_data[['alt_simulate_dp']] + 2)
    
    
    db_data = subset(db_data, known_bias < 0.6 &  known_bias > 0.4)
    cat('Reference bias filtering:', nrow(db_data), 'of ', length(known_bias) )
    
    
    new_bias = (db_data[['ref_simulate_dp']] + 1)/(db_data[['ref_simulate_dp']] + db_data[['alt_simulate_dp']] + 2)
    
    db_data[,c("p_value","FDR")]=f_binomial_test_adjust_ref_bias(db_data[,c("ref_tf_dp","alt_tf_dp")], fdr, normalize_flag = F, new_bias )
  }else{
    db_data[,c("p_value","FDR")]=f_binomial_test(db_data[,c("ref_tf_dp","alt_tf_dp")],fdr, q_value = T)
    
  }
  
  #db_data = subset(db_data, FDR < fdr | FDR > 0.25)
  db_data$ASB = db_data$FDR < fdr
  db_data$FDR = NULL
  
  db_data$fold_change = ((db_data$ref_tf_dp + 0.01)/(db_data$alt_tf_dp + 0.01) < 0.5 | (db_data$ref_tf_dp + 0.01)/(db_data$alt_tf_dp + 0.01) > 2)
  db_data$ASB[db_data$ASB]="ASB"
  db_data$ASB[db_data$ASB==FALSE | db_data$fold_change == FALSE ]="nonASB"
  db_data$fold_change = NULL
  
  #table(db_data$ASB)
  #View(db_data[,c("ref_tf_dp","alt_tf_dp","ASB","p_value")])
  #db_data$start = as.numeric(matrix(unlist(str_split(rownames(db_data),pattern = "-")), ncol=2, byrow = TRUE)[,2])
  #db_data$chr = (matrix(unlist(str_split(rownames(db_data),pattern = "-")), ncol=2, byrow = TRUE)[,1])
  db_data=f_filter_multiple_variants_in_peak(db_data)
  
  #db_data$start=NULL
  #db_data$chr=NULL
  
  train_data = db_data
  #colnames(train_data) = tolower(colnames(train_data))
  dim(train_data)
  
  
  #The p-value is only available in database cell
  if (cell == 'gm12878')
  {
    #diffpeak_method = 'PePr'
    diffpeak_cols=grep(diffpeak_method, colnames(train_data), value = T, ignore.case = T)
    replaced_cols=str_replace_all(diffpeak_cols,paste0('_', diffpeak_method),'')
    
    train_data[,replaced_cols] = train_data[,diffpeak_cols]
    #train_data$diff_pvalue=f_log_pvalue(train_data$diff_pvalue)  
  }
  
  
  homer_cols = grep('ref.*[.]homer',colnames(train_data), value = T)
  
  ref_homer_col = homer_cols[1]
  for(ref_homer_col in homer_cols){
    alt_homer_col=str_replace(ref_homer_col,pattern = 'ref_', replacement = 'alt_')
    binding_homer_col=str_replace(ref_homer_col,pattern = 'ref_', replacement = 'binding_')
    train_data[[binding_homer_col]]=(train_data[[alt_homer_col]] < pwm_fdr | train_data[[ref_homer_col]] < pwm_fdr ) 
  }
  
  #install.packages('data.table')
  #library(data.table)
  
  
  if(f_p('ref_pwm_%s', tf) %in% colnames(train_data)){
    train_data[[f_p('binding_pwm_%s', tf)]]=(train_data[[f_p('ref_pwm_%s', tf)]] < pwm_fdr | train_data[[f_p('alt_pwm_%s', tf)]] < pwm_fdr )
  }
  
  #adjust the alt ASB to the strong allele (ref)
  if (strong_allele_exchange == TRUE)
  {
    #In case of host is homo alt sites, and use this one as strong binding signal.
    if (host_alt_exchange == TRUE)
    { 
      #used for prediction of variants between individuals.
      alt_asb_row_names=rownames(train_data)  
    }else
    {
      alt_asb_row_names=rownames(subset(train_data, alt_tf_dp > ref_tf_dp))
    }
    
    
    
    ref_col=sort(grep("ref",colnames(train_data),value = TRUE))
    alt_col=grep("alt",colnames(train_data),value = TRUE)
    
    
    replaced_ref_names = (str_replace_all(ref_col,pattern = "ref","alt"))
    replaced_alt_names = (str_replace_all(alt_col,pattern = "alt","ref"))
    
    f_assert(length(setdiff(replaced_ref_names, alt_col)) == 0,"Colnames mismatch")
    
    tmp_df = train_data[alt_asb_row_names,ref_col]
    train_data[alt_asb_row_names,replaced_alt_names] = train_data[alt_asb_row_names,alt_col]
    train_data[alt_asb_row_names,replaced_ref_names] = tmp_df
    
    train_data$alt_flip = FALSE
    train_data[alt_asb_row_names, 'alt_flip'] = TRUE
    
    log_ratio_cols = grep("log_", colnames(train_data),value = TRUE)
    if (length(log_ratio_cols) > 0){
      train_data[alt_asb_row_names,log_ratio_cols] = -1 * train_data[alt_asb_row_names,log_ratio_cols]
    }
    
    
  }
  
  
  
  ref_pwm_cols=grep("ref_pwm_",colnames(train_data), value = TRUE)
  for(ref_pwm_col in ref_pwm_cols)
  {
    delta=0.00000000000001
    new_col=str_replace(string = ref_pwm_col, pattern = "ref_pwm", "pvalue_ratio")
    alt_pwm_col=str_replace(string = ref_pwm_col, pattern = "ref_pwm", "alt_pwm")
    train_data[[new_col]] = -log10((train_data[[ref_pwm_col]] + delta)/(train_data[[alt_pwm_col]] + delta))
    train_data[[ref_pwm_col]] = NULL
    train_data[[alt_pwm_col]] = NULL
  }
  
  
  if (return_opt == 'lab')
  {
    return_data = list(data = train_data, lab = target_lab)
    return (return_data)
  }else{
    return (train_data)
  }
  
  
}



t_select_highest_target_lab_and_rep <- function(){
  
  
  fdr=0.05
  test_base_dir="E:/Projects/R/data/test/t_get_epi_from_paird_database/"
  
  test_cell="gm12878"
  guest_cell = 'gm12872'
  dp_threshold=1
  fdr=0.05
  tf = 'ctcf'
  output_dir = test_base_dir
  quiet = TRUE
  cell_filter = 'gm12878'
  het_filter = 'alt'
  add_het=FALSE
  cord_name=TRUE
  strong_allele_exchange=TRUE
  raw_data=FALSE
  db_data=NULL
  quiet=TRUE
  labs = c('uw','sydh')
  target_lab = 'uw'
  rep='1'
  file_suffix ='test.database'
  
  db_file=f_create_paired_database_name(test_cell, guest_cell, tf, file_suffix)
  db_data=f_read_bed_table_remove_dup_by_keys(f_find_file(test_base_dir,db_file),key_cols = c("chr","start","cell"), quiet)
  
  
  numeric_cols=grep("strand",
                    grep("score|dp|cons|peak_.*dis|pwm|hct|value|count|depth|start|log_|minp_", colnames(db_data), value = TRUE), 
                    value=TRUE, invert=TRUE
  )
  table(db_data$het_type)
  
  for (col_name in numeric_cols)
  {
    if (any(db_data[[col_name]]=="."))
    {
      db_data[db_data[[col_name]]==".",col_name]="0"
    }
    
    db_data[[col_name]] = as.numeric(db_data[[col_name]])
  }
  
  
  #mtrace(f_select_highest_target_lab_and_rep)
  return_list=f_select_highest_target_lab_and_rep(db_data, tf ,lab_list)
  
  return_list$best[4]
}


#undebug(f_select_highest_target_lab_and_rep)

f_select_highest_target_lab_and_rep <- function(db_data, tf ,lab_list){
  
  lab_pattern= f_p( '.*%s.*_(%s)[1-9]*', tf ,paste(lab_list,collapse = '|')) 
  lab_cols=grep(lab_pattern, colnames(db_data), value = TRUE)
  
  if (length(lab_cols) == 0){
    return (NULL)
  }
  
  best_count  = -1
  best_lab = ''
  best_rep = ''
  
  lab_and_rep_stat = data.frame(lab = lab_list, rep1 = 0, rep2 = 0, rep3 =0)
  rownames(lab_and_rep_stat) = lab_list
  
  for (lab_loc in lab_list){
    
    loc_lab_cols = grep(lab_loc, lab_cols,value = TRUE)
    
    #loop over the replicates
    if (length(loc_lab_cols) == 0) next
    
    for (rep in 1:3){
      
      rep_cols = grep(f_p('.*_%s%s', lab_loc, rep), loc_lab_cols ,value = TRUE)
      
      if(length(rep_cols) != 0){
        loc_count= sum(colSums( db_data[rep_cols],na.rm = T))
        
        lab_and_rep_stat[lab_loc, f_p('rep%s',rep)] = loc_count
        
        if (loc_count > best_count){
          best_count = loc_count 
          best_lab = lab_loc
          best_rep = rep
        }
        
        
        
      }
      
    }
    
    
    
  }
  
  best_lab_reps = names(which.max(rowSums(lab_and_rep_stat[,2:4]))) 
  
  return_list=list(best=c(best_count, best_lab, best_rep, best_lab_reps), stats = lab_and_rep_stat)
  
  return (return_list)
  
}  




options(run.main=FALSE)
if (getOption('run.main', default=TRUE)) {
runTestFile('./test/t_get_epi_from_paird_database_rep.r',testFuncRegexp = '^t_.*')
}