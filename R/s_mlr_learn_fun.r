#source('./s_ggplot2_theme.R')

f_rename_by_genome_cordinate <- function(data, key_cols=c("chr","start"))
{
  data$tmp.rownames=data[[key_cols[1]]]
  for (col_name in key_cols[-1])
  {
    data$tmp.rownames=paste(data$tmp.rownames, data[[col_name]], sep="-")
  }
  rownames(data)=data$tmp.rownames
  data$tmp.rownames=NULL
  return (data)
}

f_rename_chromotin_properties <- function(feature_match){
  search_pattern = c('DNASE','ME','AC','INPUT')
  search_replacement = c('DHS','me','ac','Control')
  input_vector = toupper(feature_match)
  for(i in 1:length(search_replacement)){
    input_vector = str_replace_all(input_vector, pattern = search_pattern[i], replacement = search_replacement[i] ) 
  }
  return (input_vector)
}


f_binomial_test <- function(data, fdr=0.05, normalize_flag = FALSE, q_value = T)
{
  #loginfo("=============: fdr %s", fdr)
  
  result=as.data.frame(matrix(0,nrow=nrow(data),ncol=2))
  colnames(result)=c("p_value","p_adjust")
  flog.info('P-value normalization: %s' , normalize_flag)
  for (i in 1:nrow(data))
  {
    
    a=unlist(data[i,c("alt_tf_dp","ref_tf_dp")])
    if (normalize_flag == TRUE){
      result[i,"p_value"]=binom.test(x=a)$p.value * sum(a)^0.5 * 0.31
    }else{
      result[i,"p_value"]=binom.test(x=a)$p.value
    }
    
    
  }
  if(q_value == T){
    result$p_adjust= p.adjust(result$p_value, "BH")
  }else{
    result$p_adjust= p.adjust(result$p_value, "BH") < fdr  
  }
  
  
  return (result)
  
}

f_binomial_test_adjust_ref_bias <- function(data, fdr=0.05, normalize_flag = FALSE, known_bias)
{
  #loginfo("=============: fdr %s", fdr)
  
  result=as.data.frame(matrix(0,nrow=nrow(data),ncol=2))
  colnames(result)=c("p_value","p_adjust")
  flog.info('P-value normalization: %s' , normalize_flag)
  for (i in 1:nrow(data))
  {
    
    a=unlist(data[i,c("alt_tf_dp","ref_tf_dp")])
    if (normalize_flag == TRUE){
      result[i,"p_value"]=binom.test(x=a, p = known_bias[i])$p.value * sum(a)^0.5 * 0.31
    }else{
      result[i,"p_value"]=binom.test(x=a, p = known_bias[i] )$p.value
    }
  }
  result$p_adjust= p.adjust(result$p_value, "BH")
  return (result)
  
}

f_binomial_test2 <- function(data, fdr=0.05)
{
  #loginfo("=============: fdr %f", fdr)
  
  result=as.data.frame(matrix(0,nrow=nrow(data),ncol=2))
  colnames(result)=c("p_value","p_adjust")
  for (i in 1:nrow(data))
  {
    a=unlist(data[i,])
    result[i,"p_value"]=binom.test(x=a)$p.value
    
  }
  result$p_adjust= p.adjust(result$p_value, "BH") < fdr
  return (result)
  
}


f_filter_multiple_variants_in_peak <- function(all_data, target_lab = NA, dup_peak = FALSE)
{

    
    bed_col = f_p('peak_%s_bed_start', target_lab)
    if (bed_col %in% colnames(all_data)){
      
      lab_start_loc = paste0(all_data[['chr']], '-' , all_data[[bed_col]])
      duplicate_locs=f_duplicated(lab_start_loc)
      lab_start_loc[duplicate_locs]
      snv_loc = paste0(all_data[['chr']], '-' , all_data[['start']])
      snv_loc[duplicate_locs]
      if (dup_peak == F){
        return (all_data[!duplicate_locs,])
      }else{
        return (list(data = all_data[!duplicate_locs,], filter_beds =length(unique(lab_start_loc[duplicate_locs]))))
      }
      
      
    }else{
      
      flog.error('Missing bed_col')
      #return (NA)
      colnames(all_data)
      all_data$distance=0
      all_data$distance= c(all_data$start[-1],0) - all_data$start
      all_data$neighbour=c(rownames(all_data[-1,]),"end")
      head(all_data)
      short_distance_het = subset(all_data, subset = distance < 100 & distance > -1)
      
      neighbour_data = rbind(short_distance_het, all_data[short_distance_het$neighbour,])
      dim(all_data)
      dim(short_distance_het)
      
      
      #neighbour_data_sorted=neighbour_data[order(rownames(neighbour_data)),c("chr","start","ref_tf_dp","alt_tf_dp","ASB","best_match","neighbour")]
      
      filtered_data=all_data[!(rownames(all_data) %in% rownames(neighbour_data)),]
      
      range(abs(filtered_data$distance))
      
      #dim(neighbour_data_sorted)
      cat(dim(neighbour_data),"het sites have been filtered out for multiple var in one peak")
      filtered_data$distance = NULL
      filtered_data$neighbour = NULL
      
      
      if (dup_peak == F){
        return (filtered_data)
      }else{
        return (list(data = filtered_data, filter_beds = nrow(short_distance_het) ))
      }
      
      
    }
    

}

#debug(f_filter_multiple_variants_in_peak)

f_log_pvalue <- function(input_vector, default_value = 0)
{
  input_vector[input_vector == default_value] = 1
  return (-log10(input_vector))
}
