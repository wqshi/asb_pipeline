f_create_file_name <- function(cell,tf,suffix,lab=NULL,rep=NULL)
{
  prefix=paste0(cell,"-",tf)
  
  if (!is.null(rep))
    prefix=paste0(prefix,"-",rep)
  
  
  if (!is.null(lab))
    prefix=paste0(lab,"-",prefix)
  
  return (paste0(prefix,".",suffix))
}


f_create_file_name_for_grep <- function(cell,tf,suffix,rep=NULL)
{
  prefix=paste0(cell,"-",tf)
  
  if (!is.null(rep))
    prefix=paste0(prefix,"-",rep)
  
  
  return (paste0(prefix,".",suffix,"$"))
}



f_read_bed_table <- function(input_file, quiet=FALSE,...)
{
  data=f_read_table3(input_file, quiet)
  row.names(data)=paste0(data$chr,"-",data$start)
  
  return (data)
}

f_read_bed_table_remove_dup <- function(input_file, quiet = TRUE, ...)
{
  data=f_read_table3(input_file, quiet )
  
  if(any(lapply(data$chip_ref,FUN=str_length)!=1))
  {
    cat("filter bad bed format lines\n")
    if (quiet != TRUE)
    {
      print(data[lapply(data$chip_ref,FUN=str_length)!=1,])  
    }
    
  }
  
  
  data_subset=data[lapply(data$chip_ref,FUN=str_length)==1,]
  
  if(any(duplicated(paste0(data_subset$chr,"-",data_subset$start))))
  {
    print(data_subset[duplicated(paste0(data_subset$chr,"-",data_subset$start)),])
    data_subset=data_subset[!duplicated(paste0(data_subset$chr,"-",data_subset$start)),]
    
  }
  
  
  row.names(data_subset)=paste0(data_subset$chr,"-",data_subset$start)
  
  return (data_subset)
}


f_read_bed_table_remove_dup_by_keys <- function(input_file, key_cols ,quiet = TRUE, ...)
{
  data=f_read_table3(input_file, quiet )
  
  data$start = as.integer(data$start)
  
  data$tmp.rownames=data[[key_cols[1]]]
  for (col_name in key_cols[-1])
  {
    data$tmp.rownames=paste(data$tmp.rownames, data[[col_name]], sep="-")
  }
  
  
  #filter out the ref allele with more than 1 letter
  if ("chip_ref" %in% colnames(data))
  {
    if(any(lapply(data$chip_ref,FUN=str_length)!=1))
    {
      cat("filter bad bed format lines\n")
      print(data[lapply(data$chip_ref,FUN=str_length)!=1,])
    }
    data_subset=data[lapply(data$chip_ref,FUN=str_length)==1,]
  }else
  {
    data_subset = data
  }
  
  
  
  if(any(duplicated(data_subset$tmp.rownames)))
  {
    if(quiet == T)
    {
      cat("Duplicated rows: ", nrow(data_subset[f_duplicated(data_subset$tmp.rownames),]),'\n')
    }else
    {
      print(data_subset[f_duplicated(data_subset$tmp.rownames),])
    }
    
    data_subset=data_subset[!duplicated(data_subset$tmp.rownames),]
    
  }
  
  
  row.names(data_subset)=data_subset$tmp.rownames
  data_subset$tmp.rownames = NULL
  
  return (data_subset)
}

