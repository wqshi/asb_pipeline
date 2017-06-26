######function library############
f_assert  <- function(condition,message){
  #Description: self defined warning function for debug.
  #input:condition is the cretiria of the warning, if False, show the message and the location.
  #output:
  #By wenqiang
  #Date: 2013-02-05
  as.character(match.call()[[1]]);
  if (condition == FALSE)
  {
    print(sprintf('The assertion failed:%s, at %s',message, as.character(match.call()[[1]])));
    a=1/0
  }
}

f_ggsave<-function(filename, loc_height = 3.3 ,single = TRUE, ...){
  
  if(single == TRUE){
    loc_width = 3.3
  }else{
    loc_width = 6.6
  }
  ggsave(filename = filename, width = loc_width, height = loc_height, dpi = 310, ...)
}


f_check_na <- function(x)
{
  return (which(sapply(x, function(x)any(is.na(x)))))
}


f_as_numeric <- function(input_data, NA_value)
{
  result = as.numeric(input_data)
  result[is.na(result)] = NA_value
  return (result)
}


f_sortBycol <- function(data,index)
{
  sort_data=data[sort.list(data[,index]),]
  return (sort_data)
}

f_sort_by_col <- function(data,index, decr_flag=FALSE)
{
  sort_data=data[sort.list(data[,index],decreasing = decr_flag),]
  return (sort_data)
}

f_percentage <- function (x)
{
  percentage = x/sum(x)
  
  
  return (x/sum(x))
}


f_read_table <- function(data_prob_file,...)
{
  data_dir=file.path(getwd(),'data/hw2/')
  #data_prob_file= 'array.txt';
  #data_design_file= 'gse7191_design.txt';
  prob_data <- read.table(paste(data_dir,data_prob_file,sep=""),header=TRUE,row.names=1,...);
  
  #Data Check
  f_assert(all(is.na.data.frame(prob_data)==FALSE),"There is NaN Value in the Probdata")
  #printf()
  #range(prob_data)
  colnames(prob_data)=sub(".","_",colnames(prob_data),fixed=TRUE)
  prob_data
}



f_read_table2 <- function(file_path,quiet=FALSE,...)
{
  prob_data <- read.delim2(file_path,...)
  
  #Data Check
  f_assert(all(is.na.data.frame(prob_data)==FALSE),"There is NaN Value in the Probdata")
  
  #To remove the # for the colnames
  data_colnames=colnames(prob_data)
  data_colnames[1]=str_replace(data_colnames[1],"X.","")
  colnames(prob_data) = data_colnames
  
  
  #printf()
  #range(prob_data)
  #colnames(prob_data)=sub(".","_",colnames(prob_data),fixed=TRUE)
  if (quiet == FALSE)
  {
    cat("row number", nrow(prob_data), "\n")
    print(head(prob_data))
  }
  prob_data
}

f_read_table2_1 <- function(file_path,...)
{
  prob_data <- read.table(file_path,...)
  
  #Data Check
  f_assert(all(is.na.data.frame(prob_data)==FALSE),"There is NaN Value in the Probdata")
  
  #To remove the # for the colnames
  data_colnames=colnames(prob_data)
  data_colnames[1]=str_replace(data_colnames[1],"X.","")
  colnames(prob_data) = data_colnames
  
  cat("row number", nrow(prob_data), "\n")
  #printf()
  #range(prob_data)
  #colnames(prob_data)=sub(".","_",colnames(prob_data),fixed=TRUE)
  print(head(prob_data))
  prob_data
}



f_read_table3 <- function(file_path,quiet=FALSE,...)
{
  return (f_read_table2(file_path,quiet,header=TRUE, comment.char="",...))
}

f_read_table3_1 <- function(file_path,...)
{
  return (f_read_table2_1(file_path,header=TRUE, comment.char="",...))
}



f_check_show<-function(data)
{
  #NaN check
  f_assert(all(is.na.data.frame(data)==FALSE),"There is NaN Value in the Probdata")
  
  if (is.vector(data))
  {
    cat("\nLenght:",length(data))
    cat("\nRange:",range(data))
  }
  else
  {
    cat("\ndim:",dim(data))
  }
  
  
  cat("\nhead:\n")
  print(summary(data))
  print(head(data))
  cat("\nTail:\n")
  print(tail(data))
}

f_subset_by_col <- function(data,by_col,sel_col,threshold=1e-5, smaller = TRUE, show=TRUE)
{
  if(sel_col=='.')
  {
    sel_col=names(data)
  }
  
  if (smaller == TRUE)
  {
    subse_col=data[data[,by_col]<threshold,sel_col]
  }
  else
  {
    subse_col=data[data[,by_col]>threshold,sel_col]
  }
  if (show == TRUE)
  {
    cat("Length:",length(subse_col)) 
  }
  
  return ((subse_col))
  
}

f_subset<-function(data,by_col,value)
{
  subse_col=data[data[,by_col]==value,]
  return (subse_col)
  
}


f_venn_plot <- function(data_a,data_b,name_a,name_b,...)
{
  #f_assert(is.vector(data_a)&&is.vector(data_b),"Input to the f_venn_plot Wrong!!!!")
  #print(sprintf("list_data <- list(%s = data_a, %s = data_b)", name_a, name_b))
  eval(parse(text=sprintf("list_data <- list(%s = data_a, %s = data_b)", name_a, name_b)))
  par(mar = rep(4, 4))
  plot.new()
  venn.plot <- venn.diagram(list_data, height = 440, width = 400, ,filename = NULL, fill = c("red", "blue"), ...)
  grid.draw(venn.plot)
  venn.plot
  return (venn.plot)
}

f_venn_plot3 <- function(data_a,data_b,data_c,name_a,name_b,name_c,...)
{
  #f_assert(is.vector(data_a)&&is.vector(data_b),"Input to the f_venn_plot Wrong!!!!")
  eval(parse(text=sprintf("list_data <- list(%s = data_a, %s = data_b, %s=data_c)", name_a, name_b, name_c)))
  plot.new()
  venn.plot <- venn.diagram(list_data, filename = NULL, scaled=TRUE, fill = c("red", "blue","yellow"), ...)
  grid.draw(venn.plot)
  venn.plot
  
}




f_prepare_byCol <- function (prDat,prDes)
{ #Create gExp data with addinal factor colums from design Matrix
  
  #luckyGenes=as.vector(luckyGenes);
  f_assert(all(colnames(prDat)==row.names(prDes)),'Input Error in "f_prepare_data.R"');
  group_factors=colnames(prDes);
  subset_data_frame=data.frame(gExp=as.vector((as.matrix(prDat))),names=rep(rownames(prDat),times=ncol(prDat)));
  for (factor_name in group_factors)
  {
    eval(parse(text=sprintf("subset_data_frame=data.frame(subset_data_frame,%s=factor(rep(prDes[,factor_name],each=nrow(prDat))))",factor_name))); 
  }
  #rownames(subset_data_frame)<-as.vector(rep(as.vector(rownames(prDes)),each=length(luckyGenes)))
  subset_data_frame;
}

f_check_null <- function(str)
{
  if(is.null(str) || length(str)==0)
  {
    cat("Null check Failed",'\n')
    return (TRUE)
  }else
  {
    return (FALSE)
  }
  
}


f_check_na <- function(data)
{
  return (any(as.vector(is.na(data))))
}
##test
# prDat=data.frame(row.names=c("a","b","c"),DNA=1:3,RNA=4:6)
# prDes=data.frame(type=c("A","B"),row.names=c("DNA","RNA"))
# f_prepare_byCol(prDat,prDes)

f_plain_write <- function(data, file)
{
  cat("write file in", file, "\n")
  write.table(data,file,quote=FALSE, row.names=FALSE, col.names=FALSE,sep='\t')
  
}


f_header_write <- function(data, file)
{
  cat("write file in", file, "\n")
  write.table(data,file,quote=FALSE, row.names=FALSE, col.names=TRUE,sep='\t')
  
}


f_length <- function(data)
{
  data_len=-1
  if (class(data) == "data.frame")
  {
    data_len = nrow(data)
  }else {
    data_len = length(data)
  }
  
  return (data_len)
  
}

source('f_create_fun.r')
source('f_create_test.r')

f_logic2integer <- function(logic_vector)
{
  data=(1:length(logic_vector))[logic_vector]
  
  return (data)
  
  
}


f_merge_duplicate  <- function(.data, index){
  #Description: For a data frame, merge the duplicate rows according to one factor colum.
  #input: data should be data.frame and index should be charators
  #output:
  #By Wenqiang,  2013-06-26 
  
  data=f_sort_by_col(.data, index)
  
  duplicate_index=f_logic2integer(duplicated(data[[index]]))
  
  new_data=data[-duplicate_index,]
  rownames(new_data) = new_data[,index]
  selected_add_cols = setdiff(names(data),index)
  
  for (i in duplicate_index)
  {
    row_name = data[i, index]
    new_data[row_name,selected_add_cols] = data[i,selected_add_cols] + new_data[row_name,selected_add_cols]
  }
  
  
  return (new_data[,selected_add_cols])
}


f_assert_noDup <- function(data)
{
  f_assert(sum(duplicated(data)) == 0, "Duplicated Check Failed")
}

f_run_shell <-function(file_name, shell_cmd){
  #Description: this function gives send the shell cmd to the linux server via winscp
  #input: file_name, linux script name
  #output:
  #By Wenqiang,  2013-07-15 
  
  #Write the txt file for the winscp
  winscp_dir='./winscp/'
  txt_name=paste(file_name,".txt",sep="");
  cat(
"option batch on
option confirm off
open loire
call  chmod 777 ",shell_cmd,"
call ",shell_cmd,"
close
exit",file= paste0(winscp_dir,txt_name));
  
  #write the bat file
  bat_name=paste(file_name,".bat",sep="");
  cat("WinSCP.exe /log=./log /console /script=E:/Projects/R/winscp/",txt_name,sep="",file= paste0(winscp_dir,bat_name))
  
  system(paste0(winscp_dir,bat_name))
  
  
  
}

f_run_shell_para <-function(file_name, shell_cmd, para){
  #Description: this function gives send the shell cmd to the linux server via winscp
  #input: file_name, linux script name
  #output:
  #By Wenqiang,  2013-07-15 
  
  #Write the txt file for the winscp
  winscp_dir='./winscp/'
  txt_name=paste(file_name,".txt",sep="");
  cat(
    "option batch on
    option confirm off
    open loire
    call  chmod 777 ",shell_cmd,"
    call ",shell_cmd," ", para,"
    close
    exit",file= paste0(winscp_dir,txt_name));
  
  #write the bat file
  bat_name=paste(file_name,".bat",sep="");
  cat("WinSCP.exe /log=./log /console /script=E:/Projects/R/winscp/",txt_name,sep="",file= paste0(winscp_dir,bat_name))
  
  system(paste0(winscp_dir,bat_name))
  
  
  
}

f_read_bed <- function(file_name, show=TRUE)
{
  #Read the bed file, and assign the name of the filed, and remove the chrM items
  bed_raw=read.table(file_name, header = FALSE, stringsAsFactors=FALSE)
  
  names(bed_raw)=c("chr","start","end","name","score","strand", "thickStart", 
                       "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts")[1: ncol(bed_raw)]
  
  
  chrm_list=grep("chrM",bed_raw$chr,ignore.case=TRUE)
  
  bed_data=bed_raw
  if (length(chrm_list) > 0)
  {
    cat("Remove", length(chrm_list), "ChrM items",'\n')
    bed_data=bed_data[-grep("chrM",bed_data$chr),]
  }
  
  if (show==TRUE)
  {
    print(head(bed_raw))
  }
  
  
  return (bed_data)
}


f_add_table_to_dataframe <- function(dis, index, table_data)
{
  tmp=as.data.frame(table_data)
  rownames(tmp)=tmp$Var1
  
  f_assert( setequal(intersect(rownames(tmp),rownames(dis)),rownames(tmp)), "Input Error")
  
  dis[rownames(tmp),index]=tmp$Freq
  
  return (dis)
}

f_remove_duplicate <- function(data,id_index)
{
  return (data[!duplicated(data[,id_index]),])
}


f_find_file <- function(d_dir, file_name)
{
  
  file_true_name=grep(file_name, dir(d_dir), ignore.case=TRUE, value=TRUE)
  f_assert(length(file_true_name)==1,message=paste("f_find_file Error:",d_dir, file_name, length(file_true_name), 'found'))
  
  if (length(file_true_name)==1)
  {
    return (paste0(d_dir,"/",file_true_name))
  } else
  {
    return (NULL)
  }
  
  
}

f_test_subset <- function(data, test_size=10)
{
  if (nrow(data) == 0)
  {
    return (data)
  }
  return (data[1:min(test_size,nrow(data)),])
}

f_read_pwm_with_name <- function(tf_name)
{
  print(paste0('./data/rSNP/',toupper(tf_name),'.pwm'))
  pwm_matrix=(f_read_pwm( paste0('./data/rSNP/',toupper(tf_name),'.pwm') ))
  col_names =colnames(pwm_matrix)
  for (i in 1:ncol(pwm_matrix))
  {
    one_col = pwm_matrix[,i]
    one_col_sorted=sort(one_col,decreasing=TRUE)
    nt_names=names(which(one_col_sorted>0))
    col_names[i]=paste(nt_names,collapse="/")
  }
  colnames(pwm_matrix)=col_names
  return (pwm_matrix)
  
}

f_read_pwm <- function(path)
{
  
  pwm=read.table(path)
  pwm_true=as.matrix(pwm[,-1])
  
  row.names(pwm_true) = pwm[,1]
  return (PWM(pwm_true))
}


f_pwm_scan_both <- function(pwm,sub)
{
  #It scan the sub using the pwn, and then return the best score, only for only string
  str_len=length(sub)
  score=rep(0,str_len)
  pwm_len=ncol(pwm)
  pwm=as.matrix(pwm)
  
  score_sense=PWMscoreStartingAt(pwm=pwm,subject=sub,starting.at=1:(str_len-pwm_len+1))
  score_antis=PWMscoreStartingAt(pwm=reverseComplement(pwm),subject=sub,starting.at=1:(str_len-pwm_len+1))
  score=data.frame(score=c(score_sense,score_antis),index=c(length(score_sense):1  ,1:length(score_sense)), strand=rep(c("+","-"), each=length(score_sense)))
  idx=which.max(score$score)  
  return (score[idx,])
}


#df=c(1,2,3,4,1,2)
f_duplicated <- function(df)
{
  return (duplicated(df) | duplicated(df[length(df):1])[length(df):1])
}

f_p <- function(...)
{
  return (sprintf(...))
}



f_select_columns <- function(pattern, data, invert_flag, quiet_flag=TRUE)
{
  select_cols=grep(pattern,colnames(data),value = TRUE, invert = invert_flag)
  if (quiet_flag == FALSE)
  {
    cat("Selected colnames: ",select_cols)
  }
  return (data[, select_cols])
}

f_check_data_features <- function(data, necessary =c(""), optional = NULL )
{
  data_cols = colnames(data)
  
  
  for (necessary_pattern in necessary)
  {
    
    if (length(grep(necessary_pattern, colnames(data),value = TRUE)) == 0)
    {
      flog.error("Don't find expected pattern %s", necessary_pattern)
    }
  }
  
  return (1)
  
}


source('f_creat_markdown.R')

f_conditional_output <- function(cmd_str, output_level = 'INFO', comment='')
{
  #cmd_str = 'str(host_data)'
  #output_level = 'DEBUG'
  
  
  library('futile.logger')
  if (flog.threshold() <= output_level )
  {
    cat(comment,'\n')
    eval( parse(text = cmd_str) )
  }
}


f_conditional_output_str <- function(data, output_level = 'INFO', comment='')
{
  #cmd_str = 'str(host_data)'
  #output_level = 'DEBUG'
  
  
  library('futile.logger')
  if (flog.threshold() <= output_level )
  {
    cat(comment,'\n')
    str(data)
  }
}



f_render_markdown <- function(r_file, subname = NULL)
{
  library(knitr)
  library(rmarkdown)
  opts_chunk$set(fig.path = "./result/figures/")
  knitr::opts_chunk$set(echo=FALSE)
  knitr::opts_chunk$set(collapse=TRUE)
  knitr::opts_chunk$set(comment=NA)
  knitr::opts_chunk$set(warning = FALSE)
  
  
  file_prefix = str_replace(r_file,pattern = "[.](r|R)", replacement = '')
  output_dir = f_p('./result/%s', file_prefix)
  dir.create(output_dir, showWarnings = FALSE)
  
  render(input = r_file, output_format = 'html_document',output_dir = f_p('./result/%s', file_prefix), clean = TRUE)
  
  html_file = f_p('%s/%s.html', output_dir, file_prefix)
  
  
  
  
  para_string=f_read_file_name_from_paras(file_prefix)
  
  
  if(!is.null(para_string)){
    display_file = f_p('%s/%s.html', output_dir, para_string)
    
  }else{
    
    if (is.null(subname))
    {
      display_file = html_file
    }else
    {
      display_file = str_replace(html_file, pattern = ".html", replacement = f_p('.%s.html', subname))
    }
    
    
  }
  
  file.rename(from = html_file, to = display_file  )
  file.show(display_file)
  
}


f_built_file_name_from_paras <- function(script_prefix, para_file_name){
  write.table(file = f_p('./filename/%s.file.name.txt', script_prefix),x = para_file_name,quote = FALSE, row.names = F, col.names = F)
}

f_read_file_name_from_paras <- function(script_prefix){
  
  file_loc = f_p('./filename/%s.file.name.txt', script_prefix)
  if (file.exists(file_loc)){
    file_name=as.character(read.table(file = file_loc ))
  }else{
    file_name=NULL
  }
  
  return (file_name)
}


f_header_in_loop <- function(header, level = 3){
  
  if (level >= 3){
    cat(paste0(rep('\n', times = level + 2 ),collapse = ''))
    cat('====================================================','\n')
    
    input_len=str_length(as.character(header))
    padding_len = (50- input_len)/2
    
    left_string = paste0(rep('=', times = floor(padding_len) ),collapse = '')
    right_string = paste0(rep('=', times = ceiling(padding_len) ),collapse = '')
    
    cat(left_string, header, right_string,'\n')
    cat('====================================================','\n')    
  }

}



f_check_na_col <- function(input_df){
  
  return ( colnames(input_df)[(sapply(input_df, function(x)any(is.na(x))))] )
  
}


library(ggplot2)
axis_adjust=theme(text = element_text(size=15), axis.text.y = element_text(colour="black") ,axis.text.x = element_text(colour="black")) 
title_plot=theme(plot.title = element_text(size = 18,face="bold"),axis.title=element_text(size=14,face="bold"))



f_table_to_word <- function(input_data, file_prefix = 'example', data_dir =NULL, rowname = FALSE, browser = F){
  
  library(ReporteRs)
  library(magrittr)
  if (is.null(data_dir)){
    file_name = f_p('%s.docx', file_prefix)
  }else{
    file_name = f_p('%s/%s.docx', data_dir ,file_prefix)
    
  }
  
  
  
  docx( ) %>% 
    addFlexTable( input_data %>%
                    FlexTable( header.cell.props = cellProperties( background.color =  "#003366" ),
                               header.text.props = textBold( color = "white" ),
                               add.rownames = rowname ) %>%
                    setZebraStyle( odd = "#DDDDDD", even = "#FFFFFF" ) ) %>%
    writeDoc( file = file_name )
  
  # open the Word document
  if (browser == T){
    browseURL(file_name)  
  }
  
  
  
}


