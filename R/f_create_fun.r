f_create_fun <- function(fun_name)
{

  file_name=paste("tmp",".r",sep="");
  cat(fun_name,
" <- function(){
#Description:
#input:
#output:
#By Wenqiang, ",as.character(Sys.Date()),#"\n\n\n}\n")
"\n\n\n}\n",file=file_name);

  
  file.show(file_name)
}