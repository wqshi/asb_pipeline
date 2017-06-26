f_create_markdown  <- function(script_name){
  #Description:Give the fun_name, create the test_scripts as well.
  #input: fun_name,e.g. f_test
  #output: test file name, e.g. t_test
  #By wenqiang, 2013-02-05
  
 
  
  mark_name=str_replace(script_name,"s_","m_");
  mark_file=paste(str_split(mark_name,pattern = '[.]')[[1]][1],".Rmd",sep="");
  cat('\n---',
      "\n'",mark_name,"'",
      '\n',
      '---\n',
      '```{r}',
      "\nsource('E:/Projects/R/", script_name  ,"')\n",
      '```',
      '\ncomments\n',
      '```{r}\n',
      "#source('", script_name  ,"')\n",
      '```\n', 
      file= mark_file, sep="")
  
  file.edit(mark_file)

}



f_run_markdown <- function(mark_file, subname = NULL)
{
  library(knitr)
  opts_chunk$set(fig.path = "./result/figures/")
  output_file = paste0('./result/', str_replace(mark_file,pattern = 'Rmd', replacement = 'md'))
  knit( mark_file, output_file)
  
  if (is.null(subname))
  {
    display_file = str_replace(output_file,pattern = ".md", replacement = '.html')  
  }else
  {
    display_file = str_replace(output_file,pattern = ".md", replacement = f_p('.%s.html', subname))
  }
  
  knit2html(input = output_file, display_file)
  file.show(display_file)
}

