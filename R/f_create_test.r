f_create_test  <- function(fun_){
  #Description:Give the fun_name, create the test_scripts as well.
  #input: fun_name,e.g. f_test
  #output: test file name, e.g. t_test
  #By wenqiang, 2013-02-05
  
  
  fun_name <- deparse(substitute(fun_))
  
  test_name=str_replace(fun_name,"f_","t_");
  file_name=paste('./test/',test_name,".r",sep="");
  cat(test_name,
      " <- function(){
      #Des: The test file for ",fun_name,
      "\n#By Wenqiang, ",as.character(Sys.Date()),
      "\n\n\t
      #write the test data
      write.table(, file='./data/test/",test_name,".txt')","  
      test_data=read.table(file='./data/test/",test_name,".txt')","
      
      #debug onece
      debug(",fun_name,")
      mtrace(",fun_name,")
      
      result=",fun_name,"()
      write.table(result, file='./data/test/",test_name,".result')
      ","
      #Normal Cases
      checkEquals(",fun_name,"(), )
      checkEquals(",fun_name,"())
      #Exception Case
      checkException(",fun_name,"())  
      
}\n",
"
options(run.main=FALSE)
if (getOption('run.main', default=TRUE)) {
runTestFile('", file_name, "',testFuncRegexp = '^t_.*')
}",

file=file_name,sep="");
  file.edit(file_name)
  
}