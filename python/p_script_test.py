import p_mymodule as my
reload(my)
import os
import time
import os
import datetime as dt

class scripttest:

    test_dir = ""             # mistaken use of a class variable
    start=[]
    new_file_list=[]
    
    def __init__(self, dir_name):
        self.test_dir = dir_name
        self.start=dt.datetime.now()
        self.new_file_list=[]
    def clean_dir(self):
        my.f_grep_and_rm(self.test_dir, ".*",quiet=True)

    def new_files(self):

        self.new_file_list=[]
        for root,dirs,files in os.walk(self.test_dir):  
            for fname in files:
                path=os.path.join(root,fname)
                st=os.stat(path)    
                mtime=dt.datetime.fromtimestamp(st.st_mtime)
                if mtime>self.start:
                    self.new_file_list.append(fname)

        return self.new_file_list
        
    
    def wc_file(self, pattern=".*"):
        print "\n==========================="
        print "============WC Files========="
        print "============================="
        file_list=self.new_files()
        
        os.chdir(self.test_dir)
        for single_file in file_list:
            cmd="wc -l %s" %single_file
            print "\t".join((my.f_shell_cmd(cmd, quiet=True).replace("\n", "").split(" ")))


    def sample_file(self, pattern=".*", wc_flag=True, n=10):
        import random, sys, os
        print "\n==========================="
        print "============Sample Files========="
        print "============================="
        file_list=my.f_grep_files_from_dir(self.test_dir, pattern, path=False)
        os.chdir(self.test_dir)
        for single_file in file_list:
            print "\n=========%s============"%single_file
            if wc_flag == True:
                cmd="wc -l %s" %single_file
                print  "[File lines:]","\t".join((my.f_shell_cmd(cmd, quiet=True).replace("\n", "").split(" ")))
            print ""
            file_handle = open(single_file,"r")
            print("".join(random.sample(file_handle.readlines(), n)))
            file_handle.close()
            
    def head_file(self, pattern=".*", wc_flag=True, n=10):
        import random, sys, os
        print "\n==========================="
        print "============Sample Files========="
        print "============================="
        print self.new_files()
        file_list=my.grep_list(pattern,self.new_files())
        os.chdir(self.test_dir)
        for single_file in file_list:
            print "\n=========%s============"%single_file
            if wc_flag == True:
                cmd="wc -l %s" %single_file
                print  "[File lines:]","\t".join((my.f_shell_cmd(cmd, quiet=True).replace("\n", "").split(" ")))
            print ""
            cmd="head -n %s %s" %(n, single_file)
            my.f_shell_cmd(cmd)

if __name__ == "__main__":
    test_env=scripttest("/homed/home/shi/test/test_data/scripttest")
    print test_env.test_dir

    test_env.clean_dir()
    my.f_copy_to_dir("/homed/home/shi/test/test_data", "test.vcf", test_env.test_dir)

    test_env.new_files()
    
    test_env.wc_file()

    test_env.sample_file(".*vcf")

    test_env.head_file()

