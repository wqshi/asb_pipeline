#!/usr/bin/env python2.7
import os
import socket
#from p_project_meta_data import local_bin,wgs_dir
import re
import sys
import subprocess
import logging

server_name = socket.gethostname()
print server_name
logging.getLogger().setLevel(logging.DEBUG)
import p_mymodule as my
import argparse


if __doc__ is None:
    parser = argparse.ArgumentParser(description="Map ChIP-seq reads to the presonalized genome using nonoalign")

    parser.add_argument('--fastq_dir', help = "directory to store the fastq files.\
                                                                    If this value is empty, it will run the tests", default = '' )
    parser.add_argument('--fastq_file', help = "The compressed fastq.gz file to be mapped. \
the name of the file should follow: lab-celll-tf-rep.fastq.gz",
                                                                     default = '' )
    parser.add_argument("--mode",     help="test or unique.\
                                                                  If it's test, would only map the 1k reads from fastq files. If it's unique,\
                                                                  only keep the reads uniquely mapped to the genome",
                                                                  default="unique")
    parser.add_argument('--short_read', help = "Whether the read length is less than 35bp, e.g. 20.bp\
                                                                       Have been used to map DHS data.", default = False )

    parser.add_argument('--data_dir', help = "Copy the files to this directory to process. If not provided,\
                                                                    will use fastq_dir as default.", default = '' )
    parser.add_argument('--loc_bin', help = "The dir name for novoalign", default = '' )
    parser.add_argument('--wgs_dir', help = "The dir name for novoalign index file.", default = '' )                                                                                                                                       


                                                                 

    
    args = parser.parse_args()


    fastq_gz_file=args.fastq_file
    fastq_dir = args.fastq_dir
    mode= my.f_recieve_list_para( args.mode )
    short_read=args.short_read
    if args.data_dir == '':
        
        data_dir = '%s/novo-%s' % (fastq_dir, fastq_gz_file.split('.')[0])
    else:
        data_dir = '%s/novo-%s' % (args.data_dir, fastq_gz_file.split('.')[0])
    
    cell=my.f_parse_rep_file_name(fastq_gz_file)[1]
    desination_dir=fastq_dir
    



def f_bam_remove_dup(bam_file, data_dir, head_dir ,picard_java_lib_path):
    import p_mymodule as my
    cell_tf_name=my.f_get_prefix(bam_file)
    output_file=data_dir + "/" + cell_tf_name + ".rmdup.bam"
    input_bam= data_dir + "/" + bam_file
    cmd = "java -jar %s/MarkDuplicates.jar I=%s O=%s M=%s.duplicate_report.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true" %(picard_java_lib_path,  input_bam, output_file,  data_dir + "/"+cell_tf_name)
    print cmd
    my.f_shell_cmd(cmd)
    output_pattern = '%s.(rmdup.bam|duplicate_report.txt)' % (cell_tf_name)
    my.f_grep_and_copy(data_dir,output_pattern,head_dir)
    os.remove(output_file)
    os.remove( os.path.join( data_dir, cell_tf_name + '.duplicate_report.txt'))



def f_novo_mapping(fastq_dir, fastq_gz_file, data_dir, wgs_dir, cell, local_bin, desination_dir, mode=None, short_read=False):

    my.f_ensure_make_dir(data_dir)
    my.f_copy_to_dir(fastq_dir, fastq_gz_file, data_dir)
    my.f_unzip_targz(data_dir + "/" +fastq_gz_file)
    fastq_file=data_dir + "/" + fastq_gz_file.replace(".gz","")

    cell_nix_file=wgs_dir + "/" + cell + ".nix"
    #file_prefix = my.f_get_prefix(fastq_file);
    file_prefix = fastq_file.replace(".fastq","")# + '.'.join(mode)
    map_stats_file= file_prefix + ".stats.txt"
    map_bam_file=file_prefix + ".sam.map"

    mode_string = ''
    if 'test' in mode:
        mode_string=" -#1k "
        desination_dir = desination_dir + '/test/'
        my.f_ensure_make_dir(desination_dir)
    
    if 'unique' in mode: 
        mode_string= mode_string + " -r None"

    if mode_string == '' and mode != None:
        logging.error('Unkonwn mode: %s ' % mode)
    
    #import ipdb; ipdb.set_trace()
    match_object=re.match(".*methy.*", fastq_gz_file)
    if match_object != None:
        #-F ILMFQ
        map_cmd="%s/novoalign  -d %s -f %s -o SAM %s 2> %s > %s"% (local_bin, cell_nix_file, fastq_file, mode_string ,map_stats_file, map_bam_file )
    elif short_read == True:
        map_cmd="%s/novoalign  -d %s -f %s -l 20 -o SAM %s 2> %s > %s"% (local_bin, cell_nix_file, fastq_file, mode_string ,map_stats_file, map_bam_file )
    else:
        map_cmd="%s/novoalign -F ILMFQ -d %s -f %s -o SAM %s 2> %s > %s; echo $?"% (local_bin, cell_nix_file, fastq_file, mode_string , map_stats_file, map_bam_file )

    logging.info("Map cmd: " + map_cmd)
    first_try=my.f_shell_pipe(map_cmd)


    map_stats = ''
    if first_try == '0\n':
        map_stats = 'ILMFQ';
    else:
        logging.warning('Novo output: first try faild ' + first_try)
        
        map_cmd="%s/novoalign  -d %s -f %s -o SAM %s 2> %s > %s;echo $?"% (local_bin, cell_nix_file, fastq_file, mode_string ,map_stats_file, map_bam_file )
        logging.info("2nd Map cmd: " + map_cmd)
        sencond_try  = my.f_shell_pipe(map_cmd)
        logging.info('Sencond try: ' + sencond_try )
        if sencond_try == '0\n':
            map_stats = 'default'
        
    sort_cmd="samtools view -bS %s | samtools sort - %s"%(map_bam_file, file_prefix + ".sorted")
    my.f_shell_pipe(sort_cmd)


    bam_file = "%s.sorted.bam"%os.path.basename(file_prefix)
    #f_bam_remove_dup(bam_file, data_dir, desination_dir, picard_java_lib_path)
    
    
    
    my.f_copy_to_dir(data_dir, "%s.stats.txt"%os.path.basename(file_prefix), desination_dir)
    my.f_copy_to_dir(data_dir, "%s.sorted.bam"%os.path.basename(file_prefix), desination_dir)

    os.remove(map_stats_file)
    os.remove("%s.sorted.bam"%file_prefix)
    os.remove(fastq_file)        
    os.remove(map_bam_file)
  

    if not os.listdir(data_dir):
        logging.info('Empty dir: %s ' % data_dir)
        os.rmdir(data_dir)
    else:
       logging.warning( "Not empty dir" )
       logging.warning(os.listdir(data_dir))

    return map_stats

#Main Part
if __name__ == "__main__":

    #Run the unit tests
    if server_name == 'loire' and fastq_dir == '':
        import unittest
        import tests.test_novo_mapping as test
        
        suite = unittest.TestLoader().loadTestsFromTestCase( test.Test_novo_mapping )
        unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( suite )
##################Real Run####################
    else:
        try:
            f_novo_mapping(fastq_dir, fastq_gz_file, data_dir, wgs_dir, cell, local_bin, desination_dir, mode,  short_read=short_read)
        except:
            import traceback
            traceback.print_exc()
            logging.warning('Error happend in mapping')
            logging.warning(os.listdir(data_dir))

        import shutil
        if os.path.exists(data_dir) and mode != 'test':
            shutil.rmtree(data_dir, ignore_errors = True)


