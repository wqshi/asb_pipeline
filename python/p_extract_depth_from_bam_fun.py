
import os
import socket
import urllib
import re
import sys
import vcf
from p_mymodule import *
import p_mymodule as my
import p_pd as loc
import logging

def f_get_bam_pattern(cell,tf,file_type):
    pattern='.*%s.*%s.*(rep.)?\.%s$'%(cell, tf,file_type)
    return pattern

def f_get_reference_genome():
    server_name = f_get_server_name()
    if server_name == "loire":
        hg_file="/shared2/RAW/shi/data/wgs/hg19.fa"
    else:
        hg_file="/raid2/local/exome/data/hg19/hg19-orderMito.fasta"

    return hg_file

def f_get_files_from_dir(data_dir, cell, tf, file_type):
    file_list=os.listdir(data_dir)
    pattern = f_get_bam_pattern(cell, tf, file_type)
    match_files=grep_list(pattern, file_list)
    return [ data_dir + "/" + f for f in match_files]



def f_readdepth_from_single_vcf(vcf_file, output_file):
    import vcf
    import pandas as pd
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    chip_dp=[]
    #import ipdb; ipdb.set_trace()
    for record in vcf_reader:
        item = [record.CHROM, record.POS, record.REF]
        #print item
        if len(record.ALT) > 1:
            continue

        #print record.ALT[0]
        if record.ALT[0] is None:
            item.extend(["."])
        else:
            item.extend(record.ALT)
            
        #print record.INFO.__class__
        if "DP4" in record.INFO.keys():
            ref_dp=int(record.INFO["DP4"][0]) + int(record.INFO["DP4"][1])
            alt_dp=int(record.INFO["DP4"][2]) + int(record.INFO["DP4"][3])
            item.extend([ref_dp, alt_dp])
            #print item
            chip_dp.append(item)

    if chip_dp != []:
        chipseq_dp=pd.DataFrame(chip_dp,columns=["chr","start","chip_ref","chip_alt","chip_ref_dp","chip_alt_dp"])
        chipseq_dp.to_csv(output_file,sep="\t", index=None)
    else:
        file = open(output_file, "w")
        file.write("chr\tstart\tchip_ref\tchip_alt\tchip_ref_dp\tchip_alt_dp\n")
        file.close()


def f_extract_depth_from_bam(loc_file, bam_file, vcf_file,dp_file, mapQ = 0):
    hg19_file=f_get_reference_genome()
    cmd="~/bin/samtools mpileup -B -q %s -uf %s  -l %s %s | ~/bin/bcftools view -cg - >  %s"%( mapQ ,hg19_file, loc_file, bam_file, vcf_file)
    print cmd
    print f_shell_pipe(cmd)
    f_readdepth_from_single_vcf(vcf_file, dp_file)


def f_extract_features_on_location(loc_file, loc_tf, feature_list, cell, data_dir, extract_type ,rmdup=False, loc_cell=""):
    #Extract the dp information accoriding to the loc_flie (loc_cell and loc_tf) in bam files of feature_list, and loc_tf in the cell data
    #loc_file: the locations where to extract the binding signal
    #loc_tf / loc_cell: the source of loc_file

    #feature_list, cell, data_dir: for the bam files
    
    for tf in feature_list:

        print tf
        bam_files =[]
        if rmdup == True:
            bam_files=f_get_files_from_dir(data_dir, cell, tf, "rmdup.bam")
            rmdup_str = ".rmdup"

        if rmdup == False or bam_files==[]:
            bam_files=f_get_files_from_dir(data_dir, cell, tf, "sorted.bam")
            rmdup_str = ""
        
        if bam_files == []:
            bam_files=f_get_files_from_dir(data_dir, cell, tf, "bam")
            
        
        for bam_file in bam_files:
            print bam_file
            #Get the replicate number
            rep_object=re.match(r".*(Rep[1-9]).*",bam_file,flags=re.IGNORECASE)

            if rep_object == None:
                print "Escape the %s" % bam_file
                continue


            loc_dir = data_dir
            #my.f_get_dir_name_from_file(loc_file)
            print loc_dir
            vcf_file=f_create_file_name(data_dir=loc_dir,cell=cell,tf=tf,suffix= loc_cell + "." + loc_tf + "." + extract_type + rmdup_str + ".dp.vcf",rep=rep_object.group(1))
            
            dp_file=f_create_file_name(data_dir=loc_dir,cell=cell,tf=tf,suffix= loc_cell + "." + loc_tf + "." + extract_type + rmdup_str + ".bam.vcf.dp",rep=rep_object.group(1))
            print vcf_file, dp_file
            f_extract_depth_from_bam(loc_file, bam_file,vcf_file, dp_file)


def f_extract_features_on_location_dp(database_file, loc_tf, feature_list, cell, data_dir, rmdup=False, labs = '', guest_extract_flag = False, mapQ = 0,debug = False):
    #Extract the dp information accoriding to the loc_flie (loc_cell and loc_tf) in bam files of feature_list, and loc_tf in the cell data
    #loc_file: the locations where to extract the binding signal
    #loc_tf / loc_cell: the source of loc_file

    #feature_list, cell, data_dir: for the bam files
    #this version is for collecting all the data in to a central database file

    if debug == True:
        import ipdb; ipdb.set_trace()

    
    tf_database=loc.data_table(database_file)
    loc_file = tf_database.extract_loc(data_dir)


    error_cols = my.grep_list('.*simulate%s_dp_encode'%loc_tf, tf_database.data.columns.tolist())
    my.f_print_list(error_cols)
    #my.f_print_list(tf_database.data.columns.tolist())
    tf_database.drop_feautre(error_cols)

    error_cols = my.grep_list('.*_%s_dp_encode'%loc_tf, tf_database.data.columns.tolist())
    tf_database.drop_feautre(error_cols)
    
    #import ipdb; ipdb.set_trace()

    for tf in feature_list:

        print tf
        bam_files = []
        bam_pattern = '(%s).*%s-%s[\.-].*' % ("|".join(labs), cell, tf)
        if rmdup == True:
            bam_files = my.f_grep_files_from_dir(data_dir, bam_pattern + "rmdup.bam$")
            rmdup_str = ".rmdup"

        if rmdup == False or bam_files==[]:
            bam_files = my.f_grep_files_from_dir(data_dir, bam_pattern + "sorted.bam$")
            rmdup_str = ""
        
        if bam_files == []:
            bam_files = my.f_grep_files_from_dir(data_dir, bam_pattern + "bam$")
        
        #import ipdb; ipdb.set_trace()
        logging.info('bam files:' + data_dir + bam_pattern)
        logging.info(bam_files)
        for bam_file in bam_files:
            #print bam_file
            #Get the replicate number
            rep_object=re.match(r".*(Rep[1-9]).*",bam_file,flags=re.IGNORECASE)
            lab = re.match(r".*(%s).*%s.*"%('|'.join(labs), cell),bam_file,flags=re.IGNORECASE)

            if lab == None:
                print "Lab is missing %s in labs(%s)" % (bam_file, ' '.join(labs))
            lab = lab.group(1)
            if rep_object == None:
                print "Escape the %s" % bam_file
                continue


            loc_dir = data_dir
            #my.f_get_dir_name_from_file(loc_file)
            #print loc_dir
            #vcf_file=f_create_file_name(data_dir=loc_dir,cell=cell,tf=tf,suffix= loc_cell + "." + loc_tf + "." + extract_type + rmdup_str + ".dp.vcf",rep=rep_object.group(1))
            vcf_file=loc_dir + my.f_generate_tmp_file_name("vcf")
            
            
            #dp_file=f_create_file_name(data_dir=loc_dir,cell=cell,tf=tf,suffix= loc_cell + "." + loc_tf + "." + extract_type + rmdup_str + ".bam.vcf.dp",rep=rep_object.group(1))
            dp_file=loc_dir + my.f_generate_tmp_file_name("dp")
            
            #print vcf_file, dp_file
            f_extract_depth_from_bam(loc_file, bam_file,vcf_file, dp_file, mapQ)

            #import ipdb; ipdb.set_trace()
            rep_names={}
            if labs != '':
                for i in range(1, 8):
                    rep_names['Rep%s' % i] = '_%s%s' %( lab,  i)
            else:
                for i in range(1, 8):
                    rep_names['Rep%s' % i] = '_broad%s' % i
                

            if guest_extract_flag == True:
                rep_names['Rep1']='_guest'
            
            feature_data=tf_database.read_feature_replace_name(dp_file, ["chip_ref_dp","chip_alt_dp"],
                                                               ["ref_%s_dp%s"%(tf, rep_names[rep_object.group(1)]), "alt_%s_dp%s"% (tf, rep_names[rep_object.group(1)])])
            #feature_data=tf_database.read_feature(dp_file, tf)
            
            tf_database.merge_feature(feature_data)
       
            #os.remove(dp_file)
            #os.remove(vcf_file)

    #os.remove(loc_file)



import unittest
class TestDatabaseTable(unittest.TestCase):
    def setUp(self):
         test_name = 'unittest'
    def test_readdepth_from_single_vcf(self):
        test_dir = '/homed/home/shi/test/t_het_sites_in_narrow_peak/t_extract_depth_from_bam_fun/'
        vcf_file = test_dir + 'test.vcf'
        output_file = test_dir + 'a.dp'
        #import ipdb; ipdb.set_trace()
        f_readdepth_from_single_vcf(vcf_file, output_file)
        

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase( TestDatabaseTable )
    unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( suite )
    




