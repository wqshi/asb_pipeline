
import os
import socket
import urllib
import re
import sys
import vcf
import p_script_test as scripttest
from p_mymodule import *
import p_mymodule as my
import p_extract_depth_from_bam_fun as fun
import p_pd as loc
reload(fun)
reload(loc)
import logging
import numpy
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


#This script is the next one after the p_novo_mapping_para.py, part of the ASB pipeline.
#Based on the bam files in the directory, extract the read depth informaiton at the given locaiton.
#Then copy the dp file to the bam directory
#The important parameter to set is tf_list

#Pre:
#*target dir should contain the bam file of tf_list, and feature_list, and the *.het.loc for the tf_list.

#Processing:
#*Extract the read depth info at the given location of het.loc file from matched bam file.

#After:
#For each TF, it will have *.dp file for each feature.

#input:
#When run in the clustdell, the cell type should be given.


server_name = socket.gethostname()
logging.info(server_name)

internal_call = server_name == "loire" and len(sys.argv)<3

feature_only = False
test_ehancer = False
if internal_call:
    #Test case 1 for loire 
    #dnase_dir="/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/"
    dnase_dir="/homed/home/shi/test/t_het_sites_in_narrow_peak/full_data/tmp_depth/"
    target_dir="/homed/home/shi/test/t_het_sites_in_narrow_peak/full_data/"
    test_envi = scripttest.scripttest(target_dir)
    guest_dir = target_dir
    #loc_file   =dnase_dir + "/gm12878-ctcf.het.loc.test"
    #feature_list     =["h4k20me1", "dnase", 'simulate','chrx']
    feature_list     =["simulatectcf"]
    tf_list = ["ctcf"]
    cell = 'gm12878'
    guest_cells=[cell]
    rmdup_flag = True
    cofactor_list_raw = []
    labs = ['sydh','uw','encode']
    mapQ = 0

else:
    #loc_file="/home/shi/projects/chipseq_snp/data2/novo_ctcf/gm12878-ctcf.het.loc"

    import p_mymodule as my
    print sys.argv
    cell = sys.argv[1]
    guest_cells=[cell]
    tf_list=my.f_recieve_list_para(sys.argv[2])

    locker_file=sys.argv[3]
    mapQ = sys.argv[4]
    rmdup_flag = False
    labs = my.f_recieve_list_para( sys.argv[5] )
    target_dir = sys.argv[6]
    #target_dir="/home/shi/projects/chipseq_snp/data2/encode/%s/" % cell    
    #guest_dir = "/home/shi/projects/chipseq_snp/data2/encode/%s/" % guest_cells[0]
    guest_dir = target_dir
    dnase_dir = target_dir + "/tmp_dir/"
    #dnase_dir="/state/partition1/shi/tmp_depth/%s/" % my.f_shell_cmd('echo $JOB_ID', quiet = True).replace('\n', '')
    #my.f_scp_python_script_to_clustdell("p_extract_depth_from_bam_fun.py")
    #my.f_scp_python_script_to_clustdell("p_pd.py")
    reload(fun)
    reload(loc)
    cofactor_list_raw = []
    feature_list=[]
    #feature_list=["methy"]
    #tf_list=["inputigg","inputstd","ctcf"]
    #tf_list=["ctcf","znf143","bhlhe40","ebf1"]
    #tf_list=["znf143","ctcf","ebf1"]
    #tf_list=['brca1', 'chd2', 'elk1', 'max', 'maz', 'mxi1', 'nfya', 'nfyb', 'rad21', 'rfx5', 'smc3', 'stat3',  'tbp', 'usf2']

reload(fun)
print my.f_shell_cmd('echo $HOME', quiet = True).replace('\n','')
my.f_unique_element_in_list(guest_cells)
guest_cell = guest_cells[0]
guest_extract_flag  =  guest_cell != cell #When the host and guest cell are different, means that try to predict variation impact.
my.f_ensure_make_dir(dnase_dir)
logging.info('Node dir name: ' + dnase_dir)
cofactor_list = [ tf_name.lower() for tf_name in cofactor_list_raw ]
logging.info(cofactor_list)

reload(my)
reload(fun)
#if in clustdell, copy the tf's bam file to the local node
#import ipdb; ipdb.set_trace()

for tf in tf_list + feature_list  + cofactor_list:
    tf_bam_pattern = "*%s*bam"%tf
    print tf_bam_pattern
    f_copy_to_dir(target_dir,tf_bam_pattern,dnase_dir,"-u")

if guest_extract_flag == True:
    tf_bam_pattern = "*dnase*bam"
    print tf_bam_pattern
    f_copy_to_dir(guest_dir, tf_bam_pattern, dnase_dir,"-u")





#feature_list=["h4k20me1", "dnase"]
#target_tf_list=["ctcf"]

data_dir=dnase_dir

 #For the extraction in the same cell.


try:
    for loc_tf in tf_list:

        db_file=my.f_create_cell_pair_database_name(target_dir, cell, guest_cells, loc_tf)

        logging.info('Cofactors: %s' % ' '.join(cofactor_list)  )
        fun.f_extract_features_on_location_dp(db_file, loc_tf, cofactor_list,cell, dnase_dir,  rmdup = rmdup_flag, labs=labs, mapQ = mapQ)
        #import ipdb; ipdb.set_trace()
        fun.f_extract_features_on_location_dp(db_file, loc_tf, feature_list, cell, dnase_dir, rmdup= rmdup_flag, mapQ = mapQ, debug = False)

        if feature_only == True:
            next
        
        if "whole" in loc_tf:
            fun.f_extract_features_on_location_dp(db_file, loc_tf, [tf for tf in tf_list if tf != "whole"],        cell, dnase_dir, rmdup=True, loc_cell = cell, mapQ = mapQ)
        else:
            print "TF as feature"
            fun.f_extract_features_on_location_dp(db_file, loc_tf, [loc_tf],        cell, dnase_dir,  rmdup = rmdup_flag,  labs=labs, mapQ = mapQ, debug = False)

        if guest_extract_flag == True:
            fun.f_extract_features_on_location_dp(db_file, loc_tf, ['dnase'], guest_cell, dnase_dir, rmdup= rmdup_flag, guest_extract_flag = True, mapQ = mapQ ,debug = False
                                                  )
except:
        import traceback
        traceback.print_exc()
        logging.warning('Error happend')



#Copy the new generated files to the needed directories
if (internal_call):

    tf_database=loc.data_table(db_file)
    print tf_database.data.alt_simulatectcf_dp_broad1
    print tf_database.data.ref_simulatectcf_dp_broad1
    #print my.f_print_list(list(tf_database.data.columns))

    #Testing for the danse_guest
    #assert all( tf_database.data['alt_h4k20me1_dp'] == tf_database.data['alt_dnase_dp_guest'] )
    #assert all( tf_database.data['alt_h4k20me1_dp'] == tf_database.data['alt_h4k20me1_dp_broad1'])
    test_envi.new_files()
    test_envi.wc_file()
    test_envi.head_file('.*database')
else:
    
    import shutil
    shutil.rmtree(dnase_dir)
    logging.info('Removed the data in the node')
    my.f_write_locer_file(locker_file, server_name)
#[Mon Jul  6 16:48:54 2015] p p_run_cluster_sep.py extract-gm12878-gm12878-simulate-shi p_extract_depth_from_bam_dp.py para gm12878 gm12878 simulate simulate_q30:simulate_rmdup:simulate_mask None 20150706_164852.extract_depth.5NJ9P 0 sydh:uw:haib:uta:embl:encode 0 1
#[Mon Jul  6 16:49:25 2015] p p_run_cluster_sep.py extract-helas3-helas3-simulate-shi p_extract_depth_from_bam_dp.py para helas3 helas3 simulate simulate_q30:simulate_rmdup:simulate_mask None 20150706_164924.extract_depth.AY4OA 0 sydh:uw:haib:uta:embl:encode 0 1
