#Pre-condition:

#*download the narrowPeak file in the bed dir

#Process:
#* find every narrowPeak file, sharp peak into 100 bp.
#* Overlap 100bp peak with each cell's vcf file in the wgs_dir, and then filter the heterozygous sites (*.het.loc file)

#After:
#* for each narrowPeak file, there will be a *.het.loc file.

#Part of the ASB pipeline after the  narrowPeak download.

#Input: when submit to the clustdell, need to specifiy the cell name

#the dp version is trying to using one data structure through the whole process.

import os
import socket
import urllib
import re
import sys
import p_mymodule as my
import p_script_test as scripttest
import pandas as pd
import p_pd as dt
import logging
import pybedtools
logging.getLogger().setLevel(logging.DEBUG)
import shutil
reload(my)
def f_check_chrUn_bed(input_file, description):
    new_file = os.path.join(os.path.dirname(input_file), my.f_generate_tmp_file_name(description))
    cmd = "sed '/chrUn/d' %s > %s" %(input_file, new_file)
    cmd = "grep -P  '^chr[0-9xyXY]{1,2}\t' %s > %s" %(input_file, new_file)
    my.f_shell_cmd(cmd)
    return new_file


server_name = socket.gethostname()
print server_name

internal_call = server_name == "loire" and len(sys.argv)<3

if internal_call:
    bed_dir="/homed/home/shi/test/t_het_sites_in_narrow_peak/"
    test_envi = scripttest.scripttest(bed_dir)
    wgs_dir="/homed/home/shi/projects/wgs"
    host_cell="gmtest"
    #host_cell = 'helatest'
    guest_cells=["helatest"]
    #guest_cells=["helatest"]#use location from these files to extract signal
    tf_list=["ctcf"]
    labs=["sydh","haib","uw","uchicago","uta",'embl']
    narrowPeak_flag = False
    #pattern="*.narrowPeak"
else:
    
    print sys.argv
    cell=sys.argv[1]
    host_cell=cell
    #bed_dir="/home/shi/projects/chipseq_snp/data2/encode/%s/" % cell
    #wgs_dir="/home/shi/projects/wgs/complete_genome/"

    #pattern="*.narrowPeak"
    #pattern="*whole.narrowPeak"
    
    tf_list = my.f_recieve_list_para(sys.argv[2])
    guest_cells = my.f_recieve_list_para(sys.argv[3])
    locker_file = sys.argv[4]
    labs = my.f_recieve_list_para(sys.argv[5])
    narrowPeak_flag = True
    bed_dir = sys.argv[6]
    wgs_dir = sys.argv[7]

    #script_name = "p_pd.py"
    #my.f_scp_python_script_to_clustdell(script_name)

reload(dt)
reload(my)

logging.info('NarrowPeak mode: %s' % narrowPeak_flag)

for tf in tf_list: # for each TF of host_cell.

    #labs=["sydh","haib","uw","uchicago","uta"]
    pattern="(%s).*%s.*%s.narrowPeak$"%("|".join(labs),host_cell, tf)
    logging.info('Pattern: %s' % pattern)
    narrowPeak_list=my.f_grep_files_from_dir(bed_dir, pattern)
    logging.info('Pattern list: %s' % '\n'.join(narrowPeak_list))
    #Genertate the 101bp files.
    if tf not in ['simulate', 'subtract50k']:
        narrowPeak_list_100bp = []
        for narrowPeak_file in narrowPeak_list:


            if my.grep_file('^chrUn', narrowPeak_file) != None:
                logging.error('Detected chrUn:' + narrowPeak_file)
                tmp_narrowPeak = f_check_chrUn_bed(narrowPeak_file, 'narrowPeak')
                os.remove(narrowPeak_file)
                os.rename(tmp_narrowPeak, narrowPeak_file)
            logging.info('host narrowPeak files: ' +narrowPeak_file)        
            #narrowPeak_file=narrowPeak_list[0]
            if "Rep" in narrowPeak_file:
                print "Skip %s" % narrowPeak_file;
                continue;

            print os.path.basename(narrowPeak_file)
            output_file=os.path.dirname(narrowPeak_file) + '/' + my.f_generate_tmp_file_name('101bp.bed') 
            logging.info('Output file: ' + output_file)
            narrowPeak_list_100bp.append(output_file)

            if narrowPeak_flag == True:
                cmd="f_narrow_peak_100bp %s %s"%( narrowPeak_file, output_file)
                my.f_call_shell_fun(cmd)
            else:
                import shutil
                shutil.copy(narrowPeak_file, output_file)

                
            bed_101bp_file=narrowPeak_file + ".101bp.bed"
            if not os.path.isfile(bed_101bp_file):
                shutil.copy2(output_file, bed_101bp_file)
    else:
        narrowPeak_list_100bp = narrowPeak_list
            
    logging.info('narrowPeak_list_100bp: %s' % '\n'.join(narrowPeak_list_100bp) )
    #The merged file for host cell should be same for all the guest cell. but due to parallele problems, here renamed for guestcell
    merged_file = my.f_generate_tmp_file_name('mergePeak')
    logging.info('Mergerd narrowPeak file:' + merged_file)
    output_file = my.f_cat_multiple_files(narrowPeak_list_100bp, merged_file)



    
    output_file_checked = f_check_chrUn_bed(output_file, 'mergePeak')
    logging.debug('Before merge: %s lines' % pybedtools.BedTool(output_file_checked).count())
    merged_bed = pybedtools.BedTool(output_file_checked).sort().merge()
    logging.debug( "After merging: %s lines" % merged_bed.count())
    merged_bed.saveas(output_file)
       
    location_col_names=["chr","start","ref","alt","genotype","wgs_ref_dp","wgs_alt_dp","het_type","cell"]
    location_of_guestcells = pd.DataFrame(columns=location_col_names)
    
    for guest_cell in [host_cell] + guest_cells:
        logging.debug("guest cell %s " % guest_cell)
        vcf_file = bed_dir + my.f_generate_tmp_file_name("wgs.vcf")
        het_dp_file = bed_dir + my.f_generate_tmp_file_name(guest_cell + ".wgs.dp.het")
        alt_dp_file  = bed_dir + my.f_generate_tmp_file_name(guest_cell + ".wgs.dp.alt")
        if 'none' in guest_cell:
            continue
        #elif "gm" in guest_cell:
        #    print("Interpret as gm cell")
        #    overlab_cmd = "sed '/^[0-9]/s/^/chr/g' %s/%s.vcf | intersectBed -u -a stdin -b %s > %s" %(wgs_dir, guest_cell, output_file, vcf_file)
        else:
            overlab_cmd = "intersectBed -u -a %s/%s.vcf -b %s > %s" %(wgs_dir, guest_cell, output_file, vcf_file)

        my.f_shell_cmd(overlab_cmd)
        grep_het_cmd = "f_complete_genome_read_depth %s |  f_grep_legal_snp | sed '/^#/d' | grep -v '1/1' > %s"% (vcf_file, het_dp_file)
        grep_alt_cmd = "f_complete_genome_read_depth %s |  f_grep_legal_snp | sed '/^#/d' | grep '1[/\|]1' > %s"% (vcf_file, alt_dp_file)
        
        
        logging.info(grep_het_cmd)
        #import ipdb; ipdb.set_trace()
        my.f_shell_fun_pipe(grep_het_cmd)

        logging.info(grep_alt_cmd)
        my.f_shell_fun_pipe(grep_alt_cmd)


        vcf_header=["chr","start","ref","alt","genotype","wgs_ref_dp","wgs_alt_dp"]

        het_sites_dp=pd.io.parsers.read_csv(het_dp_file, sep="\t", parse_dates=False, header=None,  index_col=False )
        het_sites_dp.columns=vcf_header
        
        alt_sites_dp =pd.io.parsers.read_csv(alt_dp_file, sep="\t", parse_dates=False, header=None,  index_col=False )
        alt_sites_dp.columns=vcf_header
        
        het_sites_dp["het_type"]="het"
        alt_sites_dp["het_type"] = "alt"

        combined_sites=pd.concat([het_sites_dp,alt_sites_dp])
        combined_sites["cell"]=guest_cell
        #logging.debug(combined_sites.head())

        
        location_of_guestcells = pd.concat([location_of_guestcells, combined_sites])
        logging.debug(location_of_guestcells.shape)
        expected_cols=["chr","start","cell"]
        duplicated_rows = my.f_duplicated_index(location_of_guestcells, expected_cols)
        if any(duplicated_rows):
            logging.info("=======duplicated data in [%s]===========", guest_cell)
            print location_of_guestcells.ix[duplicated_rows, :]

        
        os.remove(vcf_file)
        os.remove(alt_dp_file)
        os.remove(het_dp_file)
        
    data_file= my.f_create_cell_pair_database_name(bed_dir, host_cell, guest_cells, tf)
    

    tf_db=dt.data_table(data_file, location_of_guestcells)
    tf_db.save_data()

if internal_call == False:
    my.f_write_locer_file(locker_file, server_name)
else:
    test_envi.new_files()
    test_envi.wc_file()
    test_envi.head_file()



#[Fri Sep  5 16:23:49 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 .*test.narrowPeak 1 1 1 1
#[Fri Sep  5 16:28:08 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 .*test.narrowPeak 1 1 1 1
#[Fri Sep  5 16:30:00 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Fri Sep  5 20:52:48 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Fri Sep  5 21:33:37 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Fri Sep  5 21:34:06 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Fri Sep  5 22:15:29 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Fri Sep  5 22:17:29 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Fri Sep  5 22:34:14 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Sat Sep  6 08:54:42 2014] p p_run_clustdell_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Sat Sep  6 09:09:41 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Sat Sep  6 09:32:59 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *test.narrowPeak 1 1 1 1
#[Sat Sep  6 10:25:36 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *(ctcf|bhlhe40|ebf1|znf143).narrowPeak 1 1 1 1
#[Sat Sep  6 10:26:33 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.narrowPeak 1 1 1 1
#[Mon Sep  8 10:38:16 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.narrowPeak 1 1 1 1
#[Mon Sep  8 11:09:52 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.narrowPeak 1 1 1 1
#[Mon Sep  8 11:10:10 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.narrowPeak 1 1 1 1
#[Mon Sep  8 11:11:04 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.narrowPeak 1 1 1 1
#[Mon Sep  8 11:13:22 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.narrowPeak 1 1 1 1
#[Mon Sep  8 11:34:48 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.narrowPeak 20140908_113447.het_loc 1 1 1
#[Mon Sep  8 12:05:16 2014] p p_run_cluster_sep.py het_loc p_het_sites_in_narrow_peak.py gm12878 *.(test).narrowPeak 20140908_120516.het_loc 1 1 1
#[Mon Sep  8 12:10:05 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 *.(test).narrowPeak 20140908_121004.het_loc 1 1 1
#[Mon Sep  8 12:16:56 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 test 20140908_121655.het_loc 1 1 1
#[Mon Sep  8 13:20:38 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 ctcf:bhlhe40:ebf1:znf143 20140908_132037.het_loc 1 1 1
#[Tue Sep  9 22:01:38 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 test 20140909_220137.het_loc 1 1 1
#[Tue Sep  9 22:04:36 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 test 20140909_220435.het_loc 1 1 1
#[Tue Sep  9 22:18:22 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 test 20140909_221821.het_loc 1 1 1
#[Tue Sep  9 22:29:57 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 test 20140909_222956.het_loc 1 1 1
#[Tue Sep  9 22:37:58 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 ctcf:bhlhe40:znf143:ebf1:brca1:elk1:max:nfya:nfyb:rfx5:stat3:tbp:usf2 20140909_223757.het_loc 1 1 1
#[Wed Sep 10 10:05:01 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py gm12878 test 20140910_100459.het_loc 1 1 1
#[Wed Sep 10 10:12:05 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test 20140910_101204.het_loc 1 1 1
#[Wed Sep 10 10:19:08 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test 20140910_101907.het_loc 1 1 1
#[Wed Sep 10 10:23:40 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test 20140910_102339.het_loc 1 1 1
#[Wed Sep 10 11:13:22 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test 20140910_111321.het_loc 1 1 1
#[Wed Sep 10 11:21:18 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test helas3 20140910_112117.het_loc 1 1
#[Wed Sep 10 11:31:38 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test helas3 20140910_113137.het_loc 1 1
#[Wed Sep 10 11:35:33 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test helas3 20140910_113532.het_loc 1 1
#[Wed Sep 10 11:45:15 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test helas3 20140910_114514.het_loc 1 1
#[Wed Sep 10 11:49:09 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 test helas3 20140910_114908.het_loc 1 1
#[Wed Sep 10 12:04:56 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3 ctcf:bhlhe40:znf143:ebf1:brca1:elk1:max:nfya:nfyb:rfx5:stat3:tbp:usf2 helas3 20140910_120456.het_loc 1 1
#[Wed Sep 10 12:12:54 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3test test helas3test 20140910_121251.het_loc 1 1
#[Wed Sep 10 12:13:50 2014] p p_run_cluster_sep.py test-het_loc p_het_sites_in_narrow_peak.py helas3test test helas3test 20140910_121349.het_loc 1 1
#[Wed Sep 10 12:20:31 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20140910_122030.het_loc 1 1
#[Wed Sep 10 12:33:20 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20140910_123319.het_loc 1 1
#[Wed Sep 10 12:36:44 2014] p p_run_cluster_sep.py het_loc-helas3 p_het_sites_in_narrow_peak.py helas3 ctcf:bhlhe40:znf143:ebf1:brca1:elk1:max:nfya:nfyb:rfx5:stat3:tbp:usf2 helas3 20140910_123643.het_loc 1 1
#[Wed Sep 10 13:44:20 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py testcell test testcell 20140910_134419.het_loc 1 1
#[Wed Sep 10 15:06:59 2014] p p_run_cluster_sep.py het_loc-helas3 p_het_sites_in_narrow_peak.py helas3 ctcf:znf143:brca1:elk1:max:nfya:nfyb:rfx5:stat3:tbp:usf2 helas3 20140910_150658.het_loc 1 1
#[Tue Sep 30 15:46:48 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20140930_154647.het_loc 1 1
#[Tue Sep 30 15:51:30 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20140930_155129.het_loc 1 1
#[Tue Sep 30 15:54:09 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20140930_155408.het_loc 1 1
#[Tue Sep 30 16:01:57 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20140930_160157.het_loc 1 1
#[Tue Sep 30 16:05:13 2014] p p_run_cluster_sep.py het_loc-gm12878 p_het_sites_in_narrow_peak.py gm12878 ebf1 gm12878 20140930_160512.het_loc 1 1
#[Tue Sep 30 16:21:52 2014] p p_run_cluster_sep.py het_loc-gm12878 p_het_sites_in_narrow_peak.py gm12878 ctcf:znf143:brca1:elk1:max:nfya:nfyb:rfx5:stat3:tbp:usf2 gm12878 20140930_162151.het_loc 1 1
#[Tue Sep 30 22:24:54 2014] p p_run_cluster_sep.py het_loc-gm12xxx p_het_sites_in_narrow_peak.py gm12xxx ctcf gm12864:gm12891:gm12892:gm19238:gm19239:gm19240 20140930_222453.het_loc 1 1
#[Thu Oct 16 08:57:44 2014] p p_run_cluster_sep.py het_loc-gm12xxx p_het_sites_in_narrow_peak.py gm12xxx ctcf gm12864:gm12891:gm12892:gm19238:gm19239:gm19240 20141016_085743.het_loc 1 1
#[Tue Oct 28 09:49:37 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20141028_094936.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/ 1
#[Wed Oct 29 09:30:54 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20141029_093053.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Wed Oct 29 09:33:44 2014] p p_run_cluster_sep.py het_loc-test p_het_sites_in_narrow_peak.py helas3 test helas3 20141029_093344.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Wed Oct 29 09:36:14 2014] p p_run_cluster_sep.py het_loc-gm12878 p_het_sites_in_narrow_peak.py gm12878 znf143:brca1:elk1:max:nfya:nfyb:rfx5:stat3:tbp:usf2 gm12878 20141029_093613.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/gm12878/ 1 1 1 1
#[Wed Oct 29 09:50:00 2014] p p_run_cluster_sep.py het_loc-gm12878 p_het_sites_in_narrow_peak.py gm12878 znf143:brca1:elk1:max:nfya:nfyb:rfx5:stat3:tbp:usf2 gm12878 20141029_094959.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/gm12878/ 1 1 1 1
#[Wed Oct 29 12:15:49 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak.py helas3 test helas3 20141029_121547.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Wed Oct 29 12:31:10 2014] p p_run_cluster_sep.py het_loc-gm12878-het p_het_sites_in_narrow_peak.py gm12878 bhlhe40:ctcf:znf143:ebf1:tbp:usf2 gm12878 20141029_123109.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/gm12878/ 1 1 1 1
#[Wed Oct 29 14:02:20 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak.py helas3 test helas3 20141029_140219.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov  6 14:07:07 2014] p p_run_cluster_sep.py het_loc-gm12878-het p_het_sites_in_narrow_peak.py gm12878 bhlhe40:ctcf:znf143:ebf1:tbp:usf2 gm12878 20141106_140706.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/gm12878/ 1 1 1 1
#[Thu Nov 13 17:01:13 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_170111.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:06:22 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_170621.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:08:56 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_170855.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:18:37 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_171836.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:19:33 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_171932.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:21:15 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_172115.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:22:08 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_172207.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:22:44 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_172243.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Thu Nov 13 17:25:32 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141113_172531.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 11:29:21 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_112920.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 11:35:15 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_113514.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 11:46:42 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_114642.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 11:53:39 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_115338.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 11:56:49 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_115648.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 12:00:58 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_120057.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 12:08:26 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_120826.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 12:15:40 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20141114_121539.het_loc /homed/home/shi/anthony/tfbs_pwm/rsnp/chipseq_snv/cross_cell/test/ 1 1 1 1
#[Fri Nov 14 15:20:38 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141114_152037.het_loc 1 1 1 1 1
#[Fri Nov 14 15:24:04 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141114_152403.het_loc 1 1 1 1 1
#[Fri Nov 14 15:30:20 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141114_153019.het_loc 1 1 1 1 1
#[Fri Nov 14 15:33:18 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141114_153317.het_loc 1 1 1 1 1
#[Fri Nov 14 15:55:06 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141114_155504.het_loc 1 1 1 1 1
#[Fri Nov 14 22:42:43 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141114_224243.het_loc 1 1 1 1 1
#[Fri Nov 14 22:51:24 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141114_225123.het_loc 1 1 1 1 1
#[Fri Nov 14 22:57:58 2014] p p_run_cluster_sep.py het_loc-gm12878-het p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40:ctcf:znf143:ebf1:tbp:usf2 gm12864:gm12891:gm12892:gm19238:gm19239:gm19240 20141114_225757.het_loc 1 1 1 1 1
#[Fri Nov 14 23:08:19 2014] p p_run_cluster_sep.py het_loc-gm12878-het p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864:gm12891:gm12892:gm19238:gm19239:gm19240 20141114_230818.het_loc 1 1 1 1 1
#[Sat Nov 15 10:28:43 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141115_102842.het_loc 1 1 1 1 1
#[Sat Nov 15 11:44:32 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141115_114432.het_loc 1 1 1 1 1
#[Sat Nov 15 11:46:53 2014] p p_run_cluster_sep.py het_loc-test-het p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141115_114652.het_loc 1 1 1 1 1
#[Sat Nov 15 21:24:11 2014] p p_run_cluster_sep.py het_loc-gm12878-het p_het_sites_in_narrow_peak_dp.py gm12878 tbp:usf2:znf143 helas3 20141115_212409.het_loc 1 1 1 1 1
#[Sat Nov 15 22:01:43 2014] p p_run_cluster_sep.py het_loc-helas3-het p_het_sites_in_narrow_peak_dp.py helas3 tbp:usf2:znf143 gm12878 20141115_220142.het_loc 1 1 1 1 1
#[Sun Nov 16 09:04:58 2014] p p_run_cluster_sep.py het_loc-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_090457.het_loc 1 1 1 1 1
#[Sun Nov 16 09:08:59 2014] p p_run_cluster_sep.py het_loc-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_090858.het_loc 1 1 1 1 1
#[Sun Nov 16 09:08:59 2014] p p_run_cluster_sep.py het_loc-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_090858.het_loc 1 1 1 1 1
#[Sun Nov 16 09:08:59 2014] p p_run_cluster_sep.py het_loc-gm19238-all p_het_sites_in_narrow_peak_dp.py gm19238 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_090858.het_loc 1 1 1 1 1
#[Sun Nov 16 09:08:59 2014] p p_run_cluster_sep.py het_loc-gm12892-all p_het_sites_in_narrow_peak_dp.py gm12892 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_090858.het_loc 1 1 1 1 1
#[Sun Nov 16 09:09:00 2014] p p_run_cluster_sep.py het_loc-gm19239-all p_het_sites_in_narrow_peak_dp.py gm19239 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_090858.het_loc 1 1 1 1 1
#[Sun Nov 16 09:09:00 2014] p p_run_cluster_sep.py het_loc-gm19240-all p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_090859.het_loc 1 1 1 1 1
#[Sun Nov 16 09:26:50 2014] p p_run_cluster_sep.py het_loc-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_092650.het_loc 1 1 1 1 1
#[Sun Nov 16 10:13:36 2014] p p_run_cluster_sep.py het_loc-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_101335.het_loc 1 1 1 1 1
#[Sun Nov 16 10:51:10 2014] p p_run_cluster_sep.py het_loc-gm12892-all p_het_sites_in_narrow_peak_dp.py gm12892 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_105109.het_loc 1 1 1 1 1
#[Sun Nov 16 10:51:10 2014] p p_run_cluster_sep.py het_loc-gm19238-all p_het_sites_in_narrow_peak_dp.py gm19238 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_105109.het_loc 1 1 1 1 1
#[Sun Nov 16 10:51:10 2014] p p_run_cluster_sep.py het_loc-gm19239-all p_het_sites_in_narrow_peak_dp.py gm19239 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_105109.het_loc 1 1 1 1 1
#[Sun Nov 16 10:51:10 2014] p p_run_cluster_sep.py het_loc-gm19240-all p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_105109.het_loc 1 1 1 1 1
#[Sun Nov 16 10:51:10 2014] p p_run_cluster_sep.py het_loc-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_105109.het_loc 1 1 1 1 1
#[Sun Nov 16 10:51:10 2014] p p_run_cluster_sep.py het_loc-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_105109.het_loc 1 1 1 1 1
#[Sun Nov 16 11:48:16 2014] p p_run_cluster_sep.py het_loc-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_114815LWOFC.het_loc 1 1 1 1 1
#[Sun Nov 16 11:48:16 2014] p p_run_cluster_sep.py het_loc-gm12892-all p_het_sites_in_narrow_peak_dp.py gm12892 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_114815BEAIB.het_loc 1 1 1 1 1
#[Sun Nov 16 11:48:16 2014] p p_run_cluster_sep.py het_loc-gm19239-all p_het_sites_in_narrow_peak_dp.py gm19239 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_114815MVT45.het_loc 1 1 1 1 1
#[Sun Nov 16 11:48:16 2014] p p_run_cluster_sep.py het_loc-gm19240-all p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_114815QMYK0.het_loc 1 1 1 1 1
#[Sun Nov 16 11:48:16 2014] p p_run_cluster_sep.py het_loc-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_114815DAU48.het_loc 1 1 1 1 1
#[Sun Nov 16 11:48:16 2014] p p_run_cluster_sep.py het_loc-gm19238-all p_het_sites_in_narrow_peak_dp.py gm19238 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_114815RZOKI.het_loc 1 1 1 1 1
#[Sun Nov 16 20:55:28 2014] p p_run_cluster_sep.py het_loc-gm12878-het_loc p_het_sites_in_narrow_peak_dp.py gm12878 ctcf 20141116_2055263K9DR.het_loc 1 1 1 1 1 1
#[Sun Nov 16 20:57:25 2014] p p_run_cluster_sep.py het_loc-shi-gm12878-het_loc p_het_sites_in_narrow_peak_dp.py gm12878 ctcf 20141116_205724J68UO.het_loc 1 1 1 1 1 1
#[Sun Nov 16 21:40:41 2014] p p_run_cluster_sep.py het_loc-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864:gm12891:gm12892:gm19238:gm19239:gm19240:helas3 20141116_214040TDNAB.het_loc 1 1 1 1 1
#[Sun Nov 16 21:42:34 2014] p p_run_cluster_sep.py het_loc-shi-gm19240-all p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_214233GBQS1.het_loc 1 1 1 1 1
#[Sun Nov 16 21:42:34 2014] p p_run_cluster_sep.py het_loc-shi-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_2142333XENG.het_loc 1 1 1 1 1
#[Sun Nov 16 21:42:34 2014] p p_run_cluster_sep.py het_loc-shi-gm19239-all p_het_sites_in_narrow_peak_dp.py gm19239 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_2142335X5PE.het_loc 1 1 1 1 1
#[Sun Nov 16 21:42:34 2014] p p_run_cluster_sep.py het_loc-shi-gm12892-all p_het_sites_in_narrow_peak_dp.py gm12892 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_214233XQ5D4.het_loc 1 1 1 1 1
#[Sun Nov 16 21:42:34 2014] p p_run_cluster_sep.py het_loc-shi-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_214233ZRR7P.het_loc 1 1 1 1 1
#[Sun Nov 16 21:42:34 2014] p p_run_cluster_sep.py het_loc-shi-gm19238-all p_het_sites_in_narrow_peak_dp.py gm19238 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_214233LAK6Q.het_loc 1 1 1 1 1
#[Sun Nov 16 22:18:39 2014] p p_run_cluster_sep.py het_loc-shi-gm12878-het_loc p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864:gm12891:gm12892:gm19238:gm19239:gm19240:helas3 20141116_221838EKGXB.het_loc 1 1 1 1 1
#[Sun Nov 16 22:19:53 2014] p p_run_cluster_sep.py het_loc-allq-gm12878-het_loc p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864:gm12891:gm12892:gm19238:gm19239:gm19240:helas3 20141116_221952751IH.het_loc 1 1 1 1 1
#[Sun Nov 16 22:40:37 2014] p p_run_cluster_sep.py het_loc-allq-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864:gm12891:gm12892:gm19238:gm19239:gm19240:helas3 20141116_224036FHMB2.het_loc 1 1 1 1 1
#[Sun Nov 16 23:16:54 2014] p p_run_cluster_sep.py het_loc-allq-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 tbp:usf2:znf143 helas3 20141116_231653C8UNX.het_loc 1 1 1 1 1
#[Sun Nov 16 23:17:28 2014] p p_run_cluster_sep.py het_loc-allq-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 tbp:usf2:znf143 helas3 20141116_231727NYPKP.het_loc 1 1 1 1 1
#[Sun Nov 16 23:17:41 2014] p p_run_cluster_sep.py het_loc-allq-helas3-all p_het_sites_in_narrow_peak_dp.py helas3 tbp:usf2:znf143 gm12878 20141116_231740VDHZF.het_loc 1 1 1 1 1
#[Sun Nov 16 23:30:13 2014] p p_run_cluster_sep.py het_loc-allq-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141116_233012D9WAQ.het_loc 1 1 1 1 1
#[Mon Nov 17 12:05:23 2014] p p_run_cluster_sep.py het_loc-allq-helas3-all p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20141117_120522P2QIM.het_loc 1 1 1 1 1
#[Mon Nov 17 12:16:36 2014] p p_run_cluster_sep.py het_loc-shi-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141117_121635XU87G.het_loc 1 1 1 1 1
#[Mon Nov 17 13:39:12 2014] p p_run_cluster_sep.py het_loc-shi-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141117_133911J7GIG.het_loc 1 1 1 1 1
#[Mon Nov 17 17:08:54 2014] p p_run_cluster_sep.py het_loc-shi-gm12891-het_loc p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141117_1708533DKZW.het_loc 1 1 1 1 1
#[Mon Nov 17 17:17:19 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-het_loc p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141117_17171837ADO.het_loc 1 1 1 1 1
#[Mon Nov 17 17:39:35 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-het_loc p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141117_173935EVQ6D.het_loc 1 1 1 1 1
#[Mon Nov 17 18:24:28 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-het_loc p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141117_182427WGQ08.het_loc 1 1 1 1 1
#[Mon Nov 17 18:40:46 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-het_loc p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12878:gm12891:gm19240:gm12864:gm12892 20141117_18404586S25.het_loc 1 1 1 1 1
#[Mon Nov 17 18:41:44 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-het_loc p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12892:gm12878:gm12864:gm19240 20141117_184143SW5MC.het_loc 1 1 1 1 1
#[Mon Nov 17 18:47:06 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-het_loc p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12892:gm12878:gm12864:gm19240 20141117_184705PK7WG.het_loc 1 1 1 1 1
#[Mon Nov 17 19:13:58 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12892:gm12878:gm12864:gm19240 20141117_191357KNXYE.het_loc 1 1 1 1 1
#[Mon Nov 17 19:13:58 2014] p p_run_cluster_sep.py het_loc2-shi-gm12892-all p_het_sites_in_narrow_peak_dp.py gm12892 ctcf gm19238:gm19239:gm12878:gm12891:gm12864:gm19240 20141117_191357HC015.het_loc 1 1 1 1 1
#[Mon Nov 17 19:13:58 2014] p p_run_cluster_sep.py het_loc2-shi-gm19238-all p_het_sites_in_narrow_peak_dp.py gm19238 ctcf gm19239:gm12892:gm12878:gm12891:gm12864:gm19240 20141117_191357K32M7.het_loc 1 1 1 1 1
#[Mon Nov 17 19:13:58 2014] p p_run_cluster_sep.py het_loc2-shi-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm19238:gm19239:gm12892:gm12878:gm12891:gm19240 20141117_1913578QJ2U.het_loc 1 1 1 1 1
#[Mon Nov 17 19:13:59 2014] p p_run_cluster_sep.py het_loc2-shi-gm19239-all p_het_sites_in_narrow_peak_dp.py gm19239 ctcf gm19238:gm12892:gm12878:gm12891:gm12864:gm19240 20141117_19135886ZA6.het_loc 1 1 1 1 1
#[Mon Nov 17 19:13:59 2014] p p_run_cluster_sep.py het_loc2-shi-gm19240-all p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19238:gm19239:gm12878:gm12891:gm12864:gm12892 20141117_191358KP17K.het_loc 1 1 1 1 1
#[Tue Nov 18 09:54:53 2014] p p_run_cluster_sep.py het_loc2-shi-gm12891-all p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm19238:gm19239:gm12892:gm12878:gm12864:gm19240 20141118_095452ESQL4.het_loc 1 1 1 1 1
#[Tue Nov 18 09:55:00 2014] p p_run_cluster_sep.py het_loc2-shi-gm12892-all p_het_sites_in_narrow_peak_dp.py gm12892 ctcf gm19238:gm19239:gm12878:gm12891:gm12864:gm19240 20141118_0954597ATKY.het_loc 1 1 1 1 1
#[Thu Nov 20 13:41:12 2014] p p_run_cluster_sep.py het_loc2-shi-test-all p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141120_134111het_locKDKYF 1 1 1 1 1
#[Thu Nov 20 13:43:41 2014] p p_run_cluster_sep.py het_loc2-shi-test-all p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141120_134340het_locA7Y4Z 1 1 1 1 1
#[Thu Nov 20 13:55:20 2014] p p_run_cluster_sep.py het_loc2-shi-test-all p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141120_135519het_locO8AGA 1 1 1 1 1
#[Thu Nov 20 13:58:27 2014] p p_run_cluster_sep.py het_loc2-shi-test-all p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20141120_135826het_locAP8LF 1 1 1 1 1
#[Tue Nov 25 11:29:01 2014] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864:helas3 20141125_112900.het_loc.2BGAS 1 1 1 1 1
#[Sat Jan 31 14:18:24 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-het_loc p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864:helas3 20150131_141823.het_loc.Q9R2H 1 1 1 1 1
#[Sat Jan 31 14:32:04 2015] p p_run_cluster_sep.py het_loc2-shi-test-het_loc p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20150131_143203.het_loc.YB7VN 1 1 1 1 1
#[Sat Jan 31 14:33:20 2015] p p_run_cluster_sep.py het_loc2-shi-test-het_loc p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20150131_143319.het_loc.2XVI5 1 1 1 1 1
#[Sat Jan 31 14:47:04 2015] p p_run_cluster_sep.py het_loc2-shi-test-het_loc p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20150131_144703.het_loc.VCDQ7 1 1 1 1 1
#[Sun Feb  1 11:21:03 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150201_112102.het_loc.5Q6WZ 1 1 1 1 1
#[Tue Feb  3 11:03:44 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150203_110342.het_loc.BTIY1 1 1 1 1 1
#[Thu Feb 12 17:35:04 2015] p p_run_cluster_sep.py het_loc2-shi-gm12872-all p_het_sites_in_narrow_peak_dp.py gm12872 ctcf gm12864 20150212_173503.het_loc.9DG8C 1 1 1 1 1
#[Thu Feb 12 17:37:43 2015] p p_run_cluster_sep.py het_loc2-shi-gm12872-all p_het_sites_in_narrow_peak_dp.py gm12872 ctcf gm12864 20150212_173742.het_loc.TZ3CW 1 1 1 1 1
#[Wed Feb 18 11:50:04 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150218_115003.het_loc.PD92W 1 1 1 1 1
#[Wed Feb 18 13:45:53 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 znf143 20150218_134552.het_loc.7Y6VM 1 1 1 1 1 1
#[Wed Feb 18 14:26:55 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 znf143 none 20150218_142655.het_loc.669O4 1 1 1 1 1
#[Wed Feb 18 14:30:28 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 znf143 none 20150218_143027.het_loc.TVTNS 1 1 1 1 1
#[Wed Feb 18 14:40:52 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 znf143 none 20150218_144051.het_loc.29U43 1 1 1 1 1
#[Sat Feb 21 10:32:23 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 rad21 none 20150221_103222.het_loc.3TZUH 1 1 1 1 1
#[Sat Feb 21 12:01:23 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 smc3 none 20150221_120122.het_loc.N8YO4 1 1 1 1 1
#[Wed Mar  4 16:25:59 2015] p p_run_cluster_sep.py het_loc2-shi-gm12878-all p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40:ebf1 none 20150304_162557.het_loc.HXKVP 1 1 1 1 1
#[Sat Mar  7 11:10:20 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872:gm12873:helas3:gm12864 20150307_111017.het_loc.4LV0C 1 1 1 1 1
#[Mon Mar  9 14:29:09 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872:gm12873:helas3:gm12864 20150309_142907.het_loc.FJN19 1 1 1 1 1
#[Mon Mar  9 15:17:56 2015] p p_run_cluster_sep.py het_loc2-shi-gm12873-all p_het_sites_in_narrow_peak_dp.py gm12873 ctcf gm12864 20150309_151754.het_loc.68O5Y 1 1 1 1 1
#[Mon Mar  9 16:36:21 2015] p p_run_cluster_sep.py het_loc2-shi-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf 20150309_163620.het_loc.VEUE2 1 1 1 1 1 1
#[Mon Mar  9 16:53:55 2015] p p_run_cluster_sep.py het_loc2-shi-gm12864-all p_het_sites_in_narrow_peak_dp.py gm12864 ctcf None 20150309_165354.het_loc.SAYEZ 1 1 1 1 1
#[Tue Mar 10 15:40:38 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-het_loc p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150310_154037.het_loc.MU2JP 1 1 1 1 1
#[Tue Mar 10 21:10:37 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150310_211036.het_loc.0GPMV 1 1 1 1 1
#[Tue Mar 10 21:10:54 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12873 20150310_211054.het_loc.1WUCC 1 1 1 1 1
#[Tue Mar 10 21:11:19 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm19238 20150310_211118.het_loc.VHYL6 1 1 1 1 1
#[Tue Mar 10 21:11:29 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm19239 20150310_211128.het_loc.7JWQN 1 1 1 1 1
#[Tue Mar 10 21:11:36 2015] p p_run_cluster_sep.py het_loc2-shi-ctcf-all p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm19240 20150310_211135.het_loc.MTPEI 1 1 1 1 1
#[Wed Mar 11 10:20:19 2015] p p_run_cluster_sep.py het_loc2-gm12878-all-allq p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40:ebf1 none 20150311_102018.het_loc.FVGK5 1 1 1 1 1
#[Wed Mar 11 10:31:02 2015] p p_run_cluster_sep.py het_loc2-ctcf-all-allq p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864 20150311_103101.het_loc.8NEUR 1 1 1 1 1
#[Wed Mar 11 11:01:27 2015] p p_run_cluster_sep.py het_loc2-ctcf-all-allq p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12873 20150311_110126.het_loc.4I3CC 1 1 1 1 1
#[Wed Mar 11 20:47:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-het_loc-run p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40:ebf1 none 20150311_204730.het_loc.VKG9M 1 1 1 1 1
#[Wed Mar 11 20:57:33 2015] p p_run_cluster_sep.py het_loc2-ctcf-het_loc-run p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150311_205732.het_loc.JA5YT 1 1 1 1 1
#[Wed Mar 11 21:09:50 2015] p p_run_cluster_sep.py het_loc2-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150311_210950.het_loc.WDO5E 1 1 1 1 1
#[Wed Mar 11 21:09:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12873 20150311_210958.het_loc.HBX5F 1 1 1 1 1
#[Wed Mar 11 21:10:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864 20150311_211006.het_loc.XSWD3 1 1 1 1 1
#[Fri Mar 13 10:01:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12872-het_loc-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150313_100139.het_loc.R5E0O 1 1 1 1 1
#[Fri Mar 13 10:19:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12872-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12872 20150313_101942.het_loc.V2ISU 1 1 1 1 1
#[Fri Mar 13 15:31:19 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12873-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12873 20150313_153118.het_loc.KIC61 1 1 1 1 1
#[Fri Mar 13 15:31:30 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12864-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864 20150313_153129.het_loc.QGHPD 1 1 1 1 1
#[Fri Mar 13 22:29:23 2015] p p_run_cluster_sep.py het_loc2-gm12873-gm12873-all-shi p_het_sites_in_narrow_peak_dp.py gm12873 ctcf gm12873 20150313_222923.het_loc.0F782 1 1 1 1 1
#[Sun Mar 15 10:43:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12864-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12864 20150315_104315.het_loc.DW1LZ 1 1 1 1 1
#[Mon Mar 16 10:54:42 2015] p p_run_cluster_sep.py het_loc2-gm12872-gm12872-all-shi p_het_sites_in_narrow_peak_dp.py gm12872 ctcf gm12864 20150316_105440.het_loc.FUAGC 1 1 1 1 1
#[Mon Mar 16 10:55:20 2015] p p_run_cluster_sep.py het_loc2-gm12864-gm12864-all-shi p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm12864 20150316_105519.het_loc.4F47Q 1 1 1 1 1
#[Tue Mar 17 09:55:15 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12801-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12801 20150317_095514.het_loc.U22QL 1 1 1 1 1
#[Fri Mar 20 10:01:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12878 20150320_100141.het_loc.3RTGY 1 1 1 1 1
#[Fri Mar 20 10:03:25 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12878 20150320_100324.het_loc.OEDRD 1 1 1 1 1
#[Fri Mar 20 10:03:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12891 20150320_100346.het_loc.38BJX 1 1 1 1 1
#[Fri Mar 20 10:03:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12892 20150320_100354.het_loc.XK9WY 1 1 1 1 1
#[Fri Mar 20 10:04:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm19238 20150320_100420.het_loc.CN06R 1 1 1 1 1
#[Sat Mar 21 11:08:29 2015] p p_run_cluster_sep.py het_loc2-gm12801-gm12801-all-shi p_het_sites_in_narrow_peak_dp.py gm12801 ctcf gm12801 20150321_110828.het_loc.P1G3S 1 1 1 1 1
#[Sat Mar 21 11:08:54 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-all-shi p_het_sites_in_narrow_peak_dp.py gm19238 ctcf gm19238 20150321_110853.het_loc.DWRL5 1 1 1 1 1
#[Sat Mar 21 11:09:05 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-all-shi p_het_sites_in_narrow_peak_dp.py gm19239 ctcf gm19239 20150321_110904.het_loc.STML8 1 1 1 1 1
#[Sat Mar 21 12:30:29 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-all-shi p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19240 20150321_123028.het_loc.D3B59 1 1 1 1 1
#[Tue Mar 24 11:50:18 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm19239 20150324_115017.het_loc.9AERE 1 1 1 1 1
#[Tue Mar 24 21:51:02 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-all-shi p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19240 20150324_215101.het_loc.119KF 1 1 1 1 1
#[Tue Mar 24 21:51:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm19240 20150324_215115.het_loc.EYGP8 1 1 1 1 1
#[Thu Mar 26 16:41:13 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-all-shi p_het_sites_in_narrow_peak_dp.py gm12891 ctcf gm12891 20150326_164112.het_loc.48H7F 1 1 1 1 1
#[Thu Mar 26 16:41:13 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-all-shi p_het_sites_in_narrow_peak_dp.py gm19240 ctcf gm19240 20150326_164112.het_loc.M4AES 1 1 1 1 1
#[Thu Mar 26 16:41:13 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-all-shi p_het_sites_in_narrow_peak_dp.py gm12892 ctcf gm12892 20150326_164112.het_loc.QS4DC 1 1 1 1 1
#[Thu Mar 26 16:41:14 2015] p p_run_cluster_sep.py het_loc2-gm12864-gm12864-all-shi p_het_sites_in_narrow_peak_dp.py gm12864 ctcf gm12864 20150326_164112.het_loc.ACEPR 1 1 1 1 1
#[Thu Mar 26 16:41:14 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-all-shi p_het_sites_in_narrow_peak_dp.py gm19239 ctcf gm19239 20150326_164112.het_loc.LEN2V 1 1 1 1 1
#[Thu Mar 26 16:41:14 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-all-shi p_het_sites_in_narrow_peak_dp.py gm19238 ctcf gm19238 20150326_164112.het_loc.OMVZ2 1 1 1 1 1
#[Sat Mar 28 12:00:14 2015] p p_run_cluster_sep.py het_loc2-gm12872-gm12872-all-shi p_het_sites_in_narrow_peak_dp.py gm12872 ctcf gm12864 20150328_120013.het_loc.JY6I1 1 1 1 1 1
#[Sat Mar 28 12:21:51 2015] p p_run_cluster_sep.py het_loc2-gm12872-gm12872-all-shi p_het_sites_in_narrow_peak_dp.py gm12872 ctcf gm12872 20150328_122150.het_loc.VFMEN 1 1 1 1 1
#[Fri Apr 10 22:10:11 2015] p p_run_cluster_sep.py het_loc2-test-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20150410_221010.het_loc.Z428D 1 1 1 1 1
#[Fri Apr 10 23:19:49 2015] p p_run_cluster_sep.py het_loc2-test-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20150410_231949.het_loc.D5OV6 1 1 1 1 1
#[Sat Apr 11 10:30:52 2015] p p_run_cluster_sep.py het_loc2-test-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20150411_103051.het_loc.2E47I 1 1 1 1 1
#[Sat Apr 11 11:20:10 2015] p p_run_cluster_sep.py het_loc2-test-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py helas3 test gm12878 20150411_112009.het_loc.TT7H0 1 1 1 1 1
#[Sat Apr 11 11:21:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 PU1 gm12878 20150411_112104.het_loc.CY0SQ 1 1 1 1 1
#[Sat Apr 11 11:29:49 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 PU1 gm12878 20150411_112948.het_loc.GNFX2 1 1 1 1 1
#[Sat Apr 11 11:31:43 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 PU1 gm12878 20150411_113142.het_loc.STKC2 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 11 15:26:39 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150411_152638.het_loc.T7TNH sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 12 22:22:15 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150412_222214.het_loc.8FXV1 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 12 22:32:39 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150412_223238.het_loc.8WL0G sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 12 22:32:49 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11831-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11831 20150412_223248.het_loc.I5TD3 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 12 22:33:09 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19238 20150412_223308.het_loc.246ZV sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 12 22:33:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150412_223315.het_loc.TXM25 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 12 22:33:24 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150412_223323.het_loc.IXLIS sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 11:09:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-all-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_110926.het_loc.GCJFI sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 12:16:25 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150413_121624.het_loc.UI7UF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 14:02:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150413_140201.het_loc.DZUHZ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11881-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11881 20150413_155659.het_loc.1I1FU sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11831 20150413_155659.het_loc.GYG6V sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150413_155700.het_loc.O1SAL sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11894 20150413_155659.het_loc.0557G sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_155700.het_loc.TV5I7 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11840-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11840 20150413_155700.het_loc.GYE6C sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150413_155700.het_loc.KZWCQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12043-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12043 20150413_155700.het_loc.O1L6T sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12813 20150413_155700.het_loc.4YRY9 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12891 20150413_155700.het_loc.OPADF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12776 20150413_155700.het_loc.D70QJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 15:57:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19238 20150413_155700.het_loc.201RH sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:19:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_161903.het_loc.9GB8S sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:28:12 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12892 20150413_162811.het_loc.83E6F sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:28:45 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_162844.het_loc.2CC1R sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11831 20150413_164056.het_loc.4F4H1 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11840-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11840 20150413_164056.het_loc.N1CD0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12891 20150413_164056.het_loc.HT0EV sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150413_164056.het_loc.837KR sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_164056.het_loc.C78WU sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12043-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12043 20150413_164057.het_loc.K38ZO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12776 20150413_164057.het_loc.3XARZ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11894 20150413_164057.het_loc.WI4OK sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:40:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12813 20150413_164057.het_loc.NA259 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:41:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11881-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11881 20150413_164057.het_loc.OYYDS sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:41:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150413_164057.het_loc.P2DQT sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:41:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19238 20150413_164057.het_loc.BWXS0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:41:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12892 20150413_164057.het_loc.RGY75 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19238 20150413_170455.het_loc.OUC66 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150413_170455.het_loc.CUBR4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11881-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11881 20150413_170455.het_loc.I7DNW sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11840-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11840 20150413_170455.het_loc.U8H5E sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11894 20150413_170455.het_loc.QE8AB sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12043-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12043 20150413_170455.het_loc.DXQ0Q sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11831 20150413_170455.het_loc.MTPV1 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12892 20150413_170455.het_loc.G5O6J sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_170455.het_loc.6DIIP sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12776 20150413_170455.het_loc.8BSUN sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12813 20150413_170455.het_loc.T8DCQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12891 20150413_170455.het_loc.W5C4U sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:04:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150413_170455.het_loc.CZ2YZ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 21:10:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150413_210959.het_loc.3SUCC sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 21:10:20 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 ctcf gm11830 20150413_211019.het_loc.NO2JJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 21:14:45 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150413_211444.het_loc.6KTJ1 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 21:18:56 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150413_211855.het_loc.JIGWD sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:03 2015] p p_run_cluster_sep.py het_loc2-gm11840-gm11840-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11840 pu1 gm11840 20150413_220502.het_loc.CF49U sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:03 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19238 20150413_220502.het_loc.ECW34 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150413_220502.het_loc.1RAW8 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11831 20150413_220502.het_loc.9NY8A sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm11881-gm11881-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11881 pu1 gm11881 20150413_220503.het_loc.J9UV8 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150413_220503.het_loc.YGTXA sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11840-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11840 20150413_220502.het_loc.ISLP1 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150413_220503.het_loc.DRKR5 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150413_220503.het_loc.YNNGL sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150413_220503.het_loc.N7PZ3 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm12043-gm12043-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12043 pu1 gm12043 20150413_220503.het_loc.2QZ37 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12813 20150413_220503.het_loc.G9P39 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12891 20150413_220503.het_loc.OX8DB sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:05 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150413_220502.het_loc.G4W3B sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12043-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12043 20150413_220503.het_loc.D6AHM sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150413_220503.het_loc.MXJ7G sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:05 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150413_220502.het_loc.66NIO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_220502.het_loc.GG0LQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150413_220502.het_loc.BHOOE sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150413_220503.het_loc.QWPT7 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11881-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11881 20150413_220503.het_loc.8VLJP sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12892 20150413_220503.het_loc.D0FC7 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150413_220503.het_loc.9TTYR sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12776 20150413_220503.het_loc.I2K0B sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11894 20150413_220503.het_loc.ZYPH9 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:05:07 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150413_220503.het_loc.0OHII sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:34:56 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150413_223454.het_loc.G5VZJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:34:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150413_223455.het_loc.VWDOQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:34:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12892 20150413_223455.het_loc.UT1NP sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12043-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12043 20150413_223504.het_loc.E5PJ2 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:05 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150413_223504.het_loc.KXJGF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:05 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150413_223504.het_loc.PGOT4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11894 20150413_223504.het_loc.5E979 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150413_223505.het_loc.KX9JE sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12891 20150413_223504.het_loc.TPOJP sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11831 20150413_223504.het_loc.S6R5K sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11881-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11881 20150413_223504.het_loc.11MBI sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150413_223504.het_loc.1G3GF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12776 20150413_223504.het_loc.UATXU sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12813 20150413_223505.het_loc.FXFSO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150413_223504.het_loc.WD5YA sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150413_223504.het_loc.A39S6 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150413_223504.het_loc.3ANJW sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19238 20150413_223504.het_loc.156F9 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150413_223504.het_loc.YW299 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:09 2015] p p_run_cluster_sep.py het_loc2-gm12043-gm12043-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12043 pu1 gm12043 20150413_223504.het_loc.X38WH sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:09 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150413_223504.het_loc.9X037 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:09 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150413_223504.het_loc.9TKWN sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:09 2015] p p_run_cluster_sep.py het_loc2-gm11840-gm11840-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11840 pu1 gm11840 20150413_223504.het_loc.OJF0G sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:09 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150413_223504.het_loc.IWBBC sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:09 2015] p p_run_cluster_sep.py het_loc2-gm11881-gm11881-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11881 pu1 gm11881 20150413_223504.het_loc.XDJ9M sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:35:09 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11840-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11840 20150413_223504.het_loc.RI81O sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 10:10:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11894 20150414_101018.het_loc.CO3YC sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 10:11:52 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150414_101151.het_loc.0HCRZ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 11:36:56 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150414_113655.het_loc.N0XEP sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:02:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150414_170251.het_loc.AXM33 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:24:03 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150414_172403.het_loc.GZ135 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:25:59 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150414_172558.het_loc.RH5QK sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:07 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150414_173005.het_loc.6JN8X sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:07 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150414_173006.het_loc.KGFY6 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150414_173006.het_loc.YSVNX sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150414_173006.het_loc.UULD8 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150414_173006.het_loc.3NJ1B sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150414_173006.het_loc.N02YN sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150414_173006.het_loc.4NAJK sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150414_173006.het_loc.Q5INB sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150414_173006.het_loc.U674V sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:08 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150414_173006.het_loc.4VIGG sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:19 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150414_173018.het_loc.KOMR1 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:00:05 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150414_220004.het_loc.FEA27 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:30:45 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150414_223045.het_loc.VLJUW sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:34 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150414_223233.het_loc.QBIQA sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150414_223233.het_loc.4ZF7S sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150414_223233.het_loc.L5A9J sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150414_223233.het_loc.SPEAA sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150414_223233.het_loc.VUP61 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150414_223233.het_loc.6QAYV sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150414_223233.het_loc.3KDRQ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150414_223233.het_loc.NWHHF sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150414_223233.het_loc.9ET1S sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:35 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150414_223233.het_loc.7TMMU sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:32:36 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150414_223234.het_loc.M2F2R sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:45 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150415_161143.het_loc.Z7CIA sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:45 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150415_161143.het_loc.Q7N5K sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:45 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150415_161143.het_loc.Q7TUJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:45 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150415_161143.het_loc.9FUEW sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:45 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150415_161143.het_loc.37SGF sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:45 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150415_161144.het_loc.BAP8T sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:45 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150415_161144.het_loc.NYCZ6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:11:51 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_161150.het_loc.HB2RP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:16:01 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_141559.het_loc.PO2Q2 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:13 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150415_162312.het_loc.OQ3J3 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:13 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150415_162312.het_loc.1AFVA sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:13 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150415_162312.het_loc.52MOP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:13 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150415_162312.het_loc.9K7TX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:13 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150415_162312.het_loc.LFTT4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:13 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150415_162312.het_loc.ULQEZ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:13 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150415_162312.het_loc.TC9Q5 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:23:14 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_162313.het_loc.2PZ5S sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:16:39 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_171638.het_loc.0L3W6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150415_173142.het_loc.MBLWD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150415_173142.het_loc.2AYZ5 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150415_173142.het_loc.024UJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150415_173142.het_loc.SKTD3 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150415_173142.het_loc.4TTZM sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150415_173142.het_loc.9DPWC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150415_173143.het_loc.8XPWP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150415_173142.het_loc.L2RMS sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150415_173142.het_loc.FGVOO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:31:43 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_173143.het_loc.WCT0O sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:32:45 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_173245.het_loc.TB32H sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:33:05 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150415_172905.het_loc.9Z2O6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:33:05 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150415_172904.het_loc.0OG5V sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:33:05 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150415_172904.het_loc.NV44B sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:33:05 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150415_172904.het_loc.KJF37 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:33:46 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_173346.het_loc.IPOLQ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:49:24 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_174923.het_loc.63T4L sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:50 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150415_172905.het_loc.9WTTK sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:50 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150415_173142.het_loc.8SA3F sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:50 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150415_172904.het_loc.B7Z75 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:50 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150415_172904.het_loc.WBJ3U sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:50 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150415_172904.het_loc.V0QR8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:51 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_174041.het_loc.JPW31 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150415_172904.het_loc.86RRT sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:51 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_172905.het_loc.IJ1F9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:51 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150415_172904.het_loc.NMHLU sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:51:51 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_173935.het_loc.DZ99I sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150415_210713.het_loc.IS9A9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150415_210713.het_loc.SYQ24 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150415_210713.het_loc.QVCI5 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150415_210713.het_loc.F1G38 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150415_210713.het_loc.S75KE sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150415_210713.het_loc.QD9RX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150415_210713.het_loc.2GD3P sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:15 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150415_210713.het_loc.V7V90 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:16 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150415_210714.het_loc.BCHJB sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:16 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150415_210714.het_loc.GIQFL sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:07:16 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150415_210714.het_loc.R60G4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:29:39 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150415_212938.het_loc.C27KG sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:29:40 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-myc-shi p_het_sites_in_narrow_peak_dp.py gm12892 myc gm12892 20150415_212939.het_loc.D83QW sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:29:40 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-myc-shi p_het_sites_in_narrow_peak_dp.py gm12891 myc gm12891 20150415_212939.het_loc.NC147 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:29:40 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-myc-shi p_het_sites_in_narrow_peak_dp.py gm19239 myc gm19239 20150415_212939.het_loc.5IIHC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:29:40 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-myc-shi p_het_sites_in_narrow_peak_dp.py gm19238 myc gm19238 20150415_212939.het_loc.FHD8B sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:29:40 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-myc-shi p_het_sites_in_narrow_peak_dp.py gm19240 myc gm19240 20150415_212939.het_loc.JDFLD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:45:09 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-myc-shi p_het_sites_in_narrow_peak_dp.py gm12891 myc gm12891 20150415_214509.het_loc.K8NND sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:47:27 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-myc-shi p_het_sites_in_narrow_peak_dp.py gm12891 myc gm12891 20150415_214727.het_loc.UJ3ET sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:47:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150415_214739.het_loc.H7G0G sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:14 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150418_183412.het_loc.OYXL0 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150418_183412.het_loc.2VVAN sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150418_183412.het_loc.FKWCX sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150418_183412.het_loc.QWIOB sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150418_183412.het_loc.WU33N sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150418_183413.het_loc.JUH99 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150418_183413.het_loc.RJI4L sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150418_183412.het_loc.FJLQP sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150418_183412.het_loc.8UU31 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:15 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150418_183412.het_loc.I3PK3 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:16 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150418_183413.het_loc.AC20F sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:37 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150418_183435.het_loc.FG1ID sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150418_183435.het_loc.9F2QL sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:37 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150418_183435.het_loc.JQLNL sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:37 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150418_183435.het_loc.8UR6J sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:37 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150418_183435.het_loc.NBNC3 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:38 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150418_183435.het_loc.X9EBZ sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:38 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150418_183435.het_loc.VW8GO sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:38 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150418_183435.het_loc.JNBG1 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:38 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150418_183435.het_loc.VGPY8 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:38 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150418_183436.het_loc.O5CU7 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:34:39 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150418_183436.het_loc.LJPOO sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:31 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm11830 20150418_183529.het_loc.0BRLA sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150418_183529.het_loc.747G9 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:31 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12776 20150418_183529.het_loc.CRZWY sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:31 2015] p p_run_cluster_sep.py het_loc2-gm11831-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11831 pu1 gm11831 20150418_183530.het_loc.TTMFX sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:32 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19239 pu1 gm19239 20150418_183530.het_loc.7ZYI5 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:32 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12891 20150418_183530.het_loc.S6911 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:32 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm11894 20150418_183530.het_loc.NYJCL sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:32 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150418_183530.het_loc.MRJ27 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:32 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19238 pu1 gm19238 20150418_183530.het_loc.ES9WI sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:33 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12892 20150418_183530.het_loc.2BCF6 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:33 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm19240 pu1 gm19240 20150418_183530.het_loc.UBM2Z sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11830 20150418_183544.het_loc.LO0OT sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12892 20150418_183545.het_loc.B0C1S sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11831 20150418_183545.het_loc.BNH24 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12776 20150418_183545.het_loc.CGUDG sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm11894 20150418_183545.het_loc.L6IYS sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12891 20150418_183545.het_loc.CFVFR sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19239 20150418_183545.het_loc.39NRY sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19238 20150418_183545.het_loc.0CQEE sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150418_183546.het_loc.8DKGP sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12813 20150418_183546.het_loc.4R52T sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:35:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm19240 20150418_183546.het_loc.YQCG6 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:51:46 2015] p p_run_cluster_sep.py het_loc2-gm19238-gm19238-myc-shi p_het_sites_in_narrow_peak_dp.py gm19238 myc gm19238 20150420_095145.het_loc.XTIUN sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:51:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150420_095145.het_loc.WTJK6 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:51:46 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12891-myc-shi p_het_sites_in_narrow_peak_dp.py gm12891 myc gm12891 20150420_095145.het_loc.PBZPV sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:51:46 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12892-myc-shi p_het_sites_in_narrow_peak_dp.py gm12892 myc gm12892 20150420_095145.het_loc.QRFQJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:51:46 2015] p p_run_cluster_sep.py het_loc2-gm19239-gm19239-myc-shi p_het_sites_in_narrow_peak_dp.py gm19239 myc gm19239 20150420_095145.het_loc.AKQZL sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:51:47 2015] p p_run_cluster_sep.py het_loc2-gm19240-gm19240-myc-shi p_het_sites_in_narrow_peak_dp.py gm19240 myc gm19240 20150420_095146.het_loc.SZQBR sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:52:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12891-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12891 20150420_095219.het_loc.D6ATV sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:52:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12892-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12892 20150420_095219.het_loc.CYPWJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:52:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm19238 20150420_095219.het_loc.3UAD2 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:52:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150420_095219.het_loc.2NPMZ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:52:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19239-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm19239 20150420_095219.het_loc.MH2YA sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 09:52:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm19240 20150420_095220.het_loc.NXPR8 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 11:57:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19238-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm19238 20150420_115747.het_loc.ORRU0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 20 11:57:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm19240-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm19240 20150420_115747.het_loc.INALS sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:39:30 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm11830 20150426_223929.het_loc.01P20 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:39:35 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm11831 20150426_223934.het_loc.335D6 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:39:40 2015] p p_run_cluster_sep.py het_loc2-gm11840-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11840 pu1 gm11894 20150426_223939.het_loc.D8IO3 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:39:45 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm12776 20150426_223944.het_loc.2GLJF sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:39:50 2015] p p_run_cluster_sep.py het_loc2-gm11830-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11830 pu1 gm12813 20150426_223949.het_loc.Y2ZU2 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:39:55 2015] p p_run_cluster_sep.py het_loc2-gm12776-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12776 pu1 gm12878 20150426_223954.het_loc.ZR594 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:40:00 2015] p p_run_cluster_sep.py het_loc2-gm12892-gm12891-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12892 pu1 gm12891 20150426_223959.het_loc.YB4MK sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:40:05 2015] p p_run_cluster_sep.py het_loc2-gm12891-gm12892-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12891 pu1 gm12892 20150426_224004.het_loc.YM3KV sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:40:10 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm19238-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm19238 20150426_224009.het_loc.MMBW6 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:40:15 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm19239-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm19239 20150426_224014.het_loc.B0N2A sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Apr 26 22:40:27 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm19240-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm19240 20150426_224026.het_loc.HL4QK sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 27 07:30:28 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm11830-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm11830 20150427_073027.het_loc.GJISE sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 27 07:47:15 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm12813-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm12813 20150427_074714.het_loc.BSJBO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 27 08:00:48 2015] p p_run_cluster_sep.py het_loc2-gm12813-gm11831-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12813 pu1 gm11831 20150427_080047.het_loc.ETXFR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu Apr 30 10:00:59 2015] p p_run_cluster_sep.py het_loc2-gm11894-gm12776-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11894 pu1 gm12776 20150430_100058.het_loc.Z1N7T sydh:uw:haib:uta:embl 1 1 1 1
#[Mon May 11 11:18:34 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf143-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf143 gm12878 20150511_111833.het_loc.6MUJE sydh:uw:haib:uta:embl 1 1 1 1
#[Mon May 11 12:14:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ebf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ebf1 gm12878 20150511_121426.het_loc.MHI3V sydh:uw:haib:uta:embl 1 1 1 1
#[Mon May 11 12:15:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bhlhe40-shi p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40 gm12878 20150511_121501.het_loc.U4M5G sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 13 21:36:08 2015] p p_run_cluster_sep.py het_loc2-gm11840-gm11894-pu1-shi p_het_sites_in_narrow_peak_dp.py gm11840 pu1 gm11894 20150513_213607.het_loc.TONIR sydh:uw:haib:uta:embl 1 1 1 1
#[Fri May 15 18:30:34 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150515_183033.het_loc.RTSO8 sydh:uw:haib:uta:embl 1 1 1 1
#[Fri May 15 18:30:39 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-egr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 egr1 gm12878 20150515_183038.het_loc.9BAR4 sydh:uw:haib:uta:embl 1 1 1 1
#[Fri May 15 18:30:44 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat5a-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat5a gm12878 20150515_183043.het_loc.YHK65 sydh:uw:haib:uta:embl 1 1 1 1
#[Fri May 15 18:30:49 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150515_183048.het_loc.Q3CN8 sydh:uw:haib:uta:embl 1 1 1 1
#[Fri May 15 18:30:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cdp-shi p_het_sites_in_narrow_peak_dp.py gm12878 cdp gm12878 20150515_183053.het_loc.Q57YI sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 10:15:34 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150516_101532.het_loc.0NMTF sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 10:15:38 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-egr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 egr1 gm12878 20150516_101537.het_loc.EJ7SQ sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 10:15:43 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat5a-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat5a gm12878 20150516_101542.het_loc.J75BD sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 10:15:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150516_101547.het_loc.JFUEL sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 10:15:55 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cdp-shi p_het_sites_in_narrow_peak_dp.py gm12878 cdp gm12878 20150516_101554.het_loc.GFC7Y sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 11:21:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150516_112129.het_loc.X2FA1 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 11:21:35 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-egr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 egr1 gm12878 20150516_112134.het_loc.S0IA2 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 11:21:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat5a-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat5a gm12878 20150516_112139.het_loc.5VPFA sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 11:21:45 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150516_112144.het_loc.B22C7 sydh:uw:haib:uta:embl 1 1 1 1
#[Sat May 16 11:21:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cdp-shi p_het_sites_in_narrow_peak_dp.py gm12878 cdp gm12878 20150516_112151.het_loc.YCLUR sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 10:35:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cdp-shi p_het_sites_in_narrow_peak_dp.py gm12878 cdp gm12878 20150519_103520.het_loc.BWRJ5 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 11:20:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cdp-shi p_het_sites_in_narrow_peak_dp.py gm12878 cdp gm12878 20150519_112053.het_loc.ARAE6 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 14:00:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150519_140026.het_loc.AFDIR sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 14:00:32 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd1-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd1 gm12878 20150519_140031.het_loc.LNQXD sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 14:00:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cebpb-shi p_het_sites_in_narrow_peak_dp.py gm12878 cebpb gm12878 20150519_140036.het_loc.LF45T sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 14:00:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf143-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf143 gm12878 20150519_140041.het_loc.HE1AM sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 14:00:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150519_140046.het_loc.1II5A sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 14:00:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-corest-shi p_het_sites_in_narrow_peak_dp.py gm12878 corest gm12878 20150519_140051.het_loc.FBJOC sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:22:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150519_152240.het_loc.OXQ6S sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:22:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd1-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd1 gm12878 20150519_152245.het_loc.FH73X sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:22:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cebpb-shi p_het_sites_in_narrow_peak_dp.py gm12878 cebpb gm12878 20150519_152250.het_loc.R1ZW3 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:22:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf143-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf143 gm12878 20150519_152255.het_loc.XCWS2 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150519_152300.het_loc.9RS1V sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-corest-shi p_het_sites_in_narrow_peak_dp.py gm12878 corest gm12878 20150519_152305.het_loc.59KQY sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:11 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-batf-shi p_het_sites_in_narrow_peak_dp.py gm12878 batf gm12878 20150519_152310.het_loc.FROC4 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bhlhe40-shi p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40 gm12878 20150519_152315.het_loc.NGIBS sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rfx5-shi p_het_sites_in_narrow_peak_dp.py gm12878 rfx5 gm12878 20150519_152320.het_loc.L5E1R sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf274-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf274 gm12878 20150519_152325.het_loc.78YJ8 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl11a-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl11a gm12878 20150519_152330.het_loc.YR801 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:36 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rxra-shi p_het_sites_in_narrow_peak_dp.py gm12878 rxra gm12878 20150519_152335.het_loc.QDGVH sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150519_152340.het_loc.AQWD3 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat5a-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat5a gm12878 20150519_152345.het_loc.M10K6 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ets1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ets1 gm12878 20150519_152350.het_loc.CN9GP sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:23:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pml-shi p_het_sites_in_narrow_peak_dp.py gm12878 pml gm12878 20150519_152355.het_loc.OWE03 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol3 gm12878 20150519_152400.het_loc.9ZMJ0 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5c20-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5c20 gm12878 20150519_152405.het_loc.UT7HS sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:11 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ikzf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ikzf1 gm12878 20150519_152410.het_loc.7HVUB sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150519_152415.het_loc.QOBBT sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tbp-shi p_het_sites_in_narrow_peak_dp.py gm12878 tbp gm12878 20150519_152420.het_loc.L9WFI sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat3-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat3 gm12878 20150519_152425.het_loc.LBT7R sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf3 gm12878 20150519_152430.het_loc.77WC6 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:36 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rad21-shi p_het_sites_in_narrow_peak_dp.py gm12878 rad21 gm12878 20150519_152435.het_loc.22F6E sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-egr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 egr1 gm12878 20150519_152440.het_loc.GKS47 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2a-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2a gm12878 20150519_152445.het_loc.80JIJ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2c-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2c gm12878 20150519_152450.het_loc.BI6S9 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:24:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-runx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 runx3 gm12878 20150519_152455.het_loc.H3TCR sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mxi1-shi p_het_sites_in_narrow_peak_dp.py gm12878 mxi1 gm12878 20150519_152500.het_loc.NK7CS sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrsf-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrsf gm12878 20150519_152506.het_loc.DDIE6 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:11 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-whip-shi p_het_sites_in_narrow_peak_dp.py gm12878 whip gm12878 20150519_152510.het_loc.W8TK0 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-six5-shi p_het_sites_in_narrow_peak_dp.py gm12878 six5 gm12878 20150519_152515.het_loc.0J8X7 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-input-shi p_het_sites_in_narrow_peak_dp.py gm12878 input gm12878 20150519_152521.het_loc.DWS0K sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf2 gm12878 20150519_152526.het_loc.5U8MV sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl3-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl3 gm12878 20150519_152531.het_loc.UTV9F sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:36 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tr4-shi p_het_sites_in_narrow_peak_dp.py gm12878 tr4 gm12878 20150519_152536.het_loc.DLEGZ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ctcf-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12878 20150519_152541.het_loc.IG93Y sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sin3a-shi p_het_sites_in_narrow_peak_dp.py gm12878 sin3a gm12878 20150519_152546.het_loc.5DS28 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-brca1-shi p_het_sites_in_narrow_peak_dp.py gm12878 brca1 gm12878 20150519_152551.het_loc.9V2H6 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:25:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5n19-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5n19 gm12878 20150519_152556.het_loc.W9PRJ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-yy1-shi p_het_sites_in_narrow_peak_dp.py gm12878 yy1 gm12878 20150519_152601.het_loc.JY99P sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol2 gm12878 20150519_152606.het_loc.QO361 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:12 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-taf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 taf1 gm12878 20150519_152611.het_loc.U4TWM sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-irf4-shi p_het_sites_in_narrow_peak_dp.py gm12878 irf4 gm12878 20150519_152616.het_loc.V9RKG sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:22 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sp1-shi p_het_sites_in_narrow_peak_dp.py gm12878 sp1 gm12878 20150519_152621.het_loc.KIWLE sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-smc3-shi p_het_sites_in_narrow_peak_dp.py gm12878 smc3 gm12878 20150519_152626.het_loc.H3JTI sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mta3-shi p_het_sites_in_narrow_peak_dp.py gm12878 mta3 gm12878 20150519_152631.het_loc.W5D21 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-jund-shi p_het_sites_in_narrow_peak_dp.py gm12878 jund gm12878 20150519_152636.het_loc.RPXVY sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bclaf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 bclaf1 gm12878 20150519_152641.het_loc.W9A26 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf2 gm12878 20150519_152646.het_loc.55KWD sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cmyc-shi p_het_sites_in_narrow_peak_dp.py gm12878 cmyc gm12878 20150519_152651.het_loc.KQZ7W sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:26:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfya-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfya gm12878 20150519_152656.het_loc.VYZ6W sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elk1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elk1 gm12878 20150519_152701.het_loc.DILED sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-foxm1-shi p_het_sites_in_narrow_peak_dp.py gm12878 foxm1 gm12878 20150519_152706.het_loc.KZBFG sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:13 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfe2-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfe2 gm12878 20150519_152712.het_loc.FCUMT sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf12-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf12 gm12878 20150519_152716.het_loc.CRBC5 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol24h8-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol24h8 gm12878 20150519_152721.het_loc.QM3V1 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-srf-shi p_het_sites_in_narrow_peak_dp.py gm12878 srf gm12878 20150519_152726.het_loc.2K5IS sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-dnase-shi p_het_sites_in_narrow_peak_dp.py gm12878 dnase gm12878 20150519_152731.het_loc.YUACG sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf3 gm12878 20150519_152736.het_loc.HPUGT sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrf1 gm12878 20150519_152741.het_loc.9NOWS sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfic-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfic gm12878 20150519_152746.het_loc.HB3KZ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfatc1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfatc1 gm12878 20150519_152751.het_loc.LOP7K sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:27:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-max-shi p_het_sites_in_narrow_peak_dp.py gm12878 max gm12878 20150519_152756.het_loc.RDCZ4 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-e2f4-shi p_het_sites_in_narrow_peak_dp.py gm12878 e2f4 gm12878 20150519_152801.het_loc.RM2I7 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tblr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 tblr1 gm12878 20150519_152806.het_loc.1ZI4A sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:12 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zzz3-shi p_het_sites_in_narrow_peak_dp.py gm12878 zzz3 gm12878 20150519_152811.het_loc.K84UX sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat1-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat1 gm12878 20150519_152816.het_loc.CEBQR sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:22 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfyb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfyb gm12878 20150519_152821.het_loc.C7C2U sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zeb1-shi p_het_sites_in_narrow_peak_dp.py gm12878 zeb1 gm12878 20150519_152826.het_loc.G2EZB sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:32 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfkb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfkb gm12878 20150519_152831.het_loc.USM8W sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-gabp-shi p_het_sites_in_narrow_peak_dp.py gm12878 gabp gm12878 20150519_152836.het_loc.5716E sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5 gm12878 20150519_152841.het_loc.QNKBC sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cfos-shi p_het_sites_in_narrow_peak_dp.py gm12878 cfos gm12878 20150519_152846.het_loc.XGSSY sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf1 gm12878 20150519_152851.het_loc.H335Y sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:28:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ebf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ebf1 gm12878 20150519_152856.het_loc.XHQHV sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:29:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-maz-shi p_het_sites_in_narrow_peak_dp.py gm12878 maz gm12878 20150519_152901.het_loc.9SZH7 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:29:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elf1 gm12878 20150519_152906.het_loc.J2TPZ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:29:13 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zbtb33-shi p_het_sites_in_narrow_peak_dp.py gm12878 zbtb33 gm12878 20150519_152912.het_loc.DIDVJ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:29:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pou2f2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pou2f2 gm12878 20150519_152916.het_loc.T3DCE sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:29:22 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pbx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pbx3 gm12878 20150519_152922.het_loc.LWKLG sydh:uw:haib:uta:embl 1 1 1 1
#[Tue May 19 15:29:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-p300-shi p_het_sites_in_narrow_peak_dp.py gm12878 p300 gm12878 20150519_152926.het_loc.WM2ZD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:15 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pu1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090114.het_loc.EP1XJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090120.het_loc.A3X6L sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:25 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090125.het_loc.9JKO6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:30 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090129.het_loc.8CDPY sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:35 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090134.het_loc.KICPZ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:40 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090139.het_loc.UXARA sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:45 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090144.het_loc.GD7H8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:50 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-inputstd-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090150.het_loc.PMJJ8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:01:55 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090154.het_loc.LS5F9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:00 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090200.het_loc.6KN5H sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:05 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090204.het_loc.VSV7G sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:10 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090210.het_loc.D7VTQ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:15 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090215.het_loc.ZS6Z2 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090220.het_loc.0LZRT sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:25 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090225.het_loc.337Q3 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:30 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090230.het_loc.GCO4N sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:35 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090235.het_loc.ZZ6HT sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:40 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-dnase-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090240.het_loc.6DKDE sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:45 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090245.het_loc.4DH4J sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:50 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090250.het_loc.O6BAC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:02:55 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tbp-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090255.het_loc.1G6EG sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:00 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-test-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090300.het_loc.MPA6X sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:05 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tr4-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090305.het_loc.5IQFN sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:10 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ctcf-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090310.het_loc.545TR sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:15 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brca1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090315.het_loc.9ANW6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pol2-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090320.het_loc.NDB4F sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:25 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-taf1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090325.het_loc.R02ZD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:30 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-max-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090330.het_loc.9STG6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:35 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-smc3-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090335.het_loc.X8Q01 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:40 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-maz-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090340.het_loc.ICUOV sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:45 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-jund-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090345.het_loc.U0A4Y sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-usf2-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090350.het_loc.G33HF sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:03:55 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfyb-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090355.het_loc.0NAN0 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfya-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090400.het_loc.LGMVJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090405.het_loc.WQYSD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk4-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090410.het_loc.MEIA9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf274-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090415.het_loc.CK6TB sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brg1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090420.het_loc.Y9KS6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cmyc-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090425.het_loc.6B4KZ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:30 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cjun-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090430.het_loc.JNMVK sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:35 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrf1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090435.het_loc.L7DXD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f6-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090440.het_loc.80MTW sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f4-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090445.het_loc.TY0GZ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-stat3-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090450.het_loc.CCI13 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:04:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zzz3-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090455.het_loc.8O0XT sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090500.het_loc.ENGMZ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:05 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-spt20-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090505.het_loc.WL0BI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:10 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gabp-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090510.het_loc.90G0D sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cfos-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090515.het_loc.ZA95P sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ini1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090520.het_loc.BTRGN sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf2-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090525.het_loc.S0ULM sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:31 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ap2gamma-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090530.het_loc.2DX8Q sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf1-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090536.het_loc.WN486 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-irf3-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090540.het_loc.4OGK2 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:05:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-p300-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf gm12878 20150520_090545.het_loc.2JX1C sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:22:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pu1-shi p_het_sites_in_narrow_peak_dp.py helas3 pu1 helas3 20150520_092222.het_loc.YO3F0 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:22:27 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150520_092227.het_loc.FMW4U sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:22:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150520_092232.het_loc.DSOR7 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 09:22:55 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150520_092254.het_loc.CKMXG sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:25:38 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pu1-shi p_het_sites_in_narrow_peak_dp.py helas3 pu1 helas3 20150520_102537.het_loc.B6YQA sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:25:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150520_102541.het_loc.YUKPJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:25:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150520_102546.het_loc.65XCP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:25:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150520_102551.het_loc.HNCQC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:25:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150520_102556.het_loc.FEN8T sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150520_102601.het_loc.HFACL sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150520_102606.het_loc.GN4LD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-inputstd-shi p_het_sites_in_narrow_peak_dp.py helas3 inputstd helas3 20150520_102611.het_loc.LVZWO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150520_102616.het_loc.OBORV sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150520_102621.het_loc.34D3V sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150520_102626.het_loc.P4KQH sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150520_102631.het_loc.L4VBA sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150520_102636.het_loc.9VHV1 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:42 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150520_102641.het_loc.1JW2K sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150520_102646.het_loc.67KWU sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:52 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150520_102651.het_loc.SRNZS sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:26:57 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150520_102656.het_loc.G4U2M sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-dnase-shi p_het_sites_in_narrow_peak_dp.py helas3 dnase helas3 20150520_102701.het_loc.TH9XN sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:07 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150520_102706.het_loc.2DA5R sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150520_102711.het_loc.KLD69 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:17 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tbp-shi p_het_sites_in_narrow_peak_dp.py helas3 tbp helas3 20150520_102716.het_loc.HZHXW sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-test-shi p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20150520_102721.het_loc.VIUW5 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:27 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tr4-shi p_het_sites_in_narrow_peak_dp.py helas3 tr4 helas3 20150520_102726.het_loc.PQW3H sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ctcf-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf helas3 20150520_102731.het_loc.A4KCX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brca1-shi p_het_sites_in_narrow_peak_dp.py helas3 brca1 helas3 20150520_102736.het_loc.XAEI2 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:42 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pol2-shi p_het_sites_in_narrow_peak_dp.py helas3 pol2 helas3 20150520_102741.het_loc.ELYNO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:47 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-taf1-shi p_het_sites_in_narrow_peak_dp.py helas3 taf1 helas3 20150520_102746.het_loc.630G4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:52 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-max-shi p_het_sites_in_narrow_peak_dp.py helas3 max helas3 20150520_102751.het_loc.AQ66F sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:27:57 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-smc3-shi p_het_sites_in_narrow_peak_dp.py helas3 smc3 helas3 20150520_102756.het_loc.PLFVH sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-maz-shi p_het_sites_in_narrow_peak_dp.py helas3 maz helas3 20150520_102801.het_loc.F3WSO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:07 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-jund-shi p_het_sites_in_narrow_peak_dp.py helas3 jund helas3 20150520_102806.het_loc.AK0PK sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-usf2-shi p_het_sites_in_narrow_peak_dp.py helas3 usf2 helas3 20150520_102811.het_loc.8OXBC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:17 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfyb-shi p_het_sites_in_narrow_peak_dp.py helas3 nfyb helas3 20150520_102816.het_loc.GGXAS sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfya-shi p_het_sites_in_narrow_peak_dp.py helas3 nfya helas3 20150520_102821.het_loc.02DC5 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:27 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk1-shi p_het_sites_in_narrow_peak_dp.py helas3 elk1 helas3 20150520_102826.het_loc.Q429N sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk4-shi p_het_sites_in_narrow_peak_dp.py helas3 elk4 helas3 20150520_102831.het_loc.RNU5G sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf274-shi p_het_sites_in_narrow_peak_dp.py helas3 znf274 helas3 20150520_102836.het_loc.EZL25 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:42 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brg1-shi p_het_sites_in_narrow_peak_dp.py helas3 brg1 helas3 20150520_102841.het_loc.VK1YU sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:47 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cmyc-shi p_het_sites_in_narrow_peak_dp.py helas3 cmyc helas3 20150520_102846.het_loc.HN7B5 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:52 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cjun-shi p_het_sites_in_narrow_peak_dp.py helas3 cjun helas3 20150520_102851.het_loc.E74FC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:28:57 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrf1-shi p_het_sites_in_narrow_peak_dp.py helas3 nrf1 helas3 20150520_102856.het_loc.KM62J sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f6-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f6 helas3 20150520_102901.het_loc.RC3BX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:07 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f4-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f4 helas3 20150520_102906.het_loc.EQ9A6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-stat3-shi p_het_sites_in_narrow_peak_dp.py helas3 stat3 helas3 20150520_102911.het_loc.VPTQ8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:17 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zzz3-shi p_het_sites_in_narrow_peak_dp.py helas3 zzz3 helas3 20150520_102916.het_loc.75WI9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f1 helas3 20150520_102921.het_loc.QV9NF sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:27 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-spt20-shi p_het_sites_in_narrow_peak_dp.py helas3 spt20 helas3 20150520_102927.het_loc.UO78J sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gabp-shi p_het_sites_in_narrow_peak_dp.py helas3 gabp helas3 20150520_102931.het_loc.ENKPB sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cfos-shi p_het_sites_in_narrow_peak_dp.py helas3 cfos helas3 20150520_102936.het_loc.8GFML sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:42 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ini1-shi p_het_sites_in_narrow_peak_dp.py helas3 ini1 helas3 20150520_102941.het_loc.PDBJ0 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:48 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf2-shi p_het_sites_in_narrow_peak_dp.py helas3 brf2 helas3 20150520_102947.het_loc.B2UPR sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:52 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ap2gamma-shi p_het_sites_in_narrow_peak_dp.py helas3 ap2gamma helas3 20150520_102951.het_loc.E88C8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:29:57 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf1-shi p_het_sites_in_narrow_peak_dp.py helas3 brf1 helas3 20150520_102956.het_loc.FYY3E sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:30:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-irf3-shi p_het_sites_in_narrow_peak_dp.py helas3 irf3 helas3 20150520_103001.het_loc.LWDX2 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 10:30:07 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-p300-shi p_het_sites_in_narrow_peak_dp.py helas3 p300 helas3 20150520_103007.het_loc.Q0F8N sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 15:04:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150520_150411.het_loc.2VATB sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 15:56:50 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150520_155649.het_loc.ANLSI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:02:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150520_160255.het_loc.RLUEO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:12:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150520_161255.het_loc.BB0OP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:13:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150520_161304.het_loc.7EQA7 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:13 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ap2gamma-shi p_het_sites_in_narrow_peak_dp.py helas3 ap2gamma helas3 20150520_162612.het_loc.IABEC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150520_162617.het_loc.ZP2EI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:23 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150520_162622.het_loc.L7I97 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150520_162627.het_loc.5NECZ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:33 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brca1-shi p_het_sites_in_narrow_peak_dp.py helas3 brca1 helas3 20150520_162632.het_loc.DR8FR sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:38 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf1-shi p_het_sites_in_narrow_peak_dp.py helas3 brf1 helas3 20150520_162637.het_loc.PVDQP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:43 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf2-shi p_het_sites_in_narrow_peak_dp.py helas3 brf2 helas3 20150520_162642.het_loc.YL0F9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:26:48 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brg1-shi p_het_sites_in_narrow_peak_dp.py helas3 brg1 helas3 20150520_162647.het_loc.PQ9AJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:27:31 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ap2gamma-shi p_het_sites_in_narrow_peak_dp.py helas3 ap2gamma helas3 20150520_162730.het_loc.N6QF1 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:27:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150520_162735.het_loc.0QVAX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:27:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150520_162740.het_loc.5KGSS sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:27:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150520_162745.het_loc.KLD70 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:27:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brca1-shi p_het_sites_in_narrow_peak_dp.py helas3 brca1 helas3 20150520_162750.het_loc.D9HIP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:27:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf1-shi p_het_sites_in_narrow_peak_dp.py helas3 brf1 helas3 20150520_162755.het_loc.0M5IF sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf2-shi p_het_sites_in_narrow_peak_dp.py helas3 brf2 helas3 20150520_162800.het_loc.0DBQO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brg1-shi p_het_sites_in_narrow_peak_dp.py helas3 brg1 helas3 20150520_162805.het_loc.PM6DF sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150520_162810.het_loc.FJ8ZV sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cfos-shi p_het_sites_in_narrow_peak_dp.py helas3 cfos helas3 20150520_162815.het_loc.T9WG3 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150520_162820.het_loc.JOUSH sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cjun-shi p_het_sites_in_narrow_peak_dp.py helas3 cjun helas3 20150520_162825.het_loc.YBN0J sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:31 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cmyc-shi p_het_sites_in_narrow_peak_dp.py helas3 cmyc helas3 20150520_162830.het_loc.INX29 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150520_162835.het_loc.ESI6L sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ctcf-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf helas3 20150520_162840.het_loc.BO1V1 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-dnase-shi p_het_sites_in_narrow_peak_dp.py helas3 dnase helas3 20150520_162845.het_loc.K6SEW sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f1 helas3 20150520_162850.het_loc.JOFWQ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:28:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f4-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f4 helas3 20150520_162855.het_loc.8HAQK sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f6-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f6 helas3 20150520_162900.het_loc.OCIW3 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk1-shi p_het_sites_in_narrow_peak_dp.py helas3 elk1 helas3 20150520_162905.het_loc.TEUSC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk4-shi p_het_sites_in_narrow_peak_dp.py helas3 elk4 helas3 20150520_162910.het_loc.5KCZ1 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gabp-shi p_het_sites_in_narrow_peak_dp.py helas3 gabp helas3 20150520_162915.het_loc.K7V9Z sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150520_162920.het_loc.P19SU sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ini1-shi p_het_sites_in_narrow_peak_dp.py helas3 ini1 helas3 20150520_162925.het_loc.PKQZ6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:31 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-inputstd-shi p_het_sites_in_narrow_peak_dp.py helas3 inputstd helas3 20150520_162930.het_loc.0OJ7L sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-irf3-shi p_het_sites_in_narrow_peak_dp.py helas3 irf3 helas3 20150520_162935.het_loc.AMWHH sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-jund-shi p_het_sites_in_narrow_peak_dp.py helas3 jund helas3 20150520_162940.het_loc.LCBI8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150520_162945.het_loc.EI4GB sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-max-shi p_het_sites_in_narrow_peak_dp.py helas3 max helas3 20150520_162950.het_loc.7HQZP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:29:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-maz-shi p_het_sites_in_narrow_peak_dp.py helas3 maz helas3 20150520_162955.het_loc.WV5OU sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150520_163000.het_loc.BE7FI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfya-shi p_het_sites_in_narrow_peak_dp.py helas3 nfya helas3 20150520_163005.het_loc.RF3F4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfyb-shi p_het_sites_in_narrow_peak_dp.py helas3 nfyb helas3 20150520_163010.het_loc.A30QN sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrf1-shi p_het_sites_in_narrow_peak_dp.py helas3 nrf1 helas3 20150520_163016.het_loc.W37IL sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150520_163020.het_loc.JOPSS sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pol2-shi p_het_sites_in_narrow_peak_dp.py helas3 pol2 helas3 20150520_163026.het_loc.BT54D sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150520_163031.het_loc.KCQE5 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150520_163035.het_loc.E65FS sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:43 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150520_163042.het_loc.QWP03 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:47 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150520_163046.het_loc.Z3NQ4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:53 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-smc3-shi p_het_sites_in_narrow_peak_dp.py helas3 smc3 helas3 20150520_163052.het_loc.DREO0 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:30:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-spt20-shi p_het_sites_in_narrow_peak_dp.py helas3 spt20 helas3 20150520_163057.het_loc.H2UDT sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:03 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-stat3-shi p_het_sites_in_narrow_peak_dp.py helas3 stat3 helas3 20150520_163102.het_loc.JFY1N sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:07 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-taf1-shi p_het_sites_in_narrow_peak_dp.py helas3 taf1 helas3 20150520_163106.het_loc.D47DN sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tbp-shi p_het_sites_in_narrow_peak_dp.py helas3 tbp helas3 20150520_163111.het_loc.1LD1Q sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150520_163118.het_loc.FSB64 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-test-shi p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20150520_163118.het_loc.SDBCG sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:23 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tr4-shi p_het_sites_in_narrow_peak_dp.py helas3 tr4 helas3 20150520_163122.het_loc.R83L8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-usf2-shi p_het_sites_in_narrow_peak_dp.py helas3 usf2 helas3 20150520_163127.het_loc.SI4ZE sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:33 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150520_163132.het_loc.VW9WQ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150520_163146.het_loc.091EO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zzz3-shi p_het_sites_in_narrow_peak_dp.py helas3 zzz3 helas3 20150520_163148.het_loc.UZH0D sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 16:31:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf274-shi p_het_sites_in_narrow_peak_dp.py helas3 znf274 helas3 20150520_163148.het_loc.8GULF sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150520_174610.het_loc.WANN6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150520_174613.het_loc.KLEXU sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150520_174618.het_loc.7R2NX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:24 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brca1-shi p_het_sites_in_narrow_peak_dp.py helas3 brca1 helas3 20150520_174623.het_loc.ULWGI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150520_174628.het_loc.CFQB0 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150520_174633.het_loc.FHJ7U sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150520_174638.het_loc.598BQ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ctcf-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf helas3 20150520_174644.het_loc.LU201 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-dnase-shi p_het_sites_in_narrow_peak_dp.py helas3 dnase helas3 20150520_174649.het_loc.JUX3A sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150520_174654.het_loc.PXT7E sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:46:59 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-inputstd-shi p_het_sites_in_narrow_peak_dp.py helas3 inputstd helas3 20150520_174659.het_loc.G0RJV sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150520_174704.het_loc.HQ103 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:09 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-max-shi p_het_sites_in_narrow_peak_dp.py helas3 max helas3 20150520_174709.het_loc.4KKA1 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-maz-shi p_het_sites_in_narrow_peak_dp.py helas3 maz helas3 20150520_174714.het_loc.2RY5X sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150520_174719.het_loc.ODKNC sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:24 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150520_174724.het_loc.JOM68 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pol2-shi p_het_sites_in_narrow_peak_dp.py helas3 pol2 helas3 20150520_174729.het_loc.5VPJI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150520_174734.het_loc.QBL81 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150520_174739.het_loc.D2G51 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150520_174744.het_loc.IRN41 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150520_174749.het_loc.7AJL4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-smc3-shi p_het_sites_in_narrow_peak_dp.py helas3 smc3 helas3 20150520_174754.het_loc.4773C sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:47:59 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-taf1-shi p_het_sites_in_narrow_peak_dp.py helas3 taf1 helas3 20150520_174759.het_loc.NUB1C sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:48:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tbp-shi p_het_sites_in_narrow_peak_dp.py helas3 tbp helas3 20150520_174804.het_loc.XWRRN sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:48:09 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150520_174809.het_loc.WQJUY sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:48:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-test-shi p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20150520_174814.het_loc.CHI98 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:48:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tr4-shi p_het_sites_in_narrow_peak_dp.py helas3 tr4 helas3 20150520_174819.het_loc.S2GVL sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:48:25 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150520_174824.het_loc.9WXDE sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:48:45 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150520_174844.het_loc.QGBF6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:57:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ap2gamma-shi p_het_sites_in_narrow_peak_dp.py helas3 ap2gamma helas3 20150520_175732.het_loc.G84M6 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:57:38 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf1-shi p_het_sites_in_narrow_peak_dp.py helas3 brf1 helas3 20150520_175737.het_loc.I4Q2R sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:57:43 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf2-shi p_het_sites_in_narrow_peak_dp.py helas3 brf2 helas3 20150520_175742.het_loc.MGWXS sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:57:48 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brg1-shi p_het_sites_in_narrow_peak_dp.py helas3 brg1 helas3 20150520_175747.het_loc.SS6MN sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:57:53 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cfos-shi p_het_sites_in_narrow_peak_dp.py helas3 cfos helas3 20150520_175752.het_loc.STP3L sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:57:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cjun-shi p_het_sites_in_narrow_peak_dp.py helas3 cjun helas3 20150520_175757.het_loc.ZEXI3 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:03 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cmyc-shi p_het_sites_in_narrow_peak_dp.py helas3 cmyc helas3 20150520_175802.het_loc.72IAP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:08 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f1 helas3 20150520_175807.het_loc.M2SON sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:13 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f4-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f4 helas3 20150520_175812.het_loc.GLKEP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f6-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f6 helas3 20150520_175817.het_loc.DKPDM sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:23 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk1-shi p_het_sites_in_narrow_peak_dp.py helas3 elk1 helas3 20150520_175822.het_loc.N6ZF4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk4-shi p_het_sites_in_narrow_peak_dp.py helas3 elk4 helas3 20150520_175828.het_loc.KPMJ7 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:33 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gabp-shi p_het_sites_in_narrow_peak_dp.py helas3 gabp helas3 20150520_175833.het_loc.MUO92 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:38 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ini1-shi p_het_sites_in_narrow_peak_dp.py helas3 ini1 helas3 20150520_175838.het_loc.EVZ3P sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:43 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-irf3-shi p_het_sites_in_narrow_peak_dp.py helas3 irf3 helas3 20150520_175843.het_loc.ZTVSI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:48 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-jund-shi p_het_sites_in_narrow_peak_dp.py helas3 jund helas3 20150520_175848.het_loc.Q684W sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:53 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfya-shi p_het_sites_in_narrow_peak_dp.py helas3 nfya helas3 20150520_175853.het_loc.LL0PJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:58:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfyb-shi p_het_sites_in_narrow_peak_dp.py helas3 nfyb helas3 20150520_175858.het_loc.WABIT sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:59:03 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrf1-shi p_het_sites_in_narrow_peak_dp.py helas3 nrf1 helas3 20150520_175903.het_loc.5S08D sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:59:08 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-p300-shi p_het_sites_in_narrow_peak_dp.py helas3 p300 helas3 20150520_175908.het_loc.7F8QA sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:59:13 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-spt20-shi p_het_sites_in_narrow_peak_dp.py helas3 spt20 helas3 20150520_175913.het_loc.GV6A1 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:59:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-stat3-shi p_het_sites_in_narrow_peak_dp.py helas3 stat3 helas3 20150520_175918.het_loc.IFUIK sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:59:23 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-usf2-shi p_het_sites_in_narrow_peak_dp.py helas3 usf2 helas3 20150520_175923.het_loc.K7T4J sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:59:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf274-shi p_het_sites_in_narrow_peak_dp.py helas3 znf274 helas3 20150520_175928.het_loc.A50B4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed May 20 17:59:33 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zzz3-shi p_het_sites_in_narrow_peak_dp.py helas3 zzz3 helas3 20150520_175933.het_loc.D6UVV sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150521_085618.het_loc.SMSHE sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:23 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150521_085623.het_loc.IWH6Q sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150521_085628.het_loc.LA7WU sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:33 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_085633.het_loc.D0ZPX sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-dnase-shi p_het_sites_in_narrow_peak_dp.py helas3 dnase helas3 20150521_085657.het_loc.VVN9B sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150521_085657.het_loc.J9TIU sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150521_085657.het_loc.DF2E5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150521_085657.het_loc.CPW7K sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:56:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-inputstd-shi p_het_sites_in_narrow_peak_dp.py helas3 inputstd helas3 20150521_085658.het_loc.EG9GZ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150521_085703.het_loc.ZFLY6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:08 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150521_085708.het_loc.U5OZS sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:13 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150521_085713.het_loc.LUJGV sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150521_085718.het_loc.B5YVT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:23 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150521_085723.het_loc.TU34G sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150521_085728.het_loc.Z8YYI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:33 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150521_085733.het_loc.VDLRC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150521_085738.het_loc.76IEL sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 08:57:43 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150521_085743.het_loc.68PDA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 09:23:17 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150521_092316.het_loc.YW7UF sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pu1-shi p_het_sites_in_narrow_peak_dp.py helas3 pu1 helas3 20150521_100711.het_loc.JOY87 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:17 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150521_100716.het_loc.WRL53 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150521_100721.het_loc.LJMG6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:27 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150521_100726.het_loc.BKE4N sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_100731.het_loc.8UIX3 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150521_100736.het_loc.9IN3A sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:42 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150521_100741.het_loc.8O8MS sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:47 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150521_100746.het_loc.8Q9O4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:52 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150521_100751.het_loc.YZDT9 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:07:57 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150521_100756.het_loc.I8XJ0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150521_100801.het_loc.KB96R sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:07 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150521_100806.het_loc.871U4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150521_100811.het_loc.6ZDD6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:17 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150521_100816.het_loc.4OHVF sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150521_100821.het_loc.EA2T8 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:27 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150521_100826.het_loc.EOU1M sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150521_100831.het_loc.URLI2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:08:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150521_100836.het_loc.1URXT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:12:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pu1-shi p_het_sites_in_narrow_peak_dp.py helas3 pu1 helas3 20150521_101233.het_loc.BOWJG sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:12:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150521_101238.het_loc.J01AW sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:12:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150521_101243.het_loc.5XE4H sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:12:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150521_101248.het_loc.0Q0PO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:12:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_101253.het_loc.7ED26 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:12:59 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150521_101258.het_loc.XAAIP sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150521_101303.het_loc.930AW sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:09 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150521_101308.het_loc.BF13B sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150521_101313.het_loc.MGXFV sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150521_101318.het_loc.HSNIT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:24 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150521_101323.het_loc.I6BH4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150521_101328.het_loc.HA9GQ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150521_101333.het_loc.CJ3SH sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150521_101338.het_loc.21KAQ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150521_101343.het_loc.26NAR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150521_101348.het_loc.NQU8V sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150521_101353.het_loc.FM4NN sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:13:59 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150521_101358.het_loc.EU7LF sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:18:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pu1-shi p_het_sites_in_narrow_peak_dp.py helas3 pu1 helas3 20150521_101839.het_loc.RYHG7 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:18:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150521_101844.het_loc.P7IQ6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:18:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150521_101849.het_loc.OINOK sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:18:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150521_101854.het_loc.WSX2E sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:18:59 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_101859.het_loc.7DTHB sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150521_101904.het_loc.DQL8W sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:09 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150521_101909.het_loc.28ASA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150521_101914.het_loc.J3DEA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150521_101919.het_loc.NJTNB sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:25 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150521_101924.het_loc.XYIP9 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150521_101929.het_loc.B4Z1B sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150521_101934.het_loc.Y53GK sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150521_101939.het_loc.M45W1 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150521_101944.het_loc.NAE0K sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150521_101949.het_loc.QSDT0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:19:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150521_101954.het_loc.TAFGM sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:20:00 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150521_101959.het_loc.WGO8J sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 10:20:05 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150521_102004.het_loc.9I7AM sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150521_115320.het_loc.OW6ZM sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150521_115325.het_loc.8SHT5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:31 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150521_115330.het_loc.PIZFD sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_115335.het_loc.6Y0WO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150521_115340.het_loc.E1OH5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150521_115345.het_loc.DANIA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150521_115350.het_loc.UEO7P sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:53:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150521_115355.het_loc.QFWQ6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150521_115400.het_loc.6U5JA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150521_115405.het_loc.60U69 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150521_115411.het_loc.Z72UI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150521_115416.het_loc.FYHO2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150521_115421.het_loc.V1VTM sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150521_115426.het_loc.DDT85 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:31 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150521_115431.het_loc.A0G6T sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150521_115436.het_loc.61XLD sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:54:45 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150521_115444.het_loc.ABM3R sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:59:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150521_115905.het_loc.0M2G5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:59:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150521_115910.het_loc.XY80M sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:59:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150521_115915.het_loc.VHXDI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:59:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_115920.het_loc.3QR92 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 11:59:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150521_115925.het_loc.G924T sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:10:00 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_120959.het_loc.WLQVG sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:13:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_121336.het_loc.JHL9I sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:15:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150521_121528.het_loc.GMIDP sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:15:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_121528.het_loc.GHQ78 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:17:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ap2gamma-shi p_het_sites_in_narrow_peak_dp.py helas3 ap2gamma helas3 20150521_121749.het_loc.630A7 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:17:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150521_121751.het_loc.1UXV0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:17:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150521_121753.het_loc.W5COK sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:17:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150521_121755.het_loc.QIXDI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:17:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brca1-shi p_het_sites_in_narrow_peak_dp.py helas3 brca1 helas3 20150521_121757.het_loc.X07BK sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:00 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf1-shi p_het_sites_in_narrow_peak_dp.py helas3 brf1 helas3 20150521_121759.het_loc.TNVGM sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf2-shi p_het_sites_in_narrow_peak_dp.py helas3 brf2 helas3 20150521_121801.het_loc.LZESO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brg1-shi p_het_sites_in_narrow_peak_dp.py helas3 brg1 helas3 20150521_121803.het_loc.DVFIT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150521_121805.het_loc.VNETJ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:08 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cfos-shi p_het_sites_in_narrow_peak_dp.py helas3 cfos helas3 20150521_121807.het_loc.LK49Z sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:10 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150521_121809.het_loc.X01NI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cjun-shi p_het_sites_in_narrow_peak_dp.py helas3 cjun helas3 20150521_121811.het_loc.XPDCG sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cmyc-shi p_het_sites_in_narrow_peak_dp.py helas3 cmyc helas3 20150521_121813.het_loc.CEKU4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150521_121815.het_loc.ZFDGR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ctcf-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf helas3 20150521_121817.het_loc.MPRL1 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f1 helas3 20150521_121819.het_loc.DMDJR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f4-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f4 helas3 20150521_121821.het_loc.PZ90K sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:24 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f6-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f6 helas3 20150521_121823.het_loc.WSQPF sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk1-shi p_het_sites_in_narrow_peak_dp.py helas3 elk1 helas3 20150521_121825.het_loc.J3E79 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk4-shi p_het_sites_in_narrow_peak_dp.py helas3 elk4 helas3 20150521_121827.het_loc.VZ0AL sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:30 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gabp-shi p_het_sites_in_narrow_peak_dp.py helas3 gabp helas3 20150521_121829.het_loc.6ADS8 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150521_121831.het_loc.NBLBR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ini1-shi p_het_sites_in_narrow_peak_dp.py helas3 ini1 helas3 20150521_121833.het_loc.5BSJU sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-irf3-shi p_het_sites_in_narrow_peak_dp.py helas3 irf3 helas3 20150521_121835.het_loc.53RGU sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:38 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-jund-shi p_het_sites_in_narrow_peak_dp.py helas3 jund helas3 20150521_121837.het_loc.N9U7P sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:40 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150521_121839.het_loc.O1D2O sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:42 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-max-shi p_het_sites_in_narrow_peak_dp.py helas3 max helas3 20150521_121841.het_loc.M1KH0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-maz-shi p_het_sites_in_narrow_peak_dp.py helas3 maz helas3 20150521_121843.het_loc.A2UKT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150521_121845.het_loc.EO246 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:48 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfya-shi p_het_sites_in_narrow_peak_dp.py helas3 nfya helas3 20150521_121847.het_loc.K452R sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:50 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfyb-shi p_het_sites_in_narrow_peak_dp.py helas3 nfyb helas3 20150521_121849.het_loc.XDB7N sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:52 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrf1-shi p_het_sites_in_narrow_peak_dp.py helas3 nrf1 helas3 20150521_121851.het_loc.9LANT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150521_121853.het_loc.FAN3V sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-p300-shi p_het_sites_in_narrow_peak_dp.py helas3 p300 helas3 20150521_121855.het_loc.91DPJ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:18:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pol2-shi p_het_sites_in_narrow_peak_dp.py helas3 pol2 helas3 20150521_121857.het_loc.4BKR7 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:00 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150521_121859.het_loc.TEPFV sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150521_121901.het_loc.G3ELX sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150521_121903.het_loc.IIUN2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150521_121905.het_loc.OB2QT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:08 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-smc3-shi p_het_sites_in_narrow_peak_dp.py helas3 smc3 helas3 20150521_121907.het_loc.D0ULA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:10 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-spt20-shi p_het_sites_in_narrow_peak_dp.py helas3 spt20 helas3 20150521_121909.het_loc.G7GW3 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:12 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-stat3-shi p_het_sites_in_narrow_peak_dp.py helas3 stat3 helas3 20150521_121911.het_loc.YR4NT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-taf1-shi p_het_sites_in_narrow_peak_dp.py helas3 taf1 helas3 20150521_121913.het_loc.NALDC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tbp-shi p_het_sites_in_narrow_peak_dp.py helas3 tbp helas3 20150521_121915.het_loc.APJUI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150521_121918.het_loc.SBR4O sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-test-shi p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20150521_121920.het_loc.WLPY2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tr4-shi p_het_sites_in_narrow_peak_dp.py helas3 tr4 helas3 20150521_121922.het_loc.1O94S sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:24 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-usf2-shi p_het_sites_in_narrow_peak_dp.py helas3 usf2 helas3 20150521_121924.het_loc.LK27L sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150521_121926.het_loc.9GF2Z sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150521_121928.het_loc.A5VXT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:30 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf274-shi p_het_sites_in_narrow_peak_dp.py helas3 znf274 helas3 20150521_121930.het_loc.UF3CK sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 12:19:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zzz3-shi p_het_sites_in_narrow_peak_dp.py helas3 zzz3 helas3 20150521_121932.het_loc.NEVB7 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:18 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf2 gm12878 20150521_161717.het_loc.91WCM sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:19 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf3 gm12878 20150521_161719.het_loc.GT65T sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-batf-shi p_het_sites_in_narrow_peak_dp.py gm12878 batf gm12878 20150521_161721.het_loc.DWLE6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:24 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl11a-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl11a gm12878 20150521_161723.het_loc.SWCRB sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl3-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl3 gm12878 20150521_161725.het_loc.6CC2H sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:28 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bclaf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 bclaf1 gm12878 20150521_161727.het_loc.C0GHT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:30 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bhlhe40-shi p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40 gm12878 20150521_161729.het_loc.Y0NJC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:32 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-brca1-shi p_het_sites_in_narrow_peak_dp.py gm12878 brca1 gm12878 20150521_161731.het_loc.RDWS7 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:34 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cebpb-shi p_het_sites_in_narrow_peak_dp.py gm12878 cebpb gm12878 20150521_161733.het_loc.TQX9T sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:36 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cfos-shi p_het_sites_in_narrow_peak_dp.py gm12878 cfos gm12878 20150521_161735.het_loc.QFLNG sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:38 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd1-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd1 gm12878 20150521_161737.het_loc.COHJR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150521_161739.het_loc.ZZCU4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cmyc-shi p_het_sites_in_narrow_peak_dp.py gm12878 cmyc gm12878 20150521_161741.het_loc.LMYX5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:44 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-corest-shi p_het_sites_in_narrow_peak_dp.py gm12878 corest gm12878 20150521_161743.het_loc.UH49R sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ctcf-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12878 20150521_161745.het_loc.LPEZX sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-e2f4-shi p_het_sites_in_narrow_peak_dp.py gm12878 e2f4 gm12878 20150521_161747.het_loc.WTL1D sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:50 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ebf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ebf1 gm12878 20150521_161749.het_loc.Y8J1M sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-egr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 egr1 gm12878 20150521_161751.het_loc.K8FLW sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elf1 gm12878 20150521_161753.het_loc.BAIK0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elk1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elk1 gm12878 20150521_161755.het_loc.KPJPF sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:17:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ets1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ets1 gm12878 20150521_161757.het_loc.HHV9L sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-foxm1-shi p_het_sites_in_narrow_peak_dp.py gm12878 foxm1 gm12878 20150521_161759.het_loc.NBN0C sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-gabp-shi p_het_sites_in_narrow_peak_dp.py gm12878 gabp gm12878 20150521_161801.het_loc.X1GER sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ikzf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ikzf1 gm12878 20150521_161803.het_loc.L5MHW sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-input-shi p_het_sites_in_narrow_peak_dp.py gm12878 input gm12878 20150521_161805.het_loc.WO4OJ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-irf4-shi p_het_sites_in_narrow_peak_dp.py gm12878 irf4 gm12878 20150521_161807.het_loc.ESKH1 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:10 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-jund-shi p_het_sites_in_narrow_peak_dp.py gm12878 jund gm12878 20150521_161809.het_loc.BZAP2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:12 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-max-shi p_het_sites_in_narrow_peak_dp.py gm12878 max gm12878 20150521_161811.het_loc.1BM7A sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:14 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-maz-shi p_het_sites_in_narrow_peak_dp.py gm12878 maz gm12878 20150521_161813.het_loc.JGKNP sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2a-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2a gm12878 20150521_161815.het_loc.SR4EL sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:18 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2c-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2c gm12878 20150521_161817.het_loc.W39WE sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mta3-shi p_het_sites_in_narrow_peak_dp.py gm12878 mta3 gm12878 20150521_161819.het_loc.USOGS sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:22 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mxi1-shi p_het_sites_in_narrow_peak_dp.py gm12878 mxi1 gm12878 20150521_161821.het_loc.4JQI3 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:24 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150521_161823.het_loc.TFIP3 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfatc1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfatc1 gm12878 20150521_161825.het_loc.73ZQB sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:28 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfe2-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfe2 gm12878 20150521_161827.het_loc.W8NDA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:30 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfic-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfic gm12878 20150521_161829.het_loc.C6IMI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:32 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfkb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfkb gm12878 20150521_161831.het_loc.J4FJC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:34 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfya-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfya gm12878 20150521_161833.het_loc.2E7PQ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:36 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfyb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfyb gm12878 20150521_161835.het_loc.CGH6L sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:38 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrf1 gm12878 20150521_161837.het_loc.FK8T4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrsf-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrsf gm12878 20150521_161840.het_loc.4EUXE sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-p300-shi p_het_sites_in_narrow_peak_dp.py gm12878 p300 gm12878 20150521_161842.het_loc.9W9UQ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:44 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5 gm12878 20150521_161844.het_loc.22EXR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5c20-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5c20 gm12878 20150521_161846.het_loc.TNGTL sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5n19-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5n19 gm12878 20150521_161848.het_loc.01XMO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:50 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pbx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pbx3 gm12878 20150521_161850.het_loc.5LS35 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pml-shi p_het_sites_in_narrow_peak_dp.py gm12878 pml gm12878 20150521_161852.het_loc.FLE6S sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol2 gm12878 20150521_161854.het_loc.5CU3D sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol24h8-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol24h8 gm12878 20150521_161856.het_loc.MSZ6A sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:18:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol3 gm12878 20150521_161858.het_loc.WB2N5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pou2f2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pou2f2 gm12878 20150521_161900.het_loc.BRN48 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150521_161902.het_loc.EUXQH sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rad21-shi p_het_sites_in_narrow_peak_dp.py gm12878 rad21 gm12878 20150521_161904.het_loc.IZJP4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rfx5-shi p_het_sites_in_narrow_peak_dp.py gm12878 rfx5 gm12878 20150521_161906.het_loc.3Z13D sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-runx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 runx3 gm12878 20150521_161908.het_loc.V1ZS5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:10 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rxra-shi p_het_sites_in_narrow_peak_dp.py gm12878 rxra gm12878 20150521_161910.het_loc.IJ62B sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:12 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sin3a-shi p_het_sites_in_narrow_peak_dp.py gm12878 sin3a gm12878 20150521_161912.het_loc.R6CE3 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:14 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-six5-shi p_het_sites_in_narrow_peak_dp.py gm12878 six5 gm12878 20150521_161914.het_loc.23VNM sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-smc3-shi p_het_sites_in_narrow_peak_dp.py gm12878 smc3 gm12878 20150521_161916.het_loc.GO4ZR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:18 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sp1-shi p_het_sites_in_narrow_peak_dp.py gm12878 sp1 gm12878 20150521_161918.het_loc.Z749M sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-srf-shi p_het_sites_in_narrow_peak_dp.py gm12878 srf gm12878 20150521_161920.het_loc.Y1JFK sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:23 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat1-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat1 gm12878 20150521_161922.het_loc.C0RZO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:25 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat3-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat3 gm12878 20150521_161924.het_loc.4USY9 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat5a-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat5a gm12878 20150521_161926.het_loc.Z5HUK sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:29 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-taf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 taf1 gm12878 20150521_161928.het_loc.J9ISY sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tblr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 tblr1 gm12878 20150521_161930.het_loc.GE448 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:33 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tbp-shi p_het_sites_in_narrow_peak_dp.py gm12878 tbp gm12878 20150521_161932.het_loc.2JMHA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:35 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf12-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf12 gm12878 20150521_161934.het_loc.MM3MA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf3 gm12878 20150521_161936.het_loc.WIOPR sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:39 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tr4-shi p_het_sites_in_narrow_peak_dp.py gm12878 tr4 gm12878 20150521_161938.het_loc.KM10P sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf1 gm12878 20150521_161940.het_loc.6UTOC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:43 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf2 gm12878 20150521_161942.het_loc.K2AFD sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:45 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-whip-shi p_het_sites_in_narrow_peak_dp.py gm12878 whip gm12878 20150521_161944.het_loc.GRR0J sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-yy1-shi p_het_sites_in_narrow_peak_dp.py gm12878 yy1 gm12878 20150521_161946.het_loc.OMNT4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:49 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zbtb33-shi p_het_sites_in_narrow_peak_dp.py gm12878 zbtb33 gm12878 20150521_161948.het_loc.48FCW sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zeb1-shi p_het_sites_in_narrow_peak_dp.py gm12878 zeb1 gm12878 20150521_161951.het_loc.SH9YZ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:53 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf143-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf143 gm12878 20150521_161953.het_loc.8Y2SW sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:55 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf274-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf274 gm12878 20150521_161954.het_loc.3G8JD sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 16:19:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zzz3-shi p_het_sites_in_narrow_peak_dp.py gm12878 zzz3 gm12878 20150521_161957.het_loc.7LS8H sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf2 gm12878 20150521_172139.het_loc.96WMV sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf3 gm12878 20150521_172141.het_loc.1FRFO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:44 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-batf-shi p_het_sites_in_narrow_peak_dp.py gm12878 batf gm12878 20150521_172143.het_loc.6YQEF sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl11a-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl11a gm12878 20150521_172145.het_loc.LSKPX sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl3-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl3 gm12878 20150521_172147.het_loc.P6A2T sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:50 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bclaf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 bclaf1 gm12878 20150521_172149.het_loc.TJS3P sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bhlhe40-shi p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40 gm12878 20150521_172151.het_loc.Z525V sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-brca1-shi p_het_sites_in_narrow_peak_dp.py gm12878 brca1 gm12878 20150521_172153.het_loc.6D4HE sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cebpb-shi p_het_sites_in_narrow_peak_dp.py gm12878 cebpb gm12878 20150521_172155.het_loc.VIFWA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:21:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cfos-shi p_het_sites_in_narrow_peak_dp.py gm12878 cfos gm12878 20150521_172157.het_loc.RO5O9 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd1-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd1 gm12878 20150521_172159.het_loc.CMKY6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150521_172201.het_loc.OHBM5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cmyc-shi p_het_sites_in_narrow_peak_dp.py gm12878 cmyc gm12878 20150521_172203.het_loc.DOZDI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-corest-shi p_het_sites_in_narrow_peak_dp.py gm12878 corest gm12878 20150521_172205.het_loc.DFLL9 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ctcf-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12878 20150521_172207.het_loc.UPJV3 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:10 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-e2f4-shi p_het_sites_in_narrow_peak_dp.py gm12878 e2f4 gm12878 20150521_172209.het_loc.K5YCC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:12 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ebf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ebf1 gm12878 20150521_172211.het_loc.CJ5E5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:14 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-egr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 egr1 gm12878 20150521_172213.het_loc.SWKKT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elf1 gm12878 20150521_172215.het_loc.TG75C sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:18 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elk1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elk1 gm12878 20150521_172217.het_loc.GM8R0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ets1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ets1 gm12878 20150521_172219.het_loc.239HI sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:22 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-foxm1-shi p_het_sites_in_narrow_peak_dp.py gm12878 foxm1 gm12878 20150521_172221.het_loc.1QO1C sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:24 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-gabp-shi p_het_sites_in_narrow_peak_dp.py gm12878 gabp gm12878 20150521_172223.het_loc.UPN15 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ikzf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ikzf1 gm12878 20150521_172225.het_loc.4BHFQ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:28 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-input-shi p_het_sites_in_narrow_peak_dp.py gm12878 input gm12878 20150521_172227.het_loc.R2JNL sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:30 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-irf4-shi p_het_sites_in_narrow_peak_dp.py gm12878 irf4 gm12878 20150521_172229.het_loc.AZT2T sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:32 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-jund-shi p_het_sites_in_narrow_peak_dp.py gm12878 jund gm12878 20150521_172231.het_loc.2Y2SW sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:34 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-max-shi p_het_sites_in_narrow_peak_dp.py gm12878 max gm12878 20150521_172233.het_loc.6SAE7 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:36 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-maz-shi p_het_sites_in_narrow_peak_dp.py gm12878 maz gm12878 20150521_172235.het_loc.L14BC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:38 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2a-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2a gm12878 20150521_172237.het_loc.50TLF sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2c-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2c gm12878 20150521_172239.het_loc.03SRP sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mta3-shi p_het_sites_in_narrow_peak_dp.py gm12878 mta3 gm12878 20150521_172241.het_loc.8WB7D sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:44 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mxi1-shi p_het_sites_in_narrow_peak_dp.py gm12878 mxi1 gm12878 20150521_172243.het_loc.3D6HP sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150521_172245.het_loc.5ZM1Z sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfatc1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfatc1 gm12878 20150521_172247.het_loc.862DS sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:50 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfe2-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfe2 gm12878 20150521_172249.het_loc.IXDPC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfic-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfic gm12878 20150521_172251.het_loc.IEWG4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfkb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfkb gm12878 20150521_172253.het_loc.82NDO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfya-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfya gm12878 20150521_172255.het_loc.9RDZU sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:22:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfyb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfyb gm12878 20150521_172257.het_loc.VZ1N2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrf1 gm12878 20150521_172300.het_loc.TJBN0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:03 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrsf-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrsf gm12878 20150521_172302.het_loc.PMO8C sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-p300-shi p_het_sites_in_narrow_peak_dp.py gm12878 p300 gm12878 20150521_172304.het_loc.8GQI6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5 gm12878 20150521_172306.het_loc.RJCL5 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:09 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5c20-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5c20 gm12878 20150521_172308.het_loc.PEYJY sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:11 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5n19-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5n19 gm12878 20150521_172310.het_loc.RJJCL sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:13 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pbx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pbx3 gm12878 20150521_172312.het_loc.094DZ sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:15 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pml-shi p_het_sites_in_narrow_peak_dp.py gm12878 pml gm12878 20150521_172314.het_loc.77488 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol2 gm12878 20150521_172316.het_loc.Q4X46 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:19 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol24h8-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol24h8 gm12878 20150521_172318.het_loc.4L63Q sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:21 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol3 gm12878 20150521_172320.het_loc.Q74TB sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:23 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pou2f2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pou2f2 gm12878 20150521_172322.het_loc.ZLZFB sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:25 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150521_172324.het_loc.LXTZ4 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rad21-shi p_het_sites_in_narrow_peak_dp.py gm12878 rad21 gm12878 20150521_172326.het_loc.C4UL0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:29 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rfx5-shi p_het_sites_in_narrow_peak_dp.py gm12878 rfx5 gm12878 20150521_172328.het_loc.PW8XT sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-runx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 runx3 gm12878 20150521_172330.het_loc.WF1FD sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:33 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rxra-shi p_het_sites_in_narrow_peak_dp.py gm12878 rxra gm12878 20150521_172332.het_loc.VS24L sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:35 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sin3a-shi p_het_sites_in_narrow_peak_dp.py gm12878 sin3a gm12878 20150521_172334.het_loc.MQ7S6 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-six5-shi p_het_sites_in_narrow_peak_dp.py gm12878 six5 gm12878 20150521_172336.het_loc.EEMFB sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:39 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-smc3-shi p_het_sites_in_narrow_peak_dp.py gm12878 smc3 gm12878 20150521_172338.het_loc.JH41O sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sp1-shi p_het_sites_in_narrow_peak_dp.py gm12878 sp1 gm12878 20150521_172340.het_loc.ZMAQL sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:43 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-srf-shi p_het_sites_in_narrow_peak_dp.py gm12878 srf gm12878 20150521_172342.het_loc.O67CV sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:45 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat1-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat1 gm12878 20150521_172344.het_loc.YG672 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat3-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat3 gm12878 20150521_172346.het_loc.G66Q2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:49 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat5a-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat5a gm12878 20150521_172348.het_loc.GWD3Z sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-taf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 taf1 gm12878 20150521_172350.het_loc.N9R4P sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:53 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tblr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 tblr1 gm12878 20150521_172352.het_loc.HK4IA sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:55 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tbp-shi p_het_sites_in_narrow_peak_dp.py gm12878 tbp gm12878 20150521_172354.het_loc.76F7X sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf12-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf12 gm12878 20150521_172356.het_loc.QK9FD sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:23:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf3 gm12878 20150521_172358.het_loc.B3SE0 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tr4-shi p_het_sites_in_narrow_peak_dp.py gm12878 tr4 gm12878 20150521_172400.het_loc.31JAS sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:03 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf1 gm12878 20150521_172402.het_loc.L1G75 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf2 gm12878 20150521_172404.het_loc.QY95Q sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-whip-shi p_het_sites_in_narrow_peak_dp.py gm12878 whip gm12878 20150521_172407.het_loc.1DU1J sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:09 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-yy1-shi p_het_sites_in_narrow_peak_dp.py gm12878 yy1 gm12878 20150521_172409.het_loc.1AKFE sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:11 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zbtb33-shi p_het_sites_in_narrow_peak_dp.py gm12878 zbtb33 gm12878 20150521_172411.het_loc.ACKR2 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:13 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zeb1-shi p_het_sites_in_narrow_peak_dp.py gm12878 zeb1 gm12878 20150521_172413.het_loc.91Y36 sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:15 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf143-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf143 gm12878 20150521_172415.het_loc.W3Z8Q sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf274-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf274 gm12878 20150521_172417.het_loc.RZJCO sydh:uw:haib:uta:embl 1 1 1 1
#[Thu May 21 17:24:19 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zzz3-shi p_het_sites_in_narrow_peak_dp.py gm12878 zzz3 gm12878 20150521_172419.het_loc.6W6M4 sydh:uw:haib:uta:embl 1 1 1 1
#[Sun Jun 28 15:20:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-simulate-shi p_het_sites_in_narrow_peak_dp.py gm12878 simulate gm12878 20150628_151959.het_loc.FCFJI sydh:uw:haib:uta:embl:encode 1 1 1 1
#[Sun Jun 28 15:23:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-simulate-shi p_het_sites_in_narrow_peak_dp.py gm12878 simulate gm12878 20150628_152259.het_loc.7EKCH sydh:uw:haib:uta:embl:encode 1 1 1 1
#[Sun Jun 28 15:25:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-simulate-shi p_het_sites_in_narrow_peak_dp.py gm12878 simulate gm12878 20150628_152547.het_loc.ZP9EO sydh:uw:haib:uta:embl:encode 1 1 1 1
#[Mon Jun 29 11:51:11 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-simulate-shi p_het_sites_in_narrow_peak_dp.py helas3 simulate helas3 20150629_115111.het_loc.BPT2U sydh:uw:haib:uta:embl:encode 1 1 1 1
#[Mon Jun 29 17:48:13 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ap2gamma-shi p_het_sites_in_narrow_peak_dp.py helas3 ap2gamma helas3 20150629_174812.het_loc.QS9JO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:14 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf155-shi p_het_sites_in_narrow_peak_dp.py helas3 baf155 helas3 20150629_174813.het_loc.2326X sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:15 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-baf170-shi p_het_sites_in_narrow_peak_dp.py helas3 baf170 helas3 20150629_174814.het_loc.B7ZSN sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:16 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-bdp1-shi p_het_sites_in_narrow_peak_dp.py helas3 bdp1 helas3 20150629_174815.het_loc.ENW8R sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:17 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brca1-shi p_het_sites_in_narrow_peak_dp.py helas3 brca1 helas3 20150629_174816.het_loc.Z4SS2 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:18 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf1-shi p_het_sites_in_narrow_peak_dp.py helas3 brf1 helas3 20150629_174817.het_loc.UX2F4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:19 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brf2-shi p_het_sites_in_narrow_peak_dp.py helas3 brf2 helas3 20150629_174818.het_loc.P4GKM sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:20 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-brg1-shi p_het_sites_in_narrow_peak_dp.py helas3 brg1 helas3 20150629_174819.het_loc.XCKED sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:21 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cebpb-shi p_het_sites_in_narrow_peak_dp.py helas3 cebpb helas3 20150629_174820.het_loc.AQSUD sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:22 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cfos-shi p_het_sites_in_narrow_peak_dp.py helas3 cfos helas3 20150629_174821.het_loc.C0GX9 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:23 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-chd2-shi p_het_sites_in_narrow_peak_dp.py helas3 chd2 helas3 20150629_174823.het_loc.TIU72 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:24 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cjun-shi p_het_sites_in_narrow_peak_dp.py helas3 cjun helas3 20150629_174824.het_loc.VOM1R sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:25 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-cmyc-shi p_het_sites_in_narrow_peak_dp.py helas3 cmyc helas3 20150629_174825.het_loc.V3IK6 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:26 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-corest-shi p_het_sites_in_narrow_peak_dp.py helas3 corest helas3 20150629_174826.het_loc.DAU3G sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:27 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ctcf-shi p_het_sites_in_narrow_peak_dp.py helas3 ctcf helas3 20150629_174827.het_loc.DW3BR sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:28 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f1 helas3 20150629_174828.het_loc.HKHIX sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:29 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f4-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f4 helas3 20150629_174829.het_loc.IN0OK sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:30 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-e2f6-shi p_het_sites_in_narrow_peak_dp.py helas3 e2f6 helas3 20150629_174830.het_loc.LJDFQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:31 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk1-shi p_het_sites_in_narrow_peak_dp.py helas3 elk1 helas3 20150629_174831.het_loc.DAR7Z sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:32 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-elk4-shi p_het_sites_in_narrow_peak_dp.py helas3 elk4 helas3 20150629_174832.het_loc.4DIWL sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:33 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gabp-shi p_het_sites_in_narrow_peak_dp.py helas3 gabp helas3 20150629_174833.het_loc.9UPW1 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:34 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-gtf2f1-shi p_het_sites_in_narrow_peak_dp.py helas3 gtf2f1 helas3 20150629_174834.het_loc.6OI8E sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:36 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-ini1-shi p_het_sites_in_narrow_peak_dp.py helas3 ini1 helas3 20150629_174835.het_loc.DZTY7 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:37 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-irf3-shi p_het_sites_in_narrow_peak_dp.py helas3 irf3 helas3 20150629_174836.het_loc.NHNX2 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:38 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-jund-shi p_het_sites_in_narrow_peak_dp.py helas3 jund helas3 20150629_174837.het_loc.474O6 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:39 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mafk-shi p_het_sites_in_narrow_peak_dp.py helas3 mafk helas3 20150629_174838.het_loc.GPHT4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:40 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-max-shi p_het_sites_in_narrow_peak_dp.py helas3 max helas3 20150629_174839.het_loc.PNJQO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:41 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-maz-shi p_het_sites_in_narrow_peak_dp.py helas3 maz helas3 20150629_174840.het_loc.9VO3Z sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:42 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-mxi1-shi p_het_sites_in_narrow_peak_dp.py helas3 mxi1 helas3 20150629_174841.het_loc.XGMEC sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:43 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfya-shi p_het_sites_in_narrow_peak_dp.py helas3 nfya helas3 20150629_174842.het_loc.XBKTR sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:44 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nfyb-shi p_het_sites_in_narrow_peak_dp.py helas3 nfyb helas3 20150629_174843.het_loc.TW6CZ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:45 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrf1-shi p_het_sites_in_narrow_peak_dp.py helas3 nrf1 helas3 20150629_174844.het_loc.ZZ89J sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:46 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-nrsf-shi p_het_sites_in_narrow_peak_dp.py helas3 nrsf helas3 20150629_174845.het_loc.HJOXL sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:47 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-p300-shi p_het_sites_in_narrow_peak_dp.py helas3 p300 helas3 20150629_174846.het_loc.KYOH0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:48 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-pol2-shi p_het_sites_in_narrow_peak_dp.py helas3 pol2 helas3 20150629_174847.het_loc.YX7FW sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf2 gm12878 20150629_174848.het_loc.QSXH2 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:49 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-prdm1-shi p_het_sites_in_narrow_peak_dp.py helas3 prdm1 helas3 20150629_174848.het_loc.X0GEC sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:50 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-atf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 atf3 gm12878 20150629_174849.het_loc.R45X5 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:50 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rad21-shi p_het_sites_in_narrow_peak_dp.py helas3 rad21 helas3 20150629_174849.het_loc.1BDCD sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-batf-shi p_het_sites_in_narrow_peak_dp.py gm12878 batf gm12878 20150629_174850.het_loc.TUYFM sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:51 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rfx5-shi p_het_sites_in_narrow_peak_dp.py helas3 rfx5 helas3 20150629_174850.het_loc.DE609 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl11a-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl11a gm12878 20150629_174851.het_loc.XM54Z sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:52 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-rpc155-shi p_het_sites_in_narrow_peak_dp.py helas3 rpc155 helas3 20150629_174851.het_loc.OPHOH sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:53 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bcl3-shi p_het_sites_in_narrow_peak_dp.py gm12878 bcl3 gm12878 20150629_174852.het_loc.L8JKQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:53 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-smc3-shi p_het_sites_in_narrow_peak_dp.py helas3 smc3 helas3 20150629_174852.het_loc.2F00B sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bclaf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 bclaf1 gm12878 20150629_174853.het_loc.2ZEGG sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:54 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-spt20-shi p_het_sites_in_narrow_peak_dp.py helas3 spt20 helas3 20150629_174853.het_loc.41M09 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:55 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-bhlhe40-shi p_het_sites_in_narrow_peak_dp.py gm12878 bhlhe40 gm12878 20150629_174854.het_loc.LV8VM sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:55 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-stat3-shi p_het_sites_in_narrow_peak_dp.py helas3 stat3 helas3 20150629_174855.het_loc.3N5BN sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-brca1-shi p_het_sites_in_narrow_peak_dp.py gm12878 brca1 gm12878 20150629_174855.het_loc.5VFTI sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:56 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-taf1-shi p_het_sites_in_narrow_peak_dp.py helas3 taf1 helas3 20150629_174856.het_loc.NCO1X sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cebpb-shi p_het_sites_in_narrow_peak_dp.py gm12878 cebpb gm12878 20150629_174856.het_loc.OGTKW sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:57 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tbp-shi p_het_sites_in_narrow_peak_dp.py helas3 tbp helas3 20150629_174857.het_loc.BUG7W sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cfos-shi p_het_sites_in_narrow_peak_dp.py gm12878 cfos gm12878 20150629_174857.het_loc.RPBFS sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:58 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tcf7l2-shi p_het_sites_in_narrow_peak_dp.py helas3 tcf7l2 helas3 20150629_174858.het_loc.937GQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd1-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd1 gm12878 20150629_174858.het_loc.R00GU sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:48:59 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-test-shi p_het_sites_in_narrow_peak_dp.py helas3 test helas3 20150629_174859.het_loc.IOHJ0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-chd2-shi p_het_sites_in_narrow_peak_dp.py gm12878 chd2 gm12878 20150629_174859.het_loc.MCD6Y sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:00 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-tr4-shi p_het_sites_in_narrow_peak_dp.py helas3 tr4 helas3 20150629_174900.het_loc.4P7IF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-cmyc-shi p_het_sites_in_narrow_peak_dp.py gm12878 cmyc gm12878 20150629_174900.het_loc.B33RS sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:01 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-usf2-shi p_het_sites_in_narrow_peak_dp.py helas3 usf2 helas3 20150629_174901.het_loc.GHSQO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-corest-shi p_het_sites_in_narrow_peak_dp.py gm12878 corest gm12878 20150629_174901.het_loc.3KM36 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:02 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zkscan1-shi p_het_sites_in_narrow_peak_dp.py helas3 zkscan1 helas3 20150629_174902.het_loc.S0AEG sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:03 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ctcf-shi p_het_sites_in_narrow_peak_dp.py gm12878 ctcf gm12878 20150629_174902.het_loc.5E2NT sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:03 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf143-shi p_het_sites_in_narrow_peak_dp.py helas3 znf143 helas3 20150629_174903.het_loc.FFSXE sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-e2f4-shi p_het_sites_in_narrow_peak_dp.py gm12878 e2f4 gm12878 20150629_174903.het_loc.A9H71 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:04 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-znf274-shi p_het_sites_in_narrow_peak_dp.py helas3 znf274 helas3 20150629_174904.het_loc.8ACOU sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ebf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ebf1 gm12878 20150629_174904.het_loc.RGVSY sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:06 2015] p p_run_cluster_sep.py het_loc2-helas3-helas3-zzz3-shi p_het_sites_in_narrow_peak_dp.py helas3 zzz3 helas3 20150629_174905.het_loc.UOIO3 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-egr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 egr1 gm12878 20150629_174905.het_loc.2N75M sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:07 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elf1 gm12878 20150629_174906.het_loc.Z6X9B sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-elk1-shi p_het_sites_in_narrow_peak_dp.py gm12878 elk1 gm12878 20150629_174907.het_loc.1XLL8 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:09 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ets1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ets1 gm12878 20150629_174908.het_loc.4HBYP sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:10 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-foxm1-shi p_het_sites_in_narrow_peak_dp.py gm12878 foxm1 gm12878 20150629_174909.het_loc.EHKRY sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:11 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-gabp-shi p_het_sites_in_narrow_peak_dp.py gm12878 gabp gm12878 20150629_174910.het_loc.UGKKU sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:12 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-ikzf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 ikzf1 gm12878 20150629_174911.het_loc.J1X91 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:13 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-input-shi p_het_sites_in_narrow_peak_dp.py gm12878 input gm12878 20150629_174912.het_loc.8WIBC sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:14 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-irf4-shi p_het_sites_in_narrow_peak_dp.py gm12878 irf4 gm12878 20150629_174913.het_loc.RKN02 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:15 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-jund-shi p_het_sites_in_narrow_peak_dp.py gm12878 jund gm12878 20150629_174914.het_loc.JD6BG sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:16 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-max-shi p_het_sites_in_narrow_peak_dp.py gm12878 max gm12878 20150629_174915.het_loc.KH5RY sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:17 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-maz-shi p_het_sites_in_narrow_peak_dp.py gm12878 maz gm12878 20150629_174916.het_loc.CFZDZ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:18 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2a-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2a gm12878 20150629_174917.het_loc.RQ2QY sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mef2c-shi p_het_sites_in_narrow_peak_dp.py gm12878 mef2c gm12878 20150629_174918.het_loc.RE4C3 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:20 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mta3-shi p_het_sites_in_narrow_peak_dp.py gm12878 mta3 gm12878 20150629_174919.het_loc.1E6K4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:22 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-mxi1-shi p_het_sites_in_narrow_peak_dp.py gm12878 mxi1 gm12878 20150629_174920.het_loc.WSRGS sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:23 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-myc-shi p_het_sites_in_narrow_peak_dp.py gm12878 myc gm12878 20150629_174921.het_loc.7T3PB sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:25 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfatc1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfatc1 gm12878 20150629_174922.het_loc.525WZ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:25 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfe2-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfe2 gm12878 20150629_174923.het_loc.D56RD sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfic-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfic gm12878 20150629_174924.het_loc.2DLQ3 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:26 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfkb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfkb gm12878 20150629_174925.het_loc.30LSD sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:27 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfya-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfya gm12878 20150629_174926.het_loc.FZ7YA sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:28 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nfyb-shi p_het_sites_in_narrow_peak_dp.py gm12878 nfyb gm12878 20150629_174927.het_loc.LH3CB sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:29 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrf1 gm12878 20150629_174928.het_loc.I4FS0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:30 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-nrsf-shi p_het_sites_in_narrow_peak_dp.py gm12878 nrsf gm12878 20150629_174929.het_loc.FQRE5 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:31 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-p300-shi p_het_sites_in_narrow_peak_dp.py gm12878 p300 gm12878 20150629_174930.het_loc.53H54 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:32 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5 gm12878 20150629_174931.het_loc.JKE3P sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:33 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5c20-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5c20 gm12878 20150629_174932.het_loc.GKLLQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:34 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pax5n19-shi p_het_sites_in_narrow_peak_dp.py gm12878 pax5n19 gm12878 20150629_174933.het_loc.NDDR2 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:35 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pbx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pbx3 gm12878 20150629_174934.het_loc.FCYZ8 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:36 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pml-shi p_het_sites_in_narrow_peak_dp.py gm12878 pml gm12878 20150629_174935.het_loc.XE8FJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:37 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol2 gm12878 20150629_174936.het_loc.X1VSR sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:38 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol24h8-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol24h8 gm12878 20150629_174937.het_loc.C5K1D sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:39 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pol3-shi p_het_sites_in_narrow_peak_dp.py gm12878 pol3 gm12878 20150629_174938.het_loc.LGSZK sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:40 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pou2f2-shi p_het_sites_in_narrow_peak_dp.py gm12878 pou2f2 gm12878 20150629_174939.het_loc.OZ3DT sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:41 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-pu1-shi p_het_sites_in_narrow_peak_dp.py gm12878 pu1 gm12878 20150629_174940.het_loc.5QMDH sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:42 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rad21-shi p_het_sites_in_narrow_peak_dp.py gm12878 rad21 gm12878 20150629_174942.het_loc.HQG8O sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:43 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rfx5-shi p_het_sites_in_narrow_peak_dp.py gm12878 rfx5 gm12878 20150629_174943.het_loc.2H2Q3 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:44 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-runx3-shi p_het_sites_in_narrow_peak_dp.py gm12878 runx3 gm12878 20150629_174944.het_loc.EX7IY sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:45 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-rxra-shi p_het_sites_in_narrow_peak_dp.py gm12878 rxra gm12878 20150629_174945.het_loc.GECVF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:46 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sin3a-shi p_het_sites_in_narrow_peak_dp.py gm12878 sin3a gm12878 20150629_174946.het_loc.Y87LF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:47 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-six5-shi p_het_sites_in_narrow_peak_dp.py gm12878 six5 gm12878 20150629_174947.het_loc.P1FPE sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:48 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-smc3-shi p_het_sites_in_narrow_peak_dp.py gm12878 smc3 gm12878 20150629_174948.het_loc.JE8LI sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:49 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-sp1-shi p_het_sites_in_narrow_peak_dp.py gm12878 sp1 gm12878 20150629_174949.het_loc.RPZPT sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:50 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-srf-shi p_het_sites_in_narrow_peak_dp.py gm12878 srf gm12878 20150629_174950.het_loc.TX5WM sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:51 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat1-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat1 gm12878 20150629_174951.het_loc.DNWC5 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:52 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat3-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat3 gm12878 20150629_174952.het_loc.Q2TMT sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:53 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-stat5a-shi p_het_sites_in_narrow_peak_dp.py gm12878 stat5a gm12878 20150629_174953.het_loc.D3R57 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:54 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-taf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 taf1 gm12878 20150629_174954.het_loc.E7W1Q sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:55 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tblr1-shi p_het_sites_in_narrow_peak_dp.py gm12878 tblr1 gm12878 20150629_174955.het_loc.BAPP2 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:56 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tbp-shi p_het_sites_in_narrow_peak_dp.py gm12878 tbp gm12878 20150629_174956.het_loc.CP5SY sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:57 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf12-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf12 gm12878 20150629_174957.het_loc.AXX4P sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:58 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tcf3-shi p_het_sites_in_narrow_peak_dp.py gm12878 tcf3 gm12878 20150629_174958.het_loc.5AEEU sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:49:59 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-tr4-shi p_het_sites_in_narrow_peak_dp.py gm12878 tr4 gm12878 20150629_174959.het_loc.PALU8 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:00 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf1-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf1 gm12878 20150629_175000.het_loc.74I37 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:01 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-usf2-shi p_het_sites_in_narrow_peak_dp.py gm12878 usf2 gm12878 20150629_175001.het_loc.ZNIU3 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:02 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-whip-shi p_het_sites_in_narrow_peak_dp.py gm12878 whip gm12878 20150629_175002.het_loc.6W1IO sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:03 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-yy1-shi p_het_sites_in_narrow_peak_dp.py gm12878 yy1 gm12878 20150629_175003.het_loc.U2YIT sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:04 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zbtb33-shi p_het_sites_in_narrow_peak_dp.py gm12878 zbtb33 gm12878 20150629_175004.het_loc.J23S0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:05 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zeb1-shi p_het_sites_in_narrow_peak_dp.py gm12878 zeb1 gm12878 20150629_175005.het_loc.8WJN4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:06 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf143-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf143 gm12878 20150629_175006.het_loc.6XM1E sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:08 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-znf274-shi p_het_sites_in_narrow_peak_dp.py gm12878 znf274 gm12878 20150629_175007.het_loc.5T6HA sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Jun 29 17:50:09 2015] p p_run_cluster_sep.py het_loc2-gm12878-gm12878-zzz3-shi p_het_sites_in_narrow_peak_dp.py gm12878 zzz3 gm12878 20150629_175008.het_loc.5VJV6 sydh:uw:haib:uta:embl 1 1 1 1
