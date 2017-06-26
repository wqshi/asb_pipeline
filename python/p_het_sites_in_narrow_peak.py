#Pre-condition:
#*download the narrowPeak file in the bed dir

#Process:
#* find every narrowPeak file, sharp peak into 100 bp.
#* Overlap 100bp peak with each cell's vcf file in the wgs_dir, and then filter the heterozygous sites (*.het.loc file)

#After:
#* for each narrowPeak file, there will be a *.het.loc file.

#Part of the ASB pipeline after the  narrowPeak download.

#Input: when submit to the clustdell, need to specifiy the cell name


import os
import socket
import urllib
import re
import sys
import p_mymodule as my
import p_script_test as scripttest

server_name = socket.gethostname()
print server_name
if (server_name == "loire"):
    bed_dir="/homed/home/shi/test/t_het_sites_in_narrow_peak"
    test_envi = scripttest.scripttest(bed_dir)
    wgs_dir="/homed/home/shi/projects/wgs"

    tf_list=["ctcf"]
    #pattern="*.narrowPeak"
else:
    
    print sys.argv
    cell=sys.argv[1]
    
    #pattern="*.narrowPeak"
    #pattern="*whole.narrowPeak"
    
    tf_list = my.f_recieve_list_para(sys.argv[2])
    locker_file = sys.argv[4]
    cell_list = my.f_recieve_list_para(sys.argv[3])
    syn_dir=sys.argv[5]

    #bed_dir="/home/shi/projects/chipseq_snp/data2/encode/%s/" % cell
    bed_dir=sys.argv[6]
    wgs_dir=sys.argv[7]


    
reload(my)
pattern=".*(%s).narrowPeak$"%"|".join(tf_list)
narrowPeak_list=my.f_grep_files_from_dir(bed_dir, pattern)
print narrowPeak_list
print pattern
for narrowPeak_file in narrowPeak_list:
    if "Rep" in narrowPeak_file:
        print "Skip %s" % narrowPeak_file;
        continue;
    
    print os.path.basename(narrowPeak_file)
    output_file=narrowPeak_file + ".101bp.bed"
    cmd="f_narrow_peak_100bp %s %s"%( narrowPeak_file, output_file)
    my.f_call_shell_fun(cmd)
    
    file_names=my.f_parse_file_name(os.path.basename(narrowPeak_file))
    print os.path.basename(narrowPeak_file),file_names

    file_prefix = my.f_get_prefix(narrowPeak_file)
    if "gm" in file_names[1]:

        print("Interpret as gm cell")
        overlab_cmd = "sed '/^#/d' %s/%s.vcf | sed 's/^/chr/g' | intersectBed -u -a stdin -b %s > %s" %(wgs_dir, file_names[1], output_file, file_prefix + ".wgs.vcf")
    else:
        overlab_cmd = "sed '/^#/d' %s/%s.vcf | intersectBed -u -a stdin -b %s > %s" %(wgs_dir, file_names[1], output_file, file_prefix + ".wgs.vcf")
        
    my.f_shell_cmd(overlab_cmd)
    grep_het_cmd = "f_complete_genome_read_depth %s |  f_grep_legal_snp | sed '/^#/d' | grep -v '1/1' > %s "%(file_prefix + ".wgs.vcf", file_prefix + ".het.loc" )
    grep_alt_cmd = "f_complete_genome_read_depth %s |  f_grep_legal_snp | sed '/^#/d' | grep '1[/\|]1' > %s "%(file_prefix + ".wgs.vcf", file_prefix + ".alt.loc" )
    #os.remove(file_prefix + ".wgs.vcf")
    print grep_het_cmd
    print grep_alt_cmd
    my.f_shell_fun_pipe(grep_het_cmd)
    my.f_shell_fun_pipe(grep_alt_cmd)

    

if (server_name != "loire"):
    het_pattern=my.f_create_pattern(cell_list, tf_list,".het.loc")
    alt_pattern=my.f_create_pattern(cell_list, tf_list,".alt.loc")
    my.f_grep_and_scp_to_loire(bed_dir, het_pattern,syn_dir)
    my.f_grep_and_scp_to_loire(bed_dir, alt_pattern,syn_dir)
    my.f_scp_to_loire(bed_dir, ".*101bp.bed",syn_dir)
    #my.f_scp_to_loire(bed_dir, "*all.snv.loc",syn_dir)
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
