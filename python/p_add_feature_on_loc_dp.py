
import p_mymodule as my
import re
import os.path
import sys
import p_pd as loc
import p_script_test as scripttest
import socket
import logging
import pybedtools
import pandas as pd
import p_shell_function as shell
logging.basicConfig(level=logging.DEBUG)

logging.getLogger().setLevel(logging.DEBUG)
reload(loc)


server_name = socket.gethostname()
if (server_name == "loire"):
    data_dir="/homed/home/shi/test/t_het_sites_in_narrow_peak/"
    test_envi = scripttest.scripttest(data_dir)
    #loc_file   =dnase_dir + "/gm12878-ctcf.het.loc.test"
    diff_peak_dir='/homed/home/shi/test/t_normalize_read_depth/'
    tf_list = ["ctcf"]
    cell_list=["gmtest"]
    cons_file="/homed/home/shi/projects/wgs/conservation/phyloP7way.bed"
    dsQTL_file = '/homed/home/shi/projects/wgs/dsQTL.txt'
    ctcf_QTL_file = '/homed/home/shi/projects/wgs/ctcf_qtl.csv'
    hct_dir="/homed/home/shi/anthony/chipseq_snp/python/output/"
  
    consistent_peak_file = '/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/test/p_download_dnase/consistent_ctcf.bed'
    diff_peak_file = '/homed/home/shi/anthony/tfbs_chipseq/ENCODE/dnase/test/p_download_dnase/diff_ctcf.bed'
    loop_file = '/homed/home/shi/projects/wgs/flat_loop_ends.bed'
    loop_region_file = '/homed/home/shi/projects/wgs/loop_regions.bed'
    tf_list_file='/homed/home/shi/projects/wgs/tf_list.txt'
    #cofactor_list=my.f_parse_tf_list_file(tf_list_file)
    cofactor_list=["mzf1", "arnt", "ctcf", "tal1.tcf3", "tcf3", "brca1", "klf5", "gata4",  "foxd3", "atoh1",'znf143', 'smc3']
    diff_cell_list=["helatest"]
    dipoid_file = '/homed/home/shi/projects/wgs/diploid/gmtest.snp.vcf'
    debug_flag = False
    cell = 'gmtest'
    labs=["sydh","haib","uw","uchicago","uta"]
    #cell = 'gm12878'
    #diff_cell_list=["gm12872"]
    diffpeak_method = 'pepr'
    
else:
    cell=sys.argv[1]
    diff_cell_list = my.f_recieve_list_para(sys.argv[2])
    tf_list=my.f_recieve_list_para(sys.argv[3])
    locker_file = sys.argv[4]
    labs = my.f_recieve_list_para(sys.argv[5])
    diffpeak_method = sys.argv[6]
    cofactor_list = my.f_recieve_list_para( sys.argv[7] )
    logging.info('Diffpeak_mehtod:' + diffpeak_method)
    
    cons_file="/home/shi/projects/wgs/conservation/phyloP7way.bed"
    hct_dir="/home/shi/python/output/"
    dsQTL_file = '/home/shi/projects/wgs/dsQTL.txt'
    data_dir="/home/shi/projects/chipseq_snp/data2/encode/%s/" % cell
    
    diff_peak_file = '/home/shi/projects/wgs/diff_ctcf.bed'
    tf_list_file='/home/shi/projects/wgs/tf_list.txt'
    if tf_list[0] == 'pu1':
        cofactor_list=my.f_parse_tf_list_file(tf_list_file)
    #cofactor_list=["mzf1", "yy1", "mef2a", "jun" ,"ctcf", "tal1", "tcf3", "brca1", "klf5", "gata4",  "foxd3", "atoh1", "egr2", "zfx", "egr1", "tbp",'znf143','smc3']
    #cofactor_list = cofactor_list[0:2]
    diff_peak_dir = "/home/shi/projects/chipseq_snp/data2/encode/"
    
    ctcf_QTL_file = '/home/shi/projects/wgs/ctcf_qtl.csv'
    #my.f_scp_python_script_to_clustdell("p_pd.py")
    reload(loc)
    loop_file = '/home/shi/projects/wgs/flat_loop_ends.bed'
    loop_region_file = '/home/shi/projects/wgs/loop_regions.bed'
    debug_flag = False

    dipoid_file = '/home/shi/projects/wgs/diploid/gm12878.snp.vcf'


#common path
home_dir = os.path.expanduser('~')
segway_dir = '%s/projects/wgs/segway/' % home_dir
chromhmm_dir = '%s/projects/wgs/chromHMM/' % home_dir
funseq_dir = '%s/projects/wgs/funseq/' % home_dir

def f_get_loc_hct(loc_file, feature_file, output_file):
    cmd="f_get_bed_feature_sorted %s %s > %s" % (loc_file, feature_file ,output_file)
    logging.debug(cmd)
    my.f_call_shell_fun(cmd)


def f_overlap_with_feature_bed(feature_file, loc_file, value_col, value_name, feature_extend_distanace = 0, debug= False):
    #Extend the feature_file by feature_extend_distance, and overlap with loc_file, and get the overlap data
    #loc_file: the bed extracted from the database fiile, point positions
    #value_col: the column of the desired value in the feature file, 0 based
    #value_name: the name of the value_col
    #feature_extend_distance: this is only for the point features. for region features, set it to 0.
    #return: [loc_chr, loc_start, value_col], loc_start is 1 based same as in database file
    #import ipdb; ipdb.set_trace()
    import socket
    server_name = socket.gethostname()
    if debug == True and 'loire' in server_name :
        import ipdb; ipdb.set_trace()
    feature_bed_extended= pd.io.parsers.read_csv(feature_file, header = None, sep='\t')
    bed_data = pybedtools.BedTool(loc_file)
    feature_bed_extended[1] = feature_bed_extended[1].astype(int) - feature_extend_distanace
    feature_bed_extended[2] = feature_bed_extended[2].astype(int) + feature_extend_distanace
    feature_bed_extended[3] = feature_bed_extended[value_col] # Assign the value to column 4, because after intersect, 4+ 6 colums allowed
    logging.debug(feature_bed_extended.head(3))
    
    feature_extend_bed = my.f_pd_to_bed(feature_bed_extended.ix[:,0:4].drop_duplicates())

    logging.info("Unique features %d %d"%feature_bed_extended.ix[:,0:4].drop_duplicates().shape)
    feature_regions=bed_data.intersect(feature_extend_bed,wao=True)

    
    #Here use end (2) as start in bed format. In database, start is 1 based
    tmp_data = my.f_bed_to_pd(feature_regions)
    logging.debug(tmp_data[tmp_data[4].map(str)=="."].values[1:3])

    
    feature_regions_pd=my.f_bed_to_pd(feature_regions).ix[:,[0,2, 7]] # the 7th position is for value_col
    feature_regions_pd.columns=["chr","start",value_name]
    logging.debug(feature_regions_pd[['start', value_name]].values)
    #logging.debug("=====Missing Pvalue======") #Mostly because I extend 100bp
    feature_regions_pd['start'] = feature_regions_pd['start'].astype(float)

    return feature_regions_pd




#### Main####
for loc_tf in tf_list:

    db_file=my.f_create_cell_pair_database_name(data_dir, cell, diff_cell_list, loc_tf)
    tf_database=loc.data_table(db_file)
    bed_file = tf_database.extract_bed(data_dir)
    loc_dir_file=bed_file

    #The funseq sensitive regions.
    
    sensitive_file = funseq_dir + '/funseq_sensitive_regions.bed'
    sensitive_output_file=data_dir + my.f_generate_tmp_file_name("sensitive")
    f_get_loc_hct(loc_dir_file, sensitive_file, sensitive_output_file)
    sensitive_data=tf_database.read_feature_replace_name(sensitive_output_file,["#chr","start","end","name"],["chr","bed_start","start","lineage_senstive"])
    tf_database.merge_feature(sensitive_data[["chr","start","lineage_senstive"]])
    os.remove( sensitive_output_file)
    print 'Shape of the %s %s' % sensitive_data.shape
    print 'The NULL overlap %s' % sensitive_data['lineage_senstive'].isnull().sum()
    continue

    
    #chromhmm state
    chromhmm_file = my.f_unique_element_in_list( my.f_grep_files_from_dir(chromhmm_dir, pattern = 'wgEncodeAwgSegmentationChromhmm%s.*bed' % cell ) )
    chromhmm_output_file=data_dir + my.f_generate_tmp_file_name("chromhmm")
    f_get_loc_hct(loc_dir_file, chromhmm_file, chromhmm_output_file)
    chromhmm_data=tf_database.read_feature_replace_name(chromhmm_output_file,["#chr","start","end","name"],["chr","bed_start","start","chromhmm"])
    print chromhmm_data.head()
    tf_database.merge_feature(chromhmm_data[["chr","start","chromhmm"]])
    os.remove( chromhmm_output_file)

    #continue
    

    #HCT number
    #mport ipdb; ipdb.set_trace()
    hct_list = ['hct_001','hct_0005','hct_0001']
    for hct_loc in hct_list:
        hct_file=my.f_create_file_name(data_dir='%s/%s/'%(hct_dir, hct_loc), cell=cell, tf=loc_tf, suffix="hct.bed")
        if os.path.isfile(hct_file):
            hct_output_file=data_dir + my.f_generate_tmp_file_name("loc.hct")
            logging.debug("hct_file: " + hct_file)
            if os.path.isfile(hct_file):
                f_get_loc_hct(loc_dir_file, hct_file, hct_output_file)
                feature_data=tf_database.read_feature_replace_name(hct_output_file,["#chr","start","end","name"],["chr","bed_start","start",hct_loc])
                tf_database.merge_feature(feature_data[["chr","start",hct_loc]])
                os.remove(hct_output_file)
        else:
            logging.info("Don't find hct file: %s" % hct_file )


            
    #Segway state
    segway_file = my.f_unique_element_in_list( my.f_grep_files_from_dir(segway_dir, pattern = 'segway_%s.*bed' % cell ) )
    segway_output_file=data_dir + my.f_generate_tmp_file_name("segway")
    f_get_loc_hct(loc_dir_file, segway_file, segway_output_file)
    segway_data=tf_database.read_feature_replace_name(segway_output_file,["#chr","start","end","name"],["chr","bed_start","start","segway"])
    tf_database.merge_feature(segway_data[["chr","start","segway"]])
    os.remove(segway_output_file)




    
            
    #Consistant peaks.
    consistent_peak_file = '/home/shi/projects/wgs/consistent_%s.bed'%loc_tf
    if os.path.isfile(consistent_peak_file):
        consistent_peak_pd = f_overlap_with_feature_bed(consistent_peak_file, bed_file, 3, 'overlap_count', feature_extend_distanace = 0)
        tf_database.merge_feature(consistent_peak_pd[['chr','start','overlap_count']])
    #continue
    #add the bed overlap of different labs, lab narrowPeak files related
    
    pattern="(%s).*%s.*%s.narrowPeak.101bp.bed$"%("|".join(labs), cell, loc_tf)
    peak_file_list = my.f_grep_files_from_dir(data_dir, pattern)

    
    logging.info( "Peak 101bp bed list:\n %s " % "\n".join(peak_file_list))
    for peak_file in peak_file_list:
        peak_file_pattern = ".*/(%s).*%s.*%s.narrowPeak.101bp.bed$"%("|".join(labs),cell, loc_tf)
        lab_match = re.match(peak_file_pattern, peak_file)

        if lab_match == None:
            logging.error("Missing lab in %s" % peak_file)

        lab_name = lab_match.group(1)

        print lab_match.group(1)
        lab_col = 'lab_%s'% lab_name
        labs_data = f_overlap_with_feature_bed(peak_file, bed_file, 3, lab_col, feature_extend_distanace = 0, debug = debug_flag)
        tf_database.merge_feature(labs_data[['chr','start', lab_col ]])


        #Peak max position
        peak_output_file=data_dir + my.f_generate_tmp_file_name("loc.peak_d_")
        shell.f_get_loc_peak_max( loc_dir_file, peak_file, peak_output_file)
        
        peak_dis_data=tf_database.read_feature_replace_name(peak_output_file,["#chr","start","end","name"],["chr","bed_start","start", "peak_%s_dis" % lab_name])
        logging.info('After read feature')
        tf_database.merge_feature(peak_dis_data[["chr","start","peak_%s_dis"%lab_name]])
        tf_database.show_size()



        #Peak P values, extract from narrowPeak file, not 101bp.bed
        peak_file_wide = peak_file.replace('.101bp.bed', '')
        peak_pvalue_pd = f_overlap_with_feature_bed(peak_file_wide, bed_file, 6, "peak_%s_pvalue" % lab_name , feature_extend_distanace = 0, debug = debug_flag)
        print peak_pvalue_pd.head()
        tf_database.merge_feature(peak_pvalue_pd[["chr","start","peak_%s_pvalue" % lab_name]])



    #Drop the unuse pwm columns
    ## print tf_database.data.columns
    ## pwm_tf_list = my.grep_list('^(alt|ref)_pwm_.*', tf_database.data.columns)
    ## pwm_cofactor_list = my.grep_list('^(alt|ref)_pwm_(%s)' % '|'.join(cofactor_list), tf_database.data.columns)
    ## drop_colums = list( set(pwm_tf_list) - set(pwm_cofactor_list) )

    drop_colums = my.grep_list('^peak_[u|w]_.*', tf_database.data.columns)
    logging.debug('Drop colums:' + " ".join(drop_colums))
    tf_database.drop_feautre(drop_colums)
    ## print len(drop_colums)

    #Overlap with dsQTL
    dsQTL_regions_pd = f_overlap_with_feature_bed(dsQTL_file, bed_file, 1, "dsQTL_start", feature_extend_distanace = 50)
    tf_database.merge_feature(dsQTL_regions_pd[["chr","start","dsQTL_start"]])



    #Conservation score
    #cons_file= conservation_dir + "/" + "conservation_score.bed"
    cons_output_file=data_dir + my.f_generate_tmp_file_name("loc.cons")
    f_get_loc_hct(loc_dir_file, cons_file, cons_output_file)
    cons_data=tf_database.read_feature_replace_name(cons_output_file,["#chr","start","end","name"],["chr","bed_start","start","cons"])
    tf_database.merge_feature(cons_data[["chr","start","cons"]])




    
    ############True differetial peaks##############
    for other_cell in diff_cell_list:
        print other_cell
        print tf_database.file_path
        if cell == other_cell:
            continue

        if server_name == 'loire' and debug_flag == True:
            import ipdb; ipdb.set_trace()
        #other_cell_bed_file = tf_database.extract_bed(data_dir=data_dir, filter_col = "cell", filter_val = other_cell) #This one only focus on the guest cell
        other_cell_bed_file = tf_database.extract_bed(data_dir=data_dir)
        other_cell_bed_data = pybedtools.BedTool(other_cell_bed_file)

        print other_cell_bed_data.head()
        #sys.exit()
        logging.debug(other_cell_bed_file)
        cmp_pair = "%s-%s"%(cell, other_cell)
        diff_peak_pattern = my.f_create_pattern([cmp_pair], [loc_tf], "%s.diff_peak.bed$"%diffpeak_method)
        logging.debug(diff_peak_pattern)
        diff_peak_files=my.f_grep_files_from_dir(diff_peak_dir, diff_peak_pattern)
        diff_peak_file=my.f_unique_element_in_list(diff_peak_files)
        logging.debug(diff_peak_file)

        diff_peak_data=pybedtools.BedTool(diff_peak_file)
        diff_overlap=other_cell_bed_data.intersect(diff_peak_data, wao=True)

        tmp_data = my.f_bed_to_pd(diff_overlap)
        logging.debug(tmp_data.ix[1:50])

        diff_overlap_pd=tmp_data.ix[:,[0,2,7,8,9,10]]
        logging.debug(diff_overlap_pd)


        diff_cols = [ "diff_flag" , "diff_depth_A", "diff_depth_B", 'diff_pvalue']
        diff_cols_method= ['chr', 'start'] + [ diff_col + '_' + diffpeak_method for diff_col in diff_cols]
        diff_overlap_pd.columns = diff_cols_method
        diff_overlap_pd_subset = diff_overlap_pd[diff_overlap_pd["diff_flag_%s"%diffpeak_method] != "."]
        assert all(diff_overlap_pd_subset['diff_flag_%s'%diffpeak_method] != '.')

        logging.debug('Differential data:\n')
        print diff_overlap_pd_subset.head()

        #Prepare for the mergeing with the total database
        #diff_overlap_pd_subset['cell'] = other_cell
        diff_overlap_pd_subset['start']=diff_overlap_pd_subset['start'].map(int)
        logging.debug(diff_overlap_pd_subset.ix[:,1:10])
        #my.f_print_list(list(tf_database.data.columns))
        diff_pvalue_col = 'diff_pvalue_%s' % diffpeak_method
        tf_database.merge_feature(diff_overlap_pd_subset, expected_cols = ["chr","start"], debug =debug_flag)

        #Since database file usually contains variation locations host and guest cells. This part makes sure that host_data is also updated.
        none_zero_rows = (tf_database.data[ tf_database.data['cell'] == cell][[diff_pvalue_col]].values != '.')
        host_cell_update_data=tf_database.data[ tf_database.data['cell'] == cell][none_zero_rows][['chr','start','het_type',diff_pvalue_col]].values
        logging.debug('Updated host pvalue data\n')
        print host_cell_update_data
        #assert host_cell_update_data.shape[0] > 0, 'Host data are not updated!'


            

    if cell == 'gm12878' or cell == 'gmtest':

        #This diff_peak_file is only base on the bed files of two cellls.
        #This is only for diff_count value
        diff_peak_pd = f_overlap_with_feature_bed(diff_peak_file, bed_file, 3, 'diff_count', feature_extend_distanace = 0, debug = False)
        tf_database.merge_feature(diff_peak_pd[['chr','start','diff_count']])

        #Methylation
        logging.info("======Methylation=========")
        methy_name=my.f_create_pattern( [cell], ['methy'], "bed$")
        methy_file= my.f_unique_element_in_list( my.f_grep_files_from_dir(data_dir, methy_name) )
        logging.info(methy_file)
        bed_data = pybedtools.BedTool(bed_file)
        methy_regions          = bed_data.window(methy_file, w = 100, c = True)

        methy_pd = my.f_bed_to_pd(methy_regions)
        methy_pd.columns=['chr','start','end','name','methy']
        methy_pd['start'] = methy_pd['end'].astype('int')
        tf_database.merge_feature(methy_pd[['chr','start','methy']])



        #overlap with the loops ends
        loops_pd = f_overlap_with_feature_bed(loop_file, bed_file, 3, "loop_overlap", feature_extend_distanace = 0)
        tf_database.merge_feature(loops_pd[["chr","start","loop_overlap"]])


        #overlap with the loop regions
        loops_region = f_overlap_with_feature_bed(loop_region_file, bed_file, 3, "loop_region_overlap", feature_extend_distanace = 0)
        tf_database.merge_feature(loops_region[["chr","start","loop_region_overlap"]])            



        #Overlap with the diploid data
        dipoid_pd = f_overlap_with_feature_bed(dipoid_file, bed_file, 5, "diploid", feature_extend_distanace = 0)
        print dipoid_pd.pivot_table("chr",rows="diploid",aggfunc=len)
        tf_database.merge_feature(dipoid_pd[["chr","start","diploid"]])
    #Overlap with the consistent peaks


    
    if loc_tf == 'ctcf':
        #Overlap with CTCF QTL
        #Change the csv to bed like format.
        #ctcf_QTL = pd.io.parsers.read_csv(ctcf_QTL_file, header =None, sep='\t')
        #ctcf_QTL[0] = 'chr' + ctcf_QTL[0]
        #ctcf_QTL.to_csv(ctcf_QTL_file, sep='\t', index=None, header = None)

        ctcf_QTL_regions_pd = f_overlap_with_feature_bed(ctcf_QTL_file, bed_file, 0, "ctcf_QTL", feature_extend_distanace = 0)
        tf_database.merge_feature(ctcf_QTL_regions_pd[["chr","start","ctcf_QTL"]])
        print ctcf_QTL_regions_pd.pivot_table("chr",rows="ctcf_QTL",aggfunc=len)

        
        

if (server_name != "loire"):
    #my.f_grep_and_rm(data_dir, 'tmp.*')
    #my.f_write_locer_file(locker_file, server_name)
    print 'loire'
else:

    #check the differential p-value
    #print diff_overlap_pd_subset.ix[:,'diff_pvalue_pepr'] != '.'
    print test_envi.new_files()
    test_envi.head_file()

#Transform the pwm score to pvalue.
#pwm_list =[ os.path.basename(pwm_file) for pwm_file in my.f_grep_files_from_dir(snv_dir, ".*-.*-.*pwm_all$")]
   
#[Fri Nov 14 11:35:58 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_113514.add_feature 1 1 1 1 1 1
#[Fri Nov 14 11:47:26 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_114642.add_feature 1 1 1 1 1 1
#[Fri Nov 14 11:54:23 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_115338.add_feature 1 1 1 1 1 1
#[Fri Nov 14 11:57:23 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_115648.add_feature 1 1 1 1 1 1
#[Fri Nov 14 12:01:41 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_120057.add_feature 1 1 1 1 1 1
#[Fri Nov 14 12:06:21 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_120620.add_feature 1 1 1 1 1 1
#[Fri Nov 14 12:07:03 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_120700.add_feature 1 1 1 1 1 1
#[Fri Nov 14 12:09:10 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_120826.add_feature 1 1 1 1 1 1
#[Fri Nov 14 12:16:23 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_121539.add_feature 1 1 1 1 1 1
#[Fri Nov 14 15:25:58 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_152403.add_feature 1 1 1 1 1 1
#[Fri Nov 14 15:35:41 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_153317.add_feature 1 1 1 1 1 1
#[Fri Nov 14 15:57:30 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_155504.add_feature 1 1 1 1 1 1
#[Fri Nov 14 22:45:07 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_224243.add_feature 1 1 1 1 1 1
#[Fri Nov 14 22:53:37 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141114_225123.add_feature 1 1 1 1 1 1
#[Fri Nov 14 23:21:32 2014] p p_run_cluster_sep.py add-feature-gm12878-het-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141114_230818.add_feature 1 1 1 1 1 1
#[Sat Nov 15 10:31:06 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141115_102842.add_feature 1 1 1 1 1 1
#[Sat Nov 15 11:49:26 2014] p p_run_cluster_sep.py add-feature-test-het-shi p_add_feature_on_loc_dp.py helas3 test 20141115_114652.add_feature 1 1 1 1 1 1
#[Sat Nov 15 11:52:09 2014] p p_run_cluster_sep.py add-feature-gm12878-het-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141115_115206.add_feature 1 1 1 1 1 1
#[Sat Nov 15 21:31:05 2014] p p_run_cluster_sep.py add-feature-gm12878-het-shi p_add_feature_on_loc_dp.py gm12878 tbp:usf2:znf143 20141115_212409.add_feature 1 1 1 1 1 1
#[Sat Nov 15 22:08:27 2014] p p_run_cluster_sep.py add-feature-helas3-het-shi p_add_feature_on_loc_dp.py helas3 tbp:usf2:znf143 20141115_220142.add_feature 1 1 1 1 1 1
#[Sun Nov 16 09:41:34 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141116_092650.add_feature 1 1 1 1 1 1
#[Sun Nov 16 10:13:06 2014] p p_run_cluster_sep.py add-feature-gm12891-add_feature-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141116_101305.add_feature 1 1 1 1 1 1
#[Sun Nov 16 10:28:19 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141116_101335.add_feature 1 1 1 1 1 1
#[Sun Nov 16 11:06:15 2014] p p_run_cluster_sep.py add-feature-gm12892-all-shi p_add_feature_on_loc_dp.py gm12892 ctcf 20141116_105109.add_feature 1 1 1 1 1 1
#[Sun Nov 16 11:06:15 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141116_105109.add_feature 1 1 1 1 1 1
#[Sun Nov 16 11:06:15 2014] p p_run_cluster_sep.py add-feature-gm12864-all-shi p_add_feature_on_loc_dp.py gm12864 ctcf 20141116_105109.add_feature 1 1 1 1 1 1
#[Sun Nov 16 11:06:15 2014] p p_run_cluster_sep.py add-feature-gm19239-all-shi p_add_feature_on_loc_dp.py gm19239 ctcf 20141116_105109.add_feature 1 1 1 1 1 1
#[Sun Nov 16 11:06:16 2014] p p_run_cluster_sep.py add-feature-gm19240-all-shi p_add_feature_on_loc_dp.py gm19240 ctcf 20141116_105109.add_feature 1 1 1 1 1 1
#[Sun Nov 16 11:06:16 2014] p p_run_cluster_sep.py add-feature-gm19238-all-shi p_add_feature_on_loc_dp.py gm19238 ctcf 20141116_105109.add_feature 1 1 1 1 1 1
#[Sun Nov 16 12:03:20 2014] p p_run_cluster_sep.py add-feature-gm12864-all-shi p_add_feature_on_loc_dp.py gm12864 ctcf 20141116_114815LWOFC.add_feature 1 1 1 1 1 1
#[Sun Nov 16 12:03:21 2014] p p_run_cluster_sep.py add-feature-gm12892-all-shi p_add_feature_on_loc_dp.py gm12892 ctcf 20141116_114815BEAIB.add_feature 1 1 1 1 1 1
#[Sun Nov 16 12:03:30 2014] p p_run_cluster_sep.py add-feature-gm19239-all-shi p_add_feature_on_loc_dp.py gm19239 ctcf 20141116_114815MVT45.add_feature 1 1 1 1 1 1
#[Sun Nov 16 12:18:10 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141116_114815DAU48.add_feature 1 1 1 1 1 1
#[Sun Nov 16 12:18:15 2014] p p_run_cluster_sep.py add-feature-gm19240-all-shi p_add_feature_on_loc_dp.py gm19240 ctcf 20141116_114815QMYK0.add_feature 1 1 1 1 1 1
#[Sun Nov 16 12:19:00 2014] p p_run_cluster_sep.py add-feature-gm19238-all-shi p_add_feature_on_loc_dp.py gm19238 ctcf 20141116_114815RZOKI.add_feature 1 1 1 1 1 1
#[Sun Nov 16 21:58:18 2014] p p_run_cluster_sep.py add-feature-gm19240-all-shi p_add_feature_on_loc_dp.py gm19240 ctcf 20141116_214233GBQS1.add_feature 1 1 1 1 1 1
#[Sun Nov 16 22:04:08 2014] p p_run_cluster_sep.py add-feature-gm12892-all-shi p_add_feature_on_loc_dp.py gm12892 ctcf 20141116_214233XQ5D4.add_feature 1 1 1 1 1 1
#[Sun Nov 16 22:04:18 2014] p p_run_cluster_sep.py add-feature-gm19239-all-shi p_add_feature_on_loc_dp.py gm19239 ctcf 20141116_2142335X5PE.add_feature 1 1 1 1 1 1
#[Sun Nov 16 22:13:58 2014] p p_run_cluster_sep.py add-feature-gm12864-all-shi p_add_feature_on_loc_dp.py gm12864 ctcf 20141116_214233ZRR7P.add_feature 1 1 1 1 1 1
#[Sun Nov 16 22:19:59 2014] p p_run_cluster_sep.py add-feature-gm19238-all-shi p_add_feature_on_loc_dp.py gm19238 ctcf 20141116_214233LAK6Q.add_feature 1 1 1 1 1 1
#[Sun Nov 16 22:54:21 2014] p p_run_cluster_sep.py add-feature-gm12878-all-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141116_224036FHMB2.add_feature 1 1 1 1 1 1
#[Sun Nov 16 23:24:21 2014] p p_run_cluster_sep.py add-feature-gm12878-all-shi p_add_feature_on_loc_dp.py gm12878 tbp:usf2:znf143 20141116_231727NYPKP.add_feature 1 1 1 1 1 1
#[Sun Nov 16 23:24:34 2014] p p_run_cluster_sep.py add-feature-helas3-all-shi p_add_feature_on_loc_dp.py helas3 tbp:usf2:znf143 20141116_231740VDHZF.add_feature 1 1 1 1 1 1
#[Sun Nov 16 23:45:07 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141116_233012D9WAQ.add_feature 1 1 1 1 1 1
#[Mon Nov 17 12:07:57 2014] p p_run_cluster_sep.py add-feature-helas3-all-shi p_add_feature_on_loc_dp.py helas3 ctcf 20141117_120522P2QIM.add_feature 1 1 1 1 1 1
#[Mon Nov 17 12:31:35 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141117_121635XU87G.add_feature 1 1 1 1 1 1
#[Mon Nov 17 13:54:07 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141117_133911J7GIG.add_feature 1 1 1 1 1 1
#[Mon Nov 17 19:27:42 2014] p p_run_cluster_sep.py add-feature-gm19238-all-shi p_add_feature_on_loc_dp.py gm19238 ctcf 20141117_191357K32M7.add_feature 1 1 1 1 1 1
#[Mon Nov 17 19:30:02 2014] p p_run_cluster_sep.py add-feature-gm12864-all-shi p_add_feature_on_loc_dp.py gm12864 ctcf 20141117_1913578QJ2U.add_feature 1 1 1 1 1 1
#[Mon Nov 17 19:30:03 2014] p p_run_cluster_sep.py add-feature-gm19239-all-shi p_add_feature_on_loc_dp.py gm19239 ctcf 20141117_19135886ZA6.add_feature 1 1 1 1 1 1
#[Mon Nov 17 19:41:03 2014] p p_run_cluster_sep.py add-feature-gm19240-all-shi p_add_feature_on_loc_dp.py gm19240 ctcf 20141117_191358KP17K.add_feature 1 1 1 1 1 1
#[Tue Nov 18 10:08:02 2014] p p_run_cluster_sep.py add-feature-gm12891-all-shi p_add_feature_on_loc_dp.py gm12891 ctcf 20141118_095452ESQL4.add_feature 1 1 1 1 1 1
#[Tue Nov 18 10:08:19 2014] p p_run_cluster_sep.py add-feature-gm12892-all-shi p_add_feature_on_loc_dp.py gm12892 ctcf 20141118_0954597ATKY.add_feature 1 1 1 1 1 1
#[Thu Nov 20 13:57:05 2014] p p_run_cluster_sep.py add-feature-helas3-test-shi p_add_feature_on_loc_dp.py helas3 test 20141120_135700add_featureHPAXX 1 1 1 1 1 1
#[Thu Nov 20 13:57:48 2014] p p_run_cluster_sep.py add-feature-helas3-test-shi p_add_feature_on_loc_dp.py helas3 test 20141120_135744add_featureP2VEE 1 1 1 1 1 1
#[Thu Nov 20 14:00:52 2014] p p_run_cluster_sep.py add-feature-helas3-test-shi p_add_feature_on_loc_dp.py helas3 test 20141120_135826add_featureHWO4X 1 1 1 1 1 1
#[Mon Nov 24 11:45:16 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141124_114514.add_feature.MQ80B 1 1 1 1 1 1
#[Mon Nov 24 12:00:43 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141124_120043.add_feature.6JEDI 1 1 1 1 1 1
#[Mon Nov 24 22:37:17 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141124_223716.add_feature.2QQP7 1 1 1 1 1 1
#[Mon Nov 24 22:56:12 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141124_225611.add_feature.YBBZD 1 1 1 1 1 1
#[Tue Nov 25 11:33:35 2014] p p_run_cluster_sep.py add-feature-gm12878-all-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141125_112900.add_feature.5JO9B 1 1 1 1 1 1
#[Wed Nov 26 16:34:01 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141126_163400.add_feature.PCX1P 1 1 1 1 1 1
#[Mon Dec  1 18:06:45 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141201_180643.add_feature.PK573 1 1 1 1 1 1
#[Mon Dec  1 18:12:07 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141201_181206.add_feature.1ZFAI 1 1 1 1 1 1
#[Mon Dec  1 18:18:52 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141201_181852.add_feature.P2E1T 1 1 1 1 1 1
#[Tue Dec  2 09:16:54 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_091653.add_feature.GTZ21 1 1 1 1 1 1
#[Tue Dec  2 10:15:45 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_101544.add_feature.O7QXR 1 1 1 1 1 1
#[Tue Dec  2 10:16:39 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_101638.add_feature.MBZPH 1 1 1 1 1 1
#[Tue Dec  2 10:17:22 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_101721.add_feature.FZL02 1 1 1 1 1 1
#[Tue Dec  2 10:31:07 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_103106.add_feature.74IIG 1 1 1 1 1 1
#[Tue Dec  2 10:36:36 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_103635.add_feature.3FPAY 1 1 1 1 1 1
#[Tue Dec  2 10:37:26 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_103726.add_feature.0UJSX 1 1 1 1 1 1
#[Tue Dec  2 10:41:58 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_104157.add_feature.FXN1T 1 1 1 1 1 1
#[Tue Dec  2 10:44:24 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_104424.add_feature.S90AY 1 1 1 1 1 1
#[Tue Dec  2 10:50:43 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_105042.add_feature.QREAM 1 1 1 1 1 1
#[Tue Dec  2 10:53:10 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_105309.add_feature.0SR53 1 1 1 1 1 1
#[Tue Dec  2 11:14:46 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_111445.add_feature.22XEJ 1 1 1 1 1 1
#[Tue Dec  2 11:24:22 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_112421.add_feature.Z6XF4 1 1 1 1 1 1
#[Tue Dec  2 11:26:12 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_112611.add_feature.0KMHL 1 1 1 1 1 1
#[Tue Dec  2 11:36:09 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_113608.add_feature.ALSC9 1 1 1 1 1 1
#[Tue Dec  2 11:36:46 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_113645.add_feature.C0F46 1 1 1 1 1 1
#[Tue Dec  2 14:34:19 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_143418.add_feature.4ZUWL 1 1 1 1 1 1
#[Tue Dec  2 14:37:15 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_143714.add_feature.4KO3W 1 1 1 1 1 1
#[Tue Dec  2 14:41:36 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_144135.add_feature.MOUIL 1 1 1 1 1 1
#[Tue Dec  2 14:46:15 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_144614.add_feature.P6NT4 1 1 1 1 1 1
#[Tue Dec  2 15:04:11 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_150410.add_feature.LXA6Y 1 1 1 1 1 1
#[Tue Dec  2 15:04:53 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_150453.add_feature.TMCVP 1 1 1 1 1 1
#[Tue Dec  2 15:07:47 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_150746.add_feature.EMMU4 1 1 1 1 1 1
#[Tue Dec  2 15:12:17 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_151216.add_feature.FB84U 1 1 1 1 1 1
#[Tue Dec  2 15:21:20 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_152119.add_feature.0PYOP 1 1 1 1 1 1
#[Tue Dec  2 15:26:57 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_152656.add_feature.FP9KM 1 1 1 1 1 1
#[Tue Dec  2 15:32:41 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_153240.add_feature.TRV3K 1 1 1 1 1 1
#[Tue Dec  2 15:49:08 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_154907.add_feature.J8BJK 1 1 1 1 1 1
#[Tue Dec  2 16:08:16 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_160815.add_feature.4U6AN 1 1 1 1 1 1
#[Tue Dec  2 16:10:07 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_161006.add_feature.C0A9B 1 1 1 1 1 1
#[Tue Dec  2 16:23:21 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_162320.add_feature.YA4EF 1 1 1 1 1 1
#[Tue Dec  2 16:39:23 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_163922.add_feature.9992G 1 1 1 1 1 1
#[Tue Dec  2 16:54:03 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_165402.add_feature.5WK4Y 1 1 1 1 1 1
#[Tue Dec  2 17:02:44 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_170242.add_feature.ZDR65 1 1 1 1 1 1
#[Tue Dec  2 17:06:24 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141202_170623.add_feature.TOAF5 1 1 1 1 1 1
#[Thu Dec  4 15:49:15 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141204_154915.add_feature.3H0XL 1 1 1 1 1 1
#[Thu Dec  4 16:01:52 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141204_160151.add_feature.BBUFK 1 1 1 1 1 1
#[Thu Dec  4 17:30:41 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141204_173040.add_feature.NPMMH 1 1 1 1 1 1
#[Thu Dec  4 17:35:36 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141204_173535.add_feature.SJYZ9 1 1 1 1 1 1
#[Thu Dec  4 21:00:07 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141204_210005.add_feature.ADE5K 1 1 1 1 1 1
#[Fri Dec  5 22:04:29 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141205_220428.add_feature.6XYJX 1 1 1 1 1 1
#[Fri Dec 12 16:36:51 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141212_163650.add_feature.OYC22 1 1 1 1 1 1
#[Fri Dec 12 17:10:56 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141212_171055.add_feature.RAI8N 1 1 1 1 1 1
#[Fri Dec 12 17:23:29 2014] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20141212_172328.add_feature.FOAU6 1 1 1 1 1 1
#[Thu Jan  8 14:20:57 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150108_142056.add_feature.OAL3X 1 1 1 1 1 1
#[Thu Jan  8 14:23:46 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150108_142346.add_feature.36IVV 1 1 1 1 1 1
#[Thu Jan  8 14:24:58 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150108_142457.add_feature.GAVLU 1 1 1 1 1 1
#[Wed Jan 14 12:28:18 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150114_122817.add_feature.GQDDU 1 1 1 1 1 1
#[Wed Jan 14 12:30:45 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150114_123045.add_feature.BTCVX 1 1 1 1 1 1
#[Wed Jan 14 16:58:07 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150114_165805.add_feature.6MEUP 1 1 1 1 1 1
#[Wed Jan 14 17:01:37 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150114_170136.add_feature.YMZ9W 1 1 1 1 1 1
#[Wed Jan 28 10:16:41 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150128_101640.add_feature.LHV7P 1 1 1 1 1 1
#[Wed Jan 28 10:28:17 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150128_102816.add_feature.MQKXQ 1 1 1 1 1 1
#[Wed Jan 28 10:56:10 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150128_105609.add_feature.XOXGP 1 1 1 1 1 1
#[Wed Jan 28 11:47:21 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150128_114720.add_feature.H4EAL 1 1 1 1 1 1
#[Wed Jan 28 12:08:57 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150128_120856.add_feature.WSAQ8 1 1 1 1 1 1
#[Wed Jan 28 15:31:32 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150128_153131.add_feature.970GC 1 1 1 1 1 1
#[Wed Jan 28 21:15:35 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150128_211534.add_feature.CNCQU 1 1 1 1 1 1
#[Thu Jan 29 15:15:20 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150129_151519.add_feature.96ANB 1 1 1 1 1 1
#[Thu Jan 29 15:15:45 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150129_151544.add_feature.WAQLF 1 1 1 1 1 1
#[Sat Jan 31 18:45:13 2015] p p_run_cluster_sep.py add-feature-test-add_feature-shi p_add_feature_on_loc_dp.py helas3 test 20150131_184512.add_feature.XUSBI 1 1 1 1 1 1
#[Sun Feb  1 11:25:17 2015] p p_run_cluster_sep.py add-feature-gm12878-all-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150201_112102.add_feature.Y6XET 1 1 1 1 1 1
#[Tue Feb  3 11:07:57 2015] p p_run_cluster_sep.py add-feature-ctcf-all-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150203_110342.add_feature.RIZ7C 1 1 1 1 1 1
#[Tue Feb  3 16:18:37 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150203_161836.add_feature.Y7QAH 1 1 1 1 1 1
#[Tue Feb  3 16:37:11 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150203_163710.add_feature.D5IQB 1 1 1 1 1 1
#[Wed Feb  4 15:56:16 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_155615.add_feature.M55Q7 1 1 1 1 1 1
#[Wed Feb  4 15:57:51 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_155750.add_feature.6PYO4 1 1 1 1 1 1
#[Wed Feb  4 15:59:19 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_155918.add_feature.OFDVJ 1 1 1 1 1 1
#[Wed Feb  4 15:59:59 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_155958.add_feature.FLBWT 1 1 1 1 1 1
#[Wed Feb  4 16:05:30 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_160529.add_feature.8UFX4 1 1 1 1 1 1
#[Wed Feb  4 21:45:15 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_214514.add_feature.AS3SA 1 1 1 1 1 1
#[Wed Feb  4 22:00:18 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-all p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_220017.add_feature.CYD4M 1 1 1 1 1 1
#[Wed Feb  4 22:01:32 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_220131.add_feature.6WQDV 1 1 1 1 1 1
#[Wed Feb  4 22:07:50 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_220749.add_feature.YS046 1 1 1 1 1 1
#[Wed Feb  4 22:09:57 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150204_220956.add_feature.FQW2S 1 1 1 1 1 1
#[Thu Feb  5 14:42:14 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150205_144213.add_feature.30FXE 1 1 1 1 1 1
#[Fri Feb  6 14:59:20 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150206_145919.add_feature.Z3109 1 1 1 1 1 1
#[Fri Feb  6 15:04:23 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150206_150422.add_feature.TEE7K 1 1 1 1 1 1
#[Fri Feb  6 16:14:45 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150206_161444.add_feature.L6DUE 1 1 1 1 1 1
#[Mon Feb  9 15:02:39 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150209_150238.add_feature.IP2SU 1 1 1 1 1 1
#[Mon Feb  9 15:05:54 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150209_150554.add_feature.M04G7 1 1 1 1 1 1
#[Wed Feb 11 16:56:27 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150211_165626.add_feature.SFBWA 1 1 1 1 1 1
#[Wed Feb 11 16:59:09 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150211_165909.add_feature.POC58 1 1 1 1 1 1
#[Thu Feb 12 21:41:37 2015] p p_run_cluster_sep.py add-feature-gm12872-add_feature-allq p_add_feature_on_loc_dp.py gm12872 ctcf 20150212_214137.add_feature.3W38I 1 1 1 1 1 1
#[Fri Feb 13 12:13:35 2015] p p_run_cluster_sep.py add-feature-gm12872-add_feature-allq p_add_feature_on_loc_dp.py gm12872 ctcf 20150213_121334.add_feature.VLGAH 1 1 1 1 1 1
#[Fri Feb 13 16:05:06 2015] p p_run_cluster_sep.py add-feature-gm12872-add_feature-allq p_add_feature_on_loc_dp.py gm12872 ctcf 20150213_160505.add_feature.DY3UN 1 1 1 1 1 1
#[Fri Feb 13 16:21:00 2015] p p_run_cluster_sep.py add-feature-gm12872-add_feature-allq p_add_feature_on_loc_dp.py gm12872 ctcf 20150213_162059.add_feature.LQ3Y2 1 1 1 1 1 1
#[Wed Feb 18 11:54:09 2015] p p_run_cluster_sep.py add-feature-ctcf-all-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150218_115003.add_feature.7FVC7 1 1 1 1 1 1
#[Wed Feb 18 14:42:56 2015] p p_run_cluster_sep.py add-feature-gm12878-all-allq p_add_feature_on_loc_dp.py gm12878 znf143 20150218_144051.add_feature.VF1QK 1 1 1 1 1 1
#[Sat Feb 21 10:34:16 2015] p p_run_cluster_sep.py add-feature-gm12878-all-allq p_add_feature_on_loc_dp.py gm12878 rad21 20150221_103222.add_feature.N5H73 1 1 1 1 1 1
#[Sat Feb 21 12:03:27 2015] p p_run_cluster_sep.py add-feature-gm12878-all-allq p_add_feature_on_loc_dp.py gm12878 smc3 20150221_120122.add_feature.QHG6E 1 1 1 1 1 1
#[Tue Feb 24 14:28:03 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-allq p_add_feature_on_loc_dp.py gm12878 smc3 20150224_142801.add_feature.WQBJX 1 1 1 1 1 1
#[Tue Feb 24 14:29:40 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150224_142939.add_feature.1LTE5 1 1 1 1 1 1
#[Tue Feb 24 14:30:24 2015] p p_run_cluster_sep.py add-feature-gm12878-add_feature-allq p_add_feature_on_loc_dp.py gm12878 znf143 20150224_143023.add_feature.Q9E4O 1 1 1 1 1 1
#[Tue Feb 24 17:41:12 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150224_174111.add_feature.903ZW 1 1 1 1 1 1
#[Wed Feb 25 13:36:42 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150225_133641.add_feature.L5UJL 1 1 1 1 1 1
#[Wed Feb 25 13:57:26 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150225_135725.add_feature.IO8WD 1 1 1 1 1 1
#[Thu Feb 26 10:23:26 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150226_102324.add_feature.PE011 1 1 1 1 1 1
#[Wed Mar  4 16:29:54 2015] p p_run_cluster_sep.py add-feature-gm12878-all-allq p_add_feature_on_loc_dp.py gm12878 bhlhe40:ebf1 20150304_162557.add_feature.ZWWVA 1 1 1 1 1 1
#[Mon Mar  9 14:38:29 2015] p p_run_cluster_sep.py add-feature-ctcf-all-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150309_142907.add_feature.2LIR7 1 1 1 1 1 1
#[Mon Mar  9 16:56:08 2015] p p_run_cluster_sep.py add-feature-gm12864-all-allq p_add_feature_on_loc_dp.py gm12864 ctcf 20150309_165354.add_feature.CZJFQ 1 1 1 1 1 1
#[Mon Mar  9 22:41:30 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150309_224128.add_feature.QM77B 1 1 1 1 1 1
#[Tue Mar 10 09:29:33 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150310_092931.add_feature.9GVUD 1 1 1 1 1 1
#[Tue Mar 10 10:18:48 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150310_101847.add_feature.PAFQM 1 1 1 1 1 1
#[Tue Mar 10 14:25:07 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 ctcf 20150310_142506.add_feature.RI6DN 1 1 1 1 1 1
#[Tue Mar 10 15:56:25 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150310_155624.add_feature.ILZ45 1 1 1 1 1
#[Tue Mar 10 15:58:04 2015] p p_run_cluster_sep.py add-feature-ctcf-add_feature-allq p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150310_155803.add_feature.3PSIU 1 1 1 1 1
#[Tue Mar 10 23:13:04 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-allq p_add_feature_on_loc_dp.py gm12878 gm19238 ctcf 20150310_211118.add_feature.KZCZ2 1 1 1 1 1
#[Tue Mar 10 23:13:13 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-allq p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150310_211036.add_feature.2WLCS 1 1 1 1 1
#[Tue Mar 10 23:17:00 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-allq p_add_feature_on_loc_dp.py gm12878 gm19240 ctcf 20150310_211135.add_feature.LYWXP 1 1 1 1 1
#[Wed Mar 11 10:24:16 2015] p p_run_cluster_sep.py add-feature-gm12878-ebf1-allq p_add_feature_on_loc_dp.py gm12878 none bhlhe40:ebf1 20150311_102018.add_feature.ZJX90 1 1 1 1 1
#[Wed Mar 11 11:01:48 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-allq p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150311_103101.add_feature.TT36A 1 1 1 1 1
#[Wed Mar 11 11:05:52 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-allq p_add_feature_on_loc_dp.py gm12878 gm12873 ctcf 20150311_110126.add_feature.1BSPT 1 1 1 1 1
#[Wed Mar 11 21:13:55 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150311_210950.add_feature.YSF5H 1 1 1 1 1
#[Wed Mar 11 21:14:12 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150311_211006.add_feature.IDBJO 1 1 1 1 1
#[Wed Mar 11 21:14:33 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12873 ctcf 20150311_210958.add_feature.UM2ZF 1 1 1 1 1
#[Thu Mar 12 13:58:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-run p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150312_135804.add_feature.MS2SO 1 1 1 1 1
#[Thu Mar 12 15:03:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-run p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150312_150315.add_feature.8MTB1 1 1 1 1 1
#[Thu Mar 12 15:10:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150312_151052.add_feature.XP8TC 1 1 1 1 1
#[Thu Mar 12 16:22:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150312_162233.add_feature.4NYX0 1 1 1 1 1
#[Thu Mar 12 21:35:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150312_213552.add_feature.185V9 1 1 1 1 1
#[Fri Mar 13 09:23:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150313_092344.add_feature.38KRV 1 1 1 1 1
#[Fri Mar 13 09:38:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150313_093836.add_feature.92XME 1 1 1 1 1
#[Fri Mar 13 10:07:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150313_100727.add_feature.ZLURE 1 1 1 1 1
#[Fri Mar 13 10:08:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150313_100844.add_feature.NIEXY 1 1 1 1 1
#[Fri Mar 13 10:09:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150313_100951.add_feature.SZHP8 1 1 1 1 1
#[Fri Mar 13 10:23:57 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150313_101942.add_feature.919QQ 1 1 1 1 1
#[Fri Mar 13 15:35:44 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12873 ctcf 20150313_153118.add_feature.08DBK 1 1 1 1 1
#[Fri Mar 13 15:35:46 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150313_153129.add_feature.SLNGA 1 1 1 1 1
#[Fri Mar 13 22:34:08 2015] p p_run_cluster_sep.py add-feature-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12873 gm12873 ctcf 20150313_222923.add_feature.VYKB2 1 1 1 1 1
#[Sun Mar 15 10:44:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12864-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150315_104403.add_feature.ZJZBL 1 1 1 1 1
#[Mon Mar 16 10:59:35 2015] p p_run_cluster_sep.py add-feature-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12864 gm12864 ctcf 20150316_105519.add_feature.XCXOK 1 1 1 1 1
#[Tue Mar 17 09:59:09 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12801 ctcf 20150317_095514.add_feature.LC2FV 1 1 1 1 1
#[Fri Mar 20 09:36:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150320_093638.add_feature.EIQ6O 1 1 1 1 1
#[Fri Mar 20 09:36:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12801-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12801 ctcf 20150320_093645.add_feature.0AQMK 1 1 1 1 1
#[Fri Mar 20 09:36:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm19238 ctcf 20150320_093656.add_feature.3776G 1 1 1 1 1
#[Fri Mar 20 09:37:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm19239 ctcf 20150320_093703.add_feature.34DYK 1 1 1 1 1
#[Fri Mar 20 09:37:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm19240 ctcf 20150320_093711.add_feature.3N0BF 1 1 1 1 1
#[Fri Mar 20 09:38:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm1891-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm1891 ctcf 20150320_093846.add_feature.85RAQ 1 1 1 1 1
#[Fri Mar 20 09:39:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150320_093908.add_feature.NQZXA 1 1 1 1 1
#[Fri Mar 20 09:39:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150320_093913.add_feature.L2YBS 1 1 1 1 1
#[Fri Mar 20 09:39:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12892 ctcf 20150320_093935.add_feature.Y1Y83 1 1 1 1 1
#[Fri Mar 20 10:07:20 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150320_100324.add_feature.LRFXN 1 1 1 1 1
#[Fri Mar 20 10:07:41 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150320_100346.add_feature.02X6J 1 1 1 1 1
#[Fri Mar 20 10:07:49 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12892 ctcf 20150320_100354.add_feature.OGSXZ 1 1 1 1 1
#[Fri Mar 20 10:08:25 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19238 ctcf 20150320_100420.add_feature.ND6RP 1 1 1 1 1
#[Sat Mar 21 11:12:33 2015] p p_run_cluster_sep.py add-feature-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12801 gm12801 ctcf 20150321_110828.add_feature.YGAHK 1 1 1 1 1
#[Sat Mar 21 11:12:59 2015] p p_run_cluster_sep.py add-feature-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm19238 gm19238 ctcf 20150321_110853.add_feature.8K00V 1 1 1 1 1
#[Sat Mar 21 11:13:19 2015] p p_run_cluster_sep.py add-feature-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm19239 gm19239 ctcf 20150321_110904.add_feature.L1IOS 1 1 1 1 1
#[Sat Mar 21 12:34:43 2015] p p_run_cluster_sep.py add-feature-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm19240 gm19240 ctcf 20150321_123028.add_feature.EUBX3 1 1 1 1 1
#[Tue Mar 24 09:11:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150324_091154.add_feature.ZLOC7 1 1 1 1 1
#[Tue Mar 24 09:34:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150324_093404.add_feature.PR6TM 1 1 1 1 1
#[Tue Mar 24 11:54:22 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19239 ctcf 20150324_115017.add_feature.5P87P 1 1 1 1 1
#[Tue Mar 24 21:55:07 2015] p p_run_cluster_sep.py add-feature-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm19240 gm19240 ctcf 20150324_215101.add_feature.663YQ 1 1 1 1 1
#[Tue Mar 24 21:55:10 2015] p p_run_cluster_sep.py add-feature-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19240 ctcf 20150324_215115.add_feature.ZM4DT 1 1 1 1 1
#[Thu Mar 26 16:44:58 2015] p p_run_cluster_sep.py add-feature-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12892 gm12892 ctcf 20150326_164112.add_feature.2W8OG 1 1 1 1 1
#[Thu Mar 26 16:44:58 2015] p p_run_cluster_sep.py add-feature-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12891 gm12891 ctcf 20150326_164112.add_feature.UFJMP 1 1 1 1 1
#[Thu Mar 26 16:45:19 2015] p p_run_cluster_sep.py add-feature-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm19238 gm19238 ctcf 20150326_164112.add_feature.C99FU 1 1 1 1 1
#[Thu Mar 26 16:45:19 2015] p p_run_cluster_sep.py add-feature-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm19239 gm19239 ctcf 20150326_164112.add_feature.G0TOW 1 1 1 1 1
#[Thu Mar 26 16:45:20 2015] p p_run_cluster_sep.py add-feature-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm19240 gm19240 ctcf 20150326_164112.add_feature.8EW1J 1 1 1 1 1
#[Thu Mar 26 16:45:28 2015] p p_run_cluster_sep.py add-feature-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12864 gm12864 ctcf 20150326_164112.add_feature.VVECR 1 1 1 1 1
#[Sat Mar 28 12:26:16 2015] p p_run_cluster_sep.py add-feature-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12872 gm12872 ctcf 20150328_122150.add_feature.O0KQ8 1 1 1 1 1
#[Fri Apr 10 22:12:47 2015] p p_run_cluster_sep.py add-feature-helas3-test-shi p_add_feature_on_loc_dp.py helas3 gm12878 test 20150410_221010.add_feature.8KCXG 1 1 1 1 1
#[Fri Apr 10 23:22:05 2015] p p_run_cluster_sep.py add-feature-helas3-test-shi p_add_feature_on_loc_dp.py helas3 gm12878 test 20150410_231949.add_feature.Z2R16 1 1 1 1 1
#[Sat Apr 11 10:33:18 2015] p p_run_cluster_sep.py add-feature-helas3-test-shi p_add_feature_on_loc_dp.py helas3 gm12878 test 20150411_103051.add_feature.ASVYO 1 1 1 1 1
#[Sat Apr 11 11:22:35 2015] p p_run_cluster_sep.py add-feature-helas3-test-shi p_add_feature_on_loc_dp.py helas3 gm12878 test 20150411_112009.add_feature.WC5LX 1 1 1 1 1
#[Sat Apr 11 11:35:38 2015] p p_run_cluster_sep.py add-feature-gm12878-PU1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 PU1 20150411_113142.add_feature.F893S 1 1 1 1 1
#[Sat Apr 11 15:30:44 2015] p p_run_cluster_sep.py add-feature-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150411_152638.add_feature.3JZ05 1 1 1 1 1
#[Sun Apr 12 18:24:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150412_182423.add_feature.V2UT1 1 1 1 1 1
#[Sun Apr 12 18:27:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150412_182747.add_feature.MEWR0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 11:13:33 2015] p p_run_cluster_sep.py add-feature-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150413_110926.add_feature.54QV1 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 11:25:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150413_112554.add_feature.2HCUQ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 11:42:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-add_feature-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150413_114229.add_feature.YMIIA sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 12:20:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150413_121624.add_feature.I53NX sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 14:06:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150413_140201.add_feature.KNZBD sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:01:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150413_155700.add_feature.8WXW0 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:03:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150413_155700.add_feature.PPJP4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:03:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19238 pu1 20150413_155700.add_feature.Z4IUL sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:32:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150413_162844.add_feature.G81CG sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:45:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150413_164057.add_feature.CL5TW sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:45:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150413_164057.add_feature.XM7AI sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 16:45:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12813 pu1 20150413_164057.add_feature.GIJ9J sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:08:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12892 pu1 20150413_170455.add_feature.5ALWY sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:09:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150413_170455.add_feature.0T8CJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:09:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11831 pu1 20150413_170455.add_feature.NZN87 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:09:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150413_170455.add_feature.6QIEX sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:09:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150413_170455.add_feature.2TDUJ sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:09:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150413_170455.add_feature.E9MCB sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 17:09:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150413_170455.add_feature.T3P53 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 21:14:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150413_210959.add_feature.BZQ3C sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 21:23:32 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150413_211855.add_feature.GE9FV sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:38:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12892 pu1 20150413_223455.add_feature.84ZFF sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:39:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150413_223455.add_feature.YTPQN sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:39:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12891 pu1 20150413_223504.add_feature.2S22K sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:39:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11831 pu1 20150413_223504.add_feature.UQ69G sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:39:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150413_223504.add_feature.GB9Z4 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:39:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12813 pu1 20150413_223505.add_feature.9BQ36 sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:39:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150413_223504.add_feature.HISSH sydh:uw:haib:uta:embl 1 1 1 1
#[Mon Apr 13 22:39:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19238 pu1 20150413_223504.add_feature.ELLNK sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 10:14:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150414_101018.add_feature.GMOBB sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 10:16:27 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150414_101151.add_feature.Z3EWY sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:07:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150414_170251.add_feature.XSEGL sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:30:14 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150414_172558.add_feature.G82RY sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:34:36 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150414_173006.add_feature.JDB73 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:34:54 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150414_173006.add_feature.1QY86 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:38:24 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150414_173006.add_feature.IVLEU sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:38:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150414_173006.add_feature.Z7C44 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:40:04 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150414_173006.add_feature.11M2E sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 17:42:04 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150414_173006.add_feature.HU7LG sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:35:11 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150414_223045.add_feature.QMY0N sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:36:21 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150414_223233.add_feature.JIBXQ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:36:32 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150414_223233.add_feature.V5GKZ sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:36:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150414_223233.add_feature.CH3MR sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:36:44 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150414_223233.add_feature.QT5ZR sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:36:51 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150414_223233.add_feature.8IFLT sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:36:51 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150414_223233.add_feature.W54FT sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:37:01 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150414_223233.add_feature.GIX54 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:37:02 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150414_223233.add_feature.1NAL8 sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:37:02 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150414_223233.add_feature.ORUOM sydh:uw:haib:uta:embl 1 1 1 1
#[Tue Apr 14 22:37:11 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150414_223233.add_feature.ZWC3O sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:20:16 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150415_141559.add_feature.O994Y sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150415_162312.add_feature.17FOQ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:09 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150415_162312.add_feature.36NBP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:09 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150415_162312.add_feature.TMROP sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:29 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150415_162312.add_feature.J0C85 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:30 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150415_162313.add_feature.WK5BD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:39 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150415_162312.add_feature.WH3FI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:49 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150415_162312.add_feature.4JOX4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 16:27:51 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150415_162312.add_feature.42EV7 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:20:54 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150415_171638.add_feature.88PDX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:55:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150415_172904.add_feature.3F7HO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:07 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150415_174041.add_feature.PE7W2 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:08 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150415_173935.add_feature.H387M sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:16 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150415_172905.add_feature.OIKR4 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:16 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150415_172904.add_feature.G88PJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:16 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150415_172904.add_feature.23X72 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:27 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150415_172904.add_feature.00P4N sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:27 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150415_173142.add_feature.GI996 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 17:56:28 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150415_172904.add_feature.BVBEM sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:12 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150415_210714.add_feature.DWQ2O sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:12 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150415_210713.add_feature.9PBD7 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150415_210713.add_feature.R4HZ8 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:43 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150415_210714.add_feature.VCJ2X sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:43 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150415_210714.add_feature.ESSJK sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:43 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150415_210713.add_feature.F719F sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:43 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150415_210713.add_feature.B4KZE sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:52 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150415_210713.add_feature.BZ852 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:11:52 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150415_210713.add_feature.6AIVI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:12:02 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150415_210713.add_feature.QYMCI sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:33:26 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-myc-shi p_add_feature_on_loc_dp.py gm12892 gm12892 myc 20150415_212939.add_feature.W5LP0 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:33:46 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-myc-shi p_add_feature_on_loc_dp.py gm19239 gm19239 myc 20150415_212939.add_feature.SPBUK sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:33:46 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-myc-shi p_add_feature_on_loc_dp.py gm19240 gm19240 myc 20150415_212939.add_feature.VPEP9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:33:56 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-myc-shi p_add_feature_on_loc_dp.py gm19238 gm19238 myc 20150415_212939.add_feature.OXED3 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:51:14 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-myc-shi p_add_feature_on_loc_dp.py gm12891 gm12891 myc 20150415_214727.add_feature.1WSJJ sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 21:51:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150415_214739.add_feature.FSOGA sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12891 pu1 20150415_224634.add_feature.6OEV2 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150415_224633.add_feature.FMGC9 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150415_224633.add_feature.R0H2W sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150415_224633.add_feature.7IKJO sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12813 pu1 20150415_224633.add_feature.MP4EU sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11831 pu1 20150415_224633.add_feature.VJD88 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12892 pu1 20150415_224633.add_feature.4TQDD sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150415_224633.add_feature.TCQHX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150415_224633.add_feature.CTTGX sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19238 pu1 20150415_224634.add_feature.QNPJ0 sydh:uw:haib:uta:embl 1 1 1 1
#[Wed Apr 15 22:46:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150415_224634.add_feature.D2THC sydh:uw:haib:uta:embl 1 1 1 1
#[Thu Apr 16 09:48:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11831 pu1 20150416_094759.add_feature.14ULY sydh:uw:haib:uta:embl 1 1 1 1
#[Sat Apr 18 18:39:47 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150418_183529.add_feature.46N33 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:43:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150418_183529.add_feature.ZC2FD sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:47:58 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150418_183529.add_feature.881DQ sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:48:08 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150418_183530.add_feature.L8KLV sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:48:28 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150418_183530.add_feature.K0XA5 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:48:28 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150418_183530.add_feature.QORCE sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:52:09 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150418_183530.add_feature.J6QW2 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:52:10 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150418_183530.add_feature.FRE28 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:52:28 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150418_183530.add_feature.VYX3F sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:52:29 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150418_183530.add_feature.YJ9ZW sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:56:09 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150418_183530.add_feature.AUFA0 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:56:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150418_183544.add_feature.JE3BC sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:56:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12892 pu1 20150418_183545.add_feature.UZPRV sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:56:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11831 pu1 20150418_183545.add_feature.D72HI sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 18:59:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12891 pu1 20150418_183545.add_feature.SLOQX sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 19:00:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150418_183545.add_feature.ZXKQA sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 19:00:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150418_183545.add_feature.JYDCC sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 19:00:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150418_183545.add_feature.IWUND sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 19:03:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150418_183546.add_feature.3DOOT sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 19:03:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19238 pu1 20150418_183545.add_feature.8OM1R sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 19:04:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150418_183546.add_feature.MVTE6 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sat Apr 18 19:04:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12813 pu1 20150418_183546.add_feature.LEOIS sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 20 09:55:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150420_095145.add_feature.743KE sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 09:55:43 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-myc-shi p_add_feature_on_loc_dp.py gm12891 gm12891 myc 20150420_095145.add_feature.UI3YI sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 09:55:53 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-myc-shi p_add_feature_on_loc_dp.py gm19239 gm19239 myc 20150420_095145.add_feature.WELTI sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 09:55:53 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-myc-shi p_add_feature_on_loc_dp.py gm12892 gm12892 myc 20150420_095145.add_feature.EB77K sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 09:56:03 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-myc-shi p_add_feature_on_loc_dp.py gm19238 gm19238 myc 20150420_095145.add_feature.L2KG4 sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 09:56:03 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-myc-shi p_add_feature_on_loc_dp.py gm19240 gm19240 myc 20150420_095146.add_feature.E5BL9 sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 09:56:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150420_095219.add_feature.DUD6F sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 09:56:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19239 myc 20150420_095219.add_feature.8PZR5 sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 12:01:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19240 myc 20150420_115747.add_feature.0MEZK sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 12:01:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19238 myc 20150420_115747.add_feature.R69K2 sydh:uw:haib:uta:embl macs 1 1 1
#[Mon Apr 20 13:19:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19239 myc 20150420_131901.add_feature.P22XF sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 20 13:19:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19238 myc 20150420_131914.add_feature.KH1P3 sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 20 13:19:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19240 myc 20150420_131923.add_feature.42BCJ sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:21:56 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150423_122155.add_feature.KUBLR sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:30:49 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150423_123048.add_feature.P9UM0 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:30:54 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150423_123053.add_feature.9UJF4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:30:59 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150423_123058.add_feature.O102M sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:04 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150423_123103.add_feature.JV894 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:09 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150423_123108.add_feature.WI0S8 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150423_123113.add_feature.KX8VM sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:19 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150423_123118.add_feature.V0U6B sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:25 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150423_123124.add_feature.IQ5IV sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:29 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150423_123128.add_feature.90S56 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:34 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150423_123133.add_feature.BIY3Q sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:41 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150423_123140.add_feature.I1MEP sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12813 pu1 20150423_123150.add_feature.TG8PI sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150423_123149.add_feature.GYOSQ sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150423_123150.add_feature.RYN8R sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11831 pu1 20150423_123149.add_feature.44GFB sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150423_123150.add_feature.XOWOX sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150423_123150.add_feature.HFVLR sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12891 pu1 20150423_123150.add_feature.UV3OH sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150423_123150.add_feature.6KQCG sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19238 pu1 20150423_123150.add_feature.8UBEG sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12892 pu1 20150423_123150.add_feature.HZ6YN sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 12:31:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150423_123150.add_feature.J7N28 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 23 13:45:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150423_134545.add_feature.TH51Y sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:43:46 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150424_144345.add_feature.VNKGI sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:43:49 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150424_144348.add_feature.FI06Q sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:43:54 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150424_144353.add_feature.XWQXM sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:43:59 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150424_144358.add_feature.BT5PU sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:44:04 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150424_144403.add_feature.0TVEN sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:44:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150424_144408.add_feature.PERL7 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:44:14 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150424_144413.add_feature.WDUGQ sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:44:19 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150424_144418.add_feature.X5APS sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:44:24 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150424_144423.add_feature.46QIH sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:44:29 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150424_144428.add_feature.IUU99 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:44:34 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150424_144433.add_feature.AGD51 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:45:43 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150424_144542.add_feature.1YKD2 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:45:48 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150424_144547.add_feature.1X4HE sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:45:52 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150424_144552.add_feature.AKA1E sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:45:58 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150424_144557.add_feature.T3TFK sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:46:03 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150424_144602.add_feature.94CMT sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:46:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150424_144607.add_feature.7ZUU4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:46:13 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150424_144612.add_feature.P7QVJ sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:46:18 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150424_144617.add_feature.THRM6 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:46:23 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150424_144622.add_feature.9N8TE sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:46:28 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150424_144627.add_feature.2K3ZQ sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 14:46:33 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150424_144632.add_feature.UJYQ1 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:15:52 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150424_151551.add_feature.Q0XWC sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:15:56 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150424_151556.add_feature.J3RQ5 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:01 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150424_151601.add_feature.K02OM sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:06 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150424_151606.add_feature.HHZVP sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:11 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150424_151611.add_feature.M86DF sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150424_151616.add_feature.8AR7H sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:22 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150424_151621.add_feature.UO7IY sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:27 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150424_151626.add_feature.T01WJ sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:32 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150424_151631.add_feature.IRJX8 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:16:37 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150424_151636.add_feature.XU0EG sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:23 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150424_152722.add_feature.QG2J5 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11831 pu1 20150424_152732.add_feature.OGHME sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11830 pu1 20150424_152732.add_feature.9TH9F sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150424_152732.add_feature.FNV0B sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19238 pu1 20150424_152732.add_feature.KZV4S sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12891 pu1 20150424_152732.add_feature.XX56J sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12813 pu1 20150424_152732.add_feature.D9DTY sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm11894 pu1 20150424_152732.add_feature.7BWG1 sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12892 pu1 20150424_152732.add_feature.1C05P sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19239 pu1 20150424_152732.add_feature.LMXOK sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12776 pu1 20150424_152732.add_feature.P7R6X sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 15:27:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150424_152733.add_feature.32OBA sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 16:22:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm19240 pu1 20150424_162202.add_feature.UIE5C sydh:uw:haib:uta:embl pepr 1 1 1
#[Fri Apr 24 16:22:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12891 pu1 20150424_162229.add_feature.VY2GX sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:43:51 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11831 pu1 20150426_223934.add_feature.DT0ZP sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:43:51 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12892 pu1 20150426_224004.add_feature.RCCV8 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:43:57 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12891 pu1 20150426_223959.add_feature.CLYA9 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:44:02 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12878 pu1 20150426_223954.add_feature.5IZP4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:44:16 2015] p p_run_cluster_sep.py add-feature-gm11830-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm12813 pu1 20150426_223949.add_feature.TEHHY sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:44:37 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19238 pu1 20150426_224009.add_feature.BT13H sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:44:42 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19239 pu1 20150426_224014.add_feature.6BHZ4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Sun Apr 26 22:44:53 2015] p p_run_cluster_sep.py add-feature-gm11894-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm19240 pu1 20150426_224026.add_feature.4LSK4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 27 07:35:04 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11830 pu1 20150427_073027.add_feature.QHW4D sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 27 07:51:41 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150427_074714.add_feature.FQ300 sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 27 08:05:13 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11831 pu1 20150427_080047.add_feature.I76NQ sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 27 11:27:46 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11830 pu1 20150427_112744.add_feature.YHSX4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 27 11:32:33 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11830 pu1 20150427_113232.add_feature.B7Y8C sydh:uw:haib:uta:embl pepr 1 1 1
#[Mon Apr 27 11:51:08 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150427_115106.add_feature.H0CTK sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:14:26 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11830 pu1 20150429_171425.add_feature.T4R21 sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:14:31 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11831 pu1 20150429_171430.add_feature.76AQ4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:14:36 2015] p p_run_cluster_sep.py add-feature-gm11840-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11840 gm11894 pu1 20150429_171435.add_feature.8PIOY sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:14:41 2015] p p_run_cluster_sep.py add-feature-gm11894-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm12776 pu1 20150429_171440.add_feature.KFQOK sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:14:46 2015] p p_run_cluster_sep.py add-feature-gm11830-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm12813 pu1 20150429_171445.add_feature.JJ3S4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:14:51 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12878 pu1 20150429_171450.add_feature.I1L0Z sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:14:56 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12891 pu1 20150429_171455.add_feature.QRKOH sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:15:01 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12892 pu1 20150429_171500.add_feature.P8GXY sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:15:06 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19238 pu1 20150429_171505.add_feature.D55WP sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:15:11 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19239 pu1 20150429_171510.add_feature.GU5T9 sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:15:17 2015] p p_run_cluster_sep.py add-feature-gm11894-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm19240 pu1 20150429_171516.add_feature.PZYFD sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:17:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19239 myc 20150429_171754.add_feature.D8W4B sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:17:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19238 myc 20150429_171759.add_feature.QETAV sydh:uw:haib:uta:embl pepr 1 1 1
#[Wed Apr 29 17:18:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-myc-shi p_add_feature_on_loc_dp.py gm12878 gm19240 myc 20150429_171805.add_feature.KOYUE sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:15:33 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11830 pu1 20150430_091531.add_feature.CFS85 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:15:38 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11831 pu1 20150430_091536.add_feature.UXJFV sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:15:42 2015] p p_run_cluster_sep.py add-feature-gm11840-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11840 gm11894 pu1 20150430_091541.add_feature.NO8TX sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:15:47 2015] p p_run_cluster_sep.py add-feature-gm11894-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm12776 pu1 20150430_091546.add_feature.OLTRU sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:15:52 2015] p p_run_cluster_sep.py add-feature-gm11830-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm12813 pu1 20150430_091551.add_feature.I1UN5 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:15:57 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12878 pu1 20150430_091556.add_feature.FLQ3U sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:16:02 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12891 pu1 20150430_091601.add_feature.WYQGX sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:16:08 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12892 pu1 20150430_091607.add_feature.35IYB sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:16:13 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19238 pu1 20150430_091612.add_feature.TBT0T sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:16:18 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19239 pu1 20150430_091617.add_feature.S0XD8 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 09:33:08 2015] p p_run_cluster_sep.py add-feature-gm11894-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm19240 pu1 20150430_093307.add_feature.JUUD4 sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu Apr 30 10:05:25 2015] p p_run_cluster_sep.py add-feature-gm11894-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm12776 pu1 20150430_100058.add_feature.Y49YW sydh:uw:haib:uta:embl pepr 1 1 1
#[Thu May  7 12:16:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150507_121637.add_feature.QTS2J sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:27:46 2015] p p_run_cluster_sep.py add-feature-gm12801-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12801 gm12801 ctcf 20150507_122745.add_feature.59NM2 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:27:51 2015] p p_run_cluster_sep.py add-feature-gm12864-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12864 gm12864 ctcf 20150507_122750.add_feature.1KQC1 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:27:56 2015] p p_run_cluster_sep.py add-feature-gm12872-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12872 gm12872 ctcf 20150507_122755.add_feature.MG4FE sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:28:01 2015] p p_run_cluster_sep.py add-feature-gm12873-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12873 gm12873 ctcf 20150507_122800.add_feature.P2FI2 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:28:06 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12891 gm12891 ctcf 20150507_122805.add_feature.T5STR sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:28:11 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12892 gm12892 ctcf 20150507_122810.add_feature.HJRUX sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:28:16 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm19238 gm19238 ctcf 20150507_122815.add_feature.JOYPC sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:28:21 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm19239 gm19239 ctcf 20150507_122820.add_feature.ESLZW sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 12:28:26 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm19240 gm19240 ctcf 20150507_122825.add_feature.HXRAB sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:38:48 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm19238 gm19238 ctcf 20150507_133846.add_feature.JCPEL sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:20 2015] p p_run_cluster_sep.py add-feature-gm12801-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12801 gm12801 ctcf 20150507_134119.add_feature.G53RH sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:25 2015] p p_run_cluster_sep.py add-feature-gm12864-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12864 gm12864 ctcf 20150507_134124.add_feature.PZV7M sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:30 2015] p p_run_cluster_sep.py add-feature-gm12872-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12872 gm12872 ctcf 20150507_134129.add_feature.7VWAZ sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:35 2015] p p_run_cluster_sep.py add-feature-gm12873-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12873 gm12873 ctcf 20150507_134134.add_feature.72SOF sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150507_134139.add_feature.M8L4J sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:45 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12891 gm12891 ctcf 20150507_134144.add_feature.ZDX8D sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:50 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12892 gm12892 ctcf 20150507_134149.add_feature.RVFFA sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:41:55 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm19238 gm19238 ctcf 20150507_134155.add_feature.LDHSN sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:42:00 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm19239 gm19239 ctcf 20150507_134159.add_feature.98GKG sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:42:05 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm19240 gm19240 ctcf 20150507_134204.add_feature.D4ROO sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 13:45:02 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150507_134501.add_feature.AWP44 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:08 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150507_134507.add_feature.B8IN2 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:12 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150507_134512.add_feature.6VWMX sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:17 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150507_134517.add_feature.RJ84D sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:22 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150507_134522.add_feature.3PJFF sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150507_134527.add_feature.RIQFP sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:32 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150507_134532.add_feature.LWXP2 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:38 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150507_134537.add_feature.QUM48 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:43 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150507_134542.add_feature.LP5XY sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 13:45:48 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150507_134547.add_feature.X4SMA sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:34:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi-small p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150507_153455.add_feature.T2NTL sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:09 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi-small p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150507_153708.add_feature.32XPZ sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:14 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi-small p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150507_153713.add_feature.M4OY8 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:19 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi-small p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150507_153718.add_feature.IEQYT sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:24 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi-small p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150507_153723.add_feature.VJZY4 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:29 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi-small p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150507_153728.add_feature.42R73 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:34 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi-small p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150507_153733.add_feature.XL0NQ sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:39 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi-small p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150507_153738.add_feature.BXQNO sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:44 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi-small p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150507_153743.add_feature.11BSF sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:49 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi-small p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150507_153748.add_feature.627SS sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 15:37:54 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi-small p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150507_153753.add_feature.VCQM2 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:23:29 2015] p p_run_cluster_sep.py add-feature-gm11830-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm11830 pu1 20150507_162328.add_feature.AQO4W sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:23:34 2015] p p_run_cluster_sep.py add-feature-gm11831-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm11831 gm11831 pu1 20150507_162333.add_feature.IJ526 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:23:39 2015] p p_run_cluster_sep.py add-feature-gm11894-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm11894 pu1 20150507_162338.add_feature.5UIIU sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:23:44 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12776 pu1 20150507_162343.add_feature.3V3VY sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:23:49 2015] p p_run_cluster_sep.py add-feature-gm12813-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm12813 pu1 20150507_162348.add_feature.46WMS sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:23:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150507_162353.add_feature.Z6FAP sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:23:59 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12891 pu1 20150507_162358.add_feature.ENYC0 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:24:04 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12892 pu1 20150507_162403.add_feature.6XKUB sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:24:09 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm19238 gm19238 pu1 20150507_162408.add_feature.QJVSH sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 16:24:14 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm19239 gm19239 pu1 20150507_162413.add_feature.5EXOI sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:33:31 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm19240 gm19240 pu1 20150507_173330.add_feature.C8WLI sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:33:43 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11830-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11830 pu1 20150507_173342.add_feature.OHZJ0 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:33:48 2015] p p_run_cluster_sep.py add-feature-gm12813-gm11831-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm11831 pu1 20150507_173347.add_feature.GA9XE sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:33:53 2015] p p_run_cluster_sep.py add-feature-gm11840-gm11894-pu1-shi p_add_feature_on_loc_dp.py gm11840 gm11894 pu1 20150507_173352.add_feature.ZI0TM sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:33:58 2015] p p_run_cluster_sep.py add-feature-gm11894-gm12776-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm12776 pu1 20150507_173357.add_feature.YU830 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:03 2015] p p_run_cluster_sep.py add-feature-gm11830-gm12813-pu1-shi p_add_feature_on_loc_dp.py gm11830 gm12813 pu1 20150507_173402.add_feature.G07WW sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:08 2015] p p_run_cluster_sep.py add-feature-gm12776-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12776 gm12878 pu1 20150507_173407.add_feature.NL3YI sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:13 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12891-pu1-shi p_add_feature_on_loc_dp.py gm12892 gm12891 pu1 20150507_173412.add_feature.QWM0A sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:18 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12892-pu1-shi p_add_feature_on_loc_dp.py gm12891 gm12892 pu1 20150507_173417.add_feature.5ALBR sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:23 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19238-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19238 pu1 20150507_173422.add_feature.JVHVX sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:28 2015] p p_run_cluster_sep.py add-feature-gm12813-gm19239-pu1-shi p_add_feature_on_loc_dp.py gm12813 gm19239 pu1 20150507_173427.add_feature.YQ44T sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:42 2015] p p_run_cluster_sep.py add-feature-gm11894-gm19240-pu1-shi p_add_feature_on_loc_dp.py gm11894 gm19240 pu1 20150507_173442.add_feature.BVSWZ sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May  7 17:34:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12801 ctcf 20150507_173449.add_feature.RUK0X sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:34:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150507_173454.add_feature.RBW8V sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150507_173459.add_feature.YE522 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12873 ctcf 20150507_173504.add_feature.YOYJM sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150507_173509.add_feature.4OMAE sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150507_173514.add_feature.GSAWK sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12892 ctcf 20150507_173519.add_feature.OOEXT sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19238 ctcf 20150507_173524.add_feature.5WPXS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19239 ctcf 20150507_173529.add_feature.198U2 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 17:35:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19240 ctcf 20150507_173534.add_feature.2DSGB sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12801 ctcf 20150507_221011.add_feature.BYUZB sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150507_221016.add_feature.LZUHD sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150507_221021.add_feature.UXUQ0 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12873 ctcf 20150507_221026.add_feature.VUHKY sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150507_221031.add_feature.W19SC sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150507_221036.add_feature.YL34S sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12892 ctcf 20150507_221041.add_feature.BGWRS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19238 ctcf 20150507_221046.add_feature.EW37X sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19239 ctcf 20150507_221051.add_feature.RENTC sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Thu May  7 22:10:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19240 ctcf 20150507_221056.add_feature.DIJ1W sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:06:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12801 ctcf 20150508_160638.add_feature.7U4HE sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:06:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150508_160643.add_feature.1YGX9 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:06:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150508_160648.add_feature.6NH6D sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:06:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12873 ctcf 20150508_160653.add_feature.2G23O sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:06:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150508_160658.add_feature.XE4SS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150508_160703.add_feature.EBFQ0 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12892 ctcf 20150508_160708.add_feature.I33X7 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19238 ctcf 20150508_160713.add_feature.QY72Q sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19239 ctcf 20150508_160718.add_feature.AKLVD sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19240 ctcf 20150508_160723.add_feature.WNB8K sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:29 2015] p p_run_cluster_sep.py add-feature-gm12801-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12801 gm12801 ctcf 20150508_160728.add_feature.7TD91 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:34 2015] p p_run_cluster_sep.py add-feature-gm12864-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12864 gm12864 ctcf 20150508_160733.add_feature.G88L4 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:39 2015] p p_run_cluster_sep.py add-feature-gm12872-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12872 gm12872 ctcf 20150508_160738.add_feature.14OIR sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:44 2015] p p_run_cluster_sep.py add-feature-gm12873-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12873 gm12873 ctcf 20150508_160743.add_feature.HZTEM sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150508_160748.add_feature.V4ESI sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:54 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12891 gm12891 ctcf 20150508_160753.add_feature.CZAGS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:07:59 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12892 gm12892 ctcf 20150508_160758.add_feature.RQSUG sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:08:04 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm19238 gm19238 ctcf 20150508_160804.add_feature.IDSUD sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:08:09 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm19239 gm19239 ctcf 20150508_160808.add_feature.BPLHS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:08:40 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm19240 gm19240 ctcf 20150508_160839.add_feature.PLG1P sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12801 ctcf 20150508_161005.add_feature.6C1P1 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12864 ctcf 20150508_161010.add_feature.S19R6 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12872 ctcf 20150508_161015.add_feature.0VZS1 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12873 ctcf 20150508_161020.add_feature.UDODH sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150508_161025.add_feature.H7KEV sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12891 ctcf 20150508_161030.add_feature.N3195 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12892 ctcf 20150508_161035.add_feature.4G8QP sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19238 ctcf 20150508_161040.add_feature.E0DGM sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19239 ctcf 20150508_161045.add_feature.JCNF7 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm19240 ctcf 20150508_161050.add_feature.EULWG sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:10:56 2015] p p_run_cluster_sep.py add-feature-gm12801-gm12801-ctcf-shi p_add_feature_on_loc_dp.py gm12801 gm12801 ctcf 20150508_161055.add_feature.Y5AT0 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:01 2015] p p_run_cluster_sep.py add-feature-gm12864-gm12864-ctcf-shi p_add_feature_on_loc_dp.py gm12864 gm12864 ctcf 20150508_161100.add_feature.QSY4Y sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:06 2015] p p_run_cluster_sep.py add-feature-gm12872-gm12872-ctcf-shi p_add_feature_on_loc_dp.py gm12872 gm12872 ctcf 20150508_161105.add_feature.7WKJG sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:11 2015] p p_run_cluster_sep.py add-feature-gm12873-gm12873-ctcf-shi p_add_feature_on_loc_dp.py gm12873 gm12873 ctcf 20150508_161110.add_feature.NT1VT sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150508_161115.add_feature.BOG0S sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:21 2015] p p_run_cluster_sep.py add-feature-gm12891-gm12891-ctcf-shi p_add_feature_on_loc_dp.py gm12891 gm12891 ctcf 20150508_161120.add_feature.JKVQK sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:26 2015] p p_run_cluster_sep.py add-feature-gm12892-gm12892-ctcf-shi p_add_feature_on_loc_dp.py gm12892 gm12892 ctcf 20150508_161125.add_feature.LCBVO sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:31 2015] p p_run_cluster_sep.py add-feature-gm19238-gm19238-ctcf-shi p_add_feature_on_loc_dp.py gm19238 gm19238 ctcf 20150508_161130.add_feature.PAMEK sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:36 2015] p p_run_cluster_sep.py add-feature-gm19239-gm19239-ctcf-shi p_add_feature_on_loc_dp.py gm19239 gm19239 ctcf 20150508_161135.add_feature.OS6RY sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Fri May  8 16:11:41 2015] p p_run_cluster_sep.py add-feature-gm19240-gm19240-ctcf-shi p_add_feature_on_loc_dp.py gm19240 gm19240 ctcf 20150508_161140.add_feature.XWLDH sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:MEF2A 1 1
#[Mon May 11 11:22:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150511_111833.add_feature.J9RGX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon May 11 12:13:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150511_121311.add_feature.LJ5SY sydh:uw:haib:uta:embl pepr None 1 1
#[Mon May 11 12:13:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150511_121317.add_feature.LEIAR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon May 11 12:18:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150511_121426.add_feature.JRC8T sydh:uw:haib:uta:embl pepr None 1 1
#[Mon May 11 12:18:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150511_121501.add_feature.6749D sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 15 18:34:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cdp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cdp 20150515_183053.add_feature.A9KWY sydh:uw:haib:uta:embl pepr None 1 1
#[Sat May 16 10:19:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150516_101532.add_feature.WWFJW sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Sat May 16 10:19:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150516_101537.add_feature.E7XLA sydh:uw:haib:uta:embl pepr None 1 1
#[Sat May 16 10:19:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150516_101542.add_feature.GWLUH sydh:uw:haib:uta:embl pepr None 1 1
#[Sat May 16 10:19:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cdp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cdp 20150516_101554.add_feature.VWYGV sydh:uw:haib:uta:embl pepr None 1 1
#[Sat May 16 11:25:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150516_112134.add_feature.T21LO sydh:uw:haib:uta:embl pepr None 1 1
#[Sat May 16 11:25:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150516_112129.add_feature.0GP39 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Sat May 16 11:25:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150516_112144.add_feature.VW2G1 sydh:uw:haib:uta:embl pepr None 1 1
#[Sat May 16 11:25:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150516_112139.add_feature.LOPMK sydh:uw:haib:uta:embl pepr None 1 1
#[Sat May 16 11:25:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cdp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cdp 20150516_112151.add_feature.HJJMQ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 10:39:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cdp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cdp 20150519_103520.add_feature.GZ2VW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 11:24:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cdp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cdp 20150519_112053.add_feature.2P3ED sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 11:32:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cdp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cdp 20150519_113248.add_feature.O8EAH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 14:04:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150519_140031.add_feature.2QQ5S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 14:04:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150519_140036.add_feature.B242U sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 14:04:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150519_140041.add_feature.CEHSB sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Tue May 19 14:04:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150519_140046.add_feature.VC55P sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 14:04:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150519_140051.add_feature.254N9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:09:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150519_150949.add_feature.S6AG0 sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Tue May 19 15:09:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150519_150954.add_feature.QOCFS sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Tue May 19 15:38:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150519_152551.add_feature.5ORSY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:38:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150519_152556.add_feature.QWFGZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:38:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150519_152601.add_feature.EQ66Z sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150519_152606.add_feature.5H5LY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150519_152611.add_feature.VDDGQ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150519_152621.add_feature.CKC6D sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150519_152616.add_feature.NMRYD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150519_152631.add_feature.NHLYG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150519_152636.add_feature.WJMRH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150519_152626.add_feature.PGA7R sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:41:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150519_152641.add_feature.1AK4I sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:42:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150519_152646.add_feature.IEBM3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:42:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150519_152651.add_feature.D3EGP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:44:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150519_152656.add_feature.J96OH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150519_152706.add_feature.OZ34O sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150519_152701.add_feature.LYD3D sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150519_152712.add_feature.PI8RN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150519_152716.add_feature.NUGBP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150519_152726.add_feature.NWTME sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150519_152721.add_feature.OZQGG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150519_152736.add_feature.OKBFM sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:45:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150519_152741.add_feature.04BSQ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:46:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150519_152746.add_feature.Q1QX1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:48:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150519_152751.add_feature.MBAJO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:48:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150519_152756.add_feature.W8UK6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:48:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150519_152801.add_feature.R76RP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:48:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150519_152806.add_feature.KMRCD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:49:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150519_152816.add_feature.F31YN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:49:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150519_152811.add_feature.GHUBA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:49:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150519_152821.add_feature.B7FCO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:49:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150519_152826.add_feature.EDEDF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:49:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150519_152831.add_feature.029DL sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:49:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150519_152836.add_feature.KCOE2 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:52:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150519_152841.add_feature.UYKMG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:52:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150519_152851.add_feature.YMA4R sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:52:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150519_152846.add_feature.N2B88 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:52:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150519_152856.add_feature.167IY sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Tue May 19 15:52:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150519_152901.add_feature.N88F5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:52:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150519_152906.add_feature.OOI7G sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:52:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150519_152912.add_feature.Z3384 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:53:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150519_152922.add_feature.C6QWJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:53:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150519_152916.add_feature.RE8G7 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 19 15:53:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150519_152926.add_feature.1Z8SD sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 09:23:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150520_092232.add_feature.2KWSM sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Wed May 20 09:23:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150520_092227.add_feature.FTWFO sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 09:23:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150520_092254.add_feature.HAA4Y sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:27:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150520_102626.add_feature.SIG0D sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:27:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150520_102631.add_feature.JQ4DE sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:27:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150520_102641.add_feature.GDR2J sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:27:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150520_102646.add_feature.CHZJ4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:27:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150520_102636.add_feature.0WZ5H sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:28:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150520_102656.add_feature.JT7SZ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:28:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150520_102711.add_feature.9N746 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:28:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150520_102706.add_feature.S8YE8 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:28:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150520_102721.add_feature.URBHW sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Wed May 20 10:28:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150520_102726.add_feature.7INRH sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:29:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150520_102736.add_feature.Q1WC0 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:29:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150520_102741.add_feature.YYYLB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:30:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150520_102751.add_feature.SDI3F sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:32:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150520_102756.add_feature.8ZUGE sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:33:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150520_102801.add_feature.EF047 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:35:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150520_102811.add_feature.T9ALC sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:35:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150520_102806.add_feature.KPS29 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:36:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150520_102816.add_feature.RY9HB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:41:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150520_102821.add_feature.CRETM sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 10:48:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150520_102826.add_feature.XO6XY sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:06:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150520_102831.add_feature.VICH4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:07:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150520_102836.add_feature.A36CS sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:08:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150520_102841.add_feature.O2K9D sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:09:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150520_102846.add_feature.5K417 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:09:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150520_102851.add_feature.KHG6G sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:10:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150520_102856.add_feature.UZC3R sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:10:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150520_102901.add_feature.9LQ9W sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:10:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150520_102906.add_feature.VGM7B sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:11:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150520_102916.add_feature.XIYZE sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:11:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150520_102911.add_feature.KPCP5 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:12:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150520_102921.add_feature.ISU2Q sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:13:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150520_102927.add_feature.Q17A4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:15:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150520_102931.add_feature.L2Y8R sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:16:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150520_102936.add_feature.Q1F8V sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:16:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150520_102947.add_feature.MTOI6 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:16:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150520_102941.add_feature.AFNUG sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:16:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150520_102951.add_feature.AJ9QF sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:17:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150520_102956.add_feature.NOJJO sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:17:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150520_103001.add_feature.3DZO0 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 11:17:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150520_103007.add_feature.9L0L1 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 15:05:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150520_150411.add_feature.GALZC sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Wed May 20 15:57:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150520_155649.add_feature.3BSC4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:14:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150520_161304.add_feature.XOSOY sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:30:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150520_162920.add_feature.7LL33 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:30:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150520_162915.add_feature.U6NKB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:30:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150520_162925.add_feature.38RXB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:30:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150520_162935.add_feature.6LVMM sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:30:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150520_162950.add_feature.N2FXY sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:30:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150520_162945.add_feature.LTJ6P sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:31:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150520_163005.add_feature.D0FKA sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:31:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150520_162955.add_feature.HM89Q sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:31:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150520_163010.add_feature.4LB9O sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:31:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150520_163016.add_feature.1Q4E3 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:31:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150520_163020.add_feature.K8GQA sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:32:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150520_163026.add_feature.M51LX sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:32:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150520_163031.add_feature.NWR5E sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:32:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150520_163035.add_feature.M8B0I sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:33:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150520_163042.add_feature.9ODVB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:35:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150520_163046.add_feature.OQQSS sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:36:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150520_163052.add_feature.QWODS sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:37:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150520_163057.add_feature.96TKX sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:40:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150520_163102.add_feature.XDAJU sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:41:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150520_163106.add_feature.ZCK2V sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:42:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150520_163111.add_feature.8JE5D sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:43:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150520_163118.add_feature.ZALPY sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:43:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150520_163118.add_feature.756EN sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Wed May 20 16:44:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150520_163122.add_feature.ODTLY sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:45:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150520_163127.add_feature.O40MK sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:47:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150520_163132.add_feature.RJHCQ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:48:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150520_163146.add_feature.WYR5Q sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Wed May 20 16:49:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150520_163148.add_feature.YM6BK sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 16:50:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150520_163148.add_feature.QPG70 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:47:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150520_174628.add_feature.EHGVY sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:47:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150520_174623.add_feature.UID2X sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:47:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150520_174638.add_feature.VOWDW sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:48:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150520_174704.add_feature.BFB84 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:48:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150520_174714.add_feature.ZDTNK sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:48:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150520_174709.add_feature.L352R sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:48:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150520_174719.add_feature.E2X2E sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:48:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150520_174724.add_feature.C8DUR sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:48:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150520_174729.add_feature.R1O2G sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:49:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150520_174734.add_feature.WA32P sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:49:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150520_174739.add_feature.UFE8V sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:50:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150520_174744.add_feature.L3XJ0 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:50:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150520_174749.add_feature.NAO0D sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 17:51:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150520_174754.add_feature.WQCM6 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:01:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150520_174759.add_feature.3WZSL sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:02:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150520_174804.add_feature.HL57M sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:04:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150520_174809.add_feature.8GB3D sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:05:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150520_174814.add_feature.54182 sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Wed May 20 18:12:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150520_174819.add_feature.BY0VB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:28:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150520_174824.add_feature.1B6PM sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:33:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150520_174844.add_feature.JH7DH sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Wed May 20 18:52:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150520_175737.add_feature.77Y6L sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:53:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150520_175742.add_feature.JY81W sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:53:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150520_175747.add_feature.VA4U1 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:54:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150520_175752.add_feature.OBJI1 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:54:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150520_175757.add_feature.PLGK9 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:54:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150520_175802.add_feature.BO16N sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:55:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150520_175807.add_feature.AVJ5R sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:55:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150520_175812.add_feature.PF1HB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:55:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150520_175817.add_feature.8DDW2 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:56:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150520_175828.add_feature.IGC8V sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:56:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150520_175822.add_feature.D1DRA sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:56:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150520_175833.add_feature.U01T4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:56:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150520_175838.add_feature.MPBG6 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:57:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150520_175848.add_feature.QY3OM sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:57:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150520_175843.add_feature.KCN7S sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:57:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150520_175853.add_feature.7F6SA sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:57:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150520_175858.add_feature.7UZNR sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:58:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150520_175903.add_feature.WWTOL sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:58:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150520_175908.add_feature.XSCDW sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:58:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150520_175913.add_feature.KQUA7 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:58:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150520_175918.add_feature.ZVSQ9 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:59:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150520_175923.add_feature.LZDZV sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:59:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150520_175928.add_feature.HVMNB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 20 18:59:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150520_175933.add_feature.4UYZ0 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:57:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150521_085623.add_feature.TKV4X sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:57:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150521_085628.add_feature.HF985 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:58:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150521_085708.add_feature.8VOOL sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:58:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150521_085657.add_feature.P79LK sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:58:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150521_085713.add_feature.FARI4 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:58:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150521_085718.add_feature.GR4DR sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:58:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150521_085728.add_feature.QRI5C sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:58:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150521_085723.add_feature.GC22P sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:59:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150521_085733.add_feature.Q6VWC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:59:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150521_085738.add_feature.KR7Q7 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 08:59:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150521_085743.add_feature.06MQ4 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 09:39:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150521_092316.add_feature.31A1H sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:54:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150521_115325.add_feature.CLVEP sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:54:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150521_115330.add_feature.QYRQI sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:54:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150521_115335.add_feature.9H6XW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:54:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150521_115350.add_feature.BQLF3 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:54:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150521_115345.add_feature.UDP4F sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:55:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150521_115405.add_feature.AE7UQ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:55:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150521_115400.add_feature.OVSPM sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:55:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150521_115355.add_feature.NHB2V sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:55:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150521_115411.add_feature.9AUKM sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:55:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150521_115421.add_feature.SOSA5 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:55:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150521_115426.add_feature.RSKL7 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:55:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150521_115416.add_feature.H83MX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:56:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150521_115431.add_feature.0TB3Y sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:56:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150521_115436.add_feature.2A1II sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 11:57:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150521_115444.add_feature.KJV2L sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Thu May 21 12:00:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150521_115910.add_feature.D6VSC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 12:00:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150521_115905.add_feature.101JN sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 12:00:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150521_115915.add_feature.2QZ7H sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 12:00:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150521_115925.add_feature.W37GD sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 12:00:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150521_115920.add_feature.TFBN5 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:11:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150521_131144.add_feature.P6CE1 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:11:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150521_131146.add_feature.5Q38Q sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:11:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150521_131148.add_feature.OVZBG sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:11:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150521_131150.add_feature.L5DDF sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:11:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150521_131152.add_feature.MTUBN sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:11:59 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150521_131154.add_feature.3SX5C sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150521_131156.add_feature.DC33A sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150521_131158.add_feature.8BJYJ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150521_131200.add_feature.GNNTC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150521_131202.add_feature.BU4EQ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150521_131204.add_feature.IBSG5 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150521_131207.add_feature.2OCMA sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150521_131209.add_feature.RK0JZ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150521_131211.add_feature.57YTS sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150521_131213.add_feature.B905B sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Thu May 21 13:12:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150521_131215.add_feature.DH0YL sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150521_131217.add_feature.VIE5B sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150521_131219.add_feature.IQUWN sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150521_131221.add_feature.P3K6J sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150521_131223.add_feature.WQ0RE sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:29 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150521_131225.add_feature.YNN6C sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150521_131227.add_feature.Q3ENP sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150521_131229.add_feature.6445G sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150521_131231.add_feature.A26Z2 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150521_131233.add_feature.VW7HL sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150521_131235.add_feature.HQGBA sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150521_131237.add_feature.GWMDS sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150521_131239.add_feature.2JVQX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150521_131241.add_feature.ML1OM sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150521_131243.add_feature.6HOGJ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150521_131245.add_feature.EG9YM sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150521_131247.add_feature.GDWQI sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150521_131249.add_feature.92PQ4 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150521_131251.add_feature.Y9PZ0 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:12:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150521_131253.add_feature.UK5CE sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150521_131255.add_feature.3I3CO sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150521_131257.add_feature.WZN7U sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150521_131259.add_feature.D2G5K sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150521_131301.add_feature.FF7MA sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150521_131303.add_feature.XZZWX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150521_131305.add_feature.5V8UV sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150521_131307.add_feature.CD7LI sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150521_131309.add_feature.7ENYE sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150521_131311.add_feature.UR4UV sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150521_131313.add_feature.SFUP3 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150521_131315.add_feature.8FDEF sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Thu May 21 13:13:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150521_131317.add_feature.KDJPA sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150521_131319.add_feature.OFJKV sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150521_131321.add_feature.G4KYK sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150521_131323.add_feature.7J73R sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Thu May 21 13:13:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150521_131325.add_feature.71QC5 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 13:13:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150521_131327.add_feature.PNRTH sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150521_170216.add_feature.UZUEU sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150521_170218.add_feature.2CHU0 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150521_170220.add_feature.LKP2N sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150521_170222.add_feature.Z3A4V sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150521_170224.add_feature.KZSOB sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150521_170226.add_feature.1DULW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150521_170228.add_feature.USZFM sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Thu May 21 17:02:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150521_170230.add_feature.I4UFS sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150521_170232.add_feature.7LP5C sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150521_170234.add_feature.PHV93 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150521_170236.add_feature.WI19K sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150521_170238.add_feature.P5Q7T sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150521_170240.add_feature.FEYWH sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150521_170243.add_feature.HTG4F sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150521_170244.add_feature.HFFI3 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Thu May 21 17:02:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150521_170246.add_feature.30ZI9 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150521_170248.add_feature.LRMV4 sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Thu May 21 17:02:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150521_170251.add_feature.QBXUB sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150521_170253.add_feature.H2YGE sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:02:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150521_170255.add_feature.4UJG8 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150521_170257.add_feature.60JZW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150521_170259.add_feature.438L5 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150521_170301.add_feature.X3X4H sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150521_170303.add_feature.TELHD sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150521_170305.add_feature.NUN70 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150521_170307.add_feature.6PZAT sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150521_170309.add_feature.ZQNYG sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150521_170311.add_feature.ZF9MZ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150521_170313.add_feature.N60VR sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150521_170315.add_feature.1Z9MW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150521_170317.add_feature.86OPE sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150521_170319.add_feature.WR42N sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150521_170321.add_feature.HTY5V sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150521_170323.add_feature.7Z5PU sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Thu May 21 17:03:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150521_170325.add_feature.V494U sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150521_170327.add_feature.5Y0FY sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150521_170329.add_feature.M34PC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150521_170331.add_feature.AQWSB sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150521_170333.add_feature.JT327 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150521_170335.add_feature.VFMF6 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150521_170337.add_feature.5CP4X sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150521_170339.add_feature.0IIKL sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150521_170341.add_feature.13139 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150521_170343.add_feature.JRLNP sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150521_170345.add_feature.ZH3J3 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150521_170348.add_feature.XQEF2 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150521_170349.add_feature.JVJQF sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150521_170351.add_feature.NWPDG sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:03:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150521_170353.add_feature.1WJRF sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150521_170355.add_feature.6R250 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150521_170357.add_feature.I37SX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150521_170359.add_feature.6ELRI sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150521_170401.add_feature.F8HE7 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May 21 17:04:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150521_170403.add_feature.39I5K sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150521_170405.add_feature.TMLBH sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150521_170407.add_feature.TBW8T sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150521_170409.add_feature.ZNJW6 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150521_170411.add_feature.ROLJP sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150521_170413.add_feature.7YJY0 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150521_170415.add_feature.GFJRC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150521_170417.add_feature.7XO0Y sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150521_170419.add_feature.5A0FB sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150521_170421.add_feature.R4QBT sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150521_170423.add_feature.OHMMP sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150521_170425.add_feature.YIV93 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150521_170427.add_feature.HFT41 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150521_170429.add_feature.U785R sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150521_170431.add_feature.28QZ1 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150521_170434.add_feature.4O1UW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150521_170436.add_feature.YTQQ7 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150521_170438.add_feature.TSRLY sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150521_170440.add_feature.4JEHK sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150521_170442.add_feature.SSGWQ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150521_170444.add_feature.GUMXW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150521_170446.add_feature.K00A3 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150521_170448.add_feature.0HNB4 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150521_170450.add_feature.ZGKSI sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:04:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150521_170452.add_feature.WA9H0 sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Thu May 21 17:04:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150521_170454.add_feature.A1NT2 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 17:05:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150521_170456.add_feature.4GYSS sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:55:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150521_215555.add_feature.7QKSE sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150521_215557.add_feature.RK19H sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150521_215559.add_feature.46F15 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150521_215601.add_feature.YGCSX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150521_215603.add_feature.NHXXX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150521_215605.add_feature.M8ZFR sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150521_215607.add_feature.KGWWQ sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Thu May 21 21:56:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150521_215609.add_feature.FHU7W sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150521_215611.add_feature.ZDL5H sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150521_215613.add_feature.ABL4T sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150521_215615.add_feature.XTN45 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150521_215617.add_feature.QQYUJ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150521_215619.add_feature.Z5UQF sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150521_215621.add_feature.MUEAC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150521_215623.add_feature.5TZPJ sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Thu May 21 21:56:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150521_215625.add_feature.CG887 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150521_215627.add_feature.CPZ90 sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Thu May 21 21:56:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150521_215629.add_feature.DOEUD sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150521_215631.add_feature.BTG1P sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150521_215633.add_feature.2BCJC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150521_215635.add_feature.1YK3O sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150521_215637.add_feature.X452O sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150521_215639.add_feature.QF95O sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150521_215641.add_feature.GM1X5 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150521_215643.add_feature.W1MEX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150521_215645.add_feature.F4OQL sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150521_215647.add_feature.WJQ6P sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150521_215649.add_feature.DDE6R sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150521_215651.add_feature.XIUEM sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:56:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150521_215653.add_feature.Z1NN2 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150521_215655.add_feature.6Q50D sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150521_215657.add_feature.RXZ6Y sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150521_215659.add_feature.JBONQ sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150521_215701.add_feature.7CSV0 sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Thu May 21 21:57:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150521_215703.add_feature.N6U4Y sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150521_215705.add_feature.DOJHG sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150521_215707.add_feature.NY2XR sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150521_215709.add_feature.T4H8X sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150521_215711.add_feature.UDG5W sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150521_215713.add_feature.MGP0M sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150521_215715.add_feature.RTRZX sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150521_215717.add_feature.IWLNN sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150521_215719.add_feature.19BB4 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150521_215721.add_feature.4CJ9O sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150521_215723.add_feature.095G0 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150521_215725.add_feature.KWYNR sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150521_215727.add_feature.SPWF6 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150521_215729.add_feature.2928T sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150521_215731.add_feature.UGV5I sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150521_215733.add_feature.QTOA2 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150521_215735.add_feature.2VMTL sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150521_215737.add_feature.J0Z2X sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150521_215739.add_feature.B2VP7 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Thu May 21 21:57:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150521_215741.add_feature.Y730Y sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150521_215743.add_feature.RKQCW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150521_215745.add_feature.870DW sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150521_215747.add_feature.I0HKS sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150521_215749.add_feature.MRBIR sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150521_215751.add_feature.J26ZT sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:57:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150521_215753.add_feature.XAI62 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150521_215755.add_feature.LT9J1 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150521_215757.add_feature.AFCOP sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150521_215759.add_feature.JR353 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150521_215802.add_feature.34NNI sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150521_215804.add_feature.F69CE sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150521_215806.add_feature.3S7OR sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150521_215808.add_feature.NQKEO sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150521_215810.add_feature.STLKC sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150521_215812.add_feature.7UI8Y sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150521_215814.add_feature.RF1BU sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150521_215816.add_feature.1NE9E sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150521_215818.add_feature.L9CJT sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150521_215820.add_feature.5O5II sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150521_215822.add_feature.0S2U8 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150521_215824.add_feature.HP2SK sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150521_215826.add_feature.3HXMO sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150521_215828.add_feature.A0FQU sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150521_215830.add_feature.WVEN6 sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Thu May 21 21:58:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150521_215832.add_feature.D3C41 sydh:uw:haib:uta:embl pepr None 1 1
#[Thu May 21 21:58:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150521_215834.add_feature.K1JG2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150522_145458.add_feature.1CUQC sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150522_145500.add_feature.KEV1R sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150522_145502.add_feature.EIPZ5 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150522_145504.add_feature.VQO0R sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150522_145506.add_feature.YPAGY sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150522_145508.add_feature.YSRK3 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150522_145510.add_feature.4FY6F sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150522_145512.add_feature.FZ7HO sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150522_145514.add_feature.3AU25 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150522_145516.add_feature.L8KEL sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150522_145518.add_feature.KR2VO sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150522_145520.add_feature.1CCN6 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150522_145522.add_feature.KXX5E sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:29 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150522_145524.add_feature.HYOBP sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150522_145526.add_feature.RLXN2 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Fri May 22 14:55:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150522_145528.add_feature.3KSG2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:35 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150522_145530.add_feature.5E32A sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150522_145532.add_feature.DR8M7 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150522_145534.add_feature.T3FIS sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150522_145536.add_feature.KV95Y sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150522_145538.add_feature.8IVO2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150522_145540.add_feature.WX8S0 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150522_145542.add_feature.YXIAN sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150522_145544.add_feature.J6FTF sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150522_145546.add_feature.Q09OZ sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150522_145548.add_feature.J2FU9 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150522_145551.add_feature.LXIAE sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:55:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150522_145552.add_feature.KPVUM sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150522_145554.add_feature.2NMIH sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150522_145556.add_feature.SOHJS sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150522_145558.add_feature.CSLS5 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150522_145601.add_feature.MDT04 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150522_145602.add_feature.Y212Q sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150522_145605.add_feature.BY50T sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150522_145607.add_feature.192B4 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150522_145609.add_feature.UFO1V sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:56:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150522_145611.add_feature.F649A sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150522_145735.add_feature.H4RL2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150522_145737.add_feature.GYOYT sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150522_145739.add_feature.9CZX4 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150522_145741.add_feature.RPOV9 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150522_145743.add_feature.T6259 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150522_145745.add_feature.SSETT sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150522_145747.add_feature.NHVN8 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150522_145749.add_feature.LQQP7 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150522_145751.add_feature.259M0 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:57:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150522_145753.add_feature.WHDYH sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150522_145755.add_feature.SGYXG sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150522_145757.add_feature.4HOQP sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150522_145759.add_feature.IRQFN sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150522_145801.add_feature.T8LDM sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150522_145803.add_feature.2QQ4J sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Fri May 22 14:58:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150522_145805.add_feature.FAMGG sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150522_145807.add_feature.BJ6YW sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150522_145809.add_feature.31UMM sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150522_145811.add_feature.2LTK0 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150522_145813.add_feature.CREFV sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150522_145815.add_feature.7UDTC sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150522_145817.add_feature.YSSMP sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150522_145819.add_feature.HMDAM sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150522_145821.add_feature.OG69E sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150522_145823.add_feature.KOUFU sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150522_145825.add_feature.WFCZU sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150522_145827.add_feature.3FLF2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150522_145829.add_feature.OFB5F sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150522_145831.add_feature.OSC6N sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150522_145833.add_feature.JG2A3 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150522_145835.add_feature.OYRQ5 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150522_145837.add_feature.1AAS4 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150522_145839.add_feature.SZMMI sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150522_145841.add_feature.VH9FD sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150522_145843.add_feature.KCB5H sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150522_145845.add_feature.E6NKF sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150522_145847.add_feature.88H13 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150522_145849.add_feature.PQS96 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150522_145851.add_feature.ZXFWK sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:58:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150522_145853.add_feature.D05GO sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150522_145855.add_feature.948HP sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150522_145858.add_feature.6N69D sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150522_145859.add_feature.L7J3E sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150522_145901.add_feature.VCVYB sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150522_145904.add_feature.VIZKH sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150522_145906.add_feature.9WREG sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Fri May 22 14:59:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150522_145908.add_feature.16Z7H sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150522_145910.add_feature.QNLS5 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150522_145912.add_feature.I0HPJ sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150522_145914.add_feature.7V6CA sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Fri May 22 14:59:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150522_145916.add_feature.TIKQL sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 14:59:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150522_145918.add_feature.K0BTD sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150522_150021.add_feature.7L623 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150522_150023.add_feature.HLAH2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150522_150025.add_feature.36DHM sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150522_150027.add_feature.0F9LX sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150522_150029.add_feature.Q22WU sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150522_150031.add_feature.1XLZF sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150522_150033.add_feature.ZX8O5 sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Fri May 22 15:00:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150522_150035.add_feature.NUXH0 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150522_150037.add_feature.ZZB19 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150522_150039.add_feature.C8J3O sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150522_150041.add_feature.Z38BI sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150522_150043.add_feature.RTYY2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150522_150045.add_feature.9J0GY sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150522_150047.add_feature.A608R sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150522_150049.add_feature.WOHL8 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Fri May 22 15:00:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150522_150051.add_feature.BZXLX sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:00:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150522_150053.add_feature.83YNV sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Fri May 22 15:01:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150522_150055.add_feature.OL488 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150522_150057.add_feature.4EHCO sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150522_150059.add_feature.B7DKA sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150522_150101.add_feature.XICK3 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150522_150103.add_feature.7IF41 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150522_150105.add_feature.21O9Y sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150522_150107.add_feature.D6UGR sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150522_150109.add_feature.89KUB sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150522_150111.add_feature.M1MRL sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150522_150113.add_feature.3L2AK sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150522_150115.add_feature.IW3JD sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150522_150117.add_feature.ICLNW sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150522_150119.add_feature.WTWMC sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150522_150121.add_feature.Q8WZQ sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150522_150123.add_feature.8EUSV sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150522_150125.add_feature.NH1QW sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150522_150127.add_feature.I68BL sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Fri May 22 15:01:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150522_150129.add_feature.Z5HJQ sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150522_150131.add_feature.DOXCE sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150522_150133.add_feature.HYRH0 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150522_150135.add_feature.ZIO14 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150522_150137.add_feature.3CGPL sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150522_150139.add_feature.OPHUU sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150522_150141.add_feature.HMWRT sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150522_150143.add_feature.7YRBX sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150522_150145.add_feature.2EUTB sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150522_150147.add_feature.33J0F sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150522_150149.add_feature.X7HT2 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150522_150151.add_feature.101QZ sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:01:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150522_150153.add_feature.6HXBG sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150522_150155.add_feature.CMPME sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150522_150157.add_feature.H8CTG sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150522_150159.add_feature.30B0B sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150522_150201.add_feature.SIX68 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150522_150203.add_feature.IEI7O sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150522_150206.add_feature.Q2K57 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Fri May 22 15:02:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150522_150207.add_feature.EB1NK sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150522_150209.add_feature.8YW09 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150522_150211.add_feature.8SQRO sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150522_150213.add_feature.JE81H sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150522_150216.add_feature.GIF5B sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150522_150218.add_feature.551ZX sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150522_150220.add_feature.I25F9 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150522_150222.add_feature.0WEM4 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150522_150224.add_feature.09B67 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150522_150226.add_feature.MU4HR sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150522_150228.add_feature.CB0H7 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150522_150230.add_feature.MCE6C sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150522_150232.add_feature.KDIKI sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150522_150234.add_feature.WAJTR sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150522_150236.add_feature.2KJ3T sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150522_150238.add_feature.XVRIX sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150522_150240.add_feature.V75TU sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150522_150242.add_feature.5VM5G sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150522_150244.add_feature.4G8U1 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150522_150246.add_feature.IH89I sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150522_150248.add_feature.1AAD7 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150522_150250.add_feature.LZUAY sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150522_150252.add_feature.U2ZLI sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:02:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150522_150254.add_feature.MMO6H sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:03:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150522_150256.add_feature.FHPBT sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Fri May 22 15:03:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150522_150258.add_feature.0PH85 sydh:uw:haib:uta:embl pepr None 1 1
#[Fri May 22 15:03:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150522_150300.add_feature.60RJ8 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150526_171206.add_feature.GG8UN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150526_171208.add_feature.SKSJA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150526_171210.add_feature.EM0D3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150526_171212.add_feature.PL320 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150526_171214.add_feature.OT0NN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150526_171216.add_feature.QEDYQ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150526_171218.add_feature.EOT5C sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150526_171220.add_feature.O8SZ3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150526_171222.add_feature.5MLTB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150526_171224.add_feature.V7597 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150526_171227.add_feature.JO486 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:29 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150526_171229.add_feature.KWCUZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150526_171231.add_feature.20BV5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150526_171233.add_feature.WWDZG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:35 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150526_171235.add_feature.OPWTJ sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Tue May 26 17:12:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150526_171237.add_feature.DXVW0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150526_171239.add_feature.4M2RB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150526_171241.add_feature.EYPO6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150526_171243.add_feature.W3KTH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150526_171245.add_feature.QJS79 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150526_171247.add_feature.F76HS sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150526_171249.add_feature.MRVOU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150526_171251.add_feature.WW2TW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150526_171253.add_feature.21RLL sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150526_171255.add_feature.I553F sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:12:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150526_171257.add_feature.P8POR sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150526_171259.add_feature.KU8K3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150526_171301.add_feature.DXRYP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150526_171303.add_feature.IUTHW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150526_171305.add_feature.40F56 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150526_171307.add_feature.C4XL9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150526_171309.add_feature.9U55B sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150526_171311.add_feature.1LVEH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150526_171313.add_feature.8RSDP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150526_171315.add_feature.R519Q sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150526_171317.add_feature.V9O0F sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150526_171319.add_feature.FM972 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150526_171322.add_feature.12P72 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150526_171324.add_feature.IXNGV sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150526_171326.add_feature.JU03R sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150526_171328.add_feature.D2YNG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150526_171330.add_feature.NZ2XM sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150526_171332.add_feature.B1NX4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150526_171334.add_feature.H8BSB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150526_171336.add_feature.O46PJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150526_171338.add_feature.7T3V1 sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Tue May 26 17:13:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150526_171340.add_feature.0P7J2 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150526_171342.add_feature.GAT36 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150526_171344.add_feature.IFYIP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150526_171346.add_feature.W7NVP sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Tue May 26 17:13:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150526_171348.add_feature.MOI1S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:13:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150526_171350.add_feature.BCTAG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150526_171414.add_feature.IM22P sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150526_171416.add_feature.GTSS6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150526_171418.add_feature.8HMBL sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150526_171420.add_feature.451IO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150526_171422.add_feature.N37OE sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150526_171424.add_feature.0P3U9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150526_171426.add_feature.0KMZP sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Tue May 26 17:14:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150526_171428.add_feature.O2Z1C sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150526_171430.add_feature.TYLPE sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150526_171432.add_feature.7ZCSS sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150526_171434.add_feature.ETAQD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150526_171436.add_feature.YWGVN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150526_171438.add_feature.TRU5X sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150526_171440.add_feature.YU3Z4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150526_171443.add_feature.IN12R sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Tue May 26 17:14:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150526_171445.add_feature.C4ZEC sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150526_171447.add_feature.E13XI sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Tue May 26 17:14:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150526_171449.add_feature.A43SV sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150526_171451.add_feature.F368P sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150526_171453.add_feature.8HM2R sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150526_171455.add_feature.9SAQ4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:14:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150526_171457.add_feature.3F39Q sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150526_171459.add_feature.1JU9I sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150526_171501.add_feature.HOV6X sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150526_171503.add_feature.3C9EF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150526_171505.add_feature.V52YD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150526_171507.add_feature.ZF3BT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150526_171509.add_feature.3FV0T sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150526_171511.add_feature.M2PIL sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150526_171513.add_feature.REMCS sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150526_171515.add_feature.XG49O sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150526_171517.add_feature.7MVHG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150526_171519.add_feature.B5KAT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150526_171521.add_feature.88ULK sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Tue May 26 17:15:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150526_171523.add_feature.7X9W3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150526_171525.add_feature.TYLS5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150526_171527.add_feature.BLKC9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150526_171529.add_feature.FSQTO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150526_171531.add_feature.UD7AJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150526_171533.add_feature.54IRN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150526_171535.add_feature.986N6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150526_171537.add_feature.MX781 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150526_171539.add_feature.T5JW5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150526_171542.add_feature.Q75U0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150526_171544.add_feature.AEBUX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150526_171546.add_feature.H6D49 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150526_171548.add_feature.HI4AZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150526_171550.add_feature.VNPHB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150526_171552.add_feature.PL5YP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150526_171554.add_feature.7KXTH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150526_171556.add_feature.Z4I19 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:15:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150526_171558.add_feature.FQUR9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150526_171600.add_feature.ZRMGI sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Tue May 26 17:16:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150526_171602.add_feature.ZTPN9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150526_171604.add_feature.IZHLS sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150526_171606.add_feature.R70LB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150526_171608.add_feature.3CHF1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150526_171610.add_feature.P80XD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150526_171612.add_feature.W3N31 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150526_171614.add_feature.OYW7J sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150526_171616.add_feature.NQG48 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150526_171618.add_feature.Z6LBM sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150526_171620.add_feature.M1X2S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150526_171622.add_feature.U4SNX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150526_171624.add_feature.NZAM0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150526_171626.add_feature.YLPDO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150526_171628.add_feature.O78SO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150526_171630.add_feature.70ECL sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150526_171632.add_feature.UZY4X sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150526_171634.add_feature.K66Y2 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150526_171636.add_feature.PH2ND sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150526_171639.add_feature.IDRFB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150526_171640.add_feature.FPBUN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150526_171642.add_feature.TOG1T sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150526_171644.add_feature.VRVCG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150526_171647.add_feature.OGF6T sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150526_171649.add_feature.TR7K1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150526_171651.add_feature.Y6TOM sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Tue May 26 17:16:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150526_171653.add_feature.G48SZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue May 26 17:16:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150526_171655.add_feature.YFEQI sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:20:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150527_162055.add_feature.K056V sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:20:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150527_162057.add_feature.POJHH sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150527_162059.add_feature.0FB35 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150527_162101.add_feature.UED99 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150527_162103.add_feature.HBFQ3 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150527_162105.add_feature.VWSX4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150527_162107.add_feature.TNTHM sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Wed May 27 16:21:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150527_162109.add_feature.3PHEQ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150527_162111.add_feature.4SIZS sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150527_162113.add_feature.HCQ4R sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150527_162113.add_feature.SUG8W sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150527_162115.add_feature.MT121 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150527_162116.add_feature.B82E6 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150527_162118.add_feature.5J8NF sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150527_162120.add_feature.KWNUG sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150527_162122.add_feature.2PS5U sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150527_162124.add_feature.P8MH4 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Wed May 27 16:21:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150527_162126.add_feature.3NERY sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150527_162128.add_feature.BCTJR sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Wed May 27 16:21:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150527_162130.add_feature.94V3S sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150527_162132.add_feature.ULR3Q sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150527_162134.add_feature.6CTUN sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150527_162136.add_feature.KZZNB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150527_162138.add_feature.574HP sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:21:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150527_162140.add_feature.0S5IX sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150527_164139.add_feature.TGRHI sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150527_164141.add_feature.Z21AZ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150527_164143.add_feature.7YW4J sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150527_164145.add_feature.J9EJD sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150527_164147.add_feature.V2PKQ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150527_164149.add_feature.K59UO sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150527_164151.add_feature.HCMQR sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Wed May 27 16:41:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150527_164153.add_feature.2FAHH sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150527_164155.add_feature.THQMS sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:41:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150527_164157.add_feature.IVLNQ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150527_164159.add_feature.LISHG sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150527_164201.add_feature.9YZHA sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150527_164203.add_feature.31XZQ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150527_164205.add_feature.3M0EQ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150527_164207.add_feature.5Z7OJ sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Wed May 27 16:42:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150527_164209.add_feature.GOTL2 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150527_164211.add_feature.XUNLN sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Wed May 27 16:42:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150527_164213.add_feature.AT98D sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150527_164215.add_feature.FZA7M sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150527_164217.add_feature.TXQU6 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150527_164220.add_feature.3DZ0F sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150527_164221.add_feature.9EIV1 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150527_164224.add_feature.MPNGL sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150527_164234.add_feature.UTY23 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150527_164236.add_feature.W3YP9 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150527_164238.add_feature.4VZRZ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150527_164240.add_feature.TVAGS sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150527_164242.add_feature.ZX47C sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150527_164244.add_feature.KNQEI sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150527_164246.add_feature.03DPF sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150527_164248.add_feature.93SD1 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150527_164250.add_feature.R0ONU sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150527_164252.add_feature.B3G1A sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150527_164254.add_feature.WFGTG sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Wed May 27 16:42:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150527_164256.add_feature.QN1D5 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:42:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150527_164258.add_feature.G3TL9 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150527_164300.add_feature.RG52F sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150527_164302.add_feature.9DXXI sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150527_164304.add_feature.CXHZ2 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150527_164306.add_feature.XLTA4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150527_164308.add_feature.S9CVL sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150527_164310.add_feature.27ZFE sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150527_164312.add_feature.66G5M sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150527_164314.add_feature.87X2Q sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150527_164316.add_feature.92JL6 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150527_164318.add_feature.IITMH sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150527_164320.add_feature.PT8DP sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150527_164323.add_feature.DBSLA sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150527_164325.add_feature.DC5TA sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150527_164327.add_feature.8ID4E sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150527_164329.add_feature.Q9334 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150527_164331.add_feature.W16JT sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150527_164333.add_feature.8DT48 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Wed May 27 16:43:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150527_164335.add_feature.APH1C sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150527_164337.add_feature.ENC53 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150527_164339.add_feature.XXBI6 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150527_164341.add_feature.LJ2JU sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150527_164343.add_feature.PNXKM sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150527_164345.add_feature.GM0M1 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150527_164347.add_feature.LRBCN sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150527_164349.add_feature.4HHQB sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150527_164351.add_feature.4KREM sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150527_164353.add_feature.0ERD4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150527_164355.add_feature.RC5CJ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:43:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150527_164357.add_feature.D6K36 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150527_164359.add_feature.FUZEP sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150527_164401.add_feature.HMK2R sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150527_164403.add_feature.CVWH7 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150527_164405.add_feature.UR03I sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150527_164407.add_feature.KIPI9 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150527_164409.add_feature.O22L0 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150527_164411.add_feature.XU5IL sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150527_164413.add_feature.Q6LNM sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150527_164415.add_feature.GW2E4 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150527_164417.add_feature.0IYAV sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150527_164419.add_feature.8H3IP sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150527_164421.add_feature.D6NHZ sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150527_164423.add_feature.6B4GF sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Wed May 27 16:44:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150527_164425.add_feature.WBAB0 sydh:uw:haib:uta:embl pepr None 1 1
#[Wed May 27 16:44:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150527_164427.add_feature.5VWQ6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:26:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150601_172652.add_feature.04MWJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:26:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150601_172654.add_feature.6OG4Y sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150601_172656.add_feature.0I7YA sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150601_172658.add_feature.T4AUY sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150601_172700.add_feature.9ZY0X sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150601_172702.add_feature.OQR2V sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150601_172704.add_feature.0XW8X sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150601_172706.add_feature.EW0ZC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150601_172708.add_feature.NCL9Y sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150601_172710.add_feature.QH2NV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150601_172712.add_feature.OXE9X sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150601_172714.add_feature.A4ZEL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150601_172716.add_feature.ENLU1 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150601_172718.add_feature.7MOMC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150601_172720.add_feature.CGA9H sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Mon Jun  1 17:27:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150601_172722.add_feature.MKU0P sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:29 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150601_172724.add_feature.218X4 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150601_172726.add_feature.QP261 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150601_172728.add_feature.C8C21 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:35 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150601_172730.add_feature.QSZ96 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150601_172732.add_feature.J2UPH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150601_172733.add_feature.PVV5S sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150601_172734.add_feature.J6C0Y sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150601_172735.add_feature.VLCYS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150601_172736.add_feature.H77XI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150601_172737.add_feature.PTWEM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150601_172740.add_feature.NKGRI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150601_172738.add_feature.MKRFZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150601_172739.add_feature.EJXMJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:27:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150601_172741.add_feature.14TRA sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150601_172742.add_feature.ACEDX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150601_172753.add_feature.WXQR6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150601_172743.add_feature.TRTDR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150601_172745.add_feature.FK6JG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150601_172747.add_feature.7OQ1U sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150601_172745.add_feature.C9GTU sydh:uw:haib:uta:embl pepr USF1:SP1 1 1
#[Mon Jun  1 17:28:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150601_172749.add_feature.MG13O sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150601_172755.add_feature.URLKN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150601_172749.add_feature.M5ZYT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150601_172747.add_feature.KJDED sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150601_172753.add_feature.UTUXC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150601_172801.add_feature.0HKS2 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX 1 1
#[Mon Jun  1 17:28:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150601_172807.add_feature.LUHQF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150601_172751.add_feature.4KZEL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150601_172801.add_feature.D3EOB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150601_172757.add_feature.EXY0R sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150601_172803.add_feature.GVMHI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150601_172751.add_feature.OVJZX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150601_172755.add_feature.D2GKZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150601_172757.add_feature.9W6Q4 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150601_172809.add_feature.ZXLS8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150601_172807.add_feature.O0RME sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150601_172759.add_feature.H2NPN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150601_172803.add_feature.QHLQ9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150601_172805.add_feature.52R6G sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150601_172759.add_feature.XQE5F sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150601_172811.add_feature.FJ9Z9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150601_172815.add_feature.RQP3P sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150601_172813.add_feature.CEZCL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150601_172817.add_feature.E6JQ5 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150601_172805.add_feature.S0BSU sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Mon Jun  1 17:28:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150601_172809.add_feature.E41DX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150601_172827.add_feature.P7V8W sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150601_172811.add_feature.BW7BT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150601_172829.add_feature.CQ9IT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150601_172813.add_feature.NRBAA sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150601_172815.add_feature.AM63S sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150601_172829.add_feature.LR3TP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150601_172827.add_feature.IH62G sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150601_172817.add_feature.B1DMU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150601_172823.add_feature.DZ7JO sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150601_172819.add_feature.1F51R sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150601_172821.add_feature.Y9N1A sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150601_172825.add_feature.KAHCI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150601_172819.add_feature.4K4T9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150601_172825.add_feature.ID0EC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150601_172821.add_feature.IBF2I sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150601_172823.add_feature.U249V sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Mon Jun  1 17:28:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150601_172831.add_feature.50HNT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150601_172831.add_feature.9EH3D sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Mon Jun  1 17:28:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150601_172833.add_feature.GV040 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150601_172833.add_feature.U507C sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150601_172835.add_feature.ASKVQ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150601_172835.add_feature.MPHSL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:28:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150601_172837.add_feature.PB3TU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150601_172839.add_feature.A0PII sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Mon Jun  1 17:29:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150601_172842.add_feature.8GXHP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150601_172846.add_feature.OS9TR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150601_172844.add_feature.2IK9T sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150601_172848.add_feature.B0GZG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150601_172850.add_feature.M7MOX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150601_172852.add_feature.BC7GQ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150601_172854.add_feature.OFL99 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150601_172856.add_feature.YC4T2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150601_172858.add_feature.FK5F8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150601_172900.add_feature.D0MO0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150601_172902.add_feature.WOKIX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150601_172904.add_feature.TDP2Z sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150601_172906.add_feature.K986E sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150601_172908.add_feature.1RSFF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150601_172910.add_feature.ZIR09 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150601_172912.add_feature.31OXV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150601_172914.add_feature.X457W sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150601_172916.add_feature.UANJK sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150601_172918.add_feature.ML3LB sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA 1 1
#[Mon Jun  1 17:29:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150601_172920.add_feature.SF3VU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150601_172922.add_feature.FIOPZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150601_172924.add_feature.VDMLI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150601_172926.add_feature.SAKVO sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150601_172928.add_feature.58W6Q sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150601_172930.add_feature.V5FVP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150601_172932.add_feature.V2HCN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150601_172934.add_feature.CE65B sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150601_172937.add_feature.5Z42D sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150601_172938.add_feature.H4ESD sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150601_172941.add_feature.6YL0A sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150601_172943.add_feature.0KYER sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150601_172945.add_feature.PQKOD sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150601_172947.add_feature.SSS5N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150601_172949.add_feature.MPRM7 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150601_172951.add_feature.KFA2K sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:29:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150601_172953.add_feature.ZTRDT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150601_172955.add_feature.ICEW3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150601_172957.add_feature.IF5AD sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150601_172959.add_feature.L1EAM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150601_173001.add_feature.WBOAR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150601_173003.add_feature.BKWRQ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150601_173005.add_feature.1F9ED sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150601_173007.add_feature.FLF0N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150601_173009.add_feature.V34FX sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21 1 1
#[Mon Jun  1 17:30:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150601_173011.add_feature.9N8ZE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  1 17:30:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150601_173013.add_feature.7G5X2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150608_111304.add_feature.VLT52 sydh:uw:haib:uta:embl pepr STAT3:MAX:CJUN:CMYC 1 1
#[Mon Jun  8 11:13:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150608_111306.add_feature.Q5YCL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150608_111308.add_feature.GRC7G sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150608_111310.add_feature.ZHJ1Y sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150608_111313.add_feature.NQM27 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150608_111315.add_feature.9FFSO sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150608_111317.add_feature.SYQ0N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150608_111319.add_feature.HJ495 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150608_111321.add_feature.6NEKH sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Mon Jun  8 11:13:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150608_111323.add_feature.4RJE2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150608_111325.add_feature.WSSIZ sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Mon Jun  8 11:13:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150608_111327.add_feature.DEO2S sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150608_111329.add_feature.3OA1Y sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150608_111331.add_feature.U2V6Q sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150608_111333.add_feature.2S3BS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Mon Jun  8 11:13:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150608_111335.add_feature.2H72L sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150608_111337.add_feature.R804C sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150608_111339.add_feature.HF1CB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150608_111341.add_feature.TIO8X sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150608_111343.add_feature.P40K2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150608_111345.add_feature.DLAVH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150608_111347.add_feature.M8F1B sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150608_111349.add_feature.I2PZO sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150608_111351.add_feature.P28MR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:13:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150608_111353.add_feature.BPHTI sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Mon Jun  8 11:14:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150608_111355.add_feature.OFYFM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150608_111357.add_feature.J0V72 sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Mon Jun  8 11:14:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150608_111359.add_feature.7R8QP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150608_111401.add_feature.JOELI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150608_111403.add_feature.HWWOR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150608_111405.add_feature.ENLL8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150608_111407.add_feature.W3IWI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150608_111409.add_feature.811U5 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150608_111412.add_feature.UCE98 sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Mon Jun  8 11:14:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150608_111413.add_feature.I8NNA sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Mon Jun  8 11:14:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150608_111416.add_feature.VVDRJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150608_111418.add_feature.DH0HL sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Mon Jun  8 11:14:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150608_111420.add_feature.UEQL4 sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Mon Jun  8 11:14:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150608_111422.add_feature.5Q3TB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:29 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150608_111424.add_feature.6UL2Y sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Mon Jun  8 11:14:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150608_111426.add_feature.GOH87 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150608_111428.add_feature.5QR00 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:35 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150608_111430.add_feature.2HSDG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150608_111432.add_feature.LFEU9 sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Mon Jun  8 11:14:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150608_111434.add_feature.C19JH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150608_111436.add_feature.UMP6Y sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Mon Jun  8 11:14:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150608_111438.add_feature.6YW3O sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150608_111440.add_feature.P13J4 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150608_111442.add_feature.D2974 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150608_111444.add_feature.ZK5QX sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Mon Jun  8 11:14:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150608_111446.add_feature.OJNCL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:14:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150608_111448.add_feature.QMGUR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:17:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150608_111750.add_feature.9EESV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:17:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150608_111752.add_feature.PTAAY sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:17:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150608_111754.add_feature.3WBJR sydh:uw:haib:uta:embl pepr NFKB:NFATC1:TBLR1:POU2F2:PML:TBP:FOXM1:PAX5:POL24H8:CEBPB:BCL3:ATF2:MTA3:PAX5N19:CHD2:SRF:STAT3:BCLAF1 1 1
#[Mon Jun  8 11:18:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150608_111756.add_feature.L9QWI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150608_111758.add_feature.GMVJH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150608_111800.add_feature.D5RI9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150608_111802.add_feature.7V2MM sydh:uw:haib:uta:embl pepr USF1:SP1:MAZ:CHD2:NFKB:CDP:TBP:MAX 1 1
#[Mon Jun  8 11:18:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150608_111804.add_feature.2LY5M sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150608_111806.add_feature.XR5ZX sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Mon Jun  8 11:18:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150608_111808.add_feature.90AEM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150608_111810.add_feature.0CZY8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150608_111812.add_feature.870I4 sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Mon Jun  8 11:18:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150608_111814.add_feature.PRQ4O sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150608_111816.add_feature.N40J6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150608_111818.add_feature.HJYDS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Mon Jun  8 11:18:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150608_111821.add_feature.VTSRQ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150608_111822.add_feature.78QON sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Mon Jun  8 11:18:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150608_111824.add_feature.SLM1E sydh:uw:haib:uta:embl pepr MXI1:PML 1 1
#[Mon Jun  8 11:18:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150608_111826.add_feature.4NFAW sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150608_111829.add_feature.ZI1TT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150608_111831.add_feature.EAVVT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150608_111833.add_feature.34D48 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150608_111835.add_feature.N7RGP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150608_111837.add_feature.ZIJZM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150608_111839.add_feature.QWE48 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150608_111841.add_feature.O5HH7 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150608_111843.add_feature.93LUW sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Mon Jun  8 11:18:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150608_111845.add_feature.LP82E sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Mon Jun  8 11:18:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150608_111847.add_feature.AWUHL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150608_111849.add_feature.BC2MS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150608_111851.add_feature.OXZYC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:18:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150608_111853.add_feature.6WYPU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150608_111855.add_feature.9VNSU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150608_111857.add_feature.JZ9RS sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Mon Jun  8 11:19:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150608_111859.add_feature.U6V2X sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150608_111901.add_feature.XLS8Z sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150608_111903.add_feature.VUIP5 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150608_111905.add_feature.6RZIS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150608_111907.add_feature.30T4T sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150608_111909.add_feature.SZUOE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150608_111911.add_feature.IDLH0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150608_111913.add_feature.5QQ8H sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150608_111915.add_feature.W3BFW sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Mon Jun  8 11:19:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150608_111917.add_feature.XDXB7 sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun  8 11:19:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150608_111919.add_feature.2IOPD sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun  8 11:19:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150608_111921.add_feature.T1FIP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150608_111923.add_feature.TGPUP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150608_111925.add_feature.T9RJP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150608_111927.add_feature.FXZFK sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Mon Jun  8 11:19:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150608_111929.add_feature.AI6W4 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150608_111931.add_feature.2TIW9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150608_111934.add_feature.HP6J6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150608_111936.add_feature.WCT0I sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA:CEBPB:ATF2:MTA3 1 1
#[Mon Jun  8 11:19:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150608_111938.add_feature.G8UTI sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Mon Jun  8 11:19:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150608_111940.add_feature.MQ2FE sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Mon Jun  8 11:19:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150608_111942.add_feature.R1LO1 sydh:uw:haib:uta:embl pepr PML:CHD2:YY1:POL2:ZNF143:MAX:SMC3:MAZ:TBP:ELK1:SIN3A:WHIP:ELF1:BCLAF1:RAD21:POU2F2:P300:ETS1:CMYC:GR:MXI1:NFYB:CHD1:TBLR1:ZEB1:TAF1:SP1:GABP:STAT1:FOXM1:STAT5A 1 1
#[Mon Jun  8 11:19:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150608_111944.add_feature.KBMC7 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150608_111946.add_feature.I7VRX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150608_111948.add_feature.13DCF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150608_111950.add_feature.MCHS8 sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Mon Jun  8 11:19:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150608_111952.add_feature.BTD03 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:19:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150608_111954.add_feature.C8214 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150608_111956.add_feature.76C2F sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150608_111958.add_feature.320FG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150608_112000.add_feature.OG9N0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150608_112002.add_feature.O2BUH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150608_112004.add_feature.HUZXS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150608_112006.add_feature.W8WQG sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Mon Jun  8 11:20:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150608_112008.add_feature.XIHT9 sydh:uw:haib:uta:embl pepr CHD2:MXI1:PML:ZNF384 1 1
#[Mon Jun  8 11:20:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150608_112010.add_feature.PQWCK sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150608_112012.add_feature.D8Z4Y sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150608_112014.add_feature.WONNS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150608_112016.add_feature.TID68 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150608_112018.add_feature.6M2L4 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150608_112020.add_feature.CMQAZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150608_112022.add_feature.8K7CA sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150608_112024.add_feature.XF3QC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150608_112026.add_feature.GAKVR sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Mon Jun  8 11:20:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150608_112028.add_feature.BNCXM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 11:20:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150608_112030.add_feature.T4Q57 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150608_143014.add_feature.9NSH6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150608_143016.add_feature.4O5AP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150608_143018.add_feature.O6FZV sydh:uw:haib:uta:embl pepr NFKB:NFATC1:TBLR1:POU2F2:PML:TBP:FOXM1:PAX5:POL24H8:CEBPB:BCL3:ATF2:MTA3:PAX5N19:CHD2:SRF:STAT3:BCLAF1 1 1
#[Mon Jun  8 14:30:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150608_143020.add_feature.MBB3H sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150608_143022.add_feature.UUAIJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150608_143024.add_feature.LOBKB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150608_143026.add_feature.CR7XJ sydh:uw:haib:uta:embl pepr USF1:SP1:MAZ:CHD2:NFKB:CDP:TBP:MAX 1 1
#[Mon Jun  8 14:30:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150608_143028.add_feature.E5PE9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150608_143031.add_feature.J5Y0L sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Mon Jun  8 14:30:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150608_143032.add_feature.PKS24 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150608_143034.add_feature.3L053 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150608_143037.add_feature.A2MI3 sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Mon Jun  8 14:30:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150608_143039.add_feature.HWR6O sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150608_143041.add_feature.9EV3G sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150608_143043.add_feature.2CMPO sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Mon Jun  8 14:30:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150608_143045.add_feature.84U99 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150608_143047.add_feature.OFG4D sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Mon Jun  8 14:30:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150608_143049.add_feature.PYC5O sydh:uw:haib:uta:embl pepr MXI1:PML 1 1
#[Mon Jun  8 14:30:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150608_143051.add_feature.89YSJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:30:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150608_143053.add_feature.I84EC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150608_143055.add_feature.Z23WT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150608_143057.add_feature.7Q8KM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150608_143059.add_feature.Z5BEE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150608_143101.add_feature.W93SU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150608_143103.add_feature.CNSDV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150608_143105.add_feature.36HSV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150608_143107.add_feature.OLM20 sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Mon Jun  8 14:31:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150608_143109.add_feature.U465I sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Mon Jun  8 14:31:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150608_143111.add_feature.1GRHF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150608_143113.add_feature.T6Q6T sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150608_143115.add_feature.91ZY8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150608_143117.add_feature.PE1CJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150608_143119.add_feature.E037W sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150608_143121.add_feature.FS4M3 sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Mon Jun  8 14:31:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150608_143123.add_feature.CHEX2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150608_143125.add_feature.C5IUD sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150608_143127.add_feature.6LM90 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150608_143129.add_feature.BUUSC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150608_143131.add_feature.WHTUQ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150608_143133.add_feature.Q8B24 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150608_143135.add_feature.5F5J8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150608_143137.add_feature.4C121 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150608_143139.add_feature.LX3H1 sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Mon Jun  8 14:31:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150608_143141.add_feature.8M4J3 sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun  8 14:31:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150608_143143.add_feature.SSVN2 sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun  8 14:31:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150608_143146.add_feature.WYEJB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150608_143147.add_feature.U2NKY sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150608_143149.add_feature.Z0Z12 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:31:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150608_143152.add_feature.T4DOB sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Mon Jun  8 14:31:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150608_143154.add_feature.FRDOL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150608_143156.add_feature.266WZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150608_143158.add_feature.MH2P7 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150608_143200.add_feature.BJUEG sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA:CEBPB:ATF2:MTA3 1 1
#[Mon Jun  8 14:32:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150608_143202.add_feature.HYE70 sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Mon Jun  8 14:32:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150608_143204.add_feature.TK1YJ sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Mon Jun  8 14:32:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150608_143206.add_feature.WAN9H sydh:uw:haib:uta:embl pepr PML:CHD2:YY1:POL2:ZNF143:MAX:SMC3:MAZ:TBP:ELK1:SIN3A:WHIP:ELF1:BCLAF1:RAD21:POU2F2:P300:ETS1:CMYC:GR:MXI1:NFYB:CHD1:TBLR1:ZEB1:TAF1:SP1:GABP:STAT1:FOXM1:STAT5A 1 1
#[Mon Jun  8 14:32:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150608_143208.add_feature.A8RGI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150608_143210.add_feature.T1HW8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150608_143212.add_feature.KFE7O sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150608_143214.add_feature.3ABOJ sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Mon Jun  8 14:32:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150608_143216.add_feature.D7Q0B sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150608_143218.add_feature.99KH4 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150608_143220.add_feature.URS3C sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150608_143222.add_feature.7B751 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150608_143224.add_feature.S7LVU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150608_143226.add_feature.7QEJR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150608_143228.add_feature.RGN78 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150608_143230.add_feature.82U5V sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Mon Jun  8 14:32:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150608_143232.add_feature.PAELM sydh:uw:haib:uta:embl pepr CHD2:MXI1:PML:ZNF384 1 1
#[Mon Jun  8 14:32:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150608_143234.add_feature.KJETH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150608_143236.add_feature.DEWEH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150608_143238.add_feature.YJH1V sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150608_143240.add_feature.Q6V5N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150608_143242.add_feature.0GTBW sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150608_143244.add_feature.AQXXM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150608_143246.add_feature.68JB9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150608_143248.add_feature.J3QP9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:32:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150608_143250.add_feature.SBZ8U sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Mon Jun  8 14:32:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150608_143252.add_feature.A9UHP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun  8 14:33:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150608_143254.add_feature.PX42T sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:50:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150609_115045.add_feature.H5OD0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:50:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150609_115047.add_feature.P5TPY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:50:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150609_115049.add_feature.ZSRK7 sydh:uw:haib:uta:embl pepr NFKB:NFATC1:TBLR1:POU2F2:PML:TBP:FOXM1:PAX5:POL24H8:CEBPB:BCL3:ATF2:MTA3:PAX5N19:CHD2:SRF:STAT3:BCLAF1 1 1
#[Tue Jun  9 11:50:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150609_115051.add_feature.KGCCI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:50:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150609_115053.add_feature.XFSWY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:50:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150609_115055.add_feature.6D5LJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:50:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150609_115057.add_feature.HQ6UF sydh:uw:haib:uta:embl pepr USF1:SP1:MAZ:CHD2:NFKB:CDP:TBP:MAX 1 1
#[Tue Jun  9 11:50:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150609_115059.add_feature.PET9R sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150609_115101.add_feature.6QOGD sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Tue Jun  9 11:51:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150609_115103.add_feature.XZ8WI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150609_115105.add_feature.5QPHY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150609_115107.add_feature.XVSLR sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Tue Jun  9 11:51:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150609_115109.add_feature.BHURZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150609_115111.add_feature.ZP5TM sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150609_115112.add_feature.619KT sydh:uw:haib:uta:embl pepr STAT3:MAX:CJUN:CMYC 1 1
#[Tue Jun  9 11:51:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150609_115113.add_feature.EPQ62 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Tue Jun  9 11:51:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150609_115113.add_feature.J4FVY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150609_115114.add_feature.W7IC2 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150609_115115.add_feature.7PQEU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150609_115115.add_feature.CS6BP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150609_115116.add_feature.LO7R7 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150609_115117.add_feature.FLJ91 sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Tue Jun  9 11:51:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150609_115117.add_feature.KNMX2 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150609_115118.add_feature.YHJDD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150609_115119.add_feature.6NKI3 sydh:uw:haib:uta:embl pepr MXI1:PML 1 1
#[Tue Jun  9 11:51:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150609_115119.add_feature.BZF41 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150609_115120.add_feature.9BAK1 sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Tue Jun  9 11:51:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150609_115121.add_feature.X6SYJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150609_115121.add_feature.NMPFM sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150609_115122.add_feature.CPKFN sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Tue Jun  9 11:51:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150609_115123.add_feature.V3WYP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150609_115123.add_feature.K6L07 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150609_115124.add_feature.IVJRT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150609_115125.add_feature.3J0US sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150609_115125.add_feature.QXIMS sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150609_115126.add_feature.9ELFV sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Tue Jun  9 11:51:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150609_115127.add_feature.300HK sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150609_115127.add_feature.D0GQM sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:29 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150609_115128.add_feature.0KKJQ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150609_115129.add_feature.KLKBB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150609_115129.add_feature.KN2V4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150609_115130.add_feature.GCS13 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150609_115131.add_feature.Z9TP4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150609_115131.add_feature.HRULO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150609_115132.add_feature.8PKVA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150609_115133.add_feature.93VS5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150609_115133.add_feature.YFBCD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:35 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150609_115134.add_feature.LUXL0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150609_115135.add_feature.AZU2D sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150609_115135.add_feature.HK0XC sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150609_115136.add_feature.DTSBZ sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Tue Jun  9 11:51:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150609_115137.add_feature.ESB1U sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Tue Jun  9 11:51:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150609_115137.add_feature.YBD6O sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150609_115138.add_feature.39PSG sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Tue Jun  9 11:51:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150609_115139.add_feature.JE980 sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Tue Jun  9 11:51:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150609_115140.add_feature.LG92E sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150609_115141.add_feature.79KB6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150609_115141.add_feature.K9I5A sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150609_115142.add_feature.ZAPZ4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150609_115143.add_feature.J99W4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150609_115143.add_feature.5LJK0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150609_115144.add_feature.W389P sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150609_115145.add_feature.V3MNU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150609_115145.add_feature.NQFW1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150609_115146.add_feature.7F0XD sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Tue Jun  9 11:51:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150609_115147.add_feature.QAHLR sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Tue Jun  9 11:51:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150609_115147.add_feature.T7ADE sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150609_115148.add_feature.VL0V4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150609_115149.add_feature.AUA89 sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Tue Jun  9 11:51:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150609_115149.add_feature.25C62 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150609_115150.add_feature.06VMB sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Tue Jun  9 11:51:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150609_115151.add_feature.JK78Q sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150609_115151.add_feature.V02S7 sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Tue Jun  9 11:51:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150609_115152.add_feature.PPPK0 sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Tue Jun  9 11:51:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150609_115153.add_feature.SOYMY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150609_115153.add_feature.PKKKF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150609_115154.add_feature.NT6UP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150609_115155.add_feature.TKV8F sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150609_115155.add_feature.15ZYA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150609_115156.add_feature.TJ1VA sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Tue Jun  9 11:51:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150609_115157.add_feature.DII2A sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150609_115157.add_feature.34QDU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:51:59 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150609_115158.add_feature.E8DJN sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Tue Jun  9 11:52:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150609_115159.add_feature.XGWSW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150609_115159.add_feature.I6ZGX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150609_115200.add_feature.QSYU6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150609_115201.add_feature.OFBO1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150609_115201.add_feature.45I0Y sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150609_115202.add_feature.QIW7E sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Tue Jun  9 11:52:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150609_115203.add_feature.9F5JI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150609_115203.add_feature.ZYBF6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150609_115204.add_feature.92S2M sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150609_115205.add_feature.3AGPY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150609_115207.add_feature.218CH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150609_115209.add_feature.X15A6 sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Tue Jun  9 11:52:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150609_115211.add_feature.GZJO8 sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Tue Jun  9 11:52:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150609_115213.add_feature.VXPDQ sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Tue Jun  9 11:52:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150609_115215.add_feature.BTS7O sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150609_115217.add_feature.T8IHU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150609_115219.add_feature.BP00B sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150609_115221.add_feature.18ISM sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Tue Jun  9 11:52:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150609_115223.add_feature.2807H sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150609_115225.add_feature.PQH0F sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150609_115227.add_feature.YB0QT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150609_115229.add_feature.RZ73V sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA:CEBPB:ATF2:MTA3 1 1
#[Tue Jun  9 11:52:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150609_115231.add_feature.TMFXQ sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Tue Jun  9 11:52:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150609_115233.add_feature.1BFZC sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Tue Jun  9 11:52:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150609_115235.add_feature.GHWN8 sydh:uw:haib:uta:embl pepr PML:CHD2:YY1:POL2:ZNF143:MAX:SMC3:MAZ:TBP:ELK1:SIN3A:WHIP:ELF1:BCLAF1:RAD21:POU2F2:P300:ETS1:CMYC:GR:MXI1:NFYB:CHD1:TBLR1:ZEB1:TAF1:SP1:GABP:STAT1:FOXM1:STAT5A 1 1
#[Tue Jun  9 11:52:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150609_115237.add_feature.9HMY1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150609_115240.add_feature.L6VV0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150609_115241.add_feature.6YZZY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150609_115244.add_feature.XP9ZJ sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Tue Jun  9 11:52:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150609_115246.add_feature.I528A sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150609_115248.add_feature.Q54SP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150609_115250.add_feature.FSTQE sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150609_115252.add_feature.9X1XA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150609_115254.add_feature.WS8T5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150609_115256.add_feature.RZOWE sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:52:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150609_115258.add_feature.D9RQJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150609_115300.add_feature.JHVVW sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Tue Jun  9 11:53:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150609_115302.add_feature.A9XPZ sydh:uw:haib:uta:embl pepr CHD2:MXI1:PML:ZNF384 1 1
#[Tue Jun  9 11:53:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150609_115304.add_feature.6Z1U5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150609_115306.add_feature.V2GPZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150609_115308.add_feature.0F6AF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150609_115310.add_feature.N6Z2H sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150609_115312.add_feature.03W1O sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150609_115314.add_feature.UZ2IF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150609_115316.add_feature.SWVRD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150609_115318.add_feature.IWTVI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150609_115320.add_feature.KZCFA sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Tue Jun  9 11:53:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150609_115322.add_feature.GE4B8 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 11:53:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150609_115324.add_feature.23M9Z sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150609_143525.add_feature.X5FL4 sydh:uw:haib:uta:embl pepr STAT3:MAX:CJUN:CMYC 1 1
#[Tue Jun  9 14:35:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150609_143526.add_feature.H175Y sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150609_143527.add_feature.7C834 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:29 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150609_143528.add_feature.9Z8X8 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150609_143529.add_feature.1RN4B sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:31 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150609_143530.add_feature.N0HBD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150609_143531.add_feature.GU4XK sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150609_143532.add_feature.5E874 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150609_143534.add_feature.13T4Y sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Tue Jun  9 14:35:35 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150609_143535.add_feature.LISMP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150609_143536.add_feature.3ACN4 sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Tue Jun  9 14:35:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150609_143537.add_feature.8NPGX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:38 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150609_143538.add_feature.RUT8P sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150609_143539.add_feature.GVQ0V sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150609_143540.add_feature.MBXY3 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Tue Jun  9 14:35:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150609_143541.add_feature.JCAOH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150609_143542.add_feature.GDU03 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150609_143543.add_feature.L0745 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150609_143544.add_feature.21DZJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150609_143545.add_feature.ZWMQV sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150609_143546.add_feature.EQ7H8 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150609_143547.add_feature.F2M39 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150609_143548.add_feature.HP1DR sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150609_143549.add_feature.5IULA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150609_143550.add_feature.RXGKD sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Tue Jun  9 14:35:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150609_143551.add_feature.BPIKX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150609_143552.add_feature.4WKMS sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Tue Jun  9 14:35:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150609_143553.add_feature.PGK2Q sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150609_143554.add_feature.IPK93 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150609_143555.add_feature.ZW5HU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150609_143556.add_feature.PI1KP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150609_143557.add_feature.9VBX5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:35:59 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150609_143558.add_feature.8SFM3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150609_143559.add_feature.49ON7 sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Tue Jun  9 14:36:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150609_143600.add_feature.EUBEC sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Tue Jun  9 14:36:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150609_143601.add_feature.POSG4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150609_143602.add_feature.AD8UP sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Tue Jun  9 14:36:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150609_143603.add_feature.DODK5 sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Tue Jun  9 14:36:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150609_143604.add_feature.LX468 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150609_143605.add_feature.F4PVH sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Tue Jun  9 14:36:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150609_143606.add_feature.M6QJN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150609_143607.add_feature.X5KXJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150609_143608.add_feature.I6819 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150609_143609.add_feature.BBZD1 sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Tue Jun  9 14:36:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150609_143610.add_feature.CA4YR sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:12 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150609_143611.add_feature.148W5 sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Tue Jun  9 14:36:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150609_143612.add_feature.01KFD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150609_143613.add_feature.M8NP0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150609_143614.add_feature.JSJJQ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150609_143615.add_feature.5YRXI sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Tue Jun  9 14:36:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150609_143616.add_feature.U1QDR sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150609_143616.add_feature.ITT5E sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150609_143617.add_feature.Q2T7L sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150609_143617.add_feature.LZYAI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150609_143618.add_feature.T6LFS sydh:uw:haib:uta:embl pepr NFKB:NFATC1:TBLR1:POU2F2:PML:TBP:FOXM1:PAX5:POL24H8:CEBPB:BCL3:ATF2:MTA3:PAX5N19:CHD2:SRF:STAT3:BCLAF1 1 1
#[Tue Jun  9 14:36:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150609_143619.add_feature.WYEFG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150609_143620.add_feature.RWJ14 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150609_143621.add_feature.K51T1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150609_143622.add_feature.NH5FY sydh:uw:haib:uta:embl pepr USF1:SP1:MAZ:CHD2:NFKB:CDP:TBP:MAX 1 1
#[Tue Jun  9 14:36:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150609_143623.add_feature.TWK0U sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150609_143624.add_feature.7NDA0 sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Tue Jun  9 14:36:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150609_143625.add_feature.XWEGI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150609_143626.add_feature.URGRA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150609_143627.add_feature.FHGZT sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Tue Jun  9 14:36:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150609_143628.add_feature.C1Y2M sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150609_143629.add_feature.V1WSX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150609_143630.add_feature.VRCPZ sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Tue Jun  9 14:36:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150609_143631.add_feature.ZWG5G sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150609_143632.add_feature.ZCM5Y sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Tue Jun  9 14:36:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150609_143633.add_feature.M4724 sydh:uw:haib:uta:embl pepr MXI1:PML 1 1
#[Tue Jun  9 14:36:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150609_143634.add_feature.6SSDU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150609_143635.add_feature.X0XMG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150609_143636.add_feature.69XE1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150609_143637.add_feature.CFHCB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150609_143638.add_feature.37AZ1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150609_143639.add_feature.I4BBX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150609_143640.add_feature.G3G78 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150609_143641.add_feature.62YEN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150609_143642.add_feature.BLPJP sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Tue Jun  9 14:36:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150609_143643.add_feature.M3B1I sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Tue Jun  9 14:36:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150609_143644.add_feature.10HIH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150609_143645.add_feature.A3Y96 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150609_143647.add_feature.IW4K5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150609_143647.add_feature.MC5KP sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150609_143648.add_feature.0ZV4M sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150609_143649.add_feature.22RVY sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Tue Jun  9 14:36:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150609_143650.add_feature.0ZUHB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150609_143651.add_feature.9911H sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150609_143653.add_feature.Q7NL3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150609_143654.add_feature.DBEMY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150609_143655.add_feature.XQTLH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150609_143656.add_feature.TPELN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150609_143657.add_feature.OHL7C sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150609_143658.add_feature.9GYTX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:36:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150609_143659.add_feature.HZKKK sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Tue Jun  9 14:37:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150609_143700.add_feature.MG3BM sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Tue Jun  9 14:37:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150609_143701.add_feature.JZ5ZD sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Tue Jun  9 14:37:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150609_143702.add_feature.V9UU4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150609_143703.add_feature.FZ8M5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150609_143704.add_feature.R8K9S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150609_143705.add_feature.YJ97Z sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Tue Jun  9 14:37:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150609_143706.add_feature.AH9CY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150609_143707.add_feature.5USKJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150609_143708.add_feature.2KES1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150609_143709.add_feature.KZ9F9 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA:CEBPB:ATF2:MTA3 1 1
#[Tue Jun  9 14:37:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150609_143710.add_feature.5WS0Q sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Tue Jun  9 14:37:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150609_143711.add_feature.FVFJZ sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Tue Jun  9 14:37:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150609_143712.add_feature.0UDJC sydh:uw:haib:uta:embl pepr PML:CHD2:YY1:POL2:ZNF143:MAX:SMC3:MAZ:TBP:ELK1:SIN3A:WHIP:ELF1:BCLAF1:RAD21:POU2F2:P300:ETS1:CMYC:GR:MXI1:NFYB:CHD1:TBLR1:ZEB1:TAF1:SP1:GABP:STAT1:FOXM1:STAT5A 1 1
#[Tue Jun  9 14:37:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150609_143713.add_feature.GD3J6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150609_143714.add_feature.5YPH9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150609_143715.add_feature.5OOYG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150609_143716.add_feature.DI67H sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Tue Jun  9 14:37:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150609_143717.add_feature.X0TNZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150609_143718.add_feature.K3I0I sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150609_143719.add_feature.53SGB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150609_143720.add_feature.ZXYXB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150609_143721.add_feature.TT7LF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150609_143722.add_feature.9QJX1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150609_143723.add_feature.FI6R3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150609_143724.add_feature.Y1SBQ sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Tue Jun  9 14:37:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150609_143725.add_feature.QV33P sydh:uw:haib:uta:embl pepr CHD2:MXI1:PML:ZNF384 1 1
#[Tue Jun  9 14:37:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150609_143726.add_feature.0VSWG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150609_143727.add_feature.NQH1U sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150609_143728.add_feature.OCHZH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150609_143729.add_feature.Q0XM4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150609_143730.add_feature.WP504 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150609_143731.add_feature.HBA61 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150609_143732.add_feature.LMLLY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150609_143733.add_feature.EX22G sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150609_143734.add_feature.ASNKC sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Tue Jun  9 14:37:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150609_143735.add_feature.E9QDE sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun  9 14:37:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150609_143736.add_feature.2PGJI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:49:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150629_174813.add_feature.VJV2N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:49:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150629_174814.add_feature.0XH6B sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:49:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150629_174815.add_feature.QWFK3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:49:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150629_174816.add_feature.ICGWA sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:49:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150629_174817.add_feature.K74QL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:49:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150629_174812.add_feature.6OQSD sydh:uw:haib:uta:embl pepr STAT3:MAX:CJUN:CMYC 1 1
#[Mon Jun 29 17:50:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150629_174819.add_feature.FP59L sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:50:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150629_174821.add_feature.Z6VVL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:50:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150629_174823.add_feature.6K7RN sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Mon Jun 29 17:50:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150629_174824.add_feature.0E65R sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:50:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150629_174818.add_feature.3TTK7 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:50:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150629_174820.add_feature.AMB1T sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Mon Jun 29 17:51:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150629_174830.add_feature.G79T2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:51:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150629_174825.add_feature.GALUV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:51:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150629_174826.add_feature.EZZD0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:51:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150629_174827.add_feature.78MXE sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Mon Jun 29 17:51:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150629_174828.add_feature.HHU3R sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:51:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150629_174829.add_feature.TOMJL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:52:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150629_174831.add_feature.YK69X sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:52:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150629_174832.add_feature.W3D5G sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:52:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150629_174833.add_feature.TOVLE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:52:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150629_174834.add_feature.EAISJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:52:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150629_174835.add_feature.3AZVP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:52:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150629_174836.add_feature.XETH1 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:53:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150629_174841.add_feature.3RG4R sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:53:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150629_174842.add_feature.6NSWS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:53:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150629_174837.add_feature.DMNHS sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Mon Jun 29 17:53:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150629_174838.add_feature.OYWIG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:53:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150629_174839.add_feature.9GGCJ sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Mon Jun 29 17:53:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150629_174840.add_feature.1CZZ3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:54:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150629_174843.add_feature.L5I9M sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:54:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150629_174844.add_feature.J78HH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:54:21 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150629_174845.add_feature.1E883 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:54:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150629_174846.add_feature.9LNCG sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Mon Jun 29 17:54:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150629_174847.add_feature.76BM2 sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Mon Jun 29 17:55:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150629_174848.add_feature.QWTER sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:55:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150629_174849.add_feature.05BIG sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Mon Jun 29 17:55:26 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150629_174850.add_feature.K8VCR sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Mon Jun 29 17:56:27 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150629_174851.add_feature.X8SLR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:57:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150629_174848.add_feature.T79GC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:57:28 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150629_174852.add_feature.E5T2R sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Mon Jun 29 17:58:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150629_174849.add_feature.CDF27 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:58:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150629_174850.add_feature.JAP0O sydh:uw:haib:uta:embl pepr NFKB:NFATC1:TBLR1:POU2F2:PML:TBP:FOXM1:PAX5:POL24H8:CEBPB:BCL3:ATF2:MTA3:PAX5N19:CHD2:SRF:STAT3:BCLAF1 1 1
#[Mon Jun 29 17:58:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150629_174853.add_feature.98SW4 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:59:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150629_174852.add_feature.AP8HF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:59:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150629_174855.add_feature.3RZIB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 17:59:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150629_174851.add_feature.B0P3D sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:00:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150629_174856.add_feature.SZN5E sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:00:23 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150629_174857.add_feature.B7NEE sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Mon Jun 29 18:00:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150629_174853.add_feature.T3PBP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:01:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150629_174858.add_feature.OEQO9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:02:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150629_174854.add_feature.TDGDY sydh:uw:haib:uta:embl pepr USF1:SP1:MAZ:CHD2:NFKB:CDP:TBP:MAX 1 1
#[Mon Jun 29 18:02:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150629_174855.add_feature.RGTRM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:02:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150629_174859.add_feature.JQ6G1 sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Mon Jun 29 18:03:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150629_174900.add_feature.H47YD sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:03:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150629_174856.add_feature.EGF8S sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Mon Jun 29 18:04:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150629_174901.add_feature.KX0TB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:04:13 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150629_174857.add_feature.Q288T sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:04:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150629_174858.add_feature.BQL5L sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:05:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150629_174902.add_feature.8MW8F sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:05:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150629_174859.add_feature.YFGCY sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Mon Jun 29 18:05:39 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150629_174903.add_feature.ACSY0 sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Mon Jun 29 18:06:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150629_174900.add_feature.O1ZTJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:06:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150629_174904.add_feature.T9GM2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:06:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150629_174901.add_feature.16SC9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:07:13 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150629_174905.add_feature.9DB12 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:07:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150629_174902.add_feature.PS6EO sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Mon Jun 29 18:08:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150629_174903.add_feature.8LH7D sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:09:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150629_174904.add_feature.KY01G sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Mon Jun 29 18:10:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150629_174905.add_feature.IK56N sydh:uw:haib:uta:embl pepr MXI1:PML 1 1
#[Mon Jun 29 18:10:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150629_174906.add_feature.S100E sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:10:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150629_174907.add_feature.S1U79 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:11:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150629_174908.add_feature.YNNM3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:12:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150629_174909.add_feature.808QK sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:13:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150629_174910.add_feature.8H9Q5 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:14:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150629_174911.add_feature.GT2BR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:14:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150629_174913.add_feature.FGZYB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:15:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150629_174914.add_feature.6CV77 sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Mon Jun 29 18:15:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150629_174915.add_feature.PD0VQ sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Mon Jun 29 18:17:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150629_174916.add_feature.K7JQG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:18:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150629_174917.add_feature.AIIME sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:18:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150629_174918.add_feature.K4UQX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:19:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150629_174919.add_feature.F227N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:20:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150629_174920.add_feature.QRUVU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:21:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150629_174921.add_feature.3EZ7B sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Mon Jun 29 18:22:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150629_174922.add_feature.903AU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:22:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150629_174923.add_feature.936JR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:22:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150629_174924.add_feature.BEDWD sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:24:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150629_174925.add_feature.71O5H sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:24:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150629_174926.add_feature.6UVUB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:25:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150629_174928.add_feature.J22AV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:25:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150629_174927.add_feature.Y6ZLN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:26:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150629_174929.add_feature.N11SX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:26:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150629_174930.add_feature.KSM1B sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Mon Jun 29 18:27:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150629_174931.add_feature.MP93W sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun 29 18:28:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150629_174932.add_feature.TAENA sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun 29 18:29:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150629_174933.add_feature.TMTTB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:29:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150629_174934.add_feature.QE3PU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:30:12 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150629_174935.add_feature.AGFFV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:30:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150629_174936.add_feature.BPJNP sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Mon Jun 29 18:31:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150629_174937.add_feature.ZND13 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:31:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150629_174938.add_feature.NQ7SE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:33:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150629_174939.add_feature.1972U sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:33:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150629_174940.add_feature.MASL9 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA:CEBPB:ATF2:MTA3 1 1
#[Mon Jun 29 18:33:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150629_174942.add_feature.EC549 sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Mon Jun 29 18:34:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150629_174943.add_feature.JV70D sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Mon Jun 29 18:35:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150629_174945.add_feature.WGB9Q sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:35:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150629_174944.add_feature.EVV32 sydh:uw:haib:uta:embl pepr PML:CHD2:YY1:POL2:ZNF143:MAX:SMC3:MAZ:TBP:ELK1:SIN3A:WHIP:ELF1:BCLAF1:RAD21:POU2F2:P300:ETS1:CMYC:GR:MXI1:NFYB:CHD1:TBLR1:ZEB1:TAF1:SP1:GABP:STAT1:FOXM1:STAT5A 1 1
#[Mon Jun 29 18:36:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150629_174946.add_feature.VSS4N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:37:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150629_174947.add_feature.E8ZNE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:37:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150629_174948.add_feature.YU4FF sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Mon Jun 29 18:38:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150629_174949.add_feature.KIS5E sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:39:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150629_174950.add_feature.ACYPT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:39:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150629_174951.add_feature.LQRDX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:40:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150629_174952.add_feature.ME8NK sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:41:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150629_174954.add_feature.UV5F3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:41:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150629_174953.add_feature.U8QHD sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:42:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150629_174955.add_feature.M0V0M sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:42:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150629_174956.add_feature.F2FNT sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Mon Jun 29 18:43:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150629_174957.add_feature.X5NIB sydh:uw:haib:uta:embl pepr CHD2:MXI1:PML:ZNF384 1 1
#[Mon Jun 29 18:44:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150629_174958.add_feature.AVSXP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:45:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150629_175000.add_feature.89V1R sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:45:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150629_174959.add_feature.U7OZN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:45:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150629_175001.add_feature.V8QLU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:46:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150629_175002.add_feature.W13K0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:47:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150629_175003.add_feature.33LPW sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:48:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150629_175004.add_feature.NZEL1 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:48:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150629_175006.add_feature.YNXML sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Mon Jun 29 18:48:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150629_175005.add_feature.8O5WG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:49:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150629_175007.add_feature.M2FUV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 18:50:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150629_175008.add_feature.0OJU8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150629_215800.add_feature.LZFL3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150629_215801.add_feature.FISLU sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150629_215802.add_feature.VRG1Q sydh:uw:haib:uta:embl pepr NFKB:NFATC1:TBLR1:POU2F2:PML:TBP:FOXM1:PAX5:POL24H8:CEBPB:BCL3:ATF2:MTA3:PAX5N19:CHD2:SRF:STAT3:BCLAF1 1 1
#[Mon Jun 29 21:58:11 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150629_215803.add_feature.GS9PT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150629_215804.add_feature.ZHGOX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150629_215805.add_feature.79H57 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150629_215806.add_feature.BBWWI sydh:uw:haib:uta:embl pepr USF1:SP1:MAZ:CHD2:NFKB:CDP:TBP:MAX 1 1
#[Mon Jun 29 21:58:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150629_215807.add_feature.O7V19 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150629_215808.add_feature.HJPGL sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Mon Jun 29 21:58:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150629_215809.add_feature.TFCFN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150629_215812.add_feature.4XSPO sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150629_215810.add_feature.78BAT sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150629_215813.add_feature.CSZIW sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150629_215811.add_feature.UCGJX sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Mon Jun 29 21:58:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150629_215814.add_feature.7DF3A sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Mon Jun 29 21:58:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150629_215816.add_feature.M6MG6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150629_215826.add_feature.NVIO2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:30 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150629_215826.add_feature.UB86Z sydh:uw:haib:uta:embl pepr STAT3:MAX:CJUN:CMYC 1 1
#[Mon Jun 29 21:58:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150629_215824.add_feature.JUYWI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150629_215818.add_feature.K5179 sydh:uw:haib:uta:embl pepr MXI1:PML 1 1
#[Mon Jun 29 21:58:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150629_215822.add_feature.LZK0D sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150629_215830.add_feature.IV0QZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150629_215820.add_feature.BPCZP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150629_215829.add_feature.YENPC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150629_215817.add_feature.WXUIG sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Mon Jun 29 21:58:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150629_215823.add_feature.DI6CP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:36 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150629_215827.add_feature.K4BNR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150629_215819.add_feature.HA322 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150629_215829.add_feature.RX6QF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150629_215828.add_feature.ANQMF sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Mon Jun 29 21:58:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150629_215831.add_feature.FGF04 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150629_215825.add_feature.V6XHZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150629_215827.add_feature.VUDLI sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Mon Jun 29 21:58:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150629_215839.add_feature.MQERG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150629_215834.add_feature.GTJCD sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Mon Jun 29 21:58:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150629_215834.add_feature.0Y91T sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Mon Jun 29 21:58:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150629_215838.add_feature.Y8PIV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150629_215830.add_feature.DM71E sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150629_215828.add_feature.QEWRR sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150629_215821.add_feature.PG5ZS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150629_215835.add_feature.WOGPI sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150629_215835.add_feature.BNSYV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150629_215833.add_feature.ROE7Z sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150629_215833.add_feature.Q1K9G sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150629_215838.add_feature.5QTQP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150629_215831.add_feature.IOCXA sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150629_215841.add_feature.ABIQG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150629_215836.add_feature.81VKU sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Mon Jun 29 21:58:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150629_215832.add_feature.H9S0Z sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150629_215832.add_feature.7ZHX5 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150629_215837.add_feature.PGKY9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150629_215840.add_feature.XB9IQ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150629_215839.add_feature.BT42W sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150629_215836.add_feature.EK67D sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150629_215844.add_feature.UHJ5J sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun 29 21:58:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150629_215842.add_feature.YICGS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150629_215844.add_feature.CUQFN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150629_215842.add_feature.XS9FE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150629_215837.add_feature.FXFJL sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150629_215845.add_feature.P4GGK sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Mon Jun 29 21:58:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150629_215840.add_feature.9RWOI sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Mon Jun 29 21:58:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150629_215853.add_feature.FO31G sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Mon Jun 29 21:58:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150629_215850.add_feature.NX4AM sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150629_215841.add_feature.7KXSG sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150629_215843.add_feature.ZCVB7 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150629_215843.add_feature.TLV8K sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Mon Jun 29 21:58:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150629_215845.add_feature.4FKRP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150629_215848.add_feature.1YVKH sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150629_215846.add_feature.OMSB9 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150629_215848.add_feature.FNICB sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150629_215850.add_feature.WSFDE sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Mon Jun 29 21:58:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150629_215847.add_feature.I5X3H sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:58:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150629_215856.add_feature.UCN2A sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Mon Jun 29 21:58:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150629_215846.add_feature.KAA0W sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150629_215847.add_feature.V6WXC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150629_215849.add_feature.DNA9E sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Mon Jun 29 21:59:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150629_215849.add_feature.IDGQ7 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150629_215854.add_feature.3PNWN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150629_215855.add_feature.YCPBN sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150629_215855.add_feature.NRF12 sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Mon Jun 29 21:59:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150629_215857.add_feature.EY9BZ sydh:uw:haib:uta:embl pepr PML:CHD2:YY1:POL2:ZNF143:MAX:SMC3:MAZ:TBP:ELK1:SIN3A:WHIP:ELF1:BCLAF1:RAD21:POU2F2:P300:ETS1:CMYC:GR:MXI1:NFYB:CHD1:TBLR1:ZEB1:TAF1:SP1:GABP:STAT1:FOXM1:STAT5A 1 1
#[Mon Jun 29 21:59:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150629_215852.add_feature.QVXBE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150629_215854.add_feature.EZXZB sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA:CEBPB:ATF2:MTA3 1 1
#[Mon Jun 29 21:59:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150629_215853.add_feature.FEPR0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:04 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150629_215852.add_feature.4H68Z sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150629_215858.add_feature.I6U4E sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150629_215858.add_feature.O6RQJ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150629_215859.add_feature.1SYD3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150629_215902.add_feature.EROT8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150629_215857.add_feature.KDJHF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150629_215906.add_feature.NQZJG sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Mon Jun 29 21:59:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150629_215905.add_feature.O0O0E sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150629_215901.add_feature.34Q10 sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Mon Jun 29 21:59:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150629_215859.add_feature.FFPA0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150629_215902.add_feature.AEG7G sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:11 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150629_215856.add_feature.2NIO1 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150629_215901.add_feature.PDNTO sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Mon Jun 29 21:59:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150629_215904.add_feature.F50CV sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:14 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150629_215906.add_feature.CTI1S sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:14 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150629_215904.add_feature.8L2HH sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Mon Jun 29 21:59:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150629_215900.add_feature.SU2OX sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:15 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150629_215903.add_feature.4L2Y3 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:15 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150629_215905.add_feature.9WKY2 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150629_215910.add_feature.Q3TP3 sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Mon Jun 29 21:59:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150629_215903.add_feature.CQS7Q sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Mon Jun 29 21:59:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150629_215911.add_feature.P7SAC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:18 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150629_215910.add_feature.PUQC7 sydh:uw:haib:uta:embl pepr CHD2:MXI1:PML:ZNF384 1 1
#[Mon Jun 29 21:59:18 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150629_215908.add_feature.YKUQK sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150629_215912.add_feature.9NSIO sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:19 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150629_215909.add_feature.VF2SU sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Mon Jun 29 21:59:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150629_215900.add_feature.9WGR8 sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Mon Jun 29 21:59:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150629_215907.add_feature.SKJ9Y sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150629_215911.add_feature.EETC0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150629_215908.add_feature.7NUMC sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:20 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150629_215915.add_feature.LWZE6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150629_215907.add_feature.NH1E6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150629_215916.add_feature.9EIKW sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Mon Jun 29 21:59:22 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150629_215909.add_feature.6YFU0 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150629_215914.add_feature.OG71N sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150629_215916.add_feature.IAXF8 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:24 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150629_215913.add_feature.MXRVZ sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:25 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150629_215914.add_feature.PJHDS sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150629_215917.add_feature.WT7JF sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150629_215913.add_feature.KP0TE sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150629_215918.add_feature.KRPN6 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150629_215917.add_feature.QL02O sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150629_215920.add_feature.N2U09 sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150629_215912.add_feature.Q8I1M sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Mon Jun 29 21:59:33 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150629_215918.add_feature.6TIEP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150629_215915.add_feature.MX0EP sydh:uw:haib:uta:embl pepr None 1 1
#[Mon Jun 29 21:59:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150629_215919.add_feature.TIT1U sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Mon Jun 29 21:59:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150629_215921.add_feature.HDGO8 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:20 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf2 20150630_094419.add_feature.ZWUSY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:21 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-atf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 atf3 20150630_094420.add_feature.ZA9CE sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:22 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-batf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 batf 20150630_094421.add_feature.E1QJQ sydh:uw:haib:uta:embl pepr NFKB:NFATC1:TBLR1:POU2F2:PML:TBP:FOXM1:PAX5:POL24H8:CEBPB:BCL3:ATF2:MTA3:PAX5N19:CHD2:SRF:STAT3:BCLAF1 1 1
#[Tue Jun 30 09:44:23 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl11a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl11a 20150630_094422.add_feature.Y7XSB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:24 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bcl3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bcl3 20150630_094423.add_feature.RLRVK sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:25 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bclaf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bclaf1 20150630_094424.add_feature.6XHUI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:26 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-bhlhe40-shi p_add_feature_on_loc_dp.py gm12878 gm12878 bhlhe40 20150630_094425.add_feature.1E0SG sydh:uw:haib:uta:embl pepr USF1:SP1:MAZ:CHD2:NFKB:CDP:TBP:MAX 1 1
#[Tue Jun 30 09:44:27 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-brca1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 brca1 20150630_094426.add_feature.Q1EDR sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:28 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cebpb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cebpb 20150630_094427.add_feature.5NUN3 sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Tue Jun 30 09:44:29 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cfos-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cfos 20150630_094428.add_feature.FSTR0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:30 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd1 20150630_094429.add_feature.90KJI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:31 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-chd2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 chd2 20150630_094430.add_feature.EBZW1 sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Tue Jun 30 09:44:32 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-cmyc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 cmyc 20150630_094431.add_feature.GG14T sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:33 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-corest-shi p_add_feature_on_loc_dp.py gm12878 gm12878 corest 20150630_094432.add_feature.1T3GK sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:34 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ctcf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ctcf 20150630_094433.add_feature.986OS sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Tue Jun 30 09:44:35 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-e2f4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 e2f4 20150630_094434.add_feature.B8SON sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:36 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ebf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ebf1 20150630_094435.add_feature.C3CP2 sydh:uw:haib:uta:embl pepr PU1:PAX5:TCF12:IRF4:MEF2A:P300:NFKB:BCL3 1 1
#[Tue Jun 30 09:44:37 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-egr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 egr1 20150630_094436.add_feature.JVTBP sydh:uw:haib:uta:embl pepr MXI1:PML 1 1
#[Tue Jun 30 09:44:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elf1 20150630_094437.add_feature.Z2KXF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-elk1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 elk1 20150630_094438.add_feature.APMB6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ets1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ets1 20150630_094439.add_feature.9Q1I7 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-foxm1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 foxm1 20150630_094440.add_feature.4H90I sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ap2gamma-shi p_add_feature_on_loc_dp.py helas3 helas3 ap2gamma 20150630_094440.add_feature.ZVM3Q sydh:uw:haib:uta:embl pepr STAT3:MAX:CJUN:CMYC 1 1
#[Tue Jun 30 09:44:42 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-gabp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 gabp 20150630_094441.add_feature.M3H36 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:42 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf155-shi p_add_feature_on_loc_dp.py helas3 helas3 baf155 20150630_094441.add_feature.QOLYK sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-ikzf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 ikzf1 20150630_094442.add_feature.4N4FF sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-baf170-shi p_add_feature_on_loc_dp.py helas3 helas3 baf170 20150630_094442.add_feature.X1ONW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:44 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-input-shi p_add_feature_on_loc_dp.py gm12878 gm12878 input 20150630_094443.add_feature.JICIH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-bdp1-shi p_add_feature_on_loc_dp.py helas3 helas3 bdp1 20150630_094443.add_feature.0O8AQ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-irf4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 irf4 20150630_094444.add_feature.SR6O9 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brca1-shi p_add_feature_on_loc_dp.py helas3 helas3 brca1 20150630_094444.add_feature.MWYM5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-jund-shi p_add_feature_on_loc_dp.py gm12878 gm12878 jund 20150630_094445.add_feature.Y3ZEL sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Tue Jun 30 09:44:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf1-shi p_add_feature_on_loc_dp.py helas3 helas3 brf1 20150630_094445.add_feature.7Y4DG sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:47 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-max-shi p_add_feature_on_loc_dp.py gm12878 gm12878 max 20150630_094446.add_feature.738YV sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Tue Jun 30 09:44:47 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brf2-shi p_add_feature_on_loc_dp.py helas3 helas3 brf2 20150630_094446.add_feature.OTXVR sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:48 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-maz-shi p_add_feature_on_loc_dp.py gm12878 gm12878 maz 20150630_094447.add_feature.XX8B6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:48 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-brg1-shi p_add_feature_on_loc_dp.py helas3 helas3 brg1 20150630_094448.add_feature.4RCFT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:49 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2a 20150630_094448.add_feature.9L79S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:49 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cebpb-shi p_add_feature_on_loc_dp.py helas3 helas3 cebpb 20150630_094448.add_feature.0BE9B sydh:uw:haib:uta:embl pepr CTCF:P300:SMC3:RAD21 1 1
#[Tue Jun 30 09:44:50 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mef2c-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mef2c 20150630_094449.add_feature.WM2GC sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:50 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cfos-shi p_add_feature_on_loc_dp.py helas3 helas3 cfos 20150630_094450.add_feature.YXRNN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:51 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mta3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mta3 20150630_094450.add_feature.D4XAZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:51 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-chd2-shi p_add_feature_on_loc_dp.py helas3 helas3 chd2 20150630_094451.add_feature.L2FA4 sydh:uw:haib:uta:embl pepr BRCA1 1 1
#[Tue Jun 30 09:44:52 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-mxi1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 mxi1 20150630_094451.add_feature.U1E96 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:52 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cjun-shi p_add_feature_on_loc_dp.py helas3 helas3 cjun 20150630_094452.add_feature.EAZEU sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:53 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-myc-shi p_add_feature_on_loc_dp.py gm12878 gm12878 myc 20150630_094452.add_feature.5U48S sydh:uw:haib:uta:embl pepr MAX:BRCA1:YY1:NFYB:TFAP2A:MYC-MAX 1 1
#[Tue Jun 30 09:44:53 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-cmyc-shi p_add_feature_on_loc_dp.py helas3 helas3 cmyc 20150630_094453.add_feature.HMJ5X sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:54 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfatc1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfatc1 20150630_094453.add_feature.93J7Y sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:54 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-corest-shi p_add_feature_on_loc_dp.py helas3 helas3 corest 20150630_094454.add_feature.BLWW5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:55 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ctcf-shi p_add_feature_on_loc_dp.py helas3 helas3 ctcf 20150630_094455.add_feature.YQFE1 sydh:uw:haib:uta:embl pepr ZNF143:SMC3:YY1:RAD21:XXX:POU2F2:ZNF384:RUNX3:ATF2 1 1
#[Tue Jun 30 09:44:55 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfe2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfe2 20150630_094454.add_feature.ZCRF6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:56 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfic-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfic 20150630_094455.add_feature.B5WNN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:56 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f1 20150630_094456.add_feature.DOHL6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:57 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfkb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfkb 20150630_094456.add_feature.CHPFI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:57 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f4-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f4 20150630_094457.add_feature.6KGC8 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:58 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-e2f6-shi p_add_feature_on_loc_dp.py helas3 helas3 e2f6 20150630_094458.add_feature.UAMME sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:58 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfya-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfya 20150630_094458.add_feature.B9FTA sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nfyb-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nfyb 20150630_094459.add_feature.GO1YT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:44:59 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk1-shi p_add_feature_on_loc_dp.py helas3 helas3 elk1 20150630_094459.add_feature.ZNA6Z sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:00 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrf1 20150630_094500.add_feature.4KFRN sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:00 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-elk4-shi p_add_feature_on_loc_dp.py helas3 helas3 elk4 20150630_094500.add_feature.7NVZ4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:01 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-nrsf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 nrsf 20150630_094501.add_feature.WAJQ7 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:01 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gabp-shi p_add_feature_on_loc_dp.py helas3 helas3 gabp 20150630_094501.add_feature.00XS1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:02 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-p300-shi p_add_feature_on_loc_dp.py gm12878 gm12878 p300 20150630_094502.add_feature.92Q0R sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Tue Jun 30 09:45:02 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-gtf2f1-shi p_add_feature_on_loc_dp.py helas3 helas3 gtf2f1 20150630_094502.add_feature.P0CUO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:03 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5 20150630_094503.add_feature.BBVSI sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Tue Jun 30 09:45:03 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-ini1-shi p_add_feature_on_loc_dp.py helas3 helas3 ini1 20150630_094503.add_feature.JFW4N sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:04 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5c20-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5c20 20150630_094504.add_feature.4573M sydh:uw:haib:uta:embl pepr POU2F2:MXI1:CHD2:MAX:POL24H8:PML:POL2:ELK1:WHIP 1 1
#[Tue Jun 30 09:45:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-irf3-shi p_add_feature_on_loc_dp.py helas3 helas3 irf3 20150630_094504.add_feature.7WZ95 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pax5n19-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pax5n19 20150630_094505.add_feature.F5HWB sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-jund-shi p_add_feature_on_loc_dp.py helas3 helas3 jund 20150630_094505.add_feature.DUB2F sydh:uw:haib:uta:embl pepr CMYC 1 1
#[Tue Jun 30 09:45:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pbx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pbx3 20150630_094506.add_feature.86NRX sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:07 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mafk-shi p_add_feature_on_loc_dp.py helas3 helas3 mafk 20150630_094506.add_feature.3A8P2 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:08 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pml-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pml 20150630_094507.add_feature.I4H6S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:08 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-max-shi p_add_feature_on_loc_dp.py helas3 helas3 max 20150630_094507.add_feature.ZGKQ8 sydh:uw:haib:uta:embl pepr CMYC:USF2 1 1
#[Tue Jun 30 09:45:09 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-maz-shi p_add_feature_on_loc_dp.py helas3 helas3 maz 20150630_094508.add_feature.TNJ7W sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:09 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol2 20150630_094508.add_feature.NAPLS sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Tue Jun 30 09:45:10 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol24h8-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol24h8 20150630_094509.add_feature.DVUU2 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:10 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-mxi1-shi p_add_feature_on_loc_dp.py helas3 helas3 mxi1 20150630_094509.add_feature.U45WT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:59 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrf1-shi p_add_feature_on_loc_dp.py helas3 helas3 nrf1 20150630_094557.add_feature.WV2T4 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pu1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pu1 20150630_094557.add_feature.6WIF9 sydh:uw:haib:uta:embl pepr JUN-FOS:GATA2:CEBPA:CEBPB:ATF2:MTA3 1 1
#[Tue Jun 30 09:45:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pou2f2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pou2f2 20150630_094557.add_feature.M417F sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:59 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfyb-shi p_add_feature_on_loc_dp.py helas3 helas3 nfyb 20150630_094557.add_feature.TZQDY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:59 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nfya-shi p_add_feature_on_loc_dp.py helas3 helas3 nfya 20150630_094557.add_feature.7TX3Z sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:45:59 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-pol3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 pol3 20150630_094557.add_feature.KTLWZ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-nrsf-shi p_add_feature_on_loc_dp.py helas3 helas3 nrsf 20150630_094602.add_feature.OAGVJ sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:05 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-pol2-shi p_add_feature_on_loc_dp.py helas3 helas3 pol2 20150630_094601.add_feature.L4OUJ sydh:uw:haib:uta:embl pepr TAF1:GCN5 1 1
#[Tue Jun 30 09:46:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rad21-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rad21 20150630_094601.add_feature.LKEE3 sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Tue Jun 30 09:46:05 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-runx3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 runx3 20150630_094602.add_feature.SBA3S sydh:uw:haib:uta:embl pepr PML:CHD2:YY1:POL2:ZNF143:MAX:SMC3:MAZ:TBP:ELK1:SIN3A:WHIP:ELF1:BCLAF1:RAD21:POU2F2:P300:ETS1:CMYC:GR:MXI1:NFYB:CHD1:TBLR1:ZEB1:TAF1:SP1:GABP:STAT1:FOXM1:STAT5A 1 1
#[Tue Jun 30 09:46:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-prdm1-shi p_add_feature_on_loc_dp.py helas3 helas3 prdm1 20150630_094602.add_feature.MABOT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:06 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-p300-shi p_add_feature_on_loc_dp.py helas3 helas3 p300 20150630_094601.add_feature.ZXIIZ sydh:uw:haib:uta:embl pepr CJUN:ELK4 1 1
#[Tue Jun 30 09:46:06 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rfx5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rfx5 20150630_094602.add_feature.RODC7 sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Tue Jun 30 09:46:07 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-rxra-shi p_add_feature_on_loc_dp.py gm12878 gm12878 rxra 20150630_094602.add_feature.5RHNO sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-spt20-shi p_add_feature_on_loc_dp.py helas3 helas3 spt20 20150630_094611.add_feature.312JT sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:16 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat1 20150630_094611.add_feature.1995S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:16 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rad21-shi p_add_feature_on_loc_dp.py helas3 helas3 rad21 20150630_094611.add_feature.22DNZ sydh:uw:haib:uta:embl pepr CTCF:SMC3:CFOS:CJUN:TCF7L2:P300:STAT3:CHD2:BAF155:TBP 1 1
#[Tue Jun 30 09:46:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rpc155-shi p_add_feature_on_loc_dp.py helas3 helas3 rpc155 20150630_094611.add_feature.DJV6S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-smc3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 smc3 20150630_094611.add_feature.ACHL7 sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Tue Jun 30 09:46:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-six5-shi p_add_feature_on_loc_dp.py gm12878 gm12878 six5 20150630_094611.add_feature.F5IGH sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sin3a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sin3a 20150630_094611.add_feature.WIBH1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-srf-shi p_add_feature_on_loc_dp.py gm12878 gm12878 srf 20150630_094611.add_feature.YP2EW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:17 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-sp1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 sp1 20150630_094611.add_feature.4FNNI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:17 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-rfx5-shi p_add_feature_on_loc_dp.py helas3 helas3 rfx5 20150630_094611.add_feature.Z8QQ1 sydh:uw:haib:uta:embl pepr GABP:ELK1 1 1
#[Tue Jun 30 09:46:19 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-smc3-shi p_add_feature_on_loc_dp.py helas3 helas3 smc3 20150630_094611.add_feature.BWXW4 sydh:uw:haib:uta:embl pepr CTCF:TAF1:POL2:CHD2:CJUN:MXI1:CFOS:RAD21:P300:CMYC 1 1
#[Tue Jun 30 09:46:32 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-stat3-shi p_add_feature_on_loc_dp.py helas3 helas3 stat3 20150630_094620.add_feature.ONHQR sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:34 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-taf1-shi p_add_feature_on_loc_dp.py helas3 helas3 taf1 20150630_094621.add_feature.3CNE6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tbp-shi p_add_feature_on_loc_dp.py helas3 helas3 tbp 20150630_094621.add_feature.WY37D sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Tue Jun 30 09:46:37 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tcf7l2-shi p_add_feature_on_loc_dp.py helas3 helas3 tcf7l2 20150630_094621.add_feature.UKCJ6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat3 20150630_094621.add_feature.8U3DI sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:38 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-taf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 taf1 20150630_094621.add_feature.0WP4R sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tbp-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tbp 20150630_094630.add_feature.JHUUT sydh:uw:haib:uta:embl pepr BDP1:RPC155:GCN5 1 1
#[Tue Jun 30 09:46:39 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tblr1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tblr1 20150630_094630.add_feature.C8YS5 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf12-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf12 20150630_094630.add_feature.8H3NI sydh:uw:haib:uta:embl pepr CHD2:MXI1:PML:ZNF384 1 1
#[Tue Jun 30 09:46:40 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-stat5a-shi p_add_feature_on_loc_dp.py gm12878 gm12878 stat5a 20150630_094621.add_feature.BK6A7 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:40 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-test-shi p_add_feature_on_loc_dp.py helas3 helas3 test 20150630_094630.add_feature.WUIOT sydh:uw:haib:uta:embl pepr YY1:XY 1 1
#[Tue Jun 30 09:46:41 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tcf3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tcf3 20150630_094630.add_feature.I6BZ1 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:41 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-tr4-shi p_add_feature_on_loc_dp.py helas3 helas3 tr4 20150630_094630.add_feature.NYIV6 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:43 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zkscan1-shi p_add_feature_on_loc_dp.py helas3 helas3 zkscan1 20150630_094640.add_feature.PSNFD sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:43 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-yy1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 yy1 20150630_094640.add_feature.LU17L sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:44 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-zzz3-shi p_add_feature_on_loc_dp.py helas3 helas3 zzz3 20150630_094640.add_feature.Z2Y9S sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-whip-shi p_add_feature_on_loc_dp.py gm12878 gm12878 whip 20150630_094640.add_feature.ZU2PY sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-tr4-shi p_add_feature_on_loc_dp.py gm12878 gm12878 tr4 20150630_094640.add_feature.3A79D sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:45 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-usf2-shi p_add_feature_on_loc_dp.py helas3 helas3 usf2 20150630_094640.add_feature.17EW0 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:45 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf2-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf2 20150630_094640.add_feature.BZEED sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf143-shi p_add_feature_on_loc_dp.py helas3 helas3 znf143 20150630_094640.add_feature.8F1H8 sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zbtb33-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zbtb33 20150630_094643.add_feature.MBA54 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-helas3-helas3-znf274-shi p_add_feature_on_loc_dp.py helas3 helas3 znf274 20150630_094643.add_feature.BB4IV sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zeb1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zeb1 20150630_094643.add_feature.U4NKW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-zzz3-shi p_add_feature_on_loc_dp.py gm12878 gm12878 zzz3 20150630_094643.add_feature.CMTLW sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf274-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf274 20150630_094644.add_feature.XV7U3 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-usf1-shi p_add_feature_on_loc_dp.py gm12878 gm12878 usf1 20150630_094640.add_feature.G8Y47 sydh:uw:haib:uta:embl pepr None 1 1
#[Tue Jun 30 09:46:46 2015] p p_run_cluster_sep.py add-feature-gm12878-gm12878-znf143-shi p_add_feature_on_loc_dp.py gm12878 gm12878 znf143 20150630_094643.add_feature.806SV sydh:uw:haib:uta:embl pepr CTCF:SMC3:RAD21:NRF1:NFYB:SIX5 1 1
