
import p_mymodule as my
import sys
import os
import time
script_dir="/homed/home/shi/anthony/chipseq_snp/python"


#fastq_gz_list=["gm12878-dnase-Rep1.fastq.gz", "gm12878-dnase-Rep2.fastq.gz"]

gm12878_gz_list=[
#"gm12878-dnase-Rep1.fastq.gz",
#"gm12878-dnase-Rep2.fastq.gz",
#"gm12878-dnase-Rep3.fastq.gz",
#"gm12878-dnase-Rep4.fastq.gz",
#"gm12878-dnase-Rep5.fastq.gz",
#"gm12878-h3k27ac-Rep1.fastq.gz",
#"gm12878-h3k27ac-Rep2.fastq.gz",#-
#"gm12878-h3k4me1-Rep1.fastq.gz",
#"gm12878-h3k4me1-Rep2.fastq.gz",#-
#"gm12878-h3k4me2-Rep1.fastq.gz",
#"gm12878-h3k4me2-Rep2.fastq.gz",
#"gm12878-h3k4me3-Rep1.fastq.gz",
#"gm12878-h3k4me3-Rep2.fastq.gz",
#"gm12878-h3k9ac-Rep1.fastq.gz",
#"gm12878-h3k9ac-Rep2.fastq.gz",
#"gm12878-h4k20me1-Rep1.fastq.gz",
#"gm12878-h4k20me1-Rep2.fastq.gz",#-
#"gm12878-bhlhe40-Rep1.fastq.gz",
#"gm12878-bhlhe40-Rep2.fastq.gz",
#"gm12878-ctcf-Rep1.fastq.gz",
#"gm12878-ctcf-Rep2.fastq.gz",
#"gm12878-ebf1-Rep1.fastq.gz",
#"gm12878-ebf1-Rep2.fastq.gz",
#"gm12878-stat1-Rep1.fastq.gz",
#"gm12878-stat1-Rep2.fastq.gz",
#"gm12878-znf143-Rep1.fastq.gz",
#"gm12878-znf143-Rep2.fastq.gz",
"uw-gm12878-dnase-Rep1.fastq.gz",
"uw-gm12878-dnase-Rep2.fastq.gz",
#"gm12878-test-Rep1.fastq.gz"
]

gm12878_gz_list2=[
#"haib-gm12878-taf1-Rep1.fastq.gz"
#"sydh-gm12878-brca1-Rep1.fastq.gz",
#"sydh-gm12878-chd2-Rep1.fastq.gz",
#"sydh-gm12878-elk1-Rep1.fastq.gz",
#"sydh-gm12878-max-Rep1.fastq.gz",
#"sydh-gm12878-maz-Rep1.fastq.gz",
#"sydh-gm12878-mxi1-Rep1.fastq.gz",
#"sydh-gm12878-nfya-Rep1.fastq.gz",
#"sydh-gm12878-nfyb-Rep1.fastq.gz",
#"sydh-gm12878-rad21-Rep1.fastq.gz",
#"sydh-gm12878-rfx5-Rep1.fastq.gz",
#"sydh-gm12878-smc3-Rep1.fastq.gz",
#"sydh-gm12878-stat3-Rep1.fastq.gz",
#"sydh-gm12878-tbp-Rep1.fastq.gz",
#"sydh-gm12878-usf2-Rep1.fastq.gz"
#"haib-gm12878-bclaf1-Rep1.fastq.gz",
#"sydh-gm12878-ezh2-Rep1.fastq.gz",
#"sydh-gm12878-ezh2-Rep2.fastq.gz"
]

gm12878_histone_list=[
"gm12878-h2az-Rep1.fastq.gz",
#"gm12878-h3k27ac-Rep1.fastq.gz",
"gm12878-h3k27me3-Rep1.fastq.gz",
"gm12878-h3k36me3-Rep1.fastq.gz",
#"gm12878-h3k4me1-Rep1.fastq.gz",
#"gm12878-h3k4me2-Rep1.fastq.gz",
#"gm12878-h3k4me3-Rep1.fastq.gz",
"gm12878-h3k79me2-Rep1.fastq.gz",
#"gm12878-h3k9ac-Rep1.fastq.gz",
"gm12878-h3k9me3-Rep1.fastq.gz",
#"gm12878-h4k20me1-Rep1.fastq.gz"
]
    
gm12xxx_list=[
#"haib-gm12891-taf1-Rep1.fastq.gz",
#"haib-gm12892-taf1-Rep1.fastq.gz",
#"open-gm12891--Rep1.fastq.gz",
#"open-gm12892--Rep1.fastq.gz",
#"uw-gm12864--Rep1.fastq.gz",
#"open-gm19238--Rep1.fastq.gz",
#"open-gm19239--Rep1.fastq.gz",
#"open-gm19240--Rep1.fastq.gz",
#"ut-gm19238-ctcf-Rep1.fastq.gz",
#"ut-gm19239-ctcf-Rep1.fastq.gz",
#"ut-gm19240-ctcf-Rep1.fastq.gz",
#"ut-gm12891-ctcf-Rep1.fastq.gz",
#"ut-gm12892-ctcf-Rep1.fastq.gz",
#"uw-gm12864-ctcf-Rep1.fastq.gz",
#"uta-gm19238-ctcf-Rep2.fastq.gz",
#"uta-gm19239-ctcf-Rep2.fastq.gz",
#"uta-gm19240-ctcf-Rep2.fastq.gz",
#"uta-gm12891-ctcf-Rep2.fastq.gz",
#"uta-gm12892-ctcf-Rep2.fastq.gz",
#"uw-gm12872-ctcf-Rep1.fastq.gz",
#"uw-gm12873-ctcf-Rep1.fastq.gz"
"uw-gm12801-input-Rep1.fastq.gz",
"uw-gm12801-ctcf-Rep1.fastq.gz"
]

helas3_gz_list=[
#"haib-helas3-taf1-Rep1.fastq.gz",
"helas3-h2az-Rep1.fastq.gz",
"helas3-h3k04me1-Rep1.fastq.gz",
"helas3-h3k09me3-Rep1.fastq.gz",
"helas3-h3k27ac-Rep1.fastq.gz",
"helas3-h3k27me3-Rep1.fastq.gz",
"helas3-h3k36me3-Rep1.fastq.gz",
"helas3-h3k4me2-Rep1.fastq.gz",
"helas3-h3k4me3-Rep1.fastq.gz",
"helas3-h3k79me2-Rep1.fastq.gz",
"helas3-h3k9ac-Rep1.fastq.gz",
"helas3-h4k20me1-Rep1.fastq.gz",
"helas3-methy-Rep1.fastq.gz",
"uw-helas3-dnase-Rep1.fastq.gz"]

helas3_tf_list=[
"sydh-helas3-brca1-Rep1.fastq.gz",
#"sydh-helas3-chd2-Rep1.fastq.gz",
"sydh-helas3-elk1-Rep1.fastq.gz",
"sydh-helas3-max-Rep1.fastq.gz",
#"sydh-helas3-maz-Rep1.fastq.gz",
#"sydh-helas3-mxi1-Rep1.fastq.gz",
"sydh-helas3-nfya-Rep1.fastq.gz",
"sydh-helas3-nfyb-Rep1.fastq.gz",
#"sydh-helas3-rad21-Rep1.fastq.gz",
"sydh-helas3-rfx5-Rep1.fastq.gz",
#"sydh-helas3-smc3-Rep1.fastq.gz",
"sydh-helas3-stat3-Rep1.fastq.gz",
"sydh-helas3-tbp-Rep1.fastq.gz",
"sydh-helas3-usf2-Rep1.fastq.gz",
"sydh-helas3-znf143-Rep1.fastq.gz",
"uw-helas3-ctcf-Rep1.fastq.gz"]





target_server=sys.argv[1]
cell_name=sys.argv[2]
test_flag = sys.argv[3]

print cell_name
if test_flag == 'test':
    fastq_gz_list = ['sydh-testcell-test-Rep1.fastq.gz']
elif "gm12878" in cell_name:
    fastq_gz_list=gm12878_gz_list
elif "gm12xxx" in cell_name:
    fastq_gz_list=gm12xxx_list
elif "helas3" in cell_name:
    #fastq_gz_list = helas3_gz_list
    file_list = my.f_shell_cmd( "ssh wenqiang@orcinus.westgrid.ca find /home/wenqiang/encode/helas3/ -name '*.fastq.gz'", quiet = True).split('\n')
    #ctcf_list_raw = [ os.path.basename(fastq_file) for fastq_file in my.grep_list('^(?!.*gm12xxx|.*/ut-|.*open-).*%s'%feature, file_list)]
    #ctcf_list = my.grep_list('uw-(gm12864|gm12873).*', ctcf_list_raw)

    tf_list_file='/homed/home/shi/projects/wgs/tf_list.txt'
    tf_list=my.f_parse_tf_list_file(tf_list_file)
    
    compiled_list = [ ]

    rest_list = list(set(tf_list) - set(compiled_list))
    #tf_list = ['egr1']
    map_fastq_list = [ os.path.basename(fastq_file) for fastq_file in my.grep_list('.*-helas3.*(%s)'%'|'.join(rest_list), file_list)]

    map_fastq_list.sort()
    my.f_print_list(map_fastq_list)
    print len(map_fastq_list)
    fastq_gz_list = map_fastq_list

elif 'simulate_mask' in cell_name:
    file_list = my.f_shell_cmd( "ssh shi@clustdell.cmmt.ubc.ca find /home/shi/encode/ -name 'encode-*mask*simulate*.fastq.gz'", quiet = True).split('\n')
    
    fastq_gz_list = [ os.path.basename(fastq_file) for fastq_file in file_list ]
    my.f_print_list(fastq_gz_list)
elif 'simulate' in cell_name:
    file_list1 = my.f_shell_cmd( "ssh shi@clustdell.cmmt.ubc.ca find /home/shi/encode/ -name 'encode-*simulate*.fastq.gz'", quiet = True).split('\n')
    file_list = list( set(file_list1) - set(my.grep_list('.*mask',file_list1)) - set(my.grep_list('.*simulate-',file_list1)))
    fastq_gz_list = [ os.path.basename(fastq_file) for fastq_file in file_list ]
    my.f_print_list(fastq_gz_list)
    fastq_gz_list = (my.grep_list('.*gm12864', fastq_gz_list))
elif 'ctcf' in cell_name:
    feature = 'input'
    file_list = my.f_shell_cmd( "ssh $CLUST find /home/shi/encode/gm*/ -name '*.fastq.gz'", quiet = True).split('\n')
    #ctcf_list_raw = [ os.path.basename(fastq_file) for fastq_file in my.grep_list('^(?!.*gm12xxx|.*/ut-|.*open-).*%s'%feature, file_list)]
    #ctcf_list = my.grep_list('uw-(gm12864|gm12873).*', ctcf_list_raw)

    tf_list = [ 'yy1', 'tcf3', 'egr1', 'jund', 'mef2a']
    #tf_list = ['egr1']
    tf_list = ['pu1', 'myc']
    map_fastq_list = [ os.path.basename(fastq_file) for fastq_file in my.grep_list('.*gm12878.*(%s)'%'|'.join(tf_list), file_list)]

    map_fastq_list.sort()
    my.f_print_list(map_fastq_list)
    fastq_gz_list = map_fastq_list

elif 'all' in cell_name:
    cell = 'helas3'
    fastq_list = my.f_shell_cmd( "ssh wenqiang@orcinus.westgrid.ca find /home/wenqiang/encode/%s/ -name '*%s-*.fastq.gz'"%(cell, cell), quiet = True).split('\n')
    #ctcf_list_raw = [ os.path.basename(fastq_file) for fastq_file in my.grep_list('^(?!.*gm12xxx|.*/ut-|.*open-).*%s'%feature, file_list)]
    #ctcf_list = my.grep_list('uw-(gm12864|gm12873).*', ctcf_list_raw)
    time.sleep(2)
    bam_list = my.f_shell_cmd( "ssh wenqiang@orcinus.westgrid.ca find /home/wenqiang/encode/%s/ -size +100M -name '*.sorted.bam'"% cell, quiet = True).split('\n')

    mapped_list = [ bam_file.replace('sorted.bam', 'fastq.gz') for bam_file in bam_list]
    

    rest_list = [os.path.basename(file_name) for file_name in list(set(fastq_list) - set(mapped_list))]

    rest_list.sort()
    my.f_print_list(rest_list)
    print len(rest_list)
    fastq_gz_list = rest_list
elif 'embl' in cell_name:
    
    file_list = my.f_shell_cmd( "ssh $CLUST find /home/shi/encode/gm*/ -name '*.fastq.gz'", quiet = True).split('\n')

    tf_list = ['pu1', 'myc']
    map_fastq_list = [ os.path.basename(fastq_file) for fastq_file in my.grep_list('.*-(%s)'%'|'.join(tf_list), file_list)]

    map_fastq_list.sort()
    my.f_print_list(map_fastq_list)
    fastq_gz_list = map_fastq_list

else:
    print "Unknow input"


if test_flag == 'test':
    mode = ['test']
else:
    mode = ['unique']
    #mode = None


    
for fastq_gz_file in fastq_gz_list:
    
    if "clustdell" in target_server:
        cmd="""python2.7 p_run_cluster_sep.py "novo-%s-shi" p_novo_mapping.py %s %s"""%(fastq_gz_file,fastq_gz_file, my.f_send_list_para(mode))
    else:
        cmd="""python2.7 p_run_bugaboo.py %s "novo-%s" p_novo_mapping.py %s %s"""%(target_server, fastq_gz_file.split('.')[0], fastq_gz_file, mode)
    print cmd
    
    my.f_shell_cmd(cmd)
    
    time.sleep(2)

print '%s jobs submitted' % len(fastq_gz_list)
# p p_novo_mapping_para.py clustdell gm12878 run
# p p_novo_mapping_para.py clustdell gm12878 run
# p p_novo_mapping_para.py orcinus all run
# p p_novo_mapping_para.py clustdell simulate run
