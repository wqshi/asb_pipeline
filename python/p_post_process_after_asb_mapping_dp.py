
import re
import subprocess
import os
import sys
import socket
import time
import p_mymodule as my
import p_locker as locker
reload(my)

if len(sys.argv) > 2:
    print sys.argv
    cell = sys.argv[1]
    extract_type= sys.argv[3]#step
    guest_cell = cell
    tf = sys.argv[2]

    bed_dir = sys.argv[4]
    bam_dir = sys.argv[5]
    wgs_dir = sys.argv[6]
    cnv_dir = sys.argv[7]
else:
    cell = 'gmtest'
    guest_cell = 'gmtest'
    tf = 'ctcf'
    extract_type = 'all'

    bed_dir = "/homed/home/shi/test/t_het_sites_in_narrow_peak/"
    wgs_dir="/homed/home/shi/projects/wgs"
    bam_dir="/homed/home/shi/test/t_het_sites_in_narrow_peak"
    cnv_dir='%s/cnv/'%wgs_dir
input_step=extract_type
    
diffpeak_method = 'pepr'
feature_list = []
mapQ = 30    


labs = ['sydh', 'uw', 'haib','uta','embl','encode']

cofactor_list = []




first_step=["het_loc"]
second_step=["extract_depth", "add_feature",'add_feature_light']
third_step = ["update_data"]

if input_step == "all":
    steps= ["het_loc"] + second_step + third_step
elif input_step == "first":
    steps= ['het_loc']
elif input_step == 'second':
    steps = second_step
else:
    steps= [input_step] # + ["update_data"]

#steps=["extract_depth"]   
#steps= second_step + third_step
#steps=second_step + third_step

#steps=["cross_cell"]
#steps= ["het_loc"] + second_step + third_step
#steps= second_step + ["cross_cell"] + third_step
#steps=["rmdup"]


print my.f_recieve_list_para(cell)
reload(my)

cur_time=time.strftime("%Y%m%d_%H%M%S") + my.f_id_generator(5)
server_name = socket.gethostname()

print cur_time
tf_list =  [tf]
guest_cells = [guest_cell]

if "het_loc" in steps:
    #loc_pattern = "*.(%s).narrowPeak"%"|".join(tf_list)

    het_loc_cmd="python2.7 p_het_sites_in_narrow_peak_dp.py %s %s %s %s %s %s %s" %(cell, my.f_send_list_para(tf_list), my.f_send_list_para(guest_cells), "locker.het_loc", my.f_send_list_para(labs), bed_dir, wgs_dir) 
    my.f_shell_cmd(het_loc_cmd)
else:
    print "===Skip Het Loc!==="




if "extract_depth" in steps:
    #Tried 7G, now it's 4G
    extract_cmd="python2.7 p_extract_depth_from_bam_dp.py %s %s %s %s %s %s "%\
      (cell, my.f_send_list_para(tf_list), 'extract_depth', mapQ, my.f_send_list_para(labs), bam_dir)
    print extract_cmd
    my.f_shell_cmd(extract_cmd)
else:
    print "===Skip Extract Depth!==="



if "add_feature_light" in steps:
    #small is fine, check 2G
    add_cmd="python2.7 p_add_feature_on_loc_dp_light.py  %s %s %s %s %s"%(cell, my.f_send_list_para(tf_list), my.f_send_list_para(labs), bed_dir, cnv_dir )
    my.f_shell_cmd(add_cmd)
else:
    print "=====Skip feature light=============="

#p p_post_process_after_asb_mapping_dp.py gm12878 gm12801 extract_depth


