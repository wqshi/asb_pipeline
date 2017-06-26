#This is just complementary of the p_add_feature_on_loc_dp.py. Because the old one has several time consuming jobs.
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
import vcf

logging.basicConfig(level=logging.DEBUG)

logging.getLogger().setLevel(logging.DEBUG)


def f_get_loc_peak_max(loc_file, feature_file, output_file, debug=False):

    if debug == True:
        import ipdb; ipdb.set_trace()
        
    tmp_output = os.path.join( os.path.dirname(feature_file), my.f_generate_tmp_file_name('peak_max'))
    cmd1="awk -F'\t' -v OFS='\t' '{print $1,$2,$3,int(($2+$3)/2)}' %s | sort -k1,1 -k2,2n > %s" % (feature_file, tmp_output)
    my.f_shell_cmd(cmd1)

    cmd11="echo '#chr\tstart\tend\tname' > %s " % output_file
    my.f_shell_cmd(cmd11)
    cmd2 ="sort -k1,1 -k2,2n %s | bedtools intersect -a stdin -b %s -loj -sorted | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $8}' | sed 's/\./0/g' >> %s" % ( loc_file, tmp_output, output_file)
    
    my.f_shell_cmd(cmd2)
    os.remove(tmp_output)


def f_extract_genotype_from_vcf(vcf_file, output_file, debug = False):

    if debug == True:
        import ipdb; ipdb.set_trace()
    
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    chip_dp=[]

    sample_name = vcf_reader.samples[0]
    snp_list = []

    for record in vcf_reader:
        
        if record.is_snp:
            #assert len(record.ALT) ==1
            if record.FILTER == []:
                my_filter = 'PASS'
            else:
                my_filter = record.FILTER[0]
            call=record.genotype(sample_name)

            #Information from the complete genome vcf specification
            #Dp: totoal number of reads used to call genotype
            #AD: read depth of allele1 and allele2
            #CGA_RDP: number of reads support the reference allele.
            #So the formula would be: CGA_RDP, DP - CGA_RDP.
            ref_dp = -1
            alt_dp = -1
            if hasattr(call.data, 'DP'):
                ref_dp = int(call.data.CGA_RDP)
                alt_dp = call.data.DP - ref_dp
            item = [record.CHROM, record.POS, record.REF, record.ALT[0], my_filter, ref_dp, alt_dp, record.genotype(sample_name)['GT']]
            snp_list.append(item)

    
    snp_db = pd.DataFrame(snp_list)
    
    snp_db.columns = ['chr','end','ref','alt', 'filter',  'wgs_ref_dp', 'wgs_alt_dp', 'genotype']
    snp_db['start'] = snp_db['end'] -1
    snp_db_bed = snp_db[['chr','start','end','ref','alt', 'filter',  'wgs_ref_dp', 'wgs_alt_dp', 'genotype']]
    snp_db_bed.to_csv(output_file, na_rep ='.', sep='\t', index=None, header = False)


def f_extract_genotype_from_helas3_vcf(vcf_file, output_file, debug = False):

    if debug == True:
        import ipdb; ipdb.set_trace()
    
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    chip_dp=[]

    sample_name = vcf_reader.samples[0]
    snp_list = []

    for record in vcf_reader:
        
        if record.is_snp:
            #assert len(record.ALT) ==1
            if record.FILTER == [] or record.FILTER is None:
                my_filter = 'PASS'
            else:
                my_filter = record.FILTER[0]
            call=record.genotype(sample_name)
            seq_info = record.INFO
            #Information from the complete genome vcf specification
            #Dp: totoal number of reads used to call genotype
            #AD: read depth of allele1 and allele2
            #CGA_RDP: number of reads support the reference allele.
            #So the formula would be: CGA_RDP, DP - CGA_RDP.
            ref_dp = -1
            alt_dp = -1
            if "DP4" in record.INFO.keys():
                ref_dp=int(record.INFO["DP4"][0]) + int(record.INFO["DP4"][1])
                alt_dp=int(record.INFO["DP4"][2]) + int(record.INFO["DP4"][3])
            item = [record.CHROM, record.POS, record.REF, record.ALT[0], my_filter ,ref_dp, alt_dp, record.genotype(sample_name)['GT']]
            snp_list.append(item)

    
    snp_db = pd.DataFrame(snp_list)
    
    snp_db.columns = ['chr','end','ref','alt', 'filter',  'wgs_ref_dp', 'wgs_alt_dp', 'genotype']
    snp_db['start'] = snp_db['end'] -1
    snp_db_bed = snp_db[['chr','start','end','ref','alt', 'filter',  'wgs_ref_dp', 'wgs_alt_dp', 'genotype']]
    snp_db_bed.to_csv(output_file, na_rep ='.', sep='\t', index=None, header = False)


def f_extract_genotype_from_vcf_by_shell(vcf_file, output_file):
    grep_het_cmd = "f_complete_genome_read_depth %s |  f_grep_legal_snp | sed '/^#/d' > %s"% (vcf_file, output_file)
    print grep_het_cmd
    my.f_shell_fun_pipe(grep_het_cmd)



import unittest
class TestDatabaseTable(unittest.TestCase):
    def setUp(self):
        test_name = 'unittest'
        self.vcf_file = '/homed/home/shi/projects/wgs/vcf/gmtest.chr.vcf'
    def test_extract_genotype_from_vcf(self):
        output_file = '/homed/home/shi/projects/wgs/vcf/gmtest.snp.vcf'
        f_extract_genotype_from_vcf(self.vcf_file, output_file)

    def test_extract_genotype_from_vcf_by_shell(self):
        output_file = '/homed/home/shi/projects/wgs/vcf/gmtest.snp.shell.vcf'
        f_extract_genotype_from_vcf_by_shell(self.vcf_file, output_file)

    def test_0000_helas3_genotype(self):
        vcf_file = '/homed/home/shi/projects/wgs/complete_genome/helas3.test.vcf'
        output_file = '/homed/home/shi/projects/wgs/complete_genome/helas3.snp.dp'
        f_extract_genotype_from_helas3_vcf(vcf_file, output_file, debug = True)


#my.f_sync_scripts_to_run_server()
if __name__ == "__main__":
    server_name = socket.gethostname()
    if server_name == 'loire':
        suite = unittest.TestLoader().loadTestsFromTestCase( TestDatabaseTable )
        unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( suite )
    else:
        #vcf_file = '/home/shi/projects/wgs/complete_genome/helas3.test.vcf'
        vcf_file = '/home/shi/projects/wgs/complete_genome/helas3.vcf'
        output_file = '/home/shi/projects/wgs/complete_genome/helas3.snp.dp'
        f_extract_genotype_from_helas3_vcf(vcf_file, output_file, debug = False)
















