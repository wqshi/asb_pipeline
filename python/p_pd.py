#!
import p_mymodule as my
reload(my)
import pandas as pd
import p_script_test as scripttest
from os.path import expanduser
home = expanduser("~")
import logging
from lockfile import LockFile
from lockfile import LockTimeout
import numpy as np
import sys
import os
import pybedtools
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


class data_table:
    file_path=None
    table=None

    def __init__(self, file_path, data=None):
        #import ipdb; ipdb.set_trace()
        #the start position should be 1-based.
        self.file_path=file_path
        if data is None:
            self.data=pd.io.parsers.read_csv(file_path, sep="\t")
        else:
            self.data=data
        self.data.set_index(keys=["chr","start"], inplace=True, drop=False)

    def renew_data(self, expected_cols=["chr","start"]):
        self.data=pd.io.parsers.read_csv(self.file_path, sep="\t")
        self.data.set_index(keys=expected_cols, inplace=True, drop=False)
        
    def sort(self, col_names):
        self.data.sort(col_names, ascending=True,inplace=True)
    
    def save_data(self):
        self.data.drop_duplicates().to_csv(self.file_path, index=None, sep="\t", na_rep=".")
    
    def extract_loc(self, data_dir, header=False):
        #print self.data.head()
        loc_file = data_dir + my.f_generate_tmp_file_name("loc")
        self.data[["chr","start"]].drop_duplicates().to_csv(loc_file, header=header, index=False, sep="\t")
        return loc_file

    def show_size(self, return_flag=False):
        if return_flag == False:
            print "=====>Data size: ", self.data.shape
        else:
            return "Data size: " + str(self.data.shape)
        
    
    def extract_allele(self, data_dir, header=False, extra_cols=[]):
        #print self.data.head()
        #import ipdb; ipdb.set_trace()
        allele_file = data_dir + my.f_generate_tmp_file_name("allele")
        allele_data= self.data[["chr","start"]]
        allele_data[["end"]]=self.data[["start"]] + 1
        allele_data[["ref",'alt']]=self.data[["ref","alt"]]
        allele_data[extra_cols]=self.data[extra_cols]
        allele_data.drop_duplicates(cols=["chr","start","alt"],inplace=True)
        allele_data.drop_duplicates().to_csv(allele_file, header=header, index=False, sep="\t")
        return allele_file

    def extract_peak_max(self, data_dir, filter_col=None, filter_val=None, add_info = None):

        #import ipdb; ipdb.set_trace()
        bed_file = data_dir + '/' + my.f_generate_tmp_file_name("loc.bed")
        if filter_col is not None:
            bed_data=self.data[self.data[filter_col] == filter_val][['chr','start','peak_dis']]
        else:
            bed_data=self.data[["chr","start",'peak_dis']]

        bed_data["name"]=bed_data["chr"]+"-"+bed_data["start"].map(str)
        bed_data["end"]=bed_data["peak_dis"].astype(int)
        bed_data["start"]=bed_data["peak_dis"].astype(int)-1
        
        if add_info is not None:
            bed_data["name"]=bed_data["name"]+"-"+self.data[add_info]
        
        bed_data.set_index(keys=["chr","start"], inplace=True, drop=False)
        bed_data[['chr', 'start', 'end', 'name']].drop_duplicates().to_csv(bed_file, header=False, index=False, sep="\t")
        return bed_file        
    
    
    def extract_bed(self, data_dir, filter_col=None, filter_val=None, add_info = None):
        #print self.data.head()
        #import ipdb; ipdb.set_trace()
        bed_file = data_dir + '/' + my.f_generate_tmp_file_name("loc.bed")
        if filter_col is not None:
            bed_data=self.data[self.data[filter_col] == filter_val][['chr','start']]
        else:
            bed_data=self.data[["chr","start"]]
        bed_data["end"]=bed_data["start"].astype(int)
        bed_data["start"]=bed_data["start"].astype(int)-1
        bed_data["name"]=bed_data["chr"]+"-"+bed_data["start"].map(str)
        
        if add_info is not None:
            bed_data["name"]=bed_data["name"]+"-"+self.data[add_info]
        
        bed_data.set_index(keys=["chr","start"], inplace=True, drop=False)
        bed_data.drop_duplicates().to_csv(bed_file, header=False, index=False, sep="\t")
        return bed_file
    
    def merge_feature(self, feature_data,expected_cols=["chr", "start"], debug = False):

        if debug == True:
            import ipdb; ipdb.set_trace()
        logging.debug(self.show_size(return_flag =True))
        assert set(expected_cols) < set(feature_data.columns), "Unexpected col names in feature data"
        feature_data.set_index(keys=expected_cols, inplace=True, drop=False)
        self.data.set_index(keys=expected_cols, inplace=True, drop=False)
        #print self.data.index[1:10]
        
        if (set(self.data.index) >= set(feature_data.index)) == False:
            print self.data.index
            print feature_data.index
            print set(feature_data.ix[:,"chr"]+feature_data.ix[:,"start"].map(str)) - set(self.data.ix[:,"chr"]+self.data.ix[:,"start"].map(str))
        assert set(self.data.index) >= set(feature_data.index), "Unexpected index"

        logging.debug("before lock")
        
        lock = LockFile( home + "/.lock.%s" % os.path.basename(self.file_path))
  
        
        while not lock.i_am_locking():
            try:
                lock.acquire(timeout=600)    # wait up to 60 seconds
            except LockTimeout:
                logging.debug("break the lock")
                lock.break_lock()
                lock.acquire()

        logging.debug('file is locked')
        
        self.renew_data(expected_cols)
  
        new_cols = list( set(feature_data.columns) - set(self.data.columns) )
        update_cols = list(set(feature_data.columns).intersection(self.data.columns) - set(expected_cols) - set(['ref']) )

        assert 'cell' not in update_cols, 'cell in update_cols'
        assert 'het_type' not in update_cols, 'het_type in update cols'
        
        if  update_cols !=[]:
            duplicated_rows = my.f_duplicated_index(feature_data, expected_cols)

            if any(duplicated_rows) > 0:
                logging.info("=======duplicated data===========")
                print feature_data.ix[duplicated_rows, :]

            
            update_data = feature_data.drop_duplicates(expected_cols).ix[:,update_cols]
            #print update_data.ix[0:10,:]
            self.data.update(update_data)
            print 'Updated data in when merge'
            print self.data.ix[0:10,update_cols]
            #print self.data.index[1:10]
        if new_cols !=[]:
            self.data=pd.merge(self.data, feature_data.drop_duplicates(expected_cols).ix[:,new_cols + expected_cols], on=expected_cols, how="left")

        
        self.data.set_index(keys=["chr","start"], inplace=True, drop=False)
        self.save_data()
        self.show_size(return_flag = True)
        lock.release()

    def drop_feautre(self, col_names):
        #import ipdb; ipdb.set_trace()
        self.data = self.data.drop(col_names, axis = 1)
        self.save_data()
        
    def read_feature(self, dp_file, feature_name):
        
        feature_data=pd.io.parsers.read_csv(dp_file, index_col = False, sep="\t")
        feature_data.set_index(keys=["chr","start"], inplace=True, drop=False)
        col_names = pd.DataFrame(feature_data.columns)
        col_names[col_names=="chip_ref_dp"] = "ref_%s_dp" % feature_name
        col_names[col_names=="chip_alt_dp"] = "alt_%s_dp" % feature_name
        feature_data.columns = col_names[0].tolist()
        return feature_data

    def get_cord_name(self):
        return self.data.chr + self.data.start.astype(int).astype('str')
    
    #def combine_another_TF(self, feature_data, e ):
    
    def read_feature_replace_name(self, dp_file, target_names=[], replace_names=[]):
        feature_data=pd.io.parsers.read_csv(dp_file, index_col = False, sep="\t")
        #import ipdb; ipdb.set_trace()
   
        if "chip_ref" in feature_data.columns:
   
            good_rows = np.array(feature_data['chip_alt'].str.len() == 1) & np.array(feature_data['chip_ref'].str.len() == 1)

            if any(good_rows == False):
                logging.warning('Structure variants detected:')
                print feature_data.ix[good_rows == False,:]
                
            feature_data = feature_data.ix[good_rows,:]
        
        #import ipdb; ipdb.set_trace()
        col_names = pd.DataFrame(feature_data.columns)
        assert (set(target_names) <= set(col_names[0].tolist())), "Unkown names in target names"
 
        for target_name, replace_name in zip(target_names, replace_names):
            col_names[col_names==target_name] = replace_name

        
        
        if feature_data.empty:
            feature_data=pd.DataFrame(columns=col_names[0].tolist())
        else:
            feature_data.columns = col_names[0].tolist()

        feature_data.set_index(keys=["chr","start"], inplace=True, drop=False)
        return feature_data   


    def select_cols_data(self, cols):
        return self.data.ix[:, cols].drop_duplicates()
    
    def merge_another_database(self, other_database, merged_file, expected_cols,update_cols, new_cols, suffix, debug):
        if debug == True:
            import ipdb; ipdb.set_trace()
        merged_db =  data_table(merged_file, self.select_cols_data( expected_cols + update_cols + new_cols))
        if  update_cols !=[]:
            updated_data = other_database.data.ix[:,update_cols + expected_cols].drop_duplicates()
            merged_db.data = merged_db.data.combine_first(updated_data)
    
        if new_cols !=[]:
            new_data = other_database.data.ix[:,new_cols + expected_cols].drop_duplicates()
            merged_db.data = pd.merge(merged_db.data, new_data, on = expected_cols, how="outer", suffixes = suffix)
            merged_db.data.set_index(keys=["chr","start"], inplace=True, drop=False)
        return merged_db
    
    def filter_data(self, filter_col, filter_value):
        if filter_col in self.data.columns:
            data = self.data[self.data[filter_col] == filter_value]
        else:
            data = self.data

        return data

    def add_another_database(self, other_database, expected_cols, update_cols, new_cols, suffix, debug):
        if debug == True:
            import ipdb; ipdb.set_trace()
        
        
        if  update_cols !=[]:
            updated_data = other_database.data.ix[:,update_cols + expected_cols].drop_duplicates()
            self.data = self.data.combine_first(updated_data)

        #import ipdb; ipdb.set_trace()
        new_cols_rep = new_cols[:]
        for col_name in new_cols_rep:
            if col_name + suffix in self.data.columns:
                new_cols.remove(col_name)
    
        

        if new_cols !=[]:
            tmp_data = other_database.data.ix[:, expected_cols]
            renamed_cols = [ col_name + suffix for col_name in new_cols ]
            tmp_data[renamed_cols] = other_database.data[new_cols]
            new_data = tmp_data.drop_duplicates(cols = expected_cols)
            self.data.drop_duplicates(cols = expected_cols, inplace = True)
            self.data = pd.merge(self.data, new_data, on = expected_cols, how="outer")
            self.data.set_index(keys=["chr","start"], inplace=True, drop=False)



    def write_records_fastq(self, output_dir, fastq_recoreds, prefix = ''):
        if prefix == '':
            sequence_file= os.path.join( output_dir, my.f_generate_tmp_file_name('seq') )
        else:
            sequence_file= os.path.join( output_dir, prefix )
            
        output_handle = open(sequence_file, "w")
        SeqIO.write(fastq_recoreds, output_handle, "fasta")
        output_handle.close()
        return sequence_file 
    
    def get_ref_and_alt_fastq_files_from_database(self, output_dir, hg19_file, flanking_length =30, file_prefix = None, debug = False):
        if debug == True:
            import ipdb; ipdb.set_trace()

        tf_database = self
        
        if 'pwm_ref_strand' not in tf_database.data.columns:
            tf_database.data['pwm_ref_strand'] = '+'
        
        
        allele_file=tf_database.extract_allele(output_dir, header=True, extra_cols = ['pwm_ref_strand']  )
        pssm_len = flanking_length
        allele_data=pd.io.parsers.read_csv(allele_file, sep="\t", index_col=None)


        if allele_data.shape[0] == 0:
            print "empty input"
            #return

        bed_data=allele_data[["chr","start"]]
        bed_data[["start"]]= allele_data[["start"]] - pssm_len
        bed_data[["end"]]= allele_data[["start"]]   +pssm_len -1
        bed_data["name"] = '.'
        bed_data['score'] = '.'
        bed_data[['strand']] = allele_data[['pwm_ref_strand']]

        bed_str=bed_data.to_string(header=False, index=False)
        bed_file=pybedtools.BedTool(bed_str, from_string=True)

        fasta =  pybedtools.example_filename(hg19_file)
        a = bed_file.sequence(fi=fasta)
        from Bio import SeqIO
        fasta_file= os.path.join(output_dir, my.f_generate_tmp_file_name('fasta'))

        import shutil
        shutil.copyfile(a.seqfn, fasta_file)

        allele_data['fastq'] = ''


        mutation_seq=[]
        i=-1
        mutation_pos_col=pssm_len -1

        alt_fastq_records = []
        ref_fastq_records = []
        for record in SeqIO.parse(open(fasta_file), "fasta"):
            i=i+1;
            mutation_pos=mutation_pos_col
            line=allele_data.ix[i,]

            ref_allele=line[3]
            alt_allele=line[4]
            pwm_ref_strand = line[5]

            assert ref_allele.upper() == record[mutation_pos].upper(), "Ref Allele doesn't Match'"
            record.seq.alphabet = IUPAC.unambiguous_dna
            #The following is necessary to copy
            mutation_record = SeqRecord(record, id = record.id, name = record.name, description = record.description)
            mutation_seq = record.seq.lower().tomutable()
            mutation_seq[mutation_pos]=alt_allele.upper()
            mutation_record.seq = mutation_seq.toseq()
            mutation_record.seq.alphabet = IUPAC.unambiguous_dna

            ref_record = record
            ref_seq = record.seq.lower().tomutable()
            ref_seq[mutation_pos] = ref_allele.upper()
            ref_record.seq = ref_seq.toseq()

            #if pwm_ref_strand == '-':
                #ref_record.seq = record.seq.complement()
                #mutation_record.seq = mutation_record.seq.complement()
            #print mutation_seq
            #import ipdb; ipdb.set_trace()
            
                       
            #mutation_record.id = ref_record.id
            #mutation_record.name = ref_record.name
            
            alt_fastq_records.append(mutation_record)
            ref_fastq_records.append(ref_record)
            allele_data.ix[i, 'fastq' ] = str(ref_record.seq)

        
        alt_sequence_file= self.write_records_fastq(output_dir, alt_fastq_records, prefix = file_prefix + '.alt.fastq')
        ref_sequence_file= self.write_records_fastq(output_dir, ref_fastq_records, prefix = file_prefix + '.ref.fastq')

        os.remove(allele_file)
        return [ref_sequence_file, alt_sequence_file, allele_data]


    def get_ref_and_alt_peak_fastq_files_from_database(self, output_dir, hg19_file, file_prefix = None, target_lab = '', debug = False):
        if debug == True:
            import ipdb; ipdb.set_trace()

        tf_database = self
        
        peak_start = 'peak_%s_bed_start' % target_lab
        tf_database.data[peak_start] = tf_database.data['peak_%s_dis' % target_lab] - 50
        
        print my.grep_list('lab', tf_database.data.columns.tolist())
        
        lab_data=tf_database.data[ tf_database.data['lab_%s' % target_lab] != '.' ]
        
        peak_data = lab_data.ix[:, ['chr', peak_start]]
        peak_data[peak_start] = peak_data.ix[:,peak_start].astype('float').astype('int')
        peak_data['end'] = peak_data[peak_start] + 100
        
        #old
        full_data = tf_database.data
        tf_database.data = lab_data
        allele_file=tf_database.extract_allele(output_dir, header=True)
        allele_data=pd.io.parsers.read_csv(allele_file, sep="\t", index_col=None)

        tf_database.data = full_data
        if allele_data.shape[0] == 0:
            print "empty input"

        bed_str=peak_data.to_string(header=False, index=False)
        bed_file=pybedtools.BedTool(bed_str, from_string=True)

        fasta =  pybedtools.example_filename(hg19_file)
        a = bed_file.sequence(fi=fasta)
        from Bio import SeqIO
        fasta_file= os.path.join(output_dir, my.f_generate_tmp_file_name('fasta'))

        import shutil
        shutil.copyfile(a.seqfn, fasta_file)

        allele_data['fastq'] = ''
        mutation_seq=[]
        i=-1

        #mutation_pos_col = [int(record[5]) - int(record[1])  for record in peak_regions]
        #print allele_data.ix[:,'start']
        #print peak_data.ix[:, peak_start]
        peak_data.index = allele_data.index
        mutation_pos_col = allele_data.ix[:,'start'] - peak_data.ix[:,peak_start] - 1 
        #print mutation_pos_col
        
        alt_fastq_records = []
        ref_fastq_records = []
        for record in SeqIO.parse(open(fasta_file), "fasta"):
            i=i+1;
            mutation_pos=mutation_pos_col[i]
            line=allele_data.ix[i,]

            ref_allele=line[3]
            alt_allele=line[4]
            pwm_ref_strand = line[5]
            #print  'index %s' % i
            #if i == 112:
            #    import ipdb; ipdb.set_trace()
            assert ref_allele.upper() == record[mutation_pos].upper(), "Ref Allele doesn't Match'"
            record.seq.alphabet = IUPAC.unambiguous_dna
            mutation_record = SeqRecord(record, id = record.id, name = record.name, description = record.description)
            mutation_seq = record.seq.lower().tomutable()
            mutation_seq[mutation_pos]=alt_allele.upper()
            mutation_record.seq = mutation_seq.toseq()
            mutation_record.seq.alphabet = IUPAC.unambiguous_dna

            ref_record = record
            ref_seq = record.seq.lower().tomutable()
            ref_seq[mutation_pos] = ref_allele.upper()
            ref_record.seq = ref_seq.toseq()
                        
            alt_fastq_records.append(mutation_record)
            ref_fastq_records.append(ref_record)
            allele_data.ix[i, 'fastq' ] = str(ref_record.seq)

        
        alt_sequence_file= self.write_records_fastq(output_dir, alt_fastq_records, prefix = file_prefix + '.alt.fastq')
        ref_sequence_file= self.write_records_fastq(output_dir, ref_fastq_records, prefix = file_prefix + '.ref.fastq')

        os.remove(allele_file)
        return [ref_sequence_file, alt_sequence_file, allele_data]

    def extract_snp(self, data_dir):
        #import ipdb; ipdb.set_trace()
        #Only for the het sites.
        #mainly for read simulation.
        #zero-based in the format: snp_name chr1 933790 + G A
        
        allele_file = data_dir + my.f_generate_tmp_file_name("snp")
        allele_data= self.data.ix[self.data.het_type=='het', ["chr","start"]]
        allele_data[["start"]]=allele_data[["start"]]
        allele_data["strand"]='+'
        allele_data[["ref",'alt']]=self.data.ix[self.data.het_type=='het', ["ref","alt"]]
        allele_data.drop_duplicates(cols=["chr","start","alt"],inplace=True)
        allele_data.index=range(0, allele_data.shape[0])
        allele_data.drop_duplicates().to_csv(allele_file, header=True, index=True, sep=" ")
        return allele_file
    
    def overlap_with_feature_bed(self, feature_file, loc_file, value_col, value_name, feature_extend_distanace = 0, sep='\t' ,debug= False):
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
        feature_bed_extended= pd.io.parsers.read_csv(feature_file, header = None, sep=sep)
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

    

import unittest
class TestDatabaseTable(unittest.TestCase):
    def setUp(self):
        test_dir="/homed/home/shi/test/t_het_sites_in_narrow_peak/"
        test_envi=scripttest.scripttest(test_dir)
        data_file=test_dir + "test.database"

        loc_file=test_dir + "sydh-gmtest-ctcf.het.loc"
        df1=pd.io.parsers.read_csv(loc_file, header=None, index_col = False, sep="\t")
        df1.columns=["chr","start","ref","alt","genotype","wgs_ref_dp","wgs_alt_dp"]
        df1["het"]="het"
        dt_obj=data_table(data_file, df1)

        #feature
        self.feature_name = "h4k20me1"
        #feature_file=test_dir + "/cross_cell/gmtest2-h4k20me1-Rep1.gmtest.ctcf.alt.bam.vcf.dp"
        self.feature_file=test_dir + "/h4k20me1.dp2"
        self.test_dir = test_dir
        
        #Preparing it
        self.database = dt_obj

    def test_init(self):
        self.assertTrue( self.database.data.shape[0] != 0 )

    def test_read_and_merge(self):
        #import ipdb; ipdb.set_trace()
        feature_data=self.database.read_feature(self.feature_file, self.feature_name)

        self.database.merge_feature(feature_data)
        self.assertTrue( 'alt_%s_dp' % self.feature_name in self.database.data.columns )
        self.assertTrue( 'ref_%s_dp' % self.feature_name in self.database.data.columns )

    def test_read_replace_name(self):
        #print ''
        feature_data2 = self.database.read_feature_replace_name(self.feature_file, ["chip_ref_dp","chip_alt_dp"], ["ref_%s_dp"%self.feature_name, "alt_%s_dp"%self.feature_name] )
        self.assertTrue( 'alt_%s_dp' % self.feature_name in feature_data2.columns )
        self.assertTrue( 'ref_%s_dp' % self.feature_name in feature_data2.columns )

    def test_extract(self):

        #import ipdb; ipdb.set_trace()

        feature_data=self.database.read_feature(self.feature_file, self.feature_name)
        self.database.merge_feature(feature_data)
        
        bed_file=self.database.extract_bed(self.test_dir, filter_col = "ref_h4k20me1_dp", filter_val=0)
        df1=pd.io.parsers.read_csv(bed_file, header=None, index_col = False, sep="\t")
        
        self.assertTrue( sum(self.database.data['ref_h4k20me1_dp'] == 0) == df1.shape[0] )
        #allele_file=self.database.extract_allele(self.test_dir, header=True, extra_cols = ["wgs_ref_dp"])


        
    def test_merge_tf_database(self):
        #import ipdb; ipdb.set_trace()
        
        db_file=my.f_create_file_name(self.test_dir, 'gmtest', 'ctcf',"database")
        tf_db =data_table(db_file)
        
        update_cols=my.grep_list('^(ref|alt)$', tf_db.data.columns)
        new_cols    = my.grep_list('genotype', tf_db.data.columns)
        expected_cols = ['chr', 'start']
        
        small_db = tf_db
        small_db.data = tf_db.data.ix[0:5, expected_cols + update_cols + new_cols]

        db_file2 =my.f_create_file_name(self.test_dir, 'gmtest', 'znf143',"database")
        tf_db2   = data_table(db_file2)
        tf_db2.data = tf_db2.data.ix[0:5, expected_cols + update_cols + new_cols]

        suffix = ['_tf1','_tf2']
        merged_file = self.test_dir + '/gmtest_merged_db.database'
        
        merged_db = tf_db.merge_another_database(tf_db2, merged_file, expected_cols ,update_cols, new_cols, suffix, debug = False)
        merged_db.save_data()
        self.assertTrue( set( tf_db2.get_cord_name() ) <= set(merged_db.get_cord_name()) )

        new_add_rows = tf_db2.data.index[0]
        self.assertEqual( merged_db.data.ix[new_add_rows, 'genotype_tf2'], tf_db2.data.ix[new_add_rows, 'genotype'] )

    def test_add_tf_database(self):
        #import ipdb; ipdb.set_trace()
        
        db_file=my.f_create_file_name(self.test_dir, 'gmtest', 'ctcf',"database")
        tf_db =data_table(db_file)
        
        update_cols=my.grep_list('^(ref|alt)$', tf_db.data.columns)
        new_cols    = my.grep_list('genotype', tf_db.data.columns)
        expected_cols = ['chr', 'start']
        
        small_db = tf_db
        small_db.data = tf_db.data.ix[0:5, expected_cols + update_cols + new_cols]

        db_file2 =my.f_create_file_name(self.test_dir, 'gmtest', 'znf143',"database")
        tf_db2   = data_table(db_file2)
        tf_db2.data = tf_db2.data.ix[0:5, expected_cols + update_cols + new_cols]

        suffix = ['_tf1','_tf2']
        merged_file = self.test_dir + '/gmtest_merged_db.database'
        merged_db = data_table(merged_file, data = small_db.data[expected_cols + update_cols])

        merged_db.add_another_database(tf_db, expected_cols, update_cols, new_cols, '_tf1', debug = False)
        merged_db.add_another_database(tf_db2, expected_cols, update_cols, new_cols, '_tf2', debug = False)
        
    
        
    
        merged_db.save_data()
        self.assertTrue( set( tf_db2.get_cord_name() ) <= set(merged_db.get_cord_name()) )

        new_add_rows = tf_db2.data.index[0]
        self.assertEqual( merged_db.data.ix[new_add_rows, 'genotype_tf2'], tf_db2.data.ix[new_add_rows, 'genotype'] )
       
    def test_filter_data(self):
        #import ipdb; ipdb.set_trace()
        genotype_data=self.database.filter_data('genotype','0/1')
        self.assertTrue(all(genotype_data['genotype'] == '0/1'))

        #This case, the data don't contain the column, will return all
        genotype_data2=self.database.filter_data('genotype2','0/1')
        self.assertEqual(genotype_data2.shape, self.database.data.shape)

    def test_get_fastq_sequences(self):
        hg19_file="/homed/home/shi/projects/wgs/hg19.fa"
        output_dir = "/homed/home/shi/test/t_het_sites_in_narrow_peak/"
        file_prefix = 'gmtest-tf'
        ref_file, alt_file, allele_data = self.database.get_ref_and_alt_fastq_files_from_database(output_dir ,hg19_file, flanking_length = 10, file_prefix = file_prefix, debug = False)
        
        for record in SeqIO.parse(open(alt_file), "fasta"):
            import ipdb; ipdb.set_trace()
            assert record.seq == 'ctgggatcgGtgcgggcgg', 'Sequence error' #Assert the first sequence. Check the UCSC. The >chr1:875760-875779 is 0 based.
            break
        #print alt_file

    def test_get_peak_sequences(self):
        hg19_file="/homed/home/shi/projects/wgs/hg19.fa"
        output_dir = "/homed/home/shi/test/t_het_sites_in_narrow_peak/full_data/"
        
        file_prefix = 'gmtest-tf'
        db_file=my.f_create_cell_pair_database_name(output_dir, 'gm12878', ['gm12878'], 'ctcf')
        tf_database=data_table(db_file)
        tf_database.show_size()
        #import ipdb; ipdb.set_trace()
        ref_file, alt_file, allele_data = tf_database.get_ref_and_alt_peak_fastq_files_from_database(output_dir, hg19_file, file_prefix = 'peak', target_lab = 'sydh', debug = True)
        #dt_obj.drop_feautre('alt_h4k20me1_dp')
if __name__ == "__main__":

    suite = unittest.TestLoader().loadTestsFromTestCase( TestDatabaseTable )
    unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( suite )
