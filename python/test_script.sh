project_dir=/homed/home/shi/projects/asb_pipeline/
loc_bin=$HOME/bin/novocraft/bin/
vcf_dir=$project_dir/data/raw_data/vcf2/
fastq_dir=$project_dir/data/raw_data/fastq2/
bam_dir=$fastq_dir
bed_dir=$fastq_dir
wgs_dir=$vcf_dir
cell=gm12878
tf=ctcf

#Mode option: test or unique. In the test mode, we will only map 1000 reads to test the mapping pipeline. In the unique mode, we will map all the reads. The unique mode will take 1~2 days to finish.
mode=test

#Prepare the data


ln -s $project_dir/data/ $project_dir/python/data
ln -s $project_dir/data/ $project_dir/R/data


##Step1. Prepare the data.

mkdir $vcf_dir
mkdir $fastq_dir

##Step1.1 VCF data
cd $vcf_dir

#The vcf directory should include: 
#1. hg19.fa. #The reference genome gh19.
#2. vcf file of the cell.(e.g. gm12878)

#Step1.1.1 hg19.fa: hg19 reference genome(https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz)
wget http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz --no-check-certificate -O hg19.fa.gz
gzip -d hg19.fa.gz


#Step 1.1.2 gm12878.vcf: the vcf file of the target cell. Used to retrive from (ftp://ftp2.completegenomics.com/vcf_files/Build37_2.0.0/). But the link is broken now. We provide a subset of vcf file in the github.
#wget https://github.com/wqshi/asb_pipeline/blob/master/data/raw_data/gm12878.chr22.vcf.gz?raw=true --no-check-certificate
gzip -d -c $project_dir/data/raw_data/gm12878.chr22.vcf.gz > gm12878.vcf




#Step 1.2 The ChIP-seq data
cd $fastq_dir


#The fastq directory should include:
#Step 1.2.1 Raw ChIP-seq reads for the transcription factor in the targeted cell (e.g. CTCF in GM12878).
#We can download the ChIP-seq data from the ENCODE project.

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdRawDataRep1.fastq.gz -O sydh-gm12878-ctcf-Rep1.fastq.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdRawDataRep2.fastq.gz -O sydh-gm12878-ctcf-Rep2.fastq.gz

#Step1.2.2 Get the called ChIP-seq peaks for CTCF.
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsSydhGm12878Ctcfsc15914c20UniPk.narrowPeak.gz -O sydh-gm12878-ctcf.narrowPeak.gz
gzip -d sydh-gm12878-ctcf.narrowPeak.gz

#Step 1.3 Copy number variants
cnv_dir=$project_dir/data/raw_data/cnv/
mkdir $cnv_dir
cd $cnv_dir
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibGenotype/wgEncodeHaibGenotypeGm12878RegionsRep1.bedLogR.gz
gzip -d wgEncodeHaibGenotypeGm12878RegionsRep1.bedLogR.gz 


##Step 2: map chIP-seq reads.
cd $project_dir/python/

#Create the personalized genome, about 1~2 hours.
sh s_format_vcf_gatk.sh $cell $vcf_dir > index.log 2>&1

#Map ChIP-seq reads, in test mode it will be less than 5 minustes. The normal process time woud be aroud 1~3 days 
python2.7 p_novo_mapping.py --fastq_dir $fastq_dir --fastq_file sydh-gm12878-ctcf-Rep1.fastq.gz --mode $mode --wgs_dir $wgs_dir --loc_bin $loc_bin
python2.7 p_novo_mapping.py --fastq_dir $fastq_dir --fastq_file sydh-gm12878-ctcf-Rep2.fastq.gz --mode $mode --wgs_dir $wgs_dir --loc_bin $loc_bin


#In the test case. We use existing bam file.
if [ '$mode'=='test' ];
then
   echo 'Test modep'
   $bam_dir=$bam_dir/test/
   
   cp $project_dir/data/raw_data/sydh-gm12878-ctcf-Rep1.chr22.sorted.bam $bam_dir/sydh-gm12878-ctcf-Rep1.sorted.bam
   cp $project_dir/data/raw_data/sydh-gm12878-ctcf-Rep2.chr22.sorted.bam $bam_dir/sydh-gm12878-ctcf-Rep1.sorted.bam
fi


#Extract the signal at heterozygous sites.
python2.7 p_post_process_after_asb_mapping_dp.py $cell $tf all $bed_dir $bam_dir $wgs_dir $cnv_dir

#Call ASB events.
cd ../R
Rscript3 ./s_call_asb_events.R --base_dir $bed_dir --cell $cell --tf $tf

#The final processed data files are in the $bed_dir,named as gm12878-ctcf.asb.txt.










