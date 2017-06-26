#/bin/bash
#source ~/library2/s_function.sh
f_substr_exists()
{   
    whole_str=$1
    pattern_str=$2
    echo $whole_str | grep -q -i $pattern_str
    var=$?
    return $var
}
#f_novo_readdepth()
#{
#    vcf_file=$1
#    printf "#chr\tstart\tref\talt\tgenotype\tref_depth\talt_depth\n"
#    sed '/#/d' $vcf_file | sed '/INDEL/d' |awk -F"[\t;]" '{print $1"\t"$2"\t"$4"\t"$5"\t"$12"\t"$NF}' > tmp
#    
#    awk -F"[\t]" '{print $1"\t"$2"\t"$3"\t"$4}' tmp > tmp.basic
#    awk -F"[\t]" '{print $5}' tmp | awk -F"[=,]" '{print $2+$3"\t"$4+$5}' > tmp.depth 
#    awk -F"[\t;]" '{print $6}' tmp | awk -F"[:]" '{print $1}' > tmp.genotype
#
#    paste tmp.basic tmp.genotype tmp.depth
#    #| awk -F"[\t,=]" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'
#}

f_complete_genome_read_depth()
{
    vcf_file=$1
    tmpID=$(f_generate_tmp_ID)
    cd $(dirname $vcf_file)
    #pwd
    printf "#chr\tstart\tref\talt\tgenotype\tref_depth\talt_depth\n"
    sed '/#/d' $vcf_file | sed '/INDEL/d' |awk -F"[\t;]" '{print $1"\t"$2"\t"$4"\t"$5"\t"$NF}' > tmp.$tmpID

    awk -F"[\t]" '{print $1"\t"$2"\t"$3"\t"$4}' tmp.$tmpID > tmp.basic.$tmpID
    #For the foramt, check complete genome reference
    awk -F"[\t]" '{print $5}' tmp.$tmpID | awk -F"[=:]" '{print $NF"\t"$(NF-2)-$NF}' > tmp.depth.$tmpID
    awk -F"[\t;]" '{print $5}' tmp.$tmpID | awk -F"[:]" '{print $1}' > tmp.genotype.$tmpID

    paste tmp.basic.$tmpID tmp.genotype.$tmpID tmp.depth.$tmpID
    
    rm tmp.basic.$tmpID tmp.genotype.$tmpID tmp.depth.$tmpID tmp.$tmpID
    #| awk -F"[\t,=]" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}'                                             
}

f_novo_readdepth()
{
    vcf_file=$1
    printf "#chr\tstart\tref\talt\tgenotype\tref_depth\talt_depth\n"
    sed '/#/d' $vcf_file |  sed '/INDEL/d' | 
    awk 'BEGIN{OFS = "\t" ;FS="[\t;=,:]"}
    /DP4/ {
    for (i=1; i<NF; i++) {
       if ($i == "DP4") {
           DP1=$(i+1)
           DP2=$(i+2)
           DP3=$(i+3)
           DP4=$(i+4)
       } else if ($i == "DP")
       {
          DP=$(i+1)
       }
       else if ($i == "GQ"){
         genotype=$(i+1)
         
         if (DP != DP1 + DP2 + DP3 + DP4)
         {
            print "Error"
         }
         print $1 OFS $2 OFS $4 OFS $5 OFS $(i+1) OFS DP1 + DP2 OFS DP3 + DP4; 

       }  
 }
}'
}

f_grep_legal_snp()
{
    grep -P "^#.*|chr[0-9XYxy]{1,2}\t[0-9]+\t[ATGC]\t[ATGC]\t" | grep -v "[.]" | sed "s/|/\//"
}




f_create_index()
{
#input: vcf.bz2 file
#Process: change header and create index and change the file name also
#output: new vcf gz file
    vcf_header=$1
    vcf_gz_file=$2

    vcf_file=$(echo $vcf_gz_file | sed -n 's/vcf.*[-]\(NA[0-9]*\)[-].*/\1/p'| tr "NA" "gm").vcf

    cat $vcf_header > $vcf_file
    bzip2 -c -d $vcf_gz_file | sed '/^#/d' >> $vcf_file

    bgzip -f $vcf_file
    tabix -f -p vcf $vcf_file.gz
    bgzip -c -d $vcf_file.gz >> $vcf_file
}


f_create_index2()
{
#input: vcf.bz2 file
#Process: change header and create index and change the file name also
#output: new vcf gz file
    vcf_header=$1
    vcf_gz_file=$2

    vcf_file=$(echo $vcf_gz_file | sed -n 's/.*\(NA[0-9]*\).*/\1/p'| tr "NA" "gm").vcf

    cat $vcf_header > $vcf_file
    f_unzip_file $vcf_gz_file | sed '/^#/d' >> $vcf_file

    bgzip -f $vcf_file
    tabix -f -p vcf $vcf_file.gz
    bgzip -c -d $vcf_file.gz > $vcf_file
}



f_cat_chipseq_fastq()
{
    #given the cell and tf list, cat all the fastq files together for each tf.
    chipseq_snp_dir=$1
    declare -a cell_list=("${!2}")
    declare -a tf_list=("${!3}")
    job_name=$4
    #Cat all the replicate fastaq data file into together
    echo "Merger Fasta file!!!"
    for cell in ${cell_list[@]}
    do
	for tf in ${tf_list[@]}
	do
	    comb_file=$cell-$tf.fastq
	    > $chipseq_snp_dir/$comb_file
	    for fastq_file in $(ls $chipseq_snp_dir/*.fastq | grep -i "$cell.\?$tf.*rep")
	    do
		if ( f_substr_exists $job_name "test" )
		then
		    echo $comb_file
		    #head -10000 $fastq_file >> $chipseq_snp_dir/$comb_file.test
		else
		    cat $fastq_file >> $chipseq_snp_dir/$comb_file
		fi
	    done
	done
    done
}


f_cat_chipseq_fastq_single()
{
    #given the cell and tf list, cat all the fastq files together for each tf.
    chipseq_snp_dir=$1
    cell=$2
    tf=$3
    job_name=$4
    lab=$5
    #Cat all the replicate fastaq data file into together
    echo "Merger Fasta file!!!"
    comb_file=$lab.$cell-$tf.fastq
    for fastq_file in $(ls $chipseq_snp_dir/*.fastq | grep -i "$cell.\?$tf.*rep")
    do
	if ( f_substr_exists $job_name "test" )
	then
	    echo $comb_file
	    head -10000 $fastq_file >> $chipseq_snp_dir/$comb_file.test
	else
	    cat $fastq_file >> $chipseq_snp_dir/$comb_file
	fi
    done
    
    echo "done"
}

f_rename_fastq()
{
    #given the cell and tf list, cat all the fastq files together for each tf.
    chipseq_snp_dir=$1
    cell=$2
    tf=$3
    job_name=$4
    #Cat all the replicate fastaq data file into together
    echo "Merger Fasta file!!!"
    
    for fastq_file in $(ls $chipseq_snp_dir/wg*.fastq | grep -i "$cell.\?$tf.*rep")
    do
	#rep_num=$(echo $fastq_file | sed -n 's/wg.*\(Rep[0-9]\).*/\1/p'| tr "[A-Z]" "[a-z]"| xargs -I file basename file)
	comb_file=$(f_parse_encode_name $fastq_file).fastq
	if ( f_substr_exists $job_name "test" )
	then
	    echo $comb_file
	    head -10000 $fastq_file > $chipseq_snp_dir/$comb_file.test
	else
	    mv $fastq_file $chipseq_snp_dir/$comb_file
	fi
    done
    
    echo "done"
}


f_narrow_peak_100bp()
{

    bed_file=$1
    output_file=$2
    awk -F"\t" -v OFS="\t" '{
if ( $10 ~ "[0-9]+" )
{
    print $1,$2+$10-50, $2+$10+50, "peak"
}
else
{
    print $1,($2+$3)/2-50, ($2+$3)/2+50, "center"
}


}' $bed_file | sort -k1,1 -k2,2n > $output_file

    #cut -f4 $output_file | sort | uniq > tmp
    peak_type=($( cut -f4 $output_file | sort | uniq | grep "[a-z]"))

    
    if [ ! ${#peak_type[*]} -eq 1 ] || [ ! ${peak_types[0]} == "peak" ]
    then
	echo "Bed Peak Max Problem"
    fi
}



f_narrow_peak_distance()
{

    bed_file=$1
    output_file=$2
    distance=$3
    awk -F"\t" -v OFS="\t" -v distance=$distance '{
if ( $10 ~ "[0-9]+" )
{
    print $1,$2+$10-distance, $2+$10+distance, $4, $5, $6, $7, $8, $9, $10
}
else
{
    print $1,($2+$3)/2-distance, ($2+$3)/2+distance, $4,$5, $6, $7, $8, $9, $10
}


}' $bed_file | sort -k1,1 -k2,2n > $output_file

    #cut -f4 $output_file | sort | uniq > tmp
    peak_type=($( cut -f4 $output_file | sort | uniq | grep "[a-z]"))

    
    if [ ! ${#peak_type[*]} -eq 1 ] || [ ! ${peak_types[0]} == "peak" ]
    then
	echo "Bed Peak Max Problem"
    fi
}

f_add_bed_name()
{
    awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2":"$3":"$10,$5,$6,$7,$8,$9,$10}'
}


f_parse_encode_name()
{
    encode_file_name=$(basename $1)
    if ( f_substr_exists $encode_file_name "Rep" )
    then 
	echo $encode_file_name | sed -n 's/wgEncode\([A-Z][a-z]*\)Tfbs\([A-Z][a-z0-9]*\)\([A-Z][a-z0-9]*\).*\(Rep[0-9]\).*/\1.\2-\3.\4/p' | tr '[A-Z]' '[a-z]'
    else
	echo $encode_file_name | sed -n 's/wgEncode\([A-Z][a-z]*\)Tfbs\([A-Z][a-z0-9]*\)\([A-Z][a-z0-9]*\).*/\1.\2-\3.Rep/p' | tr '[A-Z]' '[a-z]'
    fi
}

f_parse_uni_peak_name()
{
    encode_file_name=$(basename $1)
    echo $encode_file_name | sed -n 's/wgEncode.*Tfbs\([A-Z][a-z]*\)\([A-Z][a-z0-9]*\)\([A-Z][a-z0-9]*\).*/\1.\2-\3/p' | tr '[A-Z]' '[a-z]'

}



f_vcf_call()
{
    #Finishing the steps from the bam to vcf
    input_bam=$1
    vcf_header=$2
    output_vcf=$3

    samtools sort $input_bam tmp.sorted
    samtools index tmp.sorted.bam

    samtools mpileup -Bgf /raid2/local/exome/data/hg19//hg19-orderMito.fasta  tmp.sorted.bam | bcftools  view -gvc - >  tmp.input.bcf

    vcfutils.pl varFilter -Q20 -a2 tmp.input.bcf | awk '(match ($1,"##") || $6 > 10)' > tmp.output.vcf


    cat $vcf_header > $output_vcf
    sed "/#/d" tmp.output.vcf >> $output_vcf

    rm tmp.sorted.bam tmp.input.bcf tmp.output.vcf tmp.sorted.bam.bai

}


f_snv_fetch_peak_fastq()
{
    loc_file=$1
    bed_file=$2
    fa_file=$3
    output_file=$4
    
    bedtools intersect -wb -a $loc_file -b $bed_file | awk -F"\t" -v OFS="\t" '{print $5,$6,$7}' > tmp.overlap.bed

    bedtools getfasta -fi $fa_file -bed tmp.overlap.bed -fo tmp.fa

    grep "^>chr" tmp.fa | sed 's/>//g'| awk -F"[\t:-]" -v OFS="\t" '{print $1,$2+1,$3+1}' > tmp.loc
    grep -v "^>chr" tmp.fa > tmp.seq
    
    f_create_header "#chr start end name chr_bed start_bed end_bed seq" > $output_file
    
    dos2unix $loc_file 
    paste -d"\t" $loc_file tmp.loc tmp.seq | sed 's/^M//g' >> $output_file
    rm tmp.*

}

f_snv_fetch_peak_fastq2()
{
    loc_file=$1
    bed_file=$2
    fa_file=$3
    output_file=$4
    
    awk -F"\t" -v OFS="\t" '{print $1,$2-1,$3-1,$4}' $loc_file | bedtools intersect -wb -a stdin -b $bed_file | awk -F"\t" -v OFS="\t" '{print $5,$6,$7}' > tmp.overlap.bed

    bedtools getfasta -fi $fa_file -bed tmp.overlap.bed -fo tmp.fa

    grep "^>chr" tmp.fa | sed 's/>//g'| awk -F"[\t:-]" -v OFS="\t" '{print $1,$2+1,$3+1}' > tmp.loc
    grep -v "^>chr" tmp.fa > tmp.seq
    
    f_create_header "#chr start end name chr_bed start_bed end_bed seq" > $output_file
    
    dos2unix $loc_file 
    paste -d"\t" $loc_file tmp.loc tmp.seq | sed 's/^M//g' >> $output_file
    rm tmp.*

}


f_dnashape()
{    
    local fa_file=$1
    local output_file=$2
    ~/pakage/v2.5_noSim/prediction -i $fa_file -width 100

    shape_files=($fa_file.*)
    
    for shape_file in ${shape_files[*]}
    do
	type=${shape_file##$fa_file.}
	#grep "^>chr" $shape_file | sed 's/>//g'| awk -F"[\t:-]" -v OFS="\t" '{print $1,$2+20}' > tmp.loc
	grep -v "^>chr" $shape_file | awk -F"[,]" '{print $21}' > tmp.shape.$type
    done

    
    #dos2unix $loc_file 
    paste -d"\t" tmp.shape.* | sed 's/^M//g' > $output_file
    rm tmp.shape.*
    rm ${shape_files[*]}
}

f_dnashape_both_allele()
{
    allele_file=$1
    h19_file=$2
    out_file=$3
    #This function requires the format of allele_file to be:
    #chr start anthing ref_allele alt_allele

    grep -v "start" $allele_file | awk -F"\t" -v OFS="\t" '{print $1,$2-21,$2 +19,$1":"$2 }'  > tmp.overlap.bed
    grep -v "start" $allele_file | awk -F"\t" -v OFS="\t" '{print $1,$2,$2,$1"-"$2 }'  > tmp.overlap.loc


    bedtools getfasta -fi $h19_file -bed tmp.overlap.bed -fo tmp.ref.fa

    fa_file=tmp.ref.fa
    f_dnashape  $fa_file tmp.ref.shape
    
    #Generate the alt sequence, this require the 4th and 5th filed to be ref and alt allele
    python ~/library/p_vcf_mask_fasta.py $allele_file $fa_file tmp.alt.fa $PWD
    f_dnashape tmp.alt.fa tmp.alt.shape

    echo $out_file
    f_create_header "#chr start end name ref_helt ref_mgw ref_prot ref_roll alt_helt alt_mgw alt_prot alt_roll" > $out_file
    paste -d"\t" tmp.overlap.loc tmp.ref.shape tmp.alt.shape >> $out_file
    
    rm tmp.*
}


