#cp /raid2/local/exome/data/hg19/DEC2013/hg19-combined-ORDEREDPROPERLY.fa hg19.fa

cell=$1
cell_list=($1)
#cell_list=("gm12801" "gm12864" "gm12873" "gm12872")
#cell_list=("helas3")
#cell_list=("gm19238" "gm19239" "gm19240")
#cell_list=("gm12813" "gm12776" "gm11831" "gm11830" "gm11894")

#The data dir which contains the hg19.fa file, and vcf file of given cells.
#vcf_dir=./data/raw_data/vcf/
vcf_dir=$2
cd $vcf_dir
for cell in ${cell_list[@]}
do
    echo $cell
    vcftools --vcf ./$cell.vcf --remove-indels --recode --recode-INFO-all --out ./$cell.SNPs_only
    novoutil iupac -g $cell.SNPs_only.recode.vcf  hg19.fa > $cell.fa
    novoindex -t 2 $cell.nix $cell.fa
    #rm $cell.fa $cell.SNPs_only.recode.vcf
done

#check the genotype
bedtools getfasta -fi $cell.fa -bed test.bed -fo $cell.test.fa
bedtools getfasta -fi hg19.fa -bed test.bed -fo hg19.test.fa
head *.test.fa


#For the cases, alt is ".", this is fine. See the first line in the test bed
#For the error in the log file. I have checked the last homo in the chr3. And it's fine. See test_bed.


#\rm $cell.fa $cell.SNPs_only.recode.vcf  

#For GATK filtering
#vcftools --vcf input_file.vcf --remove-indels --recode --recode-INFO-all --out SNPs_only /--keep-only-indels 
#awk '( match ($1,"##") ||  $5 !~ /\[|\]|<CGA/ )' indel_only.recode.vcf > indel.vcf
#sed 's/^\([ 0-9XM]\)/chr\1/g' indel.vcf > indel2.vcf ]""))












