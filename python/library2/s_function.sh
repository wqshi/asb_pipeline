#/bin/bash
f_substr_exists()
{   
    whole_str=$1
    pattern_str=$2
    echo $whole_str | grep -q -i $pattern_str
    var=$?
    return $var
}


f_generate_tmp_ID()
{
     head -n 100 /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 15 | head -n 1
    #NEW_UUID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 15 | head -n 1)
    #echo $NEW_UUID
    #date | md5sum | fold -w 10 | head -1
}



f_get_bed_feature()
{
    loc_file=$1
    feature_file=$2 
    bedtools intersect -a $loc_file -b $feature_file -loj | awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $NF}' | sed 's/\./0/g'
}

f_get_bed_feature_sorted()
{
    loc_file=$1
    feature_file=$2
    f_create_header "#chr start end name"
    sort -k1,1 -k2,2n $loc_file | bedtools intersect -a stdin -b $feature_file -loj -sorted | awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $8}' | sed 's/\t\./\t0/g'
}

f_get_peak_max_for_snv()
{
    loc_file=$1
    feature_file=$2
    tmpID=$(f_generate_tmp_ID)
    cd $(dirname $loc_file)
    #Get the peak max position for the input file
    awk -F"\t" -v OFS="\t" '{print $1,$2,$3,int(($2+$3)/2)}' $feature_file | sort -k1,1 -k2,2n > tmp.bed.$tmpID

    f_create_header "#chr start end name"
    sort -k1,1 -k2,2n $loc_file | bedtools intersect -a stdin -b tmp.bed.$tmpID -loj -sorted | awk -F'\t' -v OFS="\t" '{print $1, $2, $3, $8}' | sed 's/\./0/g'
    \rm tmp.bed.$tmpID
}

f_file_exist()
{
    if [ ! -f $1 ];then
	echo "File Not EXIST"
	exit -1
    fi
	
}

f_sort_bed()
{
    sort -k1,1 -k2,2n $1
}

f_create_header()
{
    input=$1
    header=`echo $input | sed 's/ /\t/g'`
    printf "%s\n" "$header"
}

f_cat_sort()
{
    cat $1 $2 | sort -k1,1 -k2,2n
}

f_merge_bed_files()
{
    cell=$1
    methyl_dir=$2
    pattern=$3
    file_list=($(ls -d $methyl_dir/* | grep -i "$pattern"))
    
    f_create_header "#chr start end name"
    if [ ${#file_list[*]} -eq 2 ]
    then
        f_cat_sort ${file_list[0]} ${file_list[1]} | \
	bedtools merge -i stdin -scores sum 
    else
	echo "Error"
    fi
}


f_bed_length()
{
    bed_file=$1
    awk -F"\t" -v OFS="\t" '{sum+=$3-$2+1} END{print sum}' $bed_file
}

f_get_tf_name()
{
    #Get the tf name from the ENCODE file name
    encode_file_name=$1
    echo $encode_file_name | sed -e 's/\([A-Z][a-z]*\)/ &/g' |  awk -F'[ ]' '{print $6}'
}

f_getcell_name()
{
    #Get the tf name from the ENCODE file name
    encode_file_name=$1
    echo $encode_file_name | sed -e 's/\([A-Z][a-z0-9]*\)/ &/g' |  awk -F'[ ]' '{print $6}'
}


f_unzip_file()
{
    file=$1
    if ( echo $file | grep -q "[.]gz$")
    then
	gzip -d -c -f $file
    elif ( echo "$file" | grep -q "[.]bz2$" )
    then
	bzip2 -d -f -c $file
    elif ( echo "$file" | grep -q "[.]tgz$" )
    then
	tar -zxvf -O $file
    else
	echo "Can't find function to unzip such a file:$file"
    fi
    
}


f_DEBUG()
{
  if [ "$DEBUG" == "true" ]; then
      $@
      echo $@　　
  fi
}
