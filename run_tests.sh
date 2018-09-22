#!/bin/bash

let num_failed=0
function must_fail {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -eq $? ];  then	
	STATUS=FAILED
	let num_failed=num_failed+1
    fi
    echo $STATUS $cmd
}

function must_succeed {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -ne $? ];  then	
	STATUS=FAILED
	bash -c "$cmd"
	let num_failed=num_failed+1
    fi
    echo $STATUS $cmd
}


export PATH=$PWD/bin:$PATH
#############################################
##
##
echo "*** bam_umi_count"
#
must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --ucounts xx  -x TX"
must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --ucounts xx  -x TX"

must_succeed  " [ `./src/bam_umi_count --bam tests/test_annot.bam  --ucounts xx && grep -v % xx |wc -l |cut -f 1 -d\ ` ==  89 ]"
#must_succeed  "./src/bam_umi_count --bam tests/test_annot.bam  --ucounts test.tmp && ./tests/check_no_dups.sh test.tmp"

must_succeed  " [ `./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam  --multi_mapped --ucounts xx && grep -v % xx |wc -l |cut -f 1 -d\ ` ==  89 ]"

must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam  --ucounts lixo "

must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot3_small.bam  --ucounts lixo "

must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3_small.bam  --ucounts lixo"

must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3_small.bam -x TX --ucounts lixo"

must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3_small.bam -x TX --ucounts lixo --ignore_sample"

must_fail   "./src/bam_umi_count --min_reads 10 --bam tests/trans_small.bam  --ucounts lixo --ucounts_MM"

#must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam  --ucounts test.tmp && ./tests/check_no_dups.sh test.tmp"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_umi tests/known_umis.txt --ucounts xxu --rcounts xxr && diff -q <(sort -k1,2 xxu) <(sort -k1,2 tests/test_annot.mtx)"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx --uniq_mapped"

must_succeed  " ./src/bam_umi_count --min_reads 10 --bam tests/test_annot2.bam --ucounts xx"


must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx -x TX --max_cells 10 --max_feat 2 --feat_cell 2"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx -x TX --max_cells 20000 --max_feat 2 --feat_cell 2"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx -x TX"

must_succeed  " ./src/bam_umi_count --min_reads 4 --bam tests/test_annot2.bam --ucounts xx"

must_succeed  " ./src/bam_umi_count --min_reads 4 --bam tests/test_annot2.bam --ucounts xx --ignore_sample"

must_succeed  " [ `./src/bam_umi_count --min_reads 1 --bam tests/test_annot4.bam --ucounts xx --ignore_sample --cell_suffix '-123456789' && grep -c 123456789 xx_cols ` -eq  9 ]"


must_succeed  "[ `./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_cells tests/known_cells.txt --ucounts xx && cat xx  | wc -l ` -eq 4 ]"



must_fail "./src/bam_umi_count --min_reads 1"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam -x"

must_fail "./src/bam_umi_count --bam tests/test_annot.bam_missing"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam_missing --ucounts folder/missing_path/xx"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam_missing --ucounts xx --ureads folder/missing_path/xxx"
must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_umi tests/known_umis.txt_missing --ucounts /dev/null --dump xx "

must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_cells tests/known_cells.txt_missing --ucounts /dev/null --dump xx "

must_fail  "./src/bam_umi_count --sorted_by_cell --min_reads 1 --bam tests/test_annot.bam --known_cells tests/known_cells.txt_missing --ucounts /dev/null --dump xx "
must_succeed  "./bin/samtools sort -t CR tests/test_annot2.bam > /dev/null"
must_succeed  "./bin/samtools sort -t CR tests/test_annot2.bam | ./src/bam_umi_count --sorted_by_cell  --min_reads 4 --bam - --ucounts xx --ignore_sample"

must_succeed  "./bin/samtools sort -t CR tests/test_annot2.bam | ./src/bam_umi_count --sorted_by_cell  --min_reads 4 --bam - --ucounts xx --rcounts xy --ignore_sample"


rm -f xx* xy*

must_fail "./src/bam_umi_count"
must_succeed "./src/bam_umi_count --help"
must_succeed "./src/bam_umi_count -h"

gcov src/bam_umi_count


echo "*** fastq_trim_poly_at"
must_fail "./src/fastq_trim_poly_at"
must_succeed "./src/fastq_trim_poly_at --help"
must_fail "./src/fastq_trim_poly_at --file tests/a_1.fastq.gz"
must_fail "./src/fastq_trim_poly_at --file tests/a_1xx.fastq.gz --outfile tmp.fastq.gz"
must_fail "./src/fastq_trim_poly_at --file tests/a_1xx.fastq.gz --outfile /xxx/tmp.fastq.gz"

# no diff
must_succeed "./src/fastq_trim_poly_at --file tests/a_1.fastq.gz --outfile tmp.fastq.gz --min_poly_at_len 20 && diff <(zcat tests/a_1.fastq.gz) <(zcat tmp.fastq.gz) "

must_succeed "./src/fastq_trim_poly_at --file tests/poly_at.fastq.gz --outfile tmp.fastq.gz --min_poly_at_len 3 && diff <(zcat tests/poly_at_len3.fastq.gz) <(zcat tmp.fastq.gz) "

must_succeed "./src/fastq_trim_poly_at --file tests/poly_at.fastq.gz --outfile tmp.fastq.gz --min_poly_at_len 300 --min_len 1 && diff <(zcat tests/poly_at.fastq.gz) <(zcat tmp.fastq.gz) "

must_succeed "zcat tests/poly_at.fastq.gz | ./src/fastq_trim_poly_at --file - --outfile tmp.fastq.gz --min_poly_at_len 300 --min_len 1 && diff <(zcat tests/poly_at.fastq.gz) <(zcat tmp.fastq.gz) " 

must_succeed "diff <(zcat tests/poly_at.fastq.gz | ./src/fastq_trim_poly_at --file - --outfile -  --min_poly_at_len 300 --min_len 1|zcat ) <(zcat tests/poly_at.fastq.gz) "

gcov src/fastq_trim_poly_at


echo "*** fastq_filter_n"
must_succeed "./src/fastq_filter_n tests/test_21_2.fastq.gz > tmp && diff /dev/null tmp"
must_fail "./src/fastq_filter_n -n 100 tests/test_21_2.fastq.gz > tmp && diff -q /dev/null tmp"
must_fail "./src/fastq_filter_n tests/test_21_2.fastq.gz > tmp && diff  tests/test_21_2.fastq.gz tmp"
must_succeed "./src/fastq_filter_n tests/test_1.fastq.gz > tmp && diff -q <(zcat tests/test_1.fastq.gz) tmp"
must_fail "./src/fastq_filter_n --help"
must_fail "./src/fastq_filter_n"

gcov src/fastq_filter_n

##
echo "*** fastq_num_reads"
must_succeed "[ `./src/fastq_num_reads tests/test_21_2.fastq.gz` -eq 2 ]"
must_fail "[ `./src/fastq_num_reads tests/test_21_2.fastq.gz` -ne 2 ]"
must_succeed "[ `./src/fastq_num_reads tests/c18_10000_1.fastq.gz` -eq 10000 ]"
must_fail "[ `./src/fastq_num_reads tests/c18_10000_1.fastq.gz` -ne 10000 ]"
must_succeed "[ `./src/fastq_num_reads tests/one.fastq.gz` -eq 1 ]"
must_fail "./src/fastq_num_reads --help"
must_fail "./src/fastq_num_reads"

gcov src/fastq_num_reads

##
echo "*** fastq_truncate"
must_succeed "[ `./src/fastq_truncate tests/test_21_2.fastq.gz 1|wc -l` -eq 4 ]"
must_succeed "[ `./src/fastq_truncate tests/test_21_2.fastq.gz 0|wc -l` -eq 0 ]"
must_succeed "[ `./src/fastq_truncate tests/test_21_2.fastq.gz 2|wc -l` -eq 8 ]"
must_fail "./src/fastq_truncate tests/test_21_2.fastq.gz"
must_fail "./src/fastq_truncate --help"

gcov src/fastq_truncate

echo "*** fastq_empty"
must_succeed "./src/fastq_not_empty tests/barcode_test_1.fastq.gz" 
must_succeed "./src/fastq_not_empty tests/test_21_2.fastq.gz"
must_fail "./src/fastq_not_empty tests/test_21_2.fastq.gzaa"
must_fail "./src/fastq_not_empty /dev/null"
must_fail "./src/fastq_not_empty --help"
must_fail "./src/fastq_not_empty"

gcov src/fastq_not_empty

echo "*** fastq_info"
must_fail ./src/fastq_info 
must_fail ./src/fastq_info tests/test_e1.fastq.gz 
must_fail ./src/fastq_info tests/test_e2.fastq.gz
must_fail ./src/fastq_info tests/test_e3.fastq.gz 
must_fail ./src/fastq_info tests/test_e4.fastq.gz 
must_fail ./src/fastq_info tests/test_e5.fastq.gz 
must_fail ./src/fastq_info tests/test_e6.fastq.gz 
must_fail ./src/fastq_info tests/test_e7.fastq.gz 
must_fail ./src/fastq_info tests/test_e8.fastq.gz 
must_fail ./src/fastq_info tests/test_e9.fastq.gz
must_succeed ./src/fastq_info -r tests/test_e9.fastq.gz
must_fail ./src/fastq_info tests/test_e10.fastq.gz
must_fail ./src/fastq_info tests/test_e13.fastq.gz 
must_fail ./src/fastq_info tests/test_e14.fastq.gz 
must_fail ./src/fastq_info tests/test_e15.fastq.gz 
must_fail ./src/fastq_info tests/test_e16.fastq.gz
must_fail ./src/fastq_info -r tests/test_e10.fastq.gz
must_fail ./src/fastq_info -r tests/test_e13.fastq.gz 
must_fail ./src/fastq_info -r tests/test_e14.fastq.gz 
must_fail ./src/fastq_info -r tests/test_e15.fastq.gz 
must_fail ./src/fastq_info -r tests/test_e16.fastq.gz 
must_fail ./src/fastq_info tests/test_e17.fastq.gz
must_fail ./src/fastq_info tests/test_e19_1.fastq.gz  tests/test_e19_2.fastq.gz
must_fail ./src/fastq_info tests/test_e19_2.fastq.gz  tests/test_e19_1.fastq.gz
must_fail ./src/fastq_info tests/test_e19_1.fastq.gz tests/test_empty.fastq.gz  
must_fail ./src/fastq_info tests/test_empty.fastq.gz  tests/test_e19_1.fastq.gz
must_fail ./src/fastq_info -r -s tests/test_e19_1.fastq.gz  tests/test_e19_2.fastq.gz
must_fail ./src/fastq_info -r -s tests/test_e19_2.fastq.gz  tests/test_e19_1.fastq.gz 
must_fail ./src/fastq_info -f tests/test_dot.fastq.gz
must_fail ./src/fastq_info tests/test_empty.fastq.gz
must_fail ./src/fastq_info -r  tests/test_empty.fastq.gz
must_fail ./src/fastq_info -s -r  tests/test_empty.fastq.gz tests/test_1.fastq.gz
must_fail ./src/fastq_info -s -r  tests/test_1.fastq.gz tests/test_empty.fastq.gz 
##must_fail ./src/fastq_info tests/c18_1M_1.fastq.gz tests/c18_1M_2.fastq.gz 
##
## just checks the exit status
must_succeed ./src/fastq_info -h 
must_succeed ./src/fastq_info tests/test_dot.fastq.gz 
must_succeed 	./src/fastq_info tests/test_1.fastq.gz
must_succeed 	./src/fastq_info tests/test_30_1.fastq.gz  tests/test_30_2.fastq.gz 
must_succeed 	./src/fastq_info tests/test_2.fastq.gz 
must_succeed 	./src/fastq_info tests/test_13.fastq.gz
must_succeed 	./src/fastq_info tests/test_17.fastq.gz
must_succeed 	./src/fastq_info tests/test_21_1.fastq.gz
must_succeed 	./src/fastq_info tests/test_21_1.fastq.gz tests/test_21_2.fastq.gz
must_succeed 	./src/fastq_info -r -s tests/test_21_1.fastq.gz tests/test_21_2.fastq.gz 
must_succeed 	time -p ./src/fastq_info tests/pe_bug14.fastq.gz tests/pe_bug14.fastq.gz 
#must_succeed 	time -p ./src/fastq_info tests/c18_1M_1.fastq.gz 
#must_succeed 	time -p ./src/fastq_info tests/c18_1M_2.fastq.gz 
#must_succeed 	time -p ./src/fastq_info tests/c18_1M_1.fastq.gz 
#must_succeed 	time -p ./src/fastq_info tests/c18_1M_2.fastq.gz 
#must_succeed 	time -p ./src/fastq_info tests/c18_1M_1.fastq.gz tests/c18_1M_1.fastq.gz 
must_succeed 	time -p ./src/fastq_info tests/casava.1.8i.fastq.gz pe
must_succeed 	time -p ./src/fastq_info tests/test_solid_1.fastq.gz tests/test_solid_2.fastq.gz
must_fail 	time -p ./src/fastq_info tests/test_solid2_1.fastq.gz tests/test_solid2_2.fastq.gz
must_succeed 	time -p ./src/fastq_info tests/solexa_1.fastq.gz tests/solexa_2.fastq.gz
must_succeed 	time -p ./src/fastq_info  tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_readname_trunc_2.fastq.gz
must_succeed 	time -p ./src/fastq_info -s  tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_readname_trunc_2.fastq.gz
must_succeed 	time -p ./src/fastq_info  tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_2.fastq.gz
must_succeed 	time -p ./src/fastq_info  -r -s tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_2.fastq.gz

must_fail "./src/fastq_info --help"
##
gcov src/fastq_info


echo "*** fastq_filterpair"
must_succeed "./src/fastq_filterpair tests/test_2.fastq.gz tests/test_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz && diff <(zcat f1.fastq.gz) <(zcat tests/test_2.fastq.gz)"
must_succeed "./src/fastq_filterpair tests/a_1.fastq.gz tests/a_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz && diff <(zcat f2.fastq.gz) <(zcat tests/a_2.fastq.gz) && diff <(zcat f1.fastq.gz) <(zcat tests/a_1.fastq.gz)"
must_succeed "./src/fastq_filterpair tests/casava.1.8_2.fastq.gz tests/casava.1.8_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed "./src/fastq_filterpair tests/casava.1.8_1.fastq.gz tests/casava.1.8_1.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed ./src/fastq_filterpair tests/c18_10000_1.fastq.gz tests/c18_10000_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz
must_succeed ./src/fastq_filterpair tests/c18_10000_1.fastq.gz tests/c18_10000_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz sorted

must_fail ./src/fastq_filterpair tests/c18_10000_1.fastq.gz tests/casava.1.8_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz

must_fail ./src/fastq_filterpair tests/c18_10000_1.fastq_missing.gz tests/c18_10000_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz

must_fail ./src/fastq_filterpair tests/c18_10000_1.fastq.gz tests/c18_10000_2.fastq.gz  folder/does/not/exist/f1.fastq.gz f2.fastq.gz up.fastq.gz

must_fail "./src/fastq_filterpair --help"
#must_succeed ./src/fastq_filterpair tests/c18_1M_2.fastq.gz tests/c18_1M_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz 
#must_succeed ./src/fastq_filterpair tests/c18_1M_1.fastq.gz tests/c18_1M_1.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz 

#echo "This may take a while..."
#must_succeed ./src/fastq_filterpair tests/c18_1M_1.fastq.gz tests/c18_1M_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz
gcov src/fastq_filterpair

##
echo "*** fastq_pre_barcodes"
must_succeed ./src/fastq_pre_barcodes --index1 tests/barcode_test_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --read1 tests/barcode_test_2.fastq.gz --outfile1 test.fastq.gz

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz && diff -q  <(zcat test.fastq.gz)  <(zcat tests/pre1.fastq.gz)"

rm -f test.fastq.gz
must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz && diff -q  <(zcat test.fastq.gz)  <(zcat tests/pre1.fastq.gz)"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz && diff -q  <(zcat test.fastq.gz)  <(zcat tests/pre2.fastq.gz)"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz && diff -q  <(zcat test.fastq.gz)  <(zcat tests/pre3.fastq.gz)"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz --index2 tests/barcode_test2_1.fastq.gz --index3 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index2 --cell_offset 0 --cell_size 8 --sample_read index3 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz --index2 tests/barcode_test2_1.fastq.gz --index3 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index2 --cell_offset 0 --cell_size 8 --sample_read index3 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --read2 tests/barcode_test2_2.fastq.gz --outfile1 test_1.fastq.gz --outfile2 test_2.fastq.gz && diff -q test_1.fastq.gz test_2.fastq.gz"

# sam
must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz --index2 tests/barcode_test2_1.fastq.gz --index3 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index2 --cell_offset 0 --cell_size 8 --sample_read index3 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --read2 tests/barcode_test2_2.fastq.gz --outfile1 test_1.fastq.gz --outfile2 test_2.fastq.gz  --sam"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz --index2 tests/barcode_test2_1.fastq.gz --index3 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index2 --cell_offset 0 --cell_size 8 --sample_read index3 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz --sam"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz --sam"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz --sam"

must_succeed "./src/fastq_pre_barcodes --index1 tests/barcode_test_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --read1 tests/barcode_test_2.fastq.gz --outfile1 test.fastq.gz --sam"

must_fail "./src/fastq_pre_barcodes --interleaved read1 --read1 tests/inter.fastq.gz --index1 tests/inter.fastq.gz  --outfile1 x.gz  --umi_read index1  --umi_offset 0 --umi_size 16 --sam"

must_fail "./src/fastq_pre_barcodes --interleaved read --read1 tests/inter.fastq.gz --index1 tests/inter.fastq.gz  --outfile1 x.gz  --umi_read index1  --umi_offset 0 --umi_size 16 --sam"

must_fail "./src/fastq_pre_barcodes --interleaved read1,read2,index1 --read1 tests/inter.fastq.gz --index1 tests/inter.fastq.gz  --outfile1 x.gz  --umi_read index1  --umi_offset 0 --umi_size 16 --sam"

must_succeed "./src/fastq_pre_barcodes --interleaved index1,read1 --read1 tests/casava.1.8i.fastq.gz --index1 tests/casava.1.8i.fastq.gz  --outfile1 x1.gz  --umi_read index1  --umi_offset 0 --umi_size 16"

must_succeed "./src/fastq_pre_barcodes --interleaved read1,index1 --read1 tests/inter.fastq.gz --index1 tests/inter.fastq.gz  --outfile1 x2.gz  --umi_read index1  --umi_offset 0 --umi_size 16"

must_succeed "./src/fastq_pre_barcodes --interleaved index1,read1 --read1 tests/inter.fastq.gz --index1 tests/inter.fastq.gz  --outfile1 x1.gz  --umi_read index1  --umi_offset 0 --umi_size 16"

must_succeed "./src/fastq_pre_barcodes --interleaved read1,index1 --read1 tests/inter.fastq.gz --index1 tests/inter.fastq.gz  --outfile1 x3.gz  --umi_read index1  --umi_offset 0 --umi_size 16 --sam"

##
must_succeed "./src/fastq_split_interleaved tests/inter.fastq.gz xxx && [ -e xxx_1.fastq.gz ] && [ -e xxx_2.fastq.gz ]"

## ==x2 --interleaved read1,index1
must_succeed "./src/fastq_pre_barcodes  --read1 xxx_1.fastq.gz --index1 xxx_2.fastq.gz  --outfile1 xni.gz  --umi_read index1  --umi_offset 0 --umi_size 16 && diff xni.gz x2.gz"

## ==x1
must_succeed "./src/fastq_pre_barcodes  --read1 xxx_2.fastq.gz --index1 xxx_1.fastq.gz  --outfile1 xni.gz  --umi_read index1  --umi_offset 0 --umi_size 16 && diff xni.gz x1.gz"

rm -f xxx_*.fastq.gz x?.gz xni.gz

must_fail "./src/fastq_pre_barcodes"

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz --index2 tests/barcode_test2_1.fastq.gz --index3 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index5  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index2 --cell_offset 0 --cell_size 8 --sample_read index3 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --read2 tests/barcode_test2_2.fastq.gz --outfile1 test_1.fastq.gz --outfile2 test_2.fastq.gz "



must_fail "./src/fastq_pre_barcodes --index1 tests/test_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz "

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz --read2 tests/test1_1.fastq.gz"

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz "

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read2_offset 0 --read2_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read2 tests/barcode_test2_2.fastq.gz "

must_succeed "./src/fastq_pre_barcodes --help"

## wrapper to fastq_pre_barcodes
must_fail ./sh/fastq2bam  -b test.bam -s 10xV1a -1 tests/10xv1a_R1.fastq.gz -3 tests/10xv1a_R3.fastq.gz -2 tests/10xv1a_R2.fastq.gz 

fastq2bam -b test.bam -s 10xV1a -1 tests/10xv1a_R1.fastq.gz -2 tests/10xv1a_R3.fastq.gz -3 tests/10xv1a_R2.fastq.gz -4 tests/10xv1a_I1.fastq.gz

must_succeed ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b lixo -3 tests/tx.I2.fastq.gz
must_succeed ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b lixo2
must_fail diff <(samtools view lixo) <(samtools view lixo2)
LL="1	4	*	0	255	*	*	0	98	GAGACCATGCTCAACAGCAACATCAATGACCTGCTGATGGTGACCTACCTGGCCAATCTCACCCAGTCACAGATTGCCCTCAACGAGAAACTTGTAAA	3>BBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGBG:CFCGGGEF0GGGFGGECCFEG?FGGE>FG/;;?GFGFGG:C0	on:Z:D00408:393:CAYR8ANXX:1:1101:2629:2244@1:N:0:0	op:Z:3>BBCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGBG:CFCGGGEF0GGGFGGECCFEG?FGGE>FG/;;?GFGFGG:C0	QX:Z:AAAAACGGAT	OQ:Z:ABB=:///??	CR:Z:GTGGATTGCCTAAG	CY:Z:B@BBBEGBGGGGCF	BC:Z:ACCGAACA	QT:Z::@BBBGG/"
LLG=$(samtools view lixo)
must_succeed "[ '$LL-' == '$LLG-' ]"
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b 
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz   -b lixo2
must_fail ./sh/fastq2bam -s 10xV1i -2 tests/tx.I1.fastq.gz -b lixo2
rm -f lixo lixo2
gcov src/fastq_pre_barcodes
echo "*** bam_add_tags"

must_succeed "./src/bam_add_tags --inbam tests/trans_small.bam --outbam tmp.bam"
must_succeed "./src/bam_add_tags --inbam tests/trans_small.bam --outbam tmp.bam --tx --tx_2_gx tests/mapTrans2Gene.tsv"
must_succeed "./src/bam_add_tags --inbam tests/trans_small.bam --outbam tmp.bam --tx "
must_fail "./src/bam_add_tags --inbam tests/trans_small.bam --outbam tmp.bam --tx --tx_2_gx aaaatests/mapTrans2Gene.tsv"
rm -f tmp.bam

must_fail "./src/bam_add_tags"
must_fail "./src/bam_add_tags --inbam tests/trans_small.bam --outbam tmp.bam  --tx_2_gx tests/mapTrans2Gene.tsv"
must_fail "./src/bam_add_tags --inbam tests/trans_small.bam_missing --outbam tmp.bam  --tx --tx_2_gx tests/mapTrans2Gene.tsv"
must_fail "./src/bam_add_tags --inbam tests/trans_small.bam_missing --outbam folder/does/not/exist/tmp.bam  --tx --tx_2_gx tests/mapTrans2Gene.tsv"
must_succeed "./src/bam_add_tags --help"
must_succeed "./src/bam_add_tags -h"


gcov src/bam_add_tags

echo "*** fastq_split_interleaved"
must_succeed 	time -p ./src/fastq_split_interleaved tests/casava.1.8i.fastq.gz   out_prefix

must_fail 	time -p ./src/fastq_split_interleaved tests/casava.1.8i_e1.fastq.gz   out_prefix
must_fail 	./src/fastq_split_interleaved tests/casava.1.8i.fastq.gz a1 a2
must_fail 	./src/fastq_split_interleaved
must_fail 	./src/fastq_split_interleaved tests/one.fastq.gz out_prefix
must_fail 	./src/fastq_split_interleaved tests/test_e1.fastq.gz

rm -f out_prefix_*.fastq.gz

gcov src/fastq_split_interleaved

echo "*** fastq_validator.sh"
export PATH=$PWD/src:$PATH
must_fail ./sh/fastq_validator.sh tests/c18_10000_1.fastq.gz.bz2 tests/c18_10000_2.fastq.gz.bz2
must_fail ./sh/fastq_validator.sh tests/c18_10000_1.fastq.gz tests/c18_10000_2.fastq.gz

must_fail ./sh/fastq_validator.sh tests/SRR3587500_1.fastq.gz.bz2 tests/SRR3587500_2.fastq.gz.bz2.
must_succeed ./sh/fastq_validator.sh tests/read-I1_si-ACCGAACA_lane-001-chunk-001.fastq.gz tests/read-I2_si-ACCGAACA_lane-001-chunk-001.fastq.gz
must_succeed ./sh/fastq_validator.sh tests/read-I1_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2 tests/read-I2_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2
must_succeed ./sh/fastq_validator.sh tests/read-I1_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2 tests/read-I2_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2 

echo Failed tests: $num_failed
exit $num_failed

