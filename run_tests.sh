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
	let num_failed=num_failed+1
    fi
    echo $STATUS $cmd
}

#############################################
##
##
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
must_fail ./src/fastq_info tests/test_e10.fastq.gz
must_fail ./src/fastq_info tests/test_e13.fastq.gz 
must_fail ./src/fastq_info tests/test_e14.fastq.gz 
must_fail ./src/fastq_info tests/test_e15.fastq.gz 
must_fail ./src/fastq_info tests/test_e16.fastq.gz 
must_fail ./src/fastq_info tests/test_e17.fastq.gz
must_fail ./src/fastq_info tests/test_e19_1.fastq.gz  tests/test_e19_2.fastq.gz
must_fail ./src/fastq_info tests/test_e19_2.fastq.gz  tests/test_e19_1.fastq.gz 
must_fail ./src/fastq_info -f tests/test_dot.fastq.gz 
##must_fail ./src/fastq_info tests/c18_1M_1.fastq.gz tests/c18_1M_2.fastq.gz 
##
## just checks the exit status
must_succeed ./src/fastq_info tests/test_dot.fastq.gz 
must_succeed 	./src/fastq_info tests/test_1.fastq.gz
must_succeed 	./src/fastq_info tests/test_30_1.fastq.gz  tests/test_30_2.fastq.gz 
must_succeed 	./src/fastq_info tests/test_2.fastq.gz 
must_succeed 	./src/fastq_info tests/test_13.fastq.gz
must_succeed 	./src/fastq_info tests/test_17.fastq.gz
must_succeed 	./src/fastq_info tests/test_21_1.fastq.gz
must_succeed 	./src/fastq_info tests/test_21_1.fastq.gz tests/test_21_2.fastq.gz 
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
must_succeed 	time -p ./src/fastq_info  tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_2.fastq.gz

must_fail "./src/fastq_info --help"
##
gcov src/fastq_info


echo "*** fastq_filterpair"
must_succeed "./src/fastq_filterpair tests/test_2.fastq.gz tests/test_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz && diff <(zcat f1.fastq.gz) <(zcat tests/test_2.fastq.gz)"
must_succeed "./src/fastq_filterpair tests/a_1.fastq.gz tests/a_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz && diff <(zcat f2.fastq.gz) <(zcat tests/a_2.fastq.gz) && diff <(zcat f1.fastq.gz) <(zcat tests/a_1.fastq.gz)"
must_succeed "./src/fastq_filterpair tests/casava.1.8_2.fastq.gz tests/casava.1.8_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed "./src/fastq_filterpair tests/casava.1.8_1.fastq.gz tests/casava.1.8_1.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed ./src/fastq_filterpair tests/c18_10000_1.fastq.gz tests/c18_10000_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz

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


must_fail "./src/fastq_pre_barcodes"

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz --index2 tests/barcode_test2_1.fastq.gz --index3 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index5  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index2 --cell_offset 0 --cell_size 8 --sample_read index3 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --read2 tests/barcode_test2_2.fastq.gz --outfile1 test_1.fastq.gz --outfile2 test_2.fastq.gz "



must_fail "./src/fastq_pre_barcodes --index1 tests/test_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz "

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz --outfile1 test.fastq.gz --read2 tests/test1_1.fastq.gz"

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read1 tests/barcode_test2_2.fastq.gz "

must_fail "./src/fastq_pre_barcodes --index1 tests/barcode_test2_1.fastq.gz  --phred_encoding 33 --min_qual 1 --umi_read index1  --umi_offset 0 --umi_size 16 --read2_offset 0 --read2_size -1 --cell_read index1 --cell_offset 0 --cell_size 8 --sample_read read1 --sample_offset 0  --sample_size 4 --read2 tests/barcode_test2_2.fastq.gz "

must_succeed "./src/fastq_pre_barcodes --help"

gcov src/fastq_pre_barcodes
echo "*** bam_umi_count"
#

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --ucounts xx --ucounts_MM -x TX"
must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --ucounts xx --ucounts_tsv -x TX"

must_succeed  " [ `./src/bam_umi_count --bam tests/test_annot.bam  --ucounts /dev/stdout |wc -l |cut -f 1 -d\ ` ==  90 ]"
#must_succeed  "./src/bam_umi_count --bam tests/test_annot.bam  --ucounts test.tmp && ./tests/check_no_dups.sh test.tmp"

must_succeed  " [ `./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam  --ucounts /dev/stdout |wc -l |cut -f 1 -d\ ` ==  83 ]"

must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam  --ucounts lixo --ucounts_MM"

must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot3.bam  --ucounts lixo --ucounts_MM"

must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3.bam  --ucounts lixo --ucounts_MM"

must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3.bam -x TX --ucounts lixo --ucounts_MM"

must_fail   "./src/bam_umi_count --min_reads 10 --bam tests/trans.bam  --ucounts lixo --ucounts_MM"

#must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam  --ucounts test.tmp && ./tests/check_no_dups.sh test.tmp"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_umi tests/known_umis.txt --ucounts /dev/null --dump xx && diff -q <(sort -k 1,2 xx) <(sort -k1,2 tests/out_known_umis.txt)"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx --uniq_mapped"

must_succeed  " ./src/bam_umi_count --min_reads 10 --bam tests/test_annot2.bam --ucounts xx"


must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx -x TX"

must_succeed  " ./src/bam_umi_count --min_reads 4 --bam tests/test_annot2.bam --ucounts xx"

must_succeed  "[ `./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_cells tests/known_cells.txt --ucounts /dev/stdout --dump xx | wc -l ` -eq 3 ]"



must_fail "./src/bam_umi_count --min_reads 1"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam -x"

must_fail "./src/bam_umi_count --bam tests/test_annot.bam_missing"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam_missing --ucounts folder/missing_path/xx"
must_fail "./src/bam_umi_count --bam tests/test_annot.bam_missing --ucounts xx --dump folder/missing_path/xxx"
must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_umi tests/known_umis.txt_missing --ucounts /dev/null --dump xx "

must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_cells tests/known_cells.txt_missing --ucounts /dev/null --dump xx "


rm -f xx xy

must_fail "./src/bam_umi_count"
must_succeed "./src/bam_umi_count --help"

gcov src/bam_umi_count
echo "*** bam_add_tags"

must_succeed "./src/bam_add_tags --inbam tests/trans.bam --outbam tmp.bam"
must_succeed "./src/bam_add_tags --inbam tests/trans.bam --outbam tmp.bam --tx --tx_2_gx tests/mapTrans2Gene.tsv"
must_succeed "./src/bam_add_tags --inbam tests/trans.bam --outbam tmp.bam --tx "
must_fail "./src/bam_add_tags --inbam tests/trans.bam --outbam tmp.bam --tx --tx_2_gx aaaatests/mapTrans2Gene.tsv"
rm -f tmp.bam

must_fail "./src/bam_add_tags"
must_fail "./src/bam_add_tags --inbam tests/trans.bam --outbam tmp.bam  --tx_2_gx tests/mapTrans2Gene.tsv"
must_fail "./src/bam_add_tags --inbam tests/trans.bam_missing --outbam tmp.bam  --tx --tx_2_gx tests/mapTrans2Gene.tsv"
must_fail "./src/bam_add_tags --inbam tests/trans.bam_missing --outbam folder/does/not/exist/tmp.bam  --tx --tx_2_gx tests/mapTrans2Gene.tsv"
must_succeed "./src/bam_add_tags --help"


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

echo Failed tests: $num_failed
exit $num_failed
