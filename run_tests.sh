#!/usr/bin/env bash

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
echo "*** fastq_split_interleaved"
must_succeed 	time -p ./src/fastq_split_interleaved tests/casava.1.8i.fastq.gz   out_prefix

must_fail 	time -p ./src/fastq_split_interleaved tests/casava.1.8i_e1.fastq.gz   out_prefix
must_fail 	./src/fastq_split_interleaved tests/casava.1.8i.fastq.gz a1 a2
must_fail 	./src/fastq_split_interleaved
must_fail 	./src/fastq_split_interleaved tests/one.fastq.gz out_prefix
must_fail 	./src/fastq_split_interleaved tests/test_e1.fastq.gz
must_fail        ./src/fastq_split_interleaved 
must_fail      ./src/fastq_split_interleaved   tests/test_21_2.fastq.gz xxx
must_succeed "./src/fastq_split_interleaved tests/inter.fastq.gz tests/xxx && [ -e tests/xxx_1.fastq.gz ] && [ -e tests/xxx_2.fastq.gz ]"

##
echo "*** bam2fastq"
rm -f tmpf*.fastq*
must_fail "./src/bam2fastq"
must_fail "./src/bam2fastq -i "
must_fail "./src/bam2fastq -o "
must_fail "./src/bam2fastq --10x "
must_fail "./src/bam2fastq --bam  tests/no_qual.bam"
must_fail "./src/bam2fastq --bam  tests/missing_no_qual.bam --out tmpf"
must_succeed "./src/bam2fastq -h"
must_succeed "[ `./src/bam2fastq --bam  tests/no_qual.bam --out tmpf1 && ./src/fastq_info tmpf1.fastq.gz 2> /dev/null && zcat tmpf1.fastq.gz|wc -l|cut -f 1 -d\ ` \> 0 ] "

must_succeed " [ `./src/bam2fastq --bam  tests/test.bam --out tmpf2 && ./src/fastq_info tmpf2.fastq.gz 2> /dev/null && zcat tmpf2.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` \> 0 ]"

must_succeed " [ `./src/bam2fastq --bam  tests/test.bam --out tmpf2 -X && ./src/fastq_info tmpf2.fastq.gz 2> /dev/null && zcat tmpf2.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` \> 0 ]"

must_succeed " ./src/bam2fastq --bam  tests/test10.bam --out tmpf5  && [ ! -e tmpf5_umi.fastq.gz  ]"

must_succeed " ./src/bam2fastq --bam  tests/test10.bam --out tmpf6 -X  && [ -e tmpf6_umi.fastq.gz  ]"

must_succeed " [ `./src/bam2fastq --bam  tests/test_one_cell.bam --out tmpf3 && ./src/fastq_info tmpf3.fastq.gz 2> /dev/null && zcat tmpf3.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` == 4939 ]"

must_succeed " [ ` ./src/bam2fastq --bam  tests/test_annot.bam --out tmpf4 && ./src/fastq_info tmpf4.fastq.gz > /dev/null && zcat tmpf4.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` == 724 ]"

must_succeed " [ ` ./src/bam2fastq --bam  tests/test_annot2.bam --out tmpf5 && ./src/fastq_info tmpf5.fastq.gz > /dev/null && zcat tmpf5.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` == 57412 ]"

must_succeed " [ `./src/bam2fastq --bam  tests/trans.bam --out tmpf6 && ./src/fastq_info tmpf6.fastq.gz > /dev/null && zcat tmpf6.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` == 497978 ]"
must_succeed " [ `./src/bam2fastq --bam  tests/se.bam --out tmpf7 && ./src/fastq_info tmpf7.fastq.gz > /dev/null && zcat tmpf7.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` == 2 ]"
must_succeed " [ `./src/bam2fastq --bam  tests/pe.bam --out tmpf8 && ./src/fastq_info tmpf8_1.fastq.gz tmpf8_2.fastq.gz > /dev/null && zcat tmpf8_1.fastq.gz|grep '^@'|wc -l|cut -f 1 -d\ ` == 2 ]"
must_succeed "./src/bam2fastq --bam  tests/no_qual.bam --out tmpf1 "
#gcov src/bam2fastq

#
rm -f tmpf*.fastq*

echo "*** bam_umi_count"
#
must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --ucounts xx  -x TX --not_sorted_by_cell"

must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --ucounts xx  -x GX --not_sorted_by_cell"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam --ucounts xx  -x TX --not_sorted_by_cell"

must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam --ucounts xx  -x GX --not_sorted_by_cell"

must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam --ucounts /xx  -x GX --not_sorted_by_cell"

must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam --ucounts /xx  -x GX --not_sorted_by_cell -X"

must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam --ucounts xx  -x TX --not_sorted_by_cell --10x"
must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam --ucounts xx  -x GX --not_sorted_by_cell --10x"


#must_succeed  " [ `./src/bam_umi_count --bam tests/test_annot.bam  --ucounts xx --not_sorted_by_cell && grep -v % xx |wc -l |cut -f 1 -d\ ` ==  89 ]"
#must_succeed  "./src/bam_umi_count --bam tests/test_annot.bam  --ucounts test.tmp && ./tests/check_no_dups.sh test.tmp"

must_succeed  " [ `./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam  --multi_mapped --ucounts xx --not_sorted_by_cell && grep -v % xx |wc -l |cut -f 1 -d\ ` ==  89 ]"

#must_succeed  "./src/bam_umi_count --min_reads 1 --not_sorted_by_cell --bam tests/test_annot.bam  --ucounts tmpf "


#must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot3_small.bam  --ucounts tmpf --not_sorted_by_cell"

#must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3_small.bam  --ucounts tmpf --not_sorted_by_cell"

#must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3_small.bam -x TX --ucounts tmpf --not_sorted_by_cell"

#must_succeed  "./src/bam_umi_count --min_reads 10 --bam tests/test_annot3_small.bam -x TX --ucounts tmpf --ignore_sample --not_sorted_by_cell"

#must_fail   "./src/bam_umi_count --min_reads 10 --bam tests/trans_small.bam  --ucounts tmpf --ucounts_MM --not_sorted_by_cell"

#must_succeed  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam  --ucounts test.tmp && ./tests/check_no_dups.sh test.tmp"

#must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot.bam --known_umi tests/known_umis.txt --not_sorted_by_cell --ucounts xxu --rcounts xxr && diff -q <(sort -k1,2 xxu) <(sort -k1,2 tests/test_annot.mtx)"

#must_succeed  " ./src/bam_umi_count --min_reads 1 --not_sorted_by_cell --bam tests/test_annot2.bam --ucounts xx --not_sorted_by_cell"

#must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx --uniq_mapped --not_sorted_by_cell"

#must_succeed  " ./src/bam_umi_count --min_reads 10 --bam tests/test_annot2.bam --ucounts xx --not_sorted_by_cell"


must_fail  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx -x TX --max_cells 10 --max_feat 2 --feat_cell 2 --not_sorted_by_cell"

#must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx -x TX --max_cells 20000 --max_feat 2 --feat_cell 2 --not_sorted_by_cell"

#must_succeed  " ./src/bam_umi_count --min_reads 1 --bam tests/test_annot2.bam --ucounts xx -x TX --not_sorted_by_cell"

#must_succeed  " ./src/bam_umi_count --min_reads 4 --bam tests/test_annot2.bam --ucounts xx --not_sorted_by_cell"

#must_succeed  " ./src/bam_umi_count --min_reads 4 --bam tests/test_annot2.bam --ucounts xx --ignore_sample --not_sorted_by_cell"

must_succeed  " [ `./src/bam_umi_count --min_reads 1 --bam tests/test_annot5.bam --ucounts xx --ignore_sample --not_sorted_by_cell --cell_suffix '-123456789' && grep -c 123456789 xx_cols ` -eq  365 ]"

must_fail  "./src/bam_umi_count --min_reads 1 --bam tests/test_annot4.bam --ucounts xx --ignore_sample --not_sorted_by_cell --cell_suffix '-123456789' && grep -c 123456789 xx_cols"


must_succeed  "[ `./src/bam_umi_count --not_sorted_by_cell --min_reads 1 --bam tests/test_annot5.bam --known_cells tests/known_cells.txt --ucounts xx && cat xx  | wc -l ` -eq 4 ]"



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

#gcov src/bam_umi_count


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

#gcov src/fastq_trim_poly_at


echo "*** fastq_filter_n"
must_succeed "./src/fastq_filter_n tests/test_21_2.fastq.gz > tmp && diff /dev/null tmp"
must_fail "./src/fastq_filter_n -n 100 tests/test_21_2.fastq.gz > tmp && diff -q /dev/null tmp"
must_fail "./src/fastq_filter_n tests/test_21_2.fastq.gz > tmp && diff  tests/test_21_2.fastq.gz tmp"
must_succeed "./src/fastq_filter_n tests/test_1.fastq.gz > tmp && diff -q <(zcat tests/test_1.fastq.gz) tmp"
must_fail "./src/fastq_filter_n --help"
must_fail "./src/fastq_filter_n"

#gcov src/fastq_filter_n

##
echo "*** fastq_num_reads"
must_succeed "[ `./src/fastq_num_reads tests/test_21_2.fastq.gz` -eq 2 ]"
must_fail "[ `./src/fastq_num_reads tests/test_21_2.fastq.gz` -ne 2 ]"
must_succeed "[ `./src/fastq_num_reads tests/c18_10000_1.fastq.gz` -eq 10000 ]"
must_fail "[ `./src/fastq_num_reads tests/c18_10000_1.fastq.gz` -ne 10000 ]"
must_succeed "[ `./src/fastq_num_reads tests/one.fastq.gz` -eq 1 ]"
must_fail "./src/fastq_num_reads --help"
must_fail "./src/fastq_num_reads"

#gcov src/fastq_num_reads

##
echo "*** fastq_truncate"
must_succeed "[ `./src/fastq_truncate tests/test_21_2.fastq.gz 1|wc -l` -eq 4 ]"
must_succeed "[ `./src/fastq_truncate tests/test_21_2.fastq.gz 0|wc -l` -eq 0 ]"
must_succeed "[ `./src/fastq_truncate tests/test_21_2.fastq.gz 2|wc -l` -eq 8 ]"
must_fail "./src/fastq_truncate tests/test_21_2.fastq.gz"
must_fail "./src/fastq_truncate --help"

#gcov src/fastq_truncate

echo "*** fastq_empty"
must_succeed "./src/fastq_not_empty tests/barcode_test_1.fastq.gz" 
must_succeed "./src/fastq_not_empty tests/test_21_2.fastq.gz"
must_fail "./src/fastq_not_empty tests/test_21_2.fastq.gzaa"
must_fail "./src/fastq_not_empty /dev/null"
must_fail "./src/fastq_not_empty --help"
must_fail "./src/fastq_not_empty"

#gcov src/fastq_not_empty

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
must_fail ./src/fastq_info  tests/test_33.fastq.gz
must_succeed  ./src/fastq_info -q  tests/test_33.fastq.gz
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
must_succeed ./src/fastq_info -e tests/test_dot.fastq.gz
touch tests/empty.fastq
must_fail ./src/fastq_info tests/empty.fastq
must_succeed ./src/fastq_info -e tests/empty.fastq
must_succeed 	./src/fastq_info tests/test_1.fastq.gz
must_succeed 	"[ \"`./src/fastq_info tests/test_1.fastq.gz 2> /dev/stdout|grep 'Read length:'|cut -f 2 -d:`\" == \" 90 90 90\" ]"

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
must_fail 	time -p ./src/fastq_info  tests/casava.1.8_readname_trunc_1.err.fastq.gz tests/casava.1.8_readname_trunc_2.fastq.gz

must_fail 	time -p ./src/fastq_info   tests/casava.1.8_readname_trunc_2.fastq.gz tests/casava.1.8_readname_trunc_1.err.fastq.gz

must_fail 	time -p ./src/fastq_info  tests/casava.1.8_readname_trunc_1.err2.fastq.gz tests/casava.1.8_readname_trunc_2.fastq.gz

must_fail 	time -p ./src/fastq_info  tests/casava.1.8_readname_trunc_1.err.fastq.gz

must_succeed 	time -p ./src/fastq_info -s  tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_readname_trunc_2.fastq.gz
must_succeed 	time -p ./src/fastq_info  tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_2.fastq.gz
must_succeed 	time -p ./src/fastq_info  -r -s tests/casava.1.8_readname_trunc_1.fastq.gz tests/casava.1.8_2.fastq.gz

must_fail "./src/fastq_info --help"
##
echo "*** fastq_validator.sh"
export PATH=$PWD/src:$PATH
must_fail ./sh/fastq_validator.sh tests/c18_10000_1.fastq.gz.bz2 tests/c18_10000_2.fastq.gz.bz2
must_fail ./sh/fastq_validator.sh tests/c18_10000_1.fastq.gz tests/c18_10000_2.fastq.gz

must_fail ./sh/fastq_validator.sh tests/SRR3587500_1.fastq.gz.missing.bz2 
must_fail ./sh/fastq_validator.sh tests/a_1.fastq.err.bz2
must_fail ./sh/fastq_validator.sh tests/SRR3587500_1.fastq.gz.bz2 tests/SRR3587500_2.fastq.gz.bz2.
must_succeed ./sh/fastq_validator.sh tests/read-I1_si-ACCGAACA_lane-001-chunk-001.fastq.gz tests/read-I2_si-ACCGAACA_lane-001-chunk-001.fastq.gz
must_succeed ./sh/fastq_validator.sh tests/read-I1_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2 tests/read-I2_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2
must_succeed ./sh/fastq_validator.sh tests/read-I1_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2 tests/read-I2_si-ACCGAACA_lane-001-chunk-001.fastq.gz.bz2 

#gcov src/fastq_info


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
#gcov src/fastq_filterpair

##
echo "*** fastq_pre_barcodes"
must_succeed ./src/fastq_pre_barcodes --index1 tests/barcode_test_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --read1 tests/barcode_test_2.fastq.gz --outfile1 test.fastq.gz

must_succeed ./src/fastq_pre_barcodes --index1 tests/barcode_test_1.fastq.gz  --phred_encoding 33 --min_qual 10 --umi_read index1  --umi_offset 0 --umi_size 16 --read1_offset 0 --read1_size -1 --read1 tests/barcode_test_2.fastq.gz --outfile1 test.fastq.gz -X

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


## ==x2 --interleaved read1,index1
must_succeed "./src/fastq_pre_barcodes  --read1 tests/xxx_1.fastq.gz --index1 tests/xxx_2.fastq.gz  --outfile1 xni.gz  --umi_read index1  --umi_offset 0 --umi_size 16 && diff xni.gz x2.gz"

## ==x1
must_succeed "./src/fastq_pre_barcodes  --read1 tests/xxx_2.fastq.gz --index1 tests/xxx_1.fastq.gz  --outfile1 xni.gz  --umi_read index1  --umi_offset 0 --umi_size 16 && diff xni.gz x1.gz"

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

must_succeed fastq2bam -b test.bam -s 10xV1a -1 tests/10xv1a_R1.fastq.gz -2 tests/10xv1a_R3.fastq.gz -3 tests/10xv1a_R2.fastq.gz -4 tests/10xv1a_I1.fastq.gz

must_succeed fastq2bam -b test10.bam -s 10xV1a -1 tests/10xv1a_R1.fastq.gz -2 tests/10xv1a_R3.fastq.gz -3 tests/10xv1a_R2.fastq.gz -4 tests/10xv1a_I1.fastq.gz -X

must_succeed ./sh/fastq2bam -b test.bam -s 10xV1a -1 tests/10xv1a_R1.fastq.gz -2 tests/10xv1a_R3.fastq.gz -3 tests/10xv1a_R2.fastq.gz -4 tests/10xv1a_I1.fastq.gz -c 1 -C 2 -u 3 -U 4

must_succeed ./sh/fastq2bam -b test2.bam -s 10xV2 -1 tests/10xv1a_R1.fastq.gz -2 tests/10xv1a_R3.fastq.gz -3 tests/10xv1a_R2.fastq.gz -4 tests/10xv1a_I1.fastq.gz
must_succeed ./sh/fastq2bam -b test3.bam -s 10xV3 -1 tests/10xv1a_R1.fastq.gz -2 tests/10xv1a_R3.fastq.gz -3 tests/10xv1a_R2.fastq.gz -4 tests/10xv1a_I1.fastq.gz 
must_succeed  ./sh/fastq2bam -s drop-seq -1 tests/a_1.fastq.gz  -2 tests/a_2.fastq.gz  -b test.bam

must_fail diff <(samtools view test3.bam) <(samtools view test2.bam)

must_fail ./sh/fastq2bam -b test.bam -s 10xV1a -1 tests/10xv1a_R1.fastq.gz -2 tests/10xv1a_R3.fastq.gz -3 tests/10xv1a_R2.fastq.gz -4 tests/10xv1a_I1.fastq.gz -c 1 -C 2 -u 3 -U 
must_succeed ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b tmpf -3 tests/tx.I2.fastq.gz
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b tmpf -3 tests/tx.I2.fastq.gz -z 0 -Z 10
must_succeed ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b tmpf -3 tests/tx.I2.fastq.gz -z 0 -Z 5
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b tmpf2 -s 10
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b tmpf2 -S 10
must_fail ./sh/fastq2bam -s drop-seq -2 tests/10xv1a_R1.fastq.gz  -1 tests/10xv1a_R2.fastq.gz  -b test.bam

must_fail diff <(samtools view tmpf) <(samtools view tmpf2)

must_succeed diff <(samtools view tmpf) tests/ref1.sam

must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz -b 
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz  -2 tests/tx.I1.fastq.gz
must_fail ./sh/fastq2bam -s 10xV1i -1 tests/tx.RA.fastq.gz   -b tmpf2
must_fail ./sh/fastq2bam -s 10xV1i -2 tests/tx.I1.fastq.gz -b tmpf2
must_fail ./sh/fastq2bam -s 10xV1i -2 tests/tx.I1.fastq.gz -b 
rm -f tmpf tmpf2
#gcov src/fastq_pre_barcodes
echo "*** bam_add_tags"

must_succeed "./src/bam_add_tags --inbam tests/trans_small.bam --outbam tmp.bam"
must_succeed "./src/bam_add_tags --inbam tests/trans_small.bam --outbam tmp.10.bam --10x"
must_succeed "samtools view tmp.10.bam | grep -c 'UB:Z' > tmp.10 && samtools view tmp.bam | grep -c 'RX:Z' > tmp && diff tmp tmp.10"

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
must_fail "./src/bam_add_tags --inbam tests/trans_small.bam --outbam /tmp.bam --tx "

#gcov src/bam_add_tags



rm -f out_prefix_*.fastq.gz tmp.*.bam

#gcov src/fastq_split_interleaved

must_succeed ./src/fastq_tests
gcov src/fastq_tests
make -B -C src gcov

echo Failed tests: $num_failed
exit $num_failed

