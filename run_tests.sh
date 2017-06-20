#!/bin/bash

let num_failed=0
function must_fail {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -eq $? ];  then	
	STATUS=FAILED
	num_failed+=1
    fi
    echo $STATUS $cmd
}

function must_succeed {
    cmd=$*
    STATUS=PASSED
    bash -c "$cmd" 2> /dev/null
    if [ 0 -ne $? ];  then	
	STATUS=FAILED
	num_failed+=1
    fi
    echo $STATUS $cmd
}
#############################################
##
echo "*** fastq_filter_n"
must_succeed "./bin/fastq_filter_n tests/test_21_2.fastq.gz > tmp && diff /dev/null tmp"
must_fail "./bin/fastq_filter_n -n 100 tests/test_21_2.fastq.gz > tmp && diff /dev/null tmp"
must_fail "./bin/fastq_filter_n tests/test_21_2.fastq.gz > tmp && diff  tests/test_21_2.fastq.gz tmp"
must_succeed "./bin/fastq_filter_n tests/test_1.fastq.gz > tmp && diff <(zcat tests/test_1.fastq.gz) tmp"
##
echo "*** fastq_num_reads"
must_succeed "[ `./bin/fastq_num_reads tests/test_21_2.fastq.gz` -eq 2 ]"
must_fail "[ `./bin/fastq_num_reads tests/test_21_2.fastq.gz` -ne 2 ]"
must_succeed "[ `./bin/fastq_num_reads tests/c18_10000_1.fastq.gz` -eq 10000 ]"
must_fail "[ `./bin/fastq_num_reads tests/c18_10000_1.fastq.gz` -ne 10000 ]"
must_succeed "[ `./bin/fastq_num_reads tests/one.fastq.gz` -eq 1 ]"
##
echo "*** fastq_truncate"
must_succeed "[ `./bin/fastq_truncate tests/test_21_2.fastq.gz 1|wc -l` -eq 4 ]"
must_succeed "[ `./bin/fastq_truncate tests/test_21_2.fastq.gz 0|wc -l` -eq 0 ]"
must_succeed "[ `./bin/fastq_truncate tests/test_21_2.fastq.gz 2|wc -l` -eq 8 ]"
must_fail "./bin/fastq_truncate tests/test_21_2.fastq.gz"


echo "*** fastq_info"
must_fail ./bin/fastq_info tests/test_e1.fastq.gz 
must_fail ./bin/fastq_info tests/test_e2.fastq.gz
must_fail ./bin/fastq_info tests/test_e3.fastq.gz 
must_fail ./bin/fastq_info tests/test_e4.fastq.gz 
must_fail ./bin/fastq_info tests/test_e5.fastq.gz 
must_fail ./bin/fastq_info tests/test_e6.fastq.gz 
must_fail ./bin/fastq_info tests/test_e7.fastq.gz 
must_fail ./bin/fastq_info tests/test_e8.fastq.gz 
must_fail ./bin/fastq_info tests/test_e9.fastq.gz 
must_fail ./bin/fastq_info tests/test_e10.fastq.gz 
must_fail ./bin/fastq_info tests/test_e14.fastq.gz 
must_fail ./bin/fastq_info tests/test_e15.fastq.gz 
must_fail ./bin/fastq_info tests/test_e16.fastq.gz 
must_fail ./bin/fastq_info tests/test_e17.fastq.gz 
must_fail ./bin/fastq_info -f tests/test_dot.fastq.gz 
must_fail ./bin/fastq_info tests/c18_1M_1.fastq.gz tests/c18_1M_2.fastq.gz 
##
## just checks the exit status
must_succeed ./bin/fastq_info tests/test_dot.fastq.gz 
must_succeed 	./bin/fastq_info tests/test_1.fastq.gz 
must_succeed 	./bin/fastq_info tests/test_2.fastq.gz 
must_succeed 	./bin/fastq_info tests/test_13.fastq.gz
must_succeed 	./bin/fastq_info tests/test_17.fastq.gz
must_succeed 	./bin/fastq_info tests/test_21_1.fastq.gz
must_succeed 	./bin/fastq_info tests/test_21_1.fastq.gz tests/test_21_2.fastq.gz 
must_succeed 	time -p ./bin/fastq_info tests/pe_bug14.fastq.gz tests/pe_bug14.fastq.gz 
must_succeed 	time -p ./bin/fastq_info tests/c18_1M_1.fastq.gz 
must_succeed 	time -p ./bin/fastq_info tests/c18_1M_2.fastq.gz 
must_succeed 	time -p ./bin/fastq_info tests/c18_1M_1.fastq.gz 
must_succeed 	time -p ./bin/fastq_info tests/c18_1M_2.fastq.gz 
must_succeed 	time -p ./bin/fastq_info tests/c18_1M_1.fastq.gz tests/c18_1M_1.fastq.gz 
must_succeed 	time -p ./bin/fastq_info tests/casava.1.8i.fastq.gz pe 
must_succeed 	time -p ./bin/fastq_info tests/solexa_1.fastq.gz tests/solexa_2.fastq.gz 
##
echo "*** fastq_filterpair"
must_succeed "./bin/fastq_filterpair tests/test_2.fastq.gz tests/test_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz && diff <(zcat f1.fastq.gz) <(zcat tests/test_2.fastq.gz)"
must_succeed "./bin/fastq_filterpair tests/a_1.fastq.gz tests/a_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz && diff <(zcat f2.fastq.gz) <(zcat tests/a_2.fastq.gz) && diff <(zcat f1.fastq.gz) <(zcat tests/a_1.fastq.gz)"
must_succeed "./bin/fastq_filterpair tests/casava.1.8_2.fastq.gz tests/casava.1.8_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed "./bin/fastq_filterpair tests/casava.1.8_1.fastq.gz tests/casava.1.8_1.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz"
must_succeed ./bin/fastq_filterpair tests/c18_10000_1.fastq.gz tests/c18_10000_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz

#must_succeed ./bin/fastq_filterpair tests/c18_1M_2.fastq.gz tests/c18_1M_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz 
#must_succeed ./bin/fastq_filterpair tests/c18_1M_1.fastq.gz tests/c18_1M_1.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz 

#echo "This may take a while..."
#must_succeed ./bin/fastq_filterpair tests/c18_1M_1.fastq.gz tests/c18_1M_2.fastq.gz  f1.fastq.gz f2.fastq.gz up.fastq.gz


exit $num_failed

