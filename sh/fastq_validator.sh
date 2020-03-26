#!/bin/bash
# =========================================================
# Copyright 2012-2020,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# if not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
# wrapper to fastq_info to support BAM and bzip2

FILES=$1
PE_PARAMETER=
if [ "$2-" == "pe-" ]; then
    PE_PARAMETER=pe
else
    FILES=$*
fi

IS_PE=${#FILES[@]}

if [ "$1-" == "-" ]; then
    echo "ERROR: fastq_validator.sh file1 [file2|pe]" 
    exit 1
fi

function file_type {
    ## first check if is a bam/cram file
    samtools quickcheck $1 &> /dev/null
    if [ $? -eq 0 ]; then
	echo "bam"
    else
	## 
	x=$(file -b -i  $1|cut -f 1 -d\;|sed "s|.*/||")
	## octet-stream included to handle the bug in the file command
	## where gzip files are reported as Minix filesystem, V2, 30 char names...
	## https://bugzilla.redhat.com/show_bug.cgi?id=1014998
	if [ "$x" == "x-gzip" ] || [ "$x" == "octet-stream" ] ; then echo "gz";	   
	else
	    if [ "$x" == "x-bzip2" ]; then echo "bzip2";
	    else
		if [ "$x" == "plain" ]; then echo "fastq";
		else
		    x=$(file -b  $1| grep -c CRAM)
		    if [ $x -eq 1 ]; then echo "cram"
		    else
			echo "Unsupported file type $x" > /dev/stderr
			exit 4
		    fi		    
		fi
	    fi
	fi	
    fi
}

function file_extension {
    filename=$(basename $*)
    ext="${filename##*.}"
    if [ "$ext-" != "$filename-" ]; then
	echo $ext
    fi
}

# check extension
ext=`file_extension $1`

if [ "-$ext" == "-" ]; then
    ext=$(file_type $1)
    echo "File does not have an extension, assuming that it is '.$ext'"
fi

# Check integrity of gzip files
if [ "$ext-" == "gz-" ]; then
    for f in $FILES; do
	echo -n "Checking integrity of gzip file $f..."
	gzip -t $f
	if [ $? -eq 0 ]; then
	    echo "done."
	else
	    echo ""
	    echo "ERROR: Error in file $f: corrupted gzip file"
	    exit 1
	fi
	    
    done
    echo ""
fi

set -eT
#x
if [ "$ext-" == "bam-" ] || [ "$ext-" == "cram-" ]; then
    #hmm, this now validates bams...kind of
    # samtools version should be 1 or above
    f=$1    
    echo "BAM/CRAM file ($ext)"
    # check if the BAM contains unaligned reads
    echo "Checking for unmapped reads"
    UN=`samtools view -c -F 4 $f`
    if [ "$UN-" != "0-" ]; then
	echo "ERROR: Aligned reads found in $f." > /dev/stderr
	exit 1
    fi

    UN=`samtools view -c -f 4 $f`
    if [ "$UN-" == "0-" ]; then
	echo "ERROR: No unaligned reads found in $f." > /dev/stderr
	exit 1
    fi

    named_pipe=.`basename .$f`.pipe.fastq
    mkfifo $named_pipe
    echo "Converting BAM to fastq"
    samtools bam2fq $f > $named_pipe &
    
    FILES2PROCESS=$named_pipe
    FILES2DELETE=$named_pipe
else
    FILES2PROCESS=
    FILES2DELETE=
    for f in $FILES; do
	ext=`file_extension $f`
	if [ "-$ext" == "-" ]; then
	    ext=$(file_type $f)
	    echo "File $f does not have an extension, assuming that it is '.$ext'"
	fi
	if [ "-$ext" == "-bz2" ] || [ "-$ext" == "-bzip2" ] ; then
	    echo BZIP file
	    # check integrity
	    set +e
	    echo "Checking integrity of $f..."
	    set -o pipefail
	    tmp_file=$(mktemp  --suffix `basename .$f`.tmp.gz -p .)
	    rm -f $named_pipe
	    echo "Creating a temporary gzip version of $f as $tmp_file..."
	    bunzip2 -c $f | gzip -c > $tmp_file
	    if [ $? -ne 0 ]; then
		echo "ERROR: $f: error uncompressing bzip2 file"
		exit 2
	    fi
	    echo "Creating a temporary gzip version of $f...done."
	    echo "Checking integrity of $f...complete."
	    FILES2PROCESS="$FILES2PROCESS $tmp_file"
	    FILES2DELETE="$FILES2DELETE $tmp_file"
	else
	    FILES2PROCESS="$FILES2PROCESS $f"
	fi
    done
fi
## validating
failed=0

if [ $(echo $FILES2PROCESS|wc -w) -gt 1 ]; then
    echo "Checking each fastq file independently..."
    for f in $FILES2PROCESS; do
	echo "Checking $f..."
	set +e
	fastq_info $f
	estatus=$?
	let failed=(failed*10)+$estatus
	set -eT
	echo "Checking $f ($estatus)...done."
    done
    ##
    if [ $failed -eq 0 ]; then
	## checking both files
	## ugggly code :(
	PREV_EXT=
	for f in $FILES; do
	    ext=`file_extension $f`
	    if [ "-$ext" == "-" ]; then
		ext=$(file_type $f)
		echo "File $f does not have an extension, assuming that it is '.$ext'"
	    fi
	    if [ "-$PREV_EXT" == "-" ]; then
		PREV_EXT=$ext
	    fi
	    if [ "-$PREV_EXT" != "-$ext" ]; then
		echo "ERROR: File types differ $ext vs $PREV_EXT" > /dev/stderr
		exit 2
	    fi
	done
	echo "Checking $FILES2PROCESS"
	set +e
	fastq_info $FILES2PROCESS
	estatus=$?
	let failed=$estatus
	set -eT
    fi
else
    ## single file
    echo "Checking $FILES2PROCESS"
    set +e
    fastq_info $FILES2PROCESS $PE_PARAMETER
    estatus=$?
    let failed=$estatus
    set -eT
fi

if [ "-$FILES2DELETE" != "-" ]; then
    #echo -n "Removing named pipes...$FILES2DELETE"
    rm -f $FILES2DELETE
    #echo "done."
fi
exit $failed
