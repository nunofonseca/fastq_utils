#!/bin/bash
# =========================================================
# Copyright 2012-2017,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
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

FILES=$*

function file_extension {
    echo $*|sed -E "s/([^ \.]+)\.//g" 
}


# check extension
ext=`file_extension $1`
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

set -e
if [ "$ext-" == "bam-" ]; then
    #hmm, this now validates bams...kind of
    # samtools version should be 1 or above
    f=$1    
    echo "BAM file"
    # check if the BAM contains unaligned reads
    echo "Checking for unmapped reads"
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
	if [ "-$ext" == "-bz2" ] || [ "-$ext" == "-bzip2" ] ; then
	    echo BZIP file
	    named_pipe=.`basename .$f`.pipe.fastq
	    mkfifo $named_pipe
	    bunzip2 -k  -c $f > $named_pipe  &
	    FILES2PROCESS="$FILES2PROCESS $named_pipe"
	    FILES2DELETE="$FILES2DELETE $named_pipe"
	else
	    FILES2PROCESS="$FILES2PROCESS $f"
	fi
    done
fi
fastq_info $FILES2PROCESS

if [ "-$FILES2DELETE" != "-" ]; then
    #echo -n "Removing named pipes..."
    rm -f $FILES2DELETE
    #echo "done."
fi
exit 0
