#!/usr/bin/env python3
# Author: Francois Aguet

import argparse
import os
import struct
import subprocess
from datetime import datetime
import contextlib
import sys

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)


parser = argparse.ArgumentParser(description='Convert BAM to FASTQ using SamToFastq from Picard.')
parser.add_argument('bam_file', type=str, help='BAM file')
parser.add_argument('-p', '--prefix', type=str, default='Reads', help='Prefix for output files; usually <sample_id>')
parser.add_argument('-o', '--output_dir', default=os.getcwd(), help='Directory to which FASTQs will be written')
parser.add_argument('-m', '--memory', default='8', type=str, help='Memory, in GB')
parser.add_argument('--java', default='/opt/apps/java/1.8.0_31/bin/java', help='Path to java')
parser.add_argument('--jar', default='/opt/apps/picard-tools/2.9.4/picard.jar', help='Path to Picard jar')
parser.add_argument('--include_non_pf_reads', type=str.lower, choices=['true', 'false'], default='true', help='Sets INCLUDE_NON_PF_READS option (PF: passed filtering). SamToFastq default: false')
parser.add_argument('--include_non_primary_alignments', type=str.lower, choices=['true', 'false'], default='false', help='Sets INCLUDE_NON_PRIMARY_ALIGNMENTS option. SamToFastq default: false')
args = parser.parse_args()

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Starting SamToFastq', flush=True, file=sys.stderr)

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

with cd(args.output_dir):
    fastq1 = args.prefix+'_1.fastq.gz'
    fastq2 = args.prefix+'_2.fastq.gz'
    fastq0 = args.prefix+'_unpaired.fastq.gz'


    try:
        subprocess.check_call(args.java + ' -jar -Xmx'+args.memory+'g '+args.jar+' SamToFastq INPUT='+args.bam_file\
        +' INCLUDE_NON_PF_READS='+args.include_non_pf_reads\
        +' INCLUDE_NON_PRIMARY_ALIGNMENTS='+args.include_non_primary_alignments\
        +' VALIDATION_STRINGENCY=SILENT FASTQ=' + fastq1 + ' SECOND_END_FASTQ=' + fastq2 + ' UNPAIRED_FASTQ=' + fastq0, shell=True)
    except subprocess.CalledProcessError as e:
        print("ERROR: command failed for following reason: %s" % (e), file=sys.stderr)
        exit(1)


    # Delete unpaired reads FASTQ if empty
    with open(fastq0, 'rb') as f0:
        f0.seek(-4,2)
        if struct.unpack('<I', f0.read(4))[0]==0:  # empty file
            os.remove(fastq0)

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finished SamToFastq', flush=True, file=sys.stderr)
