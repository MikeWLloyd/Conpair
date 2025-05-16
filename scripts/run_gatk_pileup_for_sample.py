#!/usr/bin/env python3

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)

import os
import sys
import argparse
from shutil import move
from tempfile import NamedTemporaryFile

desc = "Program to run GATK Pileup on a single sample"
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('-B', '--bam', help='BAMFILE [mandatory field]', required=True)
parser.add_argument('-O', '--outfile', help='OUTPUT FILE (PILEUP) [mandatory field]', required=True)
parser.add_argument('-D', '--conpair_dir', help='CONPAIR DIR [$CONPAIR_DIR by default]')
parser.add_argument('-R', '--reference', help='REFERENCE GENOME [GRCh37 by default]')
parser.add_argument('-M', '--markers', help='MARKER FILE [GRCh37-default]')
parser.add_argument('-G', '--gatk', help='GATK JAR [$GATK by default]')
parser.add_argument('-J', '--java', help='PATH to JAVA [java by default]', default='java')
parser.add_argument('-t', '--temp_dir_java', help='temporary directory to set -Djava.io.tmpdir')
parser.add_argument('-m', '--xmx_java', help='Xmx java memory setting [default: 12g]', default='12g')
parser.add_argument(
    '--remove_chr_prefix',
    help='REMOVE CHR PREFIX FROM THE CHROMOSOME COLUMN IN THE OUTPUT FILE [false by default]',
    action='store_true'
)

args = parser.parse_args()

if not os.path.exists(args.bam):
    print(f'ERROR: Specified bamfile {args.bam} cannot be found.')
    sys.exit(1)

CONPAIR_DIR = args.conpair_dir if args.conpair_dir else os.environ.get('CONPAIR_DIR')
if not CONPAIR_DIR:
    print('ERROR: CONPAIR_DIR not specified and not set in environment.')
    sys.exit(1)

GATK = args.gatk if args.gatk else os.environ.get('GATK_JAR')
if not GATK or not os.path.exists(GATK):
    print(f'ERROR: GATK jar {GATK} cannot be found.')
    sys.exit(2)

MARKER_FILE = (
    args.markers
    if args.markers
    else os.path.join(
        CONPAIR_DIR,
        'data',
        'markers',
        'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed'
    )
)
if not os.path.exists(MARKER_FILE):
    print(f'ERROR: Marker file {MARKER_FILE} cannot be found.')
    sys.exit(2)

REFERENCE = args.reference if args.reference else os.path.join(CONPAIR_DIR, 'data', 'genomes', 'human_g1k_v37.fa')
if not os.path.exists(REFERENCE):
    print(f'ERROR: Reference genome {REFERENCE} cannot be found.')
    sys.exit(3)

if args.temp_dir_java:
    JAVA_TEMP = f"-Djava.io.tmpdir={args.temp_dir_java}"
    if not os.path.isdir(args.temp_dir_java):
        os.makedirs(args.temp_dir_java, 0o770)
else:
    JAVA_TEMP = ""

command_line = (
    f"{args.java} {JAVA_TEMP} -Xmx{args.xmx_java} -jar {GATK} -T Pileup "
    f"-R {REFERENCE} -I {args.bam} -L {MARKER_FILE} -o {args.outfile} "
    "-verbose -rf DuplicateRead --filter_reads_with_N_cigar "
    "--filter_mismatching_base_and_quals"
    "--fix_misencoded_quality_scores -fixMisencodedQuals"
)

os.system(command_line)

if args.remove_chr_prefix:
    print("Removing 'chr' prefix...")

    with NamedTemporaryFile('w', delete=False, encoding='utf-8') as tmp_source:
        with open(args.outfile, 'r', encoding='utf-8') as source_file:
            for line in source_file:
                if line.startswith("chr"):
                    tmp_source.write(line[3:])
                else:
                    tmp_source.write(line)

    move(tmp_source.name, args.outfile)
