#!/usr/bin/env python3

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)

import sys
import os
import optparse
import importlib.util
from collections import defaultdict
import numpy as np

HOMOZYGOUS_P_VALUE_THRESHOLD = 0.999

desc = """Program to estimate tumor-normal sample contamination"""
parser = optparse.OptionParser(version='%prog version 1.0 March/01/2016', description=desc)
parser.add_option('-T', '--tumor_pileup', help='TUMOR PILEUP FILE [mandatory field]', type='string', action='store')
parser.add_option('-N', '--normal_pileup', help='NORMAL PILEUP FILE [mandatory field]', type='string', action='store')
parser.add_option('-D', '--conpair_dir', help='CONPAIR DIR [default: $CONPAIR_DIR]', action='store')
parser.add_option('-M', '--markers', help='MARKER FILE [default: markers for GRCh37 from $CONPAIR_DIR/data/markers/ ]', type='string', action='store')
parser.add_option('-O', '--outfile', help='TXT OUTPUT FILE [default: stdout]', default="-", type='string', action='store')
parser.add_option('-G', '--grid', help='GRID INTERVAL [default: 0.01]', type='float', default=0.01, action='store')
parser.add_option('-Q', '--min_mapping_quality', help='MIN MAPPING QUALITY [default: 10]', default=10, type='int', action='store')

(opts, args) = parser.parse_args()

if opts.conpair_dir:
    CONPAIR_DIR = opts.conpair_dir
else:
    CONPAIR_DIR = os.environ.get('CONPAIR_DIR')
    if not CONPAIR_DIR:
        print("ERROR: CONPAIR_DIR environment variable not set and no -D option provided.")
        sys.exit(1)


def load_module(module_name, file_path):
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ContaminationModel = load_module('ContaminationModel', os.path.join(CONPAIR_DIR, 'modules', 'ContaminationModel.py'))
ContaminationMarker = load_module('ContaminationMarker', os.path.join(CONPAIR_DIR, 'modules', 'ContaminationMarker.py'))
MathOperations = load_module('MathOperations', os.path.join(CONPAIR_DIR, 'modules', 'MathOperations.py'))
Genotypes = load_module('Genotypes', os.path.join(CONPAIR_DIR, 'modules', 'Genotypes.py'))

if not opts.tumor_pileup or not opts.normal_pileup:
    parser.print_help()
    sys.exit(1)

if not os.path.exists(opts.tumor_pileup):
    print(f'ERROR: Input tumor file {opts.tumor_pileup} cannot be found.')
    sys.exit(1)

if not os.path.exists(opts.normal_pileup):
    print(f'ERROR: Input normal file {opts.normal_pileup} cannot be found.')
    sys.exit(1)

if opts.markers:
    MARKER_FILE = opts.markers
else:
    MARKER_FILE = os.path.join(
        CONPAIR_DIR,
        'data',
        'markers',
        'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt'
    )

if not os.path.exists(MARKER_FILE):
    print(f'ERROR: Marker file {MARKER_FILE} cannot be found.')
    sys.exit(2)

grid_precision = opts.grid
MMQ = opts.min_mapping_quality


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


Markers = ContaminationMarker.get_markers(MARKER_FILE)

Normal_homozygous_genotype = defaultdict(lambda: defaultdict())

checkpoints = list(drange(0.0, 1.0, grid_precision))
checkpoints.append(0.5)
Scores = ContaminationModel.create_conditional_likelihood_of_base_dict(checkpoints)
checkpoints = list(drange(0.0, 0.5, grid_precision))
checkpoints.append(0.5)

if opts.outfile != "-":
    outfile = open(opts.outfile, 'w')
else:
    outfile = sys.stdout

# PARSING THE NORMAL PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

with open(opts.normal_pileup) as file:
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = ContaminationMarker.parse_mpileup_line(line, min_map_quality=MMQ)
        try:
            marker = Markers[f"{pileup.chrom}:{pileup.pos}"]
        except KeyError:
            continue

        if not pileup.Quals[marker.ref] and not pileup.Quals[marker.alt]:
            continue

        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]

        # Compute genotype likelihoods for the current marker and pileup
        AA_likelihood, AB_likelihood, BB_likelihood = Genotypes.compute_genotype_likelihood(
            pileup.Quals[marker.ref],
            pileup.Quals[marker.alt],
            normalize=True
        )

        if AA_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {
                'genotype': marker.ref,
                'AA_likelihood': AA_likelihood,
                'AB_likelihood': AB_likelihood,
                'BB_likelihood': BB_likelihood
            }
        elif BB_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {
                'genotype': marker.alt,
                'AA_likelihood': AA_likelihood,
                'AB_likelihood': AB_likelihood,
                'BB_likelihood': BB_likelihood
            }

        p_AA, p_AB, p_BB = Genotypes.RAF2genotypeProb(RAF)
        lPAA = MathOperations.log10p(p_AA)
        lPAB = MathOperations.log10p(p_AB)
        lPBB = MathOperations.log10p(p_BB)

        priors = [lPAA*2, lPAA+lPBB, lPAA+lPAB, lPAB*2, lPAB+lPAA, lPAB+lPBB, lPBB*2, lPBB+lPAA, lPBB+lPAB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)

D = ContaminationModel.calculate_contamination_likelihood(checkpoints, Data, Scores)
ARGMAX = np.argmax(D)
cont = checkpoints[ARGMAX]

x1 = max(cont-grid_precision, 0.0)
x2 = cont
x3 = min(cont+grid_precision, 1.0)

if x2 == 0.0:
    x2 += grid_precision/100
elif x2 == 1.0:
    x2 -= grid_precision/100

# SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

optimal_val = ContaminationModel.apply_brents_algorithm(Data, Scores, x1, x2, x3)

# PRINTING THE NORMAL RESULTS

outfile.write(f"Normal sample contamination level: {round(100.0*optimal_val, 3)}%\n")


# PARSING THE TUMOR PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

with open(opts.tumor_pileup) as file:
    checkpoints = list(drange(0.0, 1.0, grid_precision))
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = ContaminationMarker.parse_mpileup_line(line, min_map_quality=MMQ)

        try:
            normal_hom_genotype = Normal_homozygous_genotype[pileup.chrom][pileup.pos]['genotype']
        except KeyError:
            continue

        try:
            marker = Markers[f"{pileup.chrom}:{pileup.pos}"]
        except KeyError:
            continue

        if not pileup.Quals[marker.ref] and not pileup.Quals[marker.alt]:
            continue

        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]

        Normal_info = Normal_homozygous_genotype[pileup.chrom][pileup.pos]
        AA_likelihood = Normal_info['AA_likelihood']
        AB_likelihood = Normal_info['AB_likelihood']
        BB_likelihood = Normal_info['BB_likelihood']
        nlPAA = MathOperations.log10p(AA_likelihood)
        nlPAB = MathOperations.log10p(AB_likelihood)
        nlPBB = MathOperations.log10p(BB_likelihood)

        p_AA, p_AB, p_BB = Genotypes.RAF2genotypeProb(RAF)
        lPAA = MathOperations.log10p(p_AA)
        lPAB = MathOperations.log10p(p_AB)
        lPBB = MathOperations.log10p(p_BB)
        priors = [lPAA+nlPAA, lPBB+nlPAA, lPAB+nlPAA, lPAB+nlPAB, lPAA+nlPAB, lPBB+nlPAB, lPBB+nlPBB, lPAA+nlPBB, lPAB+nlPBB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)

D = ContaminationModel.calculate_contamination_likelihood(checkpoints, Data, Scores)
ARGMAX = np.argmax(D)
cont = checkpoints[ARGMAX]

x1 = max(cont-grid_precision, 0.0)
x2 = cont
x3 = min(cont+grid_precision, 1.0)

if x2 == 0.0:
    x2 += grid_precision/100
elif x2 == 1.0:
    x2 -= grid_precision/100

# SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

optimal_val = ContaminationModel.apply_brents_algorithm(Data, Scores, x1, x2, x3)

# PRINTING THE TUMOR RESULTS

outfile.write(f"Tumor sample contamination level: {round(100.0 * optimal_val, 3)}%\n")

if opts.outfile != "-":
    outfile.close()
