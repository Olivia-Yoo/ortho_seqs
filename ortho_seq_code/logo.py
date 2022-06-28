import time
import sys
import os
import pandas as pd
import numpy as np
import ortho_seq_code.utils as utils
from ortho_seq_code.constants_orthoseqs import DNA_ALPHABETS, PROTEIN_ALPHABETS
import matplotlib.pyplot as plt
from ortho_seq_code.logger import Logger
import click

# constants - can decide where to put later
PROTEIN_DICT = {k: v for v, k in enumerate(PROTEIN_ALPHABETS)}
PROTEIN_DICT['n'] = 20
DNA_DICT = {k: v for v, k in enumerate(DNA_ALPHABETS)}
DNA_DICT['n'] = 4

def logo_plot(
        filename,
        molecule,
        alphbt_input,
        onefile,
):

    """Program to generate a logo plot from sequence data"""
    start_time = time.time()
    # out_dir = utils.create_dir_if_not_exists(out_dir)
    # sys.stdout = Logger(out_dir) # what does this do?
    global i # what does this do, is it necessary?
    print("") # new line for command line output

    # import sequence data
    (
        dm,
        sites,
        pop_size,
        seq,
        alphabets,
        custom_aa,
        exc
    ) = utils.get_seq_info(filename, alphbt_input, molecule, onefile)
    print(pop_size)
    print(sites)

    # TODO: generate position-weight matrix
    # first generate the position-frequency matrix (pure counts, will normalize later)
    n_options = 0
    if molecule == "protein":
        n_options = 21
        pfm = np.zeros((n_options, sites)) # for protein [amino acid (20 + n)][site of sequence]
        for i in range(pop_size):
            for j in range(sites):
                aa = seq[i][j]
                aa_idx = PROTEIN_DICT[aa]
                pfm[aa_idx][j] = pfm[aa_idx][j] + 1
    elif molecule == "dna":
        n_options = 5
        pfm = np.zeros((n_options, sites)) # for dna [base (4 + n)][site of sequence]
        print(pfm)
        for i in range(pop_size):
            for j in range(sites):
                base = seq[i][j]
                base_idx = DNA_DICT[base]
                pfm[base_idx][j] = pfm[base_idx][j] + 1
    else:
        print("Molecule type not supported.")
        sys.exit(1)
    print("PFM")
    print(pfm)

    # normalize pfm to generate position probability matrix (ppm)
    column_sums = np.sum(pfm, axis=0)
    print("column sums")
    print(column_sums)

    for i in range(sites):
        ppm = np.copy(pfm)
        ppm = np.divide(ppm, column_sums)

    print("PPM")
    print(ppm)

    # TODO: calculate information from PWM

    # TODO: construct logo plot (bits-wise and probability-wise)



# ACTUALLY IMPLEMENTING/USING THIS PART!
logo_plot("/Users/olivia.yoo/Desktop/code/ortho_seqs_work/ortho_seqs/ortho_seq_code/tests/data/nucleotide/first_order/test_seqs_2sites_dna.txt", "dna", None, True)
