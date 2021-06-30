#!/usr/bin/env python3

from Bio import SeqIO
import csv
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import shlex
import argparse

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""

def make_qc_plot(depth_pos, n_density, samplename, window=200):
    depth_df = pd.DataFrame( { 'position' : [pos[1] for pos in depth_pos], 'depth' : [dep[2] for dep in depth_pos] } )
    depth_df['depth_moving_average'] = depth_df.iloc[:,1].rolling(window=window).mean()

    n_df = pd.DataFrame( { 'position' : [pos[0] for pos in n_density], 'n_density' : [dens[1] for dens in n_density] } )

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel('Position')

    ax1.set_ylabel('Depth', color = 'g')
    ax1.set_ylim(top=10**5, bottom=1)
    ax1.set_yscale('log')
    ax1.plot(depth_df['depth_moving_average'], color = 'g')

    ax2.set_ylabel('N density', color = 'r')  
    ax2.plot(n_df['n_density'], color = 'r')
    ax2.set_ylim(top=1)

    plt.title(samplename)
    plt.savefig(samplename + '.depth.png')

def read_depth_file(bamfile):
    p = subprocess.Popen(['samtools', 'depth', '-a', '-d', '0', bamfile],
                       stdout=subprocess.PIPE)
    out, err = p.communicate()
    counter = 0

    pos_depth = []
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          pos_depth.append(ln.split("\t"))
    
    return pos_depth


def get_covered_pos(pos_depth, min_depth):
    counter = 0
    for contig, pos,depth in pos_depth:
        if int(depth) >= min_depth:
            counter = counter + 1
    
    return counter

def get_N_positions(fasta):
    n_pos =  [i for i, letter in enumerate(fasta.seq.lower()) if letter == 'n']

    return n_pos

def get_pct_N_bases(fasta):
    
    count_N = len(get_N_positions(fasta))

    pct_N_bases = count_N / len(fasta.seq) * 100

    return pct_N_bases

def get_largest_N_gap(fasta):
    n_pos = get_N_positions(fasta)

    n_pos = [0] + n_pos + [len(fasta.seq)]

    n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

    return sorted(n_gaps)[-1]

def get_ref_length(ref):
    record = SeqIO.read(ref, "fasta")
    return len(record.seq)


def sliding_window_N_density(sequence, window=10):

    sliding_window_n_density = []
    for i in range(0, len(sequence.seq), 1):
        window_mid = i + ( window / 2)
        window_seq = sequence.seq[i:i+window]
        n_count = window_seq.lower().count('n')
        n_density = n_count / window

        sliding_window_n_density.append( [ window_mid, n_density ] )

    return sliding_window_n_density

def get_num_reads(bamfile):

    st_filter = '0x900'
    command = 'samtools view -c -F{} {}'.format(st_filter, bamfile)
    what = shlex.split(command)

    return subprocess.check_output(what).decode().strip()

def assess(fasta_file, bam_file=None, ref_length=None, depth=None):

    # Unknown base calcs
    fasta = SeqIO.read(fasta_file, "fasta")

    pct_N_bases   = 0
    largest_N_gap = 0
    qc_pass       = "FALSE"

    if len(fasta.seq) != 0:

        pct_N_bases = get_pct_N_bases(fasta)
        largest_N_gap = get_largest_N_gap(fasta)

        # QC PASS / FAIL
        if largest_N_gap >= 10000 or pct_N_bases < 50.0:
            qc_pass = "TRUE"

    N_density = sliding_window_N_density(fasta)

    # The order of keys is important
    pairs = [('pct_N_bases', "{:.2f}".format(pct_N_bases)),
             ('longest_no_N_run', largest_N_gap),
             ('fasta', fasta_file),
             ('qc_pass', qc_pass)]
    
    depth_pos = None
    if bam_file != None:
       depth_pos = read_depth_file(bam_file)
       depth_covered_bases = get_covered_pos(depth_pos, depth)
       pct_covered_bases = depth_covered_bases / ref_length * 100
       # Number of aligned reads calculaton
       num_reads = get_num_reads(bam_file)

       pairs.insert(1,
           ('pct_covered_bases', "{:.2f}".format(pct_covered_bases)))
       pairs.insert(3, ('num_aligned_reads', num_reads)) 
       pairs.insert(5, ('bam', bam_file))
    else:
       # Remap key names
       pairs = [ (pair[0]+'_amd', pair[1]) for pair in pairs]
        
    return (dict(pairs), N_density, depth_pos) 
    
def go(args):
    if args.illumina:
        depth = 10
    elif args.nanopore:
        depth = 20

    ## Depth calcs
    ref_length = get_ref_length(args.ref)
    
    ## Get QC values for a pair of bam-fasta files
    (qc_values, N_density, depth_pos) = assess(
        args.fasta, args.bam, ref_length, depth);
    ## Get the keys in the order they were inserted
    column_names = list(qc_values)
    if args.ivar_md != None:
        qc_values['ivar_md'] = args.ivar_md
        column_names.insert(-1, 'ivar_md')  
    ## Prepend sample name column
    column_names.insert(0, 'sample_name')
    qc_values['sample_name'] = args.sample
    qc_line = qc_values

    ## If appropriate, get QC values for another pair of bam-fasta files
    if args.fasta_amd:
        (qc_values_amd, tmp1, tmp2) = assess(args.fasta_amd);
        if args.ivar_md != None:
            qc_values_amd['ivar_amd'] = args.ivar_amd
        ## Combine two dictionaries
        qc_line = {**qc_values, **qc_values_amd};
        ## Set correct order for the list of column names
        qc_pass_column = column_names.pop()
        column_names.extend(list(qc_values_amd))
        ## Reinstall qc pass columns as the last column
        column_names.append(qc_pass_column)

    ## Write all QC values to a CSV file
    with open(args.outfile, 'w') as csvfile:
        header = column_names
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerow(qc_line)

    make_qc_plot(depth_pos, N_density, args.sample)

def main():

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--nanopore', action='store_true',
        help='''A boolean flag describing the sequencing platform used.''')
    group.add_argument('--illumina', action='store_true',
        help='''A boolean flag describing the sequencing platform used.''')
    parser.add_argument('--outfile', required=True,
        help='''The path of the output QC summary file''')
    parser.add_argument('--sample', required=True, help='Sample name.')
    parser.add_argument('--ref', required=True,
        help='''The path of the reference FASTA file.''')
    parser.add_argument('--bam', required=True,
        help='''The path of the aligned and filtered BAM file.''')
    parser.add_argument('--fasta', required=True,
        help='''The path of a consensus fasta file produced by ivar using the
                minimum depth given by the --ivar_amd argument, required.''')
    parser.add_argument('--fasta_amd', required=False, default=None,
        help='''The path of a consensus fasta file produced by ivar using the
                minimum depth given by the --ivar_amd argument, optional.''')
    parser.add_argument('--ivar_md', required=False, default=None,
        help='''Minimum depth value used for ivar when generating the consensus
                file given by the --fasta argument, optional.''')
    parser.add_argument('--ivar_amd', required=False, default=None,
        help='''Minimum depth value used for ivar when generating the consensus
                file given by the --fasta_amd argument, optional.''') 

    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()
