####bin/mod/m5c_res_process.py
import os
import argparse
import subprocess
from Bio import SeqIO

def modify_dampened_fraction(plus_file,minus_file,tr_fa,out_file):
    transcript_seq = {}
    with open(tr_fa,'r') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            transcript_id = record.id
            sequence = str(record.seq)
            transcript_seq[transcript_id] = sequence
    transcript_site = {}
    with open(plus_file,'r') as plus:
        current_tr = None
        current_site = []
        for line in plus:
            if line.startswith('track'):
                continue
            if line.startswith('variableStep chrom='):
                if current_tr is not None:
                    transcript_site[current_tr] = current_site
                current_tr = line.strip().split(' ')[1].split('chrom=')[1] + '_F' 
                current_site = []
            elif line.split():
                position, fraction = line.split()
                position = int(position)
                base = transcript_seq[current_tr.strip('_F')][position-1]
                current_site.append((position, base,fraction))
    if current_tr is not None:
        transcript_site[current_tr] = current_site
    
    with open(minus_file,'r') as minus:
        current_tr = None
        current_site = []
        for line in minus:
            if line.startswith('track'):
                continue
            if line.startswith('variableStep'):
                if current_tr is not None:
                    transcript_site[current_tr] = current_site
                current_tr = line.strip().split(' ')[1].split('chrom=')[1] + '_R' 
                current_site = []
            elif line.split():
                position, fraction = line.split()
                position = int(position)
                base = transcript_seq[current_tr.strip('_R')][position-1]
                current_site.append((position, base,fraction))
    if current_tr is not None:
        transcript_site[current_tr] = current_site
    
    with open(out_file,'w') as out:
        out.write(f'transcript\tposition\tbase\tfraction\n')
        for transcript,info in transcript_site.items():
            for position, base,fraction in info:
                out.write(f'{transcript}\t{position}\t{base}\t{fraction}\n')


def run_m5c():
    parser = argparse.ArgumentParser(description='Run m5c analysis pipeline')
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--plus_file", required=True, help="m5c dampened fraction modified reads plus wig file")
    parser.add_argument("--minus_file", required=True, help="m5c dampened fraction modified reads minus wig file")
    parser.add_argument("--ref_tr_fa", required=True, help="Reference transcript FASTA file")
    args = parser.parse_args()
    out_dir = args.out_dir
    out_file = os.path.join(out_dir,'m5c.prediction.txt')
    plus_file = args.plus_file
    minus_file = args.minus_file
    ref_tr_fa = args.ref_tr_fa
    modify_dampened_fraction(plus_file,minus_file,ref_tr_fa,out_file)


if __name__ == '__main__':
    run_m5c()
