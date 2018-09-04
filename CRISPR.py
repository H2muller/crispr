#!/usr/bin/env python3
# Written by: Hans Müller Paul and Jacob Heldenbrand
#                           NOTES:


### Importing required libraries
import argparse
import datetime
import gffutils
import json
import re
import sys
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
from os import path

### Defining the arguments
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--fasta', required=True, dest='f', 
                    help='[required] path to input file in FASTA format'
                    )
# parser.add_argument('-g', '--gff', required=True, dest='g', 
#                     help='[required] path to input file in GFF format'
#                     )
parser.add_argument('-o', '--output', dest='o', default='output.txt',
                    help='path to output file'
                    )
parser.add_argument('-l', '--length', metavar='', dest='l', type=int, default=20,
                    help='length of the gRNA sequence, default = 20'
                    )
parser.add_argument('-L', '--flanking', metavar='', dest='L', type=int, default=500,
                    help='length of flanking region for verification, default = 500'
                    )
parser.add_argument('--cas9', action='store_true',
                    help='specifies that design will be made for the Cas9 CRISPR system'
                    )
parser.add_argument('--cpf1', action='store_true',
                    help='specifies that design will be made for the Cpf1 CRISPR system'
                    )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints each step of each iteration (for debugging)'
                    )

args = parser.parse_args()


### Defining functions
class PAM:
    def __init__(self, location, chromosome, guide_type, mainstrand=True):
        if mainstrand:
            strand = '+'
        else:
            strand = '-'
        self.strand = strand
        self.location = location
        self.chr = chromosome[1::]
        self.type = guide_type
        self.id = self.__get_id()
    
    def __get_id(self):
        return f'[{self.type}]{self.chr}.({self.location[0]}:{self.location[1]}).strand({self.strand})'


class CRISPR:
    def __init__(self, sequence, location, cutsite, chromosome, guide_type, strand):
        self.strand = strand
        self.sequence = sequence
        self.location = location
        self.cutsite = cutsite
        self.chr = chromosome
        self.type = guide_type
        self.id = self.__get_id()

    # def __flanking_region_for_sequencing(self):
    #     offset = args.l
    #     start = self-offset
    #     end = self+len(self.sgrna)+offset
    #     if start <=0:
    #         start == 0
    #     if end >= len(self.sequence):
    #         end == len(self.sequence)
    #     genomic_region = self.sequence[start:end]
    #     return genomic_region

    def __get_id(self):
        return f'[{self.type}]{self.chr}.({self.location[0]}:{self.location[1]}).strand({self.strand})'

def import_fasta_file(fasta):
    with open(fasta, 'r') as f:
        myFASTA = f.read()
        if args.verbose:
            print(f'Genome file {fasta} successfully imported')
        linecount = myFASTA.count('\n')
        if 2 * myFASTA.count('>') != linecount + 1:
            if args.verbose:
                print('formatting genome')
            from functions import formatted
            myFASTA = formatted(myFASTA)
            if args.verbose:
                print(f'Genome file {fasta} successfully formatted')
    from functions import generate_dictionary as gendict
    genome_dictionary = gendict(myFASTA)
    if args.verbose:
        print(f'The genome was successfully converted to a dictionary')
    return genome_dictionary


def import_gff_file(gff):
    with open(gff, 'r') as g:
        myGFF = g.read()
        if args.verbose:
            print(f'GFF file {gff} successfully imported')
        return myGFF
    featuredb = gffutils.create_db(myGFF, ':memory:', merge_strategy="create_unique",keep_order=True)
    if args.verbose:
        if len(featuredb) > 0:
            print('Feature database successfully generated')
    return featuredb


def find_PAM_site(target,input_sequence):
    PAM_site = [match.span() for match in re.finditer(target,input_sequence)] # added ints to start and end to match the PAM position
    return PAM_site


def get_reverse_complement(input_sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(input_sequence)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    return bases


def get_gRNA_sequence(input_sequence):
    complement = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
    RNA = list(input_sequence)
    RNA = reversed([complement.get(base, base) for base in RNA])
    RNA = ''.join(RNA)
    return RNA

"""
def get_sequence_by_type(type, featuredb, fasta):
    for parent in featuredb.features_of_type(type):
        parent_seq = parent.sequence(fasta)
        parent_seq = Seq(parent_seq, generic_dna)
        if parent.strand == '-':
            parent_seq = parent_seq.reverse_complement()
        return ('>' + parent.id + '\n' + parent_seq)


def get_features_under_transcript(type, featuredb, fasta):
    for transcript in featuredb.features_of_type('mRNA', order_by='start'): # or mRNA depending on the gff
        print('>' + transcript.id + '_' + type)
        seq_combined = ''
        for child in featuredb.children(transcript, featuretype=type, order_by='start'): # can also order by exon/intron
            child_seq = child.sequence(fasta, use_strand=False)  # use_strand is bugged in the current version
            seq_combined += child_seq
        seq_combined = Seq(seq_combined, generic_dna)
        if transcript.strand == '-':
            seq_combined = seq_combined.create_reverse_complement()
        for i in range(0, len(seq_combined), 60):
            print(seq_combined[i:i+60])
"""

def main():
    if not args.cas9 and not args.cpf1:
        sys.exit('Please select at least one CRISPR system: Cas9 or Cpf1')
    if args.verbose:
        print(args)


    ### Import genome files
    fasta_file = import_fasta_file(args.f)
    # gff_file = import_gff_file(args.g)


    ### Locate PAMs by nuclease type
    print(f'''
        Initiating PAM site detection.
        
        Please wait, this may take a while...
        ''')


    if args.cas9:
        invalid_cas9_targets = []
        for k,v in fasta_file.items():
            # + strand
            motif = re.compile(r'(?=.GG)')
            cas9_target_list = find_PAM_site(motif,v)
            for target in cas9_target_list:
                pam_location = (target[0]+1,target[0]+3)
                pam_site = PAM(pam_location,k,'cas9')
                sgrna_position = (pam_site.location[0]-(args.l+1),pam_site.location[0]-1)
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[0]:sgrna_position[1]]
                    sgrna_sequence = get_reverse_complement(sgrna_sequence)
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna = CRISPR(sgrna_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand)
                    entry = {
                            "sgRNA unique ID": sgrna.id,
                            "Guide sequence": sgrna.sequence,
                            "sgRNA position": (sgrna.location[0],sgrna.location[1]),
                            "Cut site position": sgrna.cutsite,
                            "Genomic features": 'sgrna.feature'
                            }
                    with open(args.o, mode='a') as CAS9:
                        json.dump(entry,CAS9,indent=2)
                else:
                    invalid_cas9_targets.append(pam_site)
            # - strand
            motif = re.compile(r'(?=CC.)')
            cas9_target_list2 = find_PAM_site(motif,v)
            for target in cas9_target_list2:
                pam_location = (target[0]+1,target[0]+3)
                pam_site = PAM(pam_location[::-1],k,'cas9',False)
                sgrna_position = (pam_site.location[0]+(args.l+1),pam_site.location[0]+1)
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[1]:sgrna_position[0]]
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna = CRISPR(sgrna_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand)
                    entry = {
                            "sgRNA unique ID": sgrna.id,
                            "Guide sequence": sgrna.sequence,
                            "sgRNA position": (sgrna.location[0],sgrna.location[1]),
                            "Cut site position": sgrna.cutsite,
                            "Genomic features": 'sgrna.feature'
                            }
                    with open(args.o, mode='a') as CAS9:
                        json.dump(entry,CAS9,indent=2)
                else:
                    invalid_cas9_targets.append(pam_site)
            if args.verbose:
                print (f'''
                {len(cas9_target_list + cas9_target_list2)} Cas9 PAM sites were found on {k[1::]}
                ''')

    if args.cpf1:
        invalid_cpf1_targets = []
        for k,v in fasta_file.items():
        # + strand
            motif = re.compile(r'(?=TTT.)')
            cpf1_target_list = find_PAM_site(motif,v)
            for target in cpf1_target_list:
                pam_location = (target[0]+1,target[0]+4)
                pam_site = PAM(pam_location,k,'cpf1')
                del pam_location
                sgrna_position = (pam_site.location[1]+1,pam_site.location[1]+(args.l+1))
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[0]:sgrna_position[1]]
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna = CRISPR(sgrna_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand)
                    entry = {
                            "sgRNA unique ID": sgrna.id,
                            "Guide sequence": sgrna.sequence,
                            "sgRNA position": (sgrna.location[0],sgrna.location[1]),
                            "Cut site position": sgrna.cutsite,
                            "Genomic features": 'sgrna.feature'
                            }
                    with open(args.o, mode='a') as CPF1:
                        json.dump(entry,CPF1,indent=2)
                else:
                    invalid_cpf1_targets.append(pam_site)

        # - strand
            motif = re.compile(r'(?=.AAA)')
            cpf1_target_list2 = find_PAM_site(motif,v)
            for target in cpf1_target_list2:
                pam_location = (target[0]+1,target[0]+4)
                pam_site = PAM(pam_location[::-1],k,'cpf1',False)
                del pam_location
                sgrna_position = (pam_site.location[1]-1,pam_site.location[1]-(args.l+1))
                if sgrna_position[0] >= 0 and sgrna_position[0] <= len(v) and sgrna_position[1] >= 0 and sgrna_position[1] <= len(v):
                    sgrna_sequence = v[sgrna_position[1]:sgrna_position[0]]
                    sgrna_sequence = get_gRNA_sequence(sgrna_sequence)
                    sgrna = CRISPR(sgrna_sequence,sgrna_position,sgrna_position[1]-3,pam_site.chr,pam_site.type,pam_site.strand)
                    entry = {
                            "sgRNA unique ID": sgrna.id,
                            "Guide sequence": sgrna.sequence,
                            "Cut site position": sgrna.cutsite,
                            "sgRNA chromosome": sgrna.chr,
                            "sgRNA position": (sgrna.location[0],sgrna.location[1]),
                            "sgRNA strand": sgrna.strand,
                            "Genomic features": 'sgrna.feature'
                            }
                    with open(args.o, mode='a') as CPF1:
                        json.dump(entry,CPF1,indent=2)
                else:
                    invalid_cpf1_targets.append(pam_site)
            if args.verbose:
                print (f'''
                {len(cpf1_target_list + cpf1_target_list2)} Cpf1 PAM sites were found on {k[1::]}
                ''')

    ### Gather info on each target motif



    """
    # Import files (GFF and FASTA)
    # Format files (GFF > DB; FASTA > Formatted) - done

    # on fasta:
    #     locate target (based on type: Cas9 or Cpf1) - done
    #     if args.cas9:
    #         design_cas9(location) - done
    #     if args.cpf1:
    #         design_cpf1(location)
    #     for target:
    #         nuclease type (cas9,cpf1) - done
    #         position (chr, (start,end)) - done
    #         cut position - done
    #         gRNA sequence - done
            From GFF.db:
                geneID or gene_name
                type (cds, ncs)
                feature (exon, intron, promoter)
    VERIFICATION:
    test doench score
    test number of hits with mismatches = 0,1,2,3,4
        recommended = hits-n(ploidy number) == 0 for 0-2 mismatches
    RANK:
    By  position in gene (early sequence) > conserved exon > num. off target > on target doench score
    RETURN:
        gene name / ID
        Nuclease type
        Cut site position
        guide RNA sequence
        Chr.start_pos : end_pos
        genome
        gene
        features (defined by cut site position)
    """


    ### WRITE TO OUTPUT FILE (add to head)

    # output_dict = {
    #     "Target genome": 'Sorghum Bicolor',
    #     "Total Cas9 gRNAs": len(cas9_guides),
    #     # "Total Cpf1 gRNAs": len(cpf1_guides),
    #     "Generated with": "$software name",
    #     "Date": str(datetime.datetime.now()),
    #     "Cas9 sgRNA Info": Cas9_guide_info,
    #     # "Cpf1 sgRNA Info": Cpf1_guide_info
    # }
    # with open(args.o, 'w') as f:
    #     json.dump(output_dict, f, indent=2)

### CONFIRMATION MESSAGE
    print(f'The output file has been generated at {args.o}')

if __name__ == '__main__':
    main()