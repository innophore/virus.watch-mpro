#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to process data from GISAID source files

@author: lena.parigger
"""
from Bio import SeqIO
import pandas as pd
import os
from multiprocessing import Pool
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
import sys

# reference Accession ID
reference_ID = 'EPI_ISL_402124'

# A FASTA file of reference protein sequences must be given, format: >NSP1\nSequence\n>NSP2\nSequence\n...
reference_file = 'reference_sequences.fasta'

# define path to the MSA FASTA file of SARS-CoV-2 genome sequences
path = 'path_to_msa.fasta'

# define the name of MSA FASTA file
msa_fasta = 'name_of_msa.fasta'

# SARS-CoV-2 proteins
proteins = ['NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'NSP6', 'NSP7', 'NSP8', 'NSP9', 'NSP10', 'NSP11', 'NSP12', 'NSP13',
     'NSP14', 'NSP15', 'NSP16', 'Spike', 'NS3', 'E', 'M', 'NS6', 'NS7a', 'NS7b', 'NS8', 'N', 'NS9b', 'NS9c', 'NS10']

# positions of proteins in reference genome sequence NC_045512.2
genome_positions = ['266-805', '806-2719', '2720-8554', '8555-10054', '10055-10972', '10973-11842', '11843-12091',
                    '12092-12685', '12686-13024', '13025-13441', '13442-13480', '13442-13468|13468-16236',
                    '16237-18039', '18040-19620', '19621-20658', '20659-21552', '21563-25384', '25393-26220',
                    '26245-26472', '26523-27191', '27202-27387', '27394-27759', '27756-27887', '27894-28259',
                    '28274-29533', '28284-28577', '28734-28955', '29558-29674']

# define parameters for pairwise alignment
gap_open = -10
gap_extend = -8
match_score = 5
mismatch_score = -5

# positions of proteins in MSA.fasta
msa_positions_start = []
msa_positions_stop = []

def retrieve_positions_of_proteins_in_msa():

    msa = path + msa_fasta
    for record in SeqIO.parse(msa, "fasta"):
        if record.description.split('|')[2] == reference_ID:
            reference_msa = record.seq

            for pro, pos in zip(proteins, genome_positions):
                if pos.count('|') == 0:
                    start = int(pos.split('-')[0])
                    end = int(pos.split('-')[-1])
                    count = 0
                    for ind_i, i in enumerate(reference_msa):
                        if i != '-':
                            count = count + 1
                            if count == start:
                                msa_positions_start.append(str(ind_i))
                            if count == end:
                                msa_positions_stop.append(str(ind_i + 1))
                else:
                    start1 = int(pos.split('|')[0].split('-')[0])
                    end1 = int(pos.split('|')[0].split('-')[-1])
                    start2 = int(pos.split('|')[-1].split('-')[0])
                    end2 = int(pos.split('|')[-1].split('-')[-1])
                    frameshift_start = []
                    frameshift_stop = []
                    count = 0
                    for ind_i, i in enumerate(reference_msa):
                        if i != '-':
                            count = count + 1
                            if count == start1:
                                frameshift_start.append(str(ind_i))
                            if count == end1:
                                frameshift_stop.append(str(ind_i + 1))
                            if count == start2:
                                frameshift_start.append(str(ind_i))
                            if count == end2:
                                frameshift_stop.append(str(ind_i + 1))
                    msa_positions_start.append('_'.join(frameshift_start))
                    msa_positions_stop.append('_'.join(frameshift_stop))

retrieve_positions_of_proteins_in_msa()


def retrieve_nt_sequences_from_msa():

    exists = os.path.exists(path + 'NT_sequences/')
    if not exists:
        os.mkdir(path + 'NT_sequences/')

    for protein in proteins:
        try:
            os.remove(path + 'NT_sequences/{}_nt.fasta'.format(protein))
        except FileNotFoundError:
            pass

    for record in SeqIO.parse(path + msa_fasta, "fasta"):
        for protein, start, stop in zip(proteins, msa_positions_start, msa_positions_stop):
            if start.count('_') == 0:
                open(path + 'NT_sequences/{}_nt.fasta'.
                     format(protein), 'a').write('>{}|{}\n{}\n'.format(protein, record.description, record.seq[int(start):int(stop)].replace('-', '')))
            else:
                start1 = int(start.split('_')[0])
                start2 = int(start.split('_')[-1])
                stop1 = int(stop.split('_')[0])
                stop2 = int(stop.split('_')[-1])
                open(path + 'NT_sequences/{}_nt.fasta'.format(protein), 'a').write('>{}|{}\n{}\n'.format(protein, record.description,
                                                record.seq[int(start1):int(stop1)].replace('-', '') +
                                                record.seq[int(start2):int(stop2)].replace('-', '')))

retrieve_nt_sequences_from_msa()


def translate_and_split_at_stop():

    exists = os.path.exists(path + 'ALL/')
    if not exists:
        os.mkdir(path + 'ALL/')

    global translate

    def translate(protein):

        infile = path + 'NT_sequences/{}_nt.fasta'.format(protein)
        outfile = path + 'ALL/{}.fasta'.format(protein)

        try:
            os.remove(outfile)
        except FileNotFoundError:
            pass

        out = open(outfile, 'a')
        for record in SeqIO.parse(infile, "fasta"):
            trans = record.seq.translate(1)
            if trans.count('*') != 0:
                out.write('>{}\n{}*\n'.format(record.description, record.seq.translate(1).split('*')[0].replace('J', 'X').replace('B', 'X').replace('Z', 'X')))
            else:
                out.write('>{}\n{}\n'.format(record.description, record.seq.translate(1).replace('J', 'X').replace('B', 'X').replace('Z', 'X')))
        out.close()

    if __name__ == '__main__':
        with Pool(8) as p:
            p.map(translate, proteins)

translate_and_split_at_stop()


def process_protein_sequences():

    for protein in proteins:

        print('Working on protein {}'.format(protein))

        def drop_duplicate_sequences():
            print('FASTA files of unique sequences are produced...')

            exists = os.path.exists(path + 'UNIQUE/')
            if not exists:
                os.mkdir(path + 'UNIQUE/')

            outfile = '{}_unique.fasta'.format(protein)

            try:
                os.remove('%sUNIQUE/%s' % (path, outfile))
            except FileNotFoundError:
                pass

            infile = open(path + 'ALL/{}.fasta'.format(protein), 'r').readlines()

            df = pd.DataFrame(list(infile[1::2]), columns=['Sequence'])
            df.drop_duplicates(subset=['Sequence'], inplace=True)
            df_list = df['Sequence'].tolist()
            nrs = list(range(1, len(df_list) + 1))

            out2 = open('%sUNIQUE/%s' % (path, outfile), 'a')
            for item, nr in zip(df_list, nrs):
                out2.write('>%s\n%s' % (nr, item))
            out2.close()

        drop_duplicate_sequences()


        def align_unique_sequences():
            print('Unique sequences are aligned to reference sequences...')

            exists = os.path.exists(path + 'ALIGNMENTS/')
            if not exists:
                os.mkdir(path + 'ALIGNMENTS/')

            references = open(path + 'reference_sequences.fasta', 'r').readlines()
            reference_header = [i[1:-1] for i in references[0::2]]
            reference_sequence = [i[:-1] for i in references[1::2]]

            outfile = path + 'ALIGNMENTS/{}_pairwise.txt'.format(protein)
            infile = open(path + 'UNIQUE/{}_unique.fasta'.format(protein), 'r').readlines()

            refseq2 = []
            for x, y in zip(reference_header, reference_sequence):
                if x == protein:
                    refseq2.append(y)
            refseq = refseq2[0]

            try:
                os.remove(outfile)
            except FileNotFoundError:
                pass

            file_seq = [i[:-1] if i[-1] == '\n' else i for i in infile[1::2]]  # write the sequences into a list

            out = open(outfile, 'a')  # open the outfile

            for seq in file_seq:
                seq2 = Seq('{}'.format(seq))
                alignments = pairwise2.align.globalms(refseq, seq2, match_score, mismatch_score,
                                                      gap_open, gap_extend, one_alignment_only=True)
                try:
                    alignment = alignments[0]
                    out.write('%s\n' % (format_alignment(*alignment)))
                except IndexError:
                    pass

        align_unique_sequences()


        def recover_aligned_sequences():
            print('The aligned sequences are written into a FASTA file...')

            exists = os.path.exists(path + 'RECOVERED/')
            if not exists:
                os.mkdir(path + 'RECOVERED/')

            outfile = path + 'RECOVERED/{}_recovered.fasta'.format(protein)
            infile = open(path + 'ALIGNMENTS/{}_pairwise.txt'.format(protein), 'r').readlines()

            try:
                os.remove(outfile)
            except FileNotFoundError:
                pass

            sequence = infile[2::5]

            unique_id = list(range(1, len(sequence) + 1))  # new IDs are given (not the same as in UNIQUE)
            cured_sequences = []
            for se in sequence:
                if se[-1] == '\n':
                    cured_sequences.append(se.replace('-', '')[:-1])
                else:
                    cured_sequences.append(se.replace('-', ''))

            out = open(outfile, 'a')
            for ids, seq in zip(unique_id, cured_sequences):
                out.write('>%s\n%s\n' % (ids, seq))
            out.close()

        recover_aligned_sequences()


        def add_header_to_alignment():
            print('The new IDs referring to the recovered sequences are added to the alignment file...')

            exists = os.path.exists(path + 'ALIGNMENTS/')
            if not exists:
                os.mkdir(path + 'ALIGNMENTS/')

            in_alignment = open(path + 'ALIGNMENTS/{}_pairwise.txt'.format(protein), 'r').readlines()
            in_recovered = open(path + 'RECOVERED/{}_recovered.fasta'.format(protein), 'r').readlines()
            out_alignment = path + 'ALIGNMENTS/{}_pairwise_with_headers.txt'.format(protein)

            try:
                os.remove(out_alignment)
            except FileNotFoundError:
                pass

            recovered_header = [int(i[1:-1]) for i in in_recovered[0::2]]
            recovered_sequences = [i[:-1] if i[-1] == '\n' else i for i in in_recovered[1::2]]
            df_recovered = pd.DataFrame(list(zip(recovered_header, recovered_sequences)), columns=['header', 'sequences'])

            alignment_refseq = in_alignment[0::5]
            alignment_align = in_alignment[1::5]
            alignment_sequ = [i[:-1] for i in in_alignment[2::5]]
            alignment_score = in_alignment[3::5]
            empty = in_alignment[4::5]
            sequences = [i.replace('-', '') for i in alignment_sequ]

            df_alignment = pd.DataFrame(list(zip(alignment_refseq, alignment_align, alignment_sequ, alignment_score, empty, sequences)), columns=['refseq', 'align', 'sequ', 'score', 'empty', 'sequences'])

            df_merged = df_recovered.merge(df_alignment, how='left', on='sequences')
            df_merged.sort_values('header', inplace=True)
            print('These numbers should be the same, otherwise not all aligned sequences were found in the recovered '
                  'sequences: {} - {} - {}'.format(len(df_recovered), len(df_alignment), len(df_merged)))
            out_header = df_merged['header'].tolist()
            out_refseq = df_merged['refseq'].tolist()
            out_align = df_merged['align'].tolist()
            out_sequences = df_merged['sequ'].tolist()
            out_score = df_merged['score'].tolist()
            out_empty = df_merged['empty'].tolist()

            out = open(out_alignment, 'a')
            for he, re, al, se, sc, em in zip(out_header, out_refseq, out_align, out_sequences, out_score, out_empty):
                out.write('{}\n{}{}{}\n{}{}'.format(he,re,al,se,sc,em))
            out.close()

        add_header_to_alignment()


        def search_in_all_sequences():
            print('The recovered sequences are merged with the sequences in ALL/ to retrieve their abundance...')

            exists = os.path.exists(path + 'ALLSEQ/')
            if not exists:
                os.mkdir(path + 'ALLSEQ/')

            infile_recovered = path + 'RECOVERED/{}_recovered.fasta'.format(protein)
            out_allseq = path + 'ALLSEQ/{}_allseq.txt'.format(protein)

            in_all = path + 'ALL/{}.fasta'.format(protein)

            # Make a dataframe of all sequences
            # To speed up the merge process, sequences are devided by the occurrence of X
            headers_with_x = []
            sequences_with_x = []
            headers_without_x = []
            sequences_without_x = []
            for he, se in zip(open(in_all, 'r').readlines()[0::2], open(in_all, 'r').readlines()[1::2]):
                if se.count('X') == 0:
                    headers_without_x.append(he[:-1])
                    if se[-1] == '\n':
                        sequences_without_x.append(se[:-1])
                    else:
                        sequences_without_x.append(se)
                else:
                    headers_with_x.append(he[:-1])
                    if se[-1] == '\n':
                        sequences_with_x.append(se[:-1])
                    else:
                        sequences_with_x.append(se)

            df_all_without_x = pd.DataFrame(list(zip(headers_without_x, sequences_without_x)), columns=['Header', 'Sequences'])
            print('Length of dataframe with all sequences without X: ' + str(len(df_all_without_x)))
            df_all_without_x.drop_duplicates('Header', inplace=True)
            print('when dropped: {}'.format(len(df_all_without_x)))
            df_all_with_x = pd.DataFrame(list(zip(headers_with_x, sequences_with_x)), columns=['Header', 'Sequences'])
            print('Length of dataframe with all sequences with X: ' + str(len(df_all_with_x)))
            df_all_with_x.drop_duplicates('Header', inplace=True)
            print('when dropped: {}'.format(len(df_all_with_x)))

            # Make dataframe of unique sequences
            # To speed up the merge process, sequences are devided by the occurrence of X
            unique_sequences_with_x = []
            unique_headers_with_x = []
            unique_sequences_without_x = []
            unique_headers_without_x = []
            for he_rec, se_rec in zip(open(infile_recovered, 'r').readlines()[0::2], open(infile_recovered, 'r').readlines()[1::2]):
                if se_rec.count('X') == 0:
                    unique_headers_without_x.append(he_rec[1:-1])
                    if se_rec[-1] == '\n':
                        unique_sequences_without_x.append(se_rec[:-1])
                    else:
                        unique_sequences_without_x.append(se_rec)
                else:
                    unique_headers_with_x.append(he_rec[1:-1])
                    if se_rec[-1] == '\n':
                        unique_sequences_with_x.append(se_rec[:-1])
                    else:
                        unique_sequences_with_x.append(se_rec)

            df_unique_with_x = pd.DataFrame(list(zip(unique_headers_with_x, unique_sequences_with_x)),
                                     columns=['ID', 'Sequences'])
            print('Length of dataframe with unique sequences with x: ' + str(len(df_unique_with_x)))
            df_unique_with_x.drop_duplicates('ID', inplace=True)
            print('when dropped: {}'.format(len(df_unique_with_x)))

            df_unique_without_x = pd.DataFrame(list(zip(unique_headers_without_x, unique_sequences_without_x)),
                                     columns=['ID', 'Sequences'])
            print('Length of dataframe with unique sequences without x: ' + str(len(df_unique_without_x)))
            df_unique_without_x.drop_duplicates('ID', inplace=True)
            print('when dropped: {}'.format(len(df_unique_without_x)))

            # Merge the two dataframes and drop duplicates (there should not be any duplicates but just in case)
            dfmerged_with_X = df_all_with_x.merge(df_unique_with_x, how='left', on=['Sequences'])
            dfmerged_without_X = df_all_without_x.merge(df_unique_without_x, how='left', on=['Sequences'])

            try:
                os.remove(out_allseq)
            except FileNotFoundError:
                pass

            df_cured1 = dfmerged_without_X.append(dfmerged_with_X)
            df_cured = df_cured1[df_cured1['ID'].notna()]
            df_astype = df_cured.astype({'ID': 'int'})
            df_astype.sort_values('ID', inplace=True)
            list_sort = df_astype.values.tolist()
            out = open(out_allseq, 'a')
            for entry in list_sort:
                out.write('>{}_{}\n'.format(entry[2], entry[0][1:]))
            out.close()

            print('***search is checked***')
            df_astype.drop_duplicates('ID', inplace=True)
            ids_dropped = df_astype['ID'].tolist()
            print('These must be the same, otherwise not all unique sequences were found: {} - {}'
                  .format(len(ids_dropped), len(unique_headers_with_x) + len(unique_headers_without_x)))

        search_in_all_sequences()


        def create_the_mutation_list():
            print('Mutations are parsed and written into a csv file...')

            exists = os.path.exists(path + 'MUTLISTS/')
            if not exists:
                os.mkdir(path + 'MUTLISTS/')

            list_fullalignment = open(path + 'ALIGNMENTS/{}_pairwise_with_headers.txt'.format(protein), 'r').readlines()
            infile_allseq = open(path + 'ALLSEQ/{}_allseq.txt'.format(protein), 'r').readlines()
            outfile_complete = path + 'MUTLISTS/{}_complete_mutation_list.csv'.format(protein)

            # Write the aligned sequences to a list
            headers = []
            for he in list_fullalignment[0::6]:
                headers.append(int(he[:-1]))
            refseq = list_fullalignment[1::6]
            alignm = list_fullalignment[2::6]
            sequence = []
            for se in list_fullalignment[3::6]:
                if se[-1] == '\n':
                    sequence.append(se.replace('-', '')[:-1])
                else:
                    sequence.append(se.replace('-', ''))

            mut_nr = []
            ex_nr = []
            ins_nr = []
            del_nr = []
            mutations = []
            protein_size = []
            for se1 in sequence:

                if se1[-1] == '*':
                    protein_size.append(len(se1) - 1)
                else:
                    protein_size.append(len(se1))
            print('*** wrote list ***')

            # Write information to dataframe
            position = list(range(1, 10000))
            for re, al, se in zip(refseq, alignm, list_fullalignment[3::6]):
                ex_nr1 = []
                ins_nr1 = []
                del_nr1 = []
                mutations1 = []
                # Set count for the "real" index (insertions lead to '-' in reference sequence, which alters the position number)
                count = 0
                for character1, character2, character3, pos in zip(re, al, se, position):
                    # Insertions
                    if character1 == '-':
                        if character3 == '\n':
                            pass
                        else:
                            count = count - 1
                            ins_nr1.append('.')
                            mutations1.append('+%i%s' % (pos + count, character3))
                    # Deletions
                    elif character3 == '-':
                        if character1 == '\n':
                            pass
                        else:
                            del_nr1.append('.')
                            mutations1.append('%s%i-' % (character1, pos + count))
                    # Amino-acid exchanges
                    elif character2 == '.':
                        if character3 == '\n':
                            pass
                        else:
                            ex_nr1.append('.')
                            mutations1.append('%s%i%s' % (character1, pos + count, character3))

                ex_nr.append(len(ex_nr1))
                ins_nr.append(len(ins_nr1))
                del_nr.append(len(del_nr1))
                mutations.append(', '.join(mutations1))
                mut_nr.append(len(ex_nr1) + len(ins_nr1) + len(del_nr1))

            print('*** mutations successfully parsed ***')

            # make a list of alignment scores
            scores = [float(i[:-1].replace(' ', '').replace('Score=', '')) if i[-1] == '\n' else float(
                i.replace(' ', '').replace('Score=', '')) for i in list_fullalignment[4::6]]

            # make a list of the number of X in the sequences
            X_count = []
            abundance_dummy = []
            for se in sequence:
                X_count.append(se.count('X'))
                abundance_dummy.append(0)

            df_partly = pd.DataFrame(list(zip(headers, protein_size, X_count, scores, mut_nr, ex_nr, ins_nr, del_nr,
                         mutations, sequence, abundance_dummy)),
                columns=['ID', 'Protein length', 'Number_of_X', 'Alignment_score', 'Mutation number', 'Amino acid exchanges',
                         'Insertions', 'Deletions', 'Mutations', 'Sequence', 'Abundance'])

            # find the abundance of unique sequence in allseq.txt
            id_all = []
            header_all = []
            for entry in infile_allseq:
                id_all.append(int(entry[1:].split('_')[0]))
                header_all.append(entry[entry.index('_') + 1:])

            df_allseq = pd.DataFrame(list(zip(id_all, header_all)), columns=['ID', 'Header'])
            abundance = df_allseq['ID'].value_counts().rename_axis('ID').reset_index(name='Abundance')

            print('*** abundance of unique sequences found ***')

            # make a list of alignment scores
            scores = [float(i[:-1].replace(' ', '').replace('Score=', '')) if i[-1] == '\n' else float(
                i.replace(' ', '').replace('Score=', '')) for i in list_fullalignment[4::6]]

            # make a list of the number of X in the sequences
            X_count = []
            for se in sequence:
                X_count.append(se.count('X'))

            df_partly = pd.DataFrame(
                list(zip(headers, protein_size, X_count, scores, mut_nr, ex_nr, ins_nr, del_nr,
                         mutations, sequence)),
                columns=['ID', 'Protein length', 'Number_of_X', 'Alignment_score', 'Mutation number',
                         'Amino acid exchanges',
                         'Insertions', 'Deletions', 'Mutations', 'Sequence'])

            if len(abundance) == len(df_partly):
                df_partly_merged1 = df_partly.merge(abundance, how='left', on='ID')
            else:
                print(
                    'THERE IS SOMETHING WRONG - length of abundance = {} does not equal length of df_partly = {}'.format(
                        len(abundance), len(df_partly)))
                sys.exit(1)

            df_partly_merged = df_partly_merged1[df_partly_merged1['Abundance'].notna()]

            df_partly_merged.to_csv(outfile_complete)

        create_the_mutation_list()


        def find_high_quality_amino_acid_exchanges():
            print('High quality amino acid exchanges are retrieved from all mutations...')

            references = open(path + 'reference_sequences.fasta', 'r').readlines()
            reference_header = [i[1:-1] for i in references[0::2]]
            reference_sequence = [i[:-1] for i in references[1::2]]

            refseq2 = []
            for x, y in zip(reference_header, reference_sequence):
                if x == protein:
                    refseq2.append(y)
            refseqwt = refseq2[0]

            df_complete = pd.read_csv(path + 'MUTLISTS/{}_complete_mutation_list.csv'.format(protein))
            df_complete_out = path + 'MUTLISTS/{}_complete_mutation_list.csv'.format(protein)

            # find high quality exchanges
            mut_strings1 = df_complete['Mutations'].tolist()
            mut_strings = []
            float_positions = []

            for mu_st in mut_strings1:
                if isinstance(mu_st, float) is True:
                    mut_strings.append('')
                    float_positions.append('')
                else:
                    altered_ins = []
                    float_pos1 = []  # to find the mutated positions later by number
                    single_mu = mu_st.split(', ')
                    for m_s in single_mu:
                        if m_s[0] == '+':
                            new_string = m_s[0] + str(int(m_s[1:-1]) + 0.6) + m_s[-1]
                            altered_ins.append(new_string)
                            float_pos1.append(float(m_s[1:-1] + '.6'))
                        else:
                            altered_ins.append(m_s)
                            float_pos1.append(float(m_s[1:-1]))
                    mut_strings.append(', '.join(altered_ins))
                    float_positions.append(float_pos1)
            high_qual_muts = []
            high_qual_muts.append('')
            for mut_string, float_pos in zip(mut_strings[1:], float_positions[1:]):
                high_qual_mutations = []
                if isinstance(mut_string, float) is True:
                    high_qual_mutations.append('')
                else:
                    mut_string_list = mut_string.split(', ')
                    indices = []
                    count = 0
                    for item1, item2 in zip(mut_string_list, mut_string_list[1:]):
                        count = count + 1
                        if item2[1:-1].endswith('.6') and item1[1:-1].endswith('.6'):
                            if int(float(item2[1:-1])) - int(float(item1[1:-1])) > 0:
                                indices.append(count)
                        else:
                            if int(float(item2[1:-1])) - int(float(item1[1:-1])) > 1:
                                indices.append(count)

                    def partition(mut, indi):
                        return [mut[i:j] for i, j in zip([0] + indi, indi + [None])]

                    mut_list = partition(mut_string_list, indices)
                    joined_mutlist = [', '.join(i) for i in mut_list]
                    # if insertions or deletions are involved, mutation cannot certainly be assigned to a position
                    if len(joined_mutlist) == 1:
                        if joined_mutlist[0].count('-') == 0 and joined_mutlist[0].count('+') == 0 and joined_mutlist[
                            0].count('X') == 0:
                            high_qual_mutations.append(joined_mutlist[0])

                    # if two separate mutated regions, they could influence each other
                    elif len(joined_mutlist) == 2:
                        if joined_mutlist[0].count('-') == 0 and joined_mutlist[0].count('+') == 0 and joined_mutlist[
                            0].count('X') == 0:
                            if int(round(float(mut_list[1][0][1:-1]))) - int(float(mut_list[0][-1][1:-1])) <= 3:
                                if int(float(mut_list[0][0][
                                             1:-1])) > 3:
                                    positions2 = list(
                                        range(int(round(float(mut_list[0][-1][1:-1]) + 0.6)),
                                              int(round(float(mut_list[0][-1][1:-1]) + 0.6)) + 10))
                                    nr = []
                                    insdels = []
                                    count_in = 0
                                    for index_pos, pos in enumerate(
                                            positions2):
                                        if index_pos >= len(positions2) - count_in:
                                            break

                                        count_ins = float_pos.count(float(pos) + 0.6)
                                        if count_ins >= len(positions2) - index_pos:
                                            count_real_ins = len(positions2) - index_pos
                                        else:
                                            count_real_ins = count_ins

                                        count_in = count_in + count_real_ins
                                        nr.append(count_real_ins)
                                        insdels.append(count_real_ins)
                                        nr.append(float_pos.count(float(pos)))
                                        try:
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                        except IndexError:
                                            pass
                                    if insdels.count(0) == len(insdels):
                                        high_qual_mutations.append(joined_mutlist[0])
                                    else:
                                        if sum(nr) <= 0.5 * len(positions2) and sum(nr) <= 0.5 * (
                                                len(refseqwt) - (int(float(mut_list[0][-1][1:-1]) + 0.6) - 1)):
                                            high_qual_mutations.append(joined_mutlist[0])

                                else:
                                    if joined_mutlist[-1].count('-') == 0 and joined_mutlist[
                                        -1].count('+') == 0:
                                        high_qual_mutations.append(joined_mutlist[0])
                            else:
                                positions2 = list(
                                    range(int(round(float(mut_list[0][-1][1:-1]) + 0.6)),
                                          int(round(float(mut_list[0][-1][1:-1]) + 0.6)) + 10))

                                nr = []
                                insdels = []
                                count_in = 0
                                for index_pos, pos in enumerate(
                                        positions2):
                                    if index_pos >= len(positions2) - count_in:
                                        break

                                    count_ins = float_pos.count(float(pos) + 0.6)
                                    if count_ins >= len(positions2) - index_pos:
                                        count_real_ins = len(positions2) - index_pos
                                    else:
                                        count_real_ins = count_ins

                                    count_in = count_in + count_real_ins
                                    nr.append(count_real_ins)
                                    insdels.append(count_real_ins)
                                    nr.append(float_pos.count(float(pos)))
                                    try:
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                    except IndexError:
                                        pass
                                if insdels.count(0) == len(insdels):
                                    high_qual_mutations.append(joined_mutlist[0])
                                else:
                                    if sum(nr) <= 0.5 * len(positions2) and sum(nr) <= 0.5 * (
                                            len(refseqwt) - (int(float(mut_list[0][-1][1:-1]) + 0.6) - 1)):
                                        high_qual_mutations.append(joined_mutlist[0])

                        if joined_mutlist[1].count('-') == 0 and joined_mutlist[1].count('+') == 0 and joined_mutlist[
                            1].count('X') == 0:
                            if int(round(float(mut_list[1][0][1:-1]))) - int(float(mut_list[0][-1][1:-1])) <= 3:
                                if int(float(mut_list[1][-1][1:-1])) < len(
                                        refseqwt) - 2:
                                    positions1 = list(
                                        range(int(round(float(mut_list[1][0][1:-1]))) - 10,
                                              int(round(float(mut_list[1][0][1:-1])))))[::-1]
                                    nr = []
                                    insdels = []
                                    count_in = 0
                                    for index_pos, pos in enumerate(
                                            positions1):
                                        if index_pos >= len(positions1) - count_in:
                                            break

                                        count_ins = float_pos.count(float(pos) + 0.6)
                                        if count_ins >= len(positions1) - index_pos:
                                            count_real_ins = len(positions1) - index_pos
                                        else:
                                            count_real_ins = count_ins

                                        count_in = count_in + count_real_ins
                                        nr.append(count_real_ins)
                                        insdels.append(count_real_ins)
                                        nr.append(float_pos.count(float(pos)))
                                        try:
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                        except IndexError:
                                            pass
                                    if insdels.count(0) == len(insdels):
                                        high_qual_mutations.append(joined_mutlist[1])
                                    else:
                                        if sum(nr) <= 0.5 * len(positions1) and sum(nr) <= 0.5 * int(
                                                round(float(mut_list[1][0][1:-1]))):
                                            high_qual_mutations.append(joined_mutlist[1])


                                else:
                                    if joined_mutlist[0].count('-') == 0 and joined_mutlist[0].count('+') == 0 and \
                                            joined_mutlist[0].count('X') == 0:
                                        high_qual_mutations.append(joined_mutlist[1])
                            else:
                                positions1 = list(
                                    range(int(round(float(mut_list[1][0][1:-1]))) - 10,
                                          int(round(float(mut_list[1][0][1:-1])))))[::-1]

                                nr = []
                                insdels = []
                                count_in = 0
                                for index_pos, pos in enumerate(
                                        positions1):
                                    if index_pos >= len(positions1) - count_in:
                                        break

                                    count_ins = float_pos.count(float(pos) + 0.6)
                                    if count_ins >= len(positions1) - index_pos:
                                        count_real_ins = len(positions1) - index_pos
                                    else:
                                        count_real_ins = count_ins

                                    count_in = count_in + count_real_ins
                                    nr.append(count_real_ins)
                                    insdels.append(count_real_ins)
                                    nr.append(float_pos.count(float(pos)))
                                    try:
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                    except IndexError:
                                        pass
                                if insdels.count(0) == len(insdels):
                                    high_qual_mutations.append(joined_mutlist[1])
                                else:
                                    if sum(nr) <= 0.5 * len(positions1) and sum(nr) <= 0.5 * int(
                                            round(float(mut_list[1][0][1:-1]))):
                                        high_qual_mutations.append(joined_mutlist[1])

                    # if at least three separated regions
                    else:
                        if joined_mutlist[0].count('-') == 0 and joined_mutlist[0].count(
                                '+') == 0 and joined_mutlist[0].count(
                            'X') == 0:
                            if int(float(mut_list[0][0][1:-1])) > 3:
                                positions2 = list(range(int(round(float(mut_list[0][-1][1:-1]) + 0.6)),
                                          int(round(float(mut_list[0][-1][1:-1]) + 0.6)) + 10))
                                nr = []
                                insdels = []
                                count_in = 0
                                for index_pos, pos in enumerate(
                                        positions2):
                                    if index_pos >= len(positions2) - count_in:
                                        break

                                    count_ins = float_pos.count(float(pos) + 0.6)
                                    if count_ins >= len(positions2) - index_pos:
                                        count_real_ins = len(positions2) - index_pos
                                    else:
                                        count_real_ins = count_ins

                                    count_in = count_in + count_real_ins
                                    nr.append(count_real_ins)
                                    insdels.append(count_real_ins)
                                    nr.append(float_pos.count(float(pos)))
                                    try:
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                    except IndexError:
                                        pass
                                if insdels.count(0) == len(insdels):
                                    high_qual_mutations.append(joined_mutlist[0])
                                else:
                                    if sum(nr) <= 0.5 * len(positions2) and sum(nr) <= 0.5 * (
                                            len(refseqwt) - (int(float(mut_list[0][-1][1:-1]) + 0.6) - 1)):
                                        high_qual_mutations.append(joined_mutlist[0])

                            else:
                                if int(round(float(mut_list[1][0][1:-1]))) - int(float(mut_list[0][-1][
                                                                                       1:-1])) <= 3:
                                    if joined_mutlist[1].count('-') == 0 and joined_mutlist[1].count(
                                            '+') == 0 and joined_mutlist[1].count(
                                        'X') == 0:
                                        positions2 = list(
                                            range(int(round(float(mut_list[0][-1][1:-1]) + 0.6)),
                                                  int(round(float(mut_list[0][-1][1:-1]) + 0.6)) + 10))
                                        nr = []
                                        insdels = []
                                        count_in = 0
                                        for index_pos, pos in enumerate(
                                                positions2):
                                            if index_pos >= len(positions2) - count_in:
                                                break

                                            count_ins = float_pos.count(float(pos) + 0.6)
                                            if count_ins >= len(positions2) - index_pos:
                                                count_real_ins = len(positions2) - index_pos
                                            else:
                                                count_real_ins = count_ins

                                            count_in = count_in + count_real_ins
                                            nr.append(count_real_ins)
                                            insdels.append(count_real_ins)
                                            nr.append(float_pos.count(float(pos)))
                                            try:
                                                insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                                insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                            except IndexError:
                                                pass
                                        if insdels.count(0) == len(insdels):
                                            high_qual_mutations.append(joined_mutlist[0])
                                        else:
                                            if sum(nr) <= 0.5 * len(positions2) and sum(nr) <= 0.5 * (
                                                    len(refseqwt) - (int(float(mut_list[0][-1][1:-1]) + 0.6) - 1)):
                                                high_qual_mutations.append(joined_mutlist[0])

                                else:
                                    positions2 = list(
                                        range(int(round(float(mut_list[0][-1][1:-1]) + 0.6)),
                                              int(round(float(mut_list[0][-1][1:-1]) + 0.6)) + 10))
                                    nr = []
                                    insdels = []
                                    count_in = 0
                                    for index_pos, pos in enumerate(
                                            positions2):
                                        if index_pos >= len(positions2) - count_in:
                                            break

                                        count_ins = float_pos.count(float(pos) + 0.6)
                                        if count_ins >= len(positions2) - index_pos:
                                            count_real_ins = len(positions2) - index_pos
                                        else:
                                            count_real_ins = count_ins

                                        count_in = count_in + count_real_ins
                                        nr.append(count_real_ins)
                                        insdels.append(count_real_ins)
                                        nr.append(float_pos.count(float(pos)))
                                        try:
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                        except IndexError:
                                            pass
                                    if insdels.count(0) == len(insdels):
                                        high_qual_mutations.append(joined_mutlist[0])
                                    else:
                                        if sum(nr) <= 0.5 * len(positions2) and sum(nr) <= 0.5 * (
                                                len(refseqwt) - (int(float(mut_list[0][-1][1:-1]) + 0.6) - 1)):
                                            high_qual_mutations.append(joined_mutlist[0])

                        for jmut1, jmut2, jmut3, mut1, mut2, mut3 in zip(joined_mutlist, joined_mutlist[1:],
                                                                         joined_mutlist[2:],
                                                                         mut_list, mut_list[1:], mut_list[2:]):

                            if jmut2.count('-') == 0 and jmut2.count('+') == 0 and jmut2.count('X') == 0:
                                positions1 = list(
                                    range(int(round(float(mut2[0][1:-1]))) - 10,
                                          int(round(float(mut2[0][1:-1])))))[::-1]
                                positions2 = list(
                                    range(int(round(float(mut2[-1][1:-1]) + 0.6)),
                                          int(round(float(mut2[-1][1:-1]) + 0.6)) + 10))

                                if int(round(float(mut2[0][1:-1]))) - int(float(mut1[-1][1:-1])) <= 3 and int(
                                        round(float(mut3[0][1:-1]))) - int(float(
                                    mut2[-1][1:-1])) <= 3:
                                    if (jmut1.count('-') == 0 and jmut1.count('+') == 0 and jmut1.count('X') == 0) or \
                                            (jmut3.count('-') == 0 and jmut3.count('+') == 0 and jmut3.count('X') == 0):

                                        nr1 = []
                                        nr2 = []
                                        insdels1 = []
                                        insdels2 = []
                                        done = []
                                        count_in = 0
                                        for index_pos, pos in enumerate(
                                                positions2):
                                            if index_pos >= len(positions2) - count_in:
                                                break

                                            count_ins = float_pos.count(float(pos) + 0.6)
                                            if count_ins >= len(positions2) - index_pos:
                                                count_real_ins = len(positions2) - index_pos
                                            else:
                                                count_real_ins = count_ins

                                            count_in = count_in + count_real_ins
                                            nr2.append(count_real_ins)
                                            insdels2.append(count_real_ins)
                                            nr2.append(float_pos.count(float(pos)))
                                            try:
                                                insdels2.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                                insdels2.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                            except IndexError:
                                                pass
                                        if insdels2.count(0) == len(insdels2):
                                            done.append('.')
                                        else:
                                            if sum(nr2) <= 0.5 * len(positions2) and sum(nr2) <= 0.5 * (
                                                    len(refseqwt) - (int(float(mut2[-1][1:-1]) + 0.6) - 1)):
                                                done.append('.')
                                        if len(done) == 1:
                                            count_in = 0
                                            for index_pos, pos in enumerate(
                                                    positions1):
                                                if index_pos >= len(positions1) - count_in:
                                                    break

                                                count_ins = float_pos.count(float(pos) + 0.6)
                                                if count_ins >= len(positions1) - index_pos:
                                                    count_real_ins = len(positions1) - index_pos
                                                else:
                                                    count_real_ins = count_ins

                                                count_in = count_in + count_real_ins
                                                nr1.append(count_real_ins)
                                                insdels1.append(count_real_ins)
                                                nr1.append(float_pos.count(float(pos)))
                                                try:
                                                    insdels1.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                                    insdels1.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                                except IndexError:
                                                    pass
                                            if insdels1.count(0) == len(insdels1):
                                                high_qual_mutations.append(jmut2)
                                            else:
                                                if sum(nr1) <= 0.5 * len(positions1) and sum(nr1) <= 0.5 * int(
                                                        round(float(mut2[0][1:-1]))):
                                                    high_qual_mutations.append(jmut2)

                                else:
                                    nr1 = []
                                    nr2 = []
                                    insdels1 = []
                                    insdels2 = []
                                    done1 = []
                                    count_in = 0
                                    for index_pos, pos in enumerate(
                                            positions2):
                                        if index_pos >= len(positions2) - count_in:
                                            break

                                        count_ins = float_pos.count(float(pos) + 0.6)
                                        if count_ins >= len(positions2) - index_pos:
                                            count_real_ins = len(positions2) - index_pos
                                        else:
                                            count_real_ins = count_ins

                                        count_in = count_in + count_real_ins
                                        nr2.append(count_real_ins)
                                        insdels2.append(count_real_ins)
                                        nr2.append(float_pos.count(float(pos)))
                                        try:
                                            insdels2.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                            insdels2.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                        except IndexError:
                                            pass
                                    if insdels2.count(0) == len(insdels2):
                                        done1.append('.')
                                    else:
                                        if sum(nr2) <= 0.5 * len(positions2) and sum(nr2) <= 0.5 * (
                                                len(refseqwt) - (int(float(mut2[-1][1:-1]) + 0.6) - 1)):
                                            done1.append('.')
                                    if len(done1) == 1:
                                        count_in = 0
                                        for index_pos, pos in enumerate(
                                                positions1):
                                            if index_pos >= len(positions1) - count_in:
                                                break

                                            count_ins = float_pos.count(float(pos) + 0.6)
                                            if count_ins >= len(positions1) - index_pos:
                                                count_real_ins = len(positions1) - index_pos
                                            else:
                                                count_real_ins = count_ins

                                            count_in = count_in + count_real_ins
                                            nr1.append(count_real_ins)
                                            insdels1.append(count_real_ins)
                                            nr1.append(float_pos.count(float(pos)))
                                            try:
                                                insdels1.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                                insdels1.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                            except IndexError:
                                                pass
                                        if insdels1.count(0) == len(insdels1):
                                            high_qual_mutations.append(jmut2)
                                        else:
                                            if sum(nr1) <= 0.5 * len(positions1) and sum(nr1) <= 0.5 * int(
                                                    round(float(mut2[0][1:-1]))):
                                                high_qual_mutations.append(jmut2)

                        if joined_mutlist[-1].count('-') == 0 and joined_mutlist[-1].count('+') == 0 and joined_mutlist[
                            -1].count('X') == 0:
                            if int(float(mut_list[-1][-1][1:-1])) < len(
                                    refseqwt) - 2:
                                positions1 = list(range(int(round(float(mut_list[-1][0][1:-1]))) - 10,
                                          int(round(float(mut_list[-1][0][1:-1])))))[::-1]
                                nr = []
                                insdels = []
                                count_in = 0
                                for index_pos, pos in enumerate(positions1):
                                    if index_pos >= len(positions1) - count_in:
                                        break

                                    count_ins = float_pos.count(float(pos) + 0.6)
                                    if count_ins >= len(positions1) - index_pos:
                                        count_real_ins = len(positions1) - index_pos
                                    else:
                                        count_real_ins = count_ins

                                    count_in = count_in + count_real_ins
                                    nr.append(count_real_ins)
                                    insdels.append(count_real_ins)
                                    nr.append(float_pos.count(float(pos)))

                                    try:
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                        insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                    except IndexError:
                                        pass

                                if insdels.count(0) == len(insdels):
                                    high_qual_mutations.append(joined_mutlist[-1])
                                else:
                                    if sum(nr) <= 0.5 * len(positions1) and sum(nr) <= 0.5 * int(
                                            round(float(mut_list[-1][0][1:-1]))):
                                        high_qual_mutations.append(joined_mutlist[-1])

                            else:
                                if int(round(float(mut_list[-1][0][1:-1]))) - int(float(
                                        mut_list[-2][-1][1:-1])) <= 3:
                                    if joined_mutlist[-2].count('-') == 0 and joined_mutlist[-2].count(
                                            '+') == 0 and joined_mutlist[-2].count('X') == 0:
                                        positions1 = list(
                                            range(int(round(float(mut_list[-1][0][1:-1]))) - 10,
                                                  int(round(float(mut_list[-1][0][1:-1])))))[::-1]
                                        nr = []
                                        insdels = []
                                        count_in = 0
                                        for index_pos, pos in enumerate(
                                                positions1):
                                            if index_pos >= len(positions1) - count_in:
                                                break

                                            count_ins = float_pos.count(float(pos) + 0.6)
                                            if count_ins >= len(positions1) - index_pos:
                                                count_real_ins = len(positions1) - index_pos
                                            else:
                                                count_real_ins = count_ins

                                            count_in = count_in + count_real_ins
                                            nr.append(count_real_ins)
                                            insdels.append(count_real_ins)
                                            nr.append(float_pos.count(float(pos)))
                                            try:
                                                insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                                insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                            except IndexError:
                                                pass
                                        if insdels.count(0) == len(insdels):
                                            high_qual_mutations.append(joined_mutlist[-1])
                                        else:
                                            if sum(nr) <= 0.5 * len(positions1) and sum(nr) <= 0.5 * int(
                                                    round(float(mut_list[-1][0][1:-1]))):
                                                high_qual_mutations.append(joined_mutlist[-1])
                                else:
                                    positions1 = list(
                                        range(int(round(float(mut_list[-1][0][1:-1]))) - 10,
                                              int(round(float(mut_list[-1][0][1:-1])))))[::-1]
                                    nr = []
                                    insdels = []
                                    count_in = 0
                                    for index_pos, pos in enumerate(
                                            positions1):
                                        if index_pos >= len(positions1) - count_in:
                                            break

                                        count_ins = float_pos.count(float(pos) + 0.6)
                                        if count_ins >= len(positions1) - index_pos:
                                            count_real_ins = len(positions1) - index_pos
                                        else:
                                            count_real_ins = count_ins

                                        count_in = count_in + count_real_ins
                                        nr.append(count_real_ins)
                                        insdels.append(count_real_ins)
                                        nr.append(float_pos.count(float(pos)))
                                        try:
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + '-'))
                                            insdels.append(mut_string.count(refseqwt[pos - 1] + str(pos) + 'X'))
                                        except IndexError:
                                            pass
                                    if insdels.count(0) == len(insdels):
                                        high_qual_mutations.append(joined_mutlist[-1])
                                    else:
                                        if sum(nr) <= 0.5 * len(positions1) and sum(nr) <= 0.5 * int(
                                                round(float(mut_list[-1][0][1:-1]))):
                                            high_qual_mutations.append(joined_mutlist[-1])

                # remove mutations to 'X' and '*' from list
                final_mutations = []
                for final_mut in ', '.join(high_qual_mutations).split(', '):
                    if final_mut == '':
                        final_mutations.append(final_mut)
                    elif final_mut[-1] != 'X' and final_mut[-1] != '*' and final_mut[0] != '*':
                        final_mutations.append(final_mut)

                for new_m in final_mutations:
                    if isinstance(new_m, float) is False:
                        if new_m.count('*') != 0:
                            print('Caution: There is still an asterisk')

                high_qual_muts.append(', '.join(final_mutations))

            df_complete['high_quality_exchanges'] = high_qual_muts

            df_complete.to_csv(df_complete_out)


        find_high_quality_amino_acid_exchanges()

process_protein_sequences()

