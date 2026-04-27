#!/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import gzip
import multiprocessing
from tkinter import END
import pandas as pd
import numpy as np
from pandarallel import pandarallel
import concurrent.futures
import multiprocessing
from sklearn.metrics import pairwise_distances, accuracy_score
import seaborn as sns
import matplotlib.pyplot as plt
import datetime
import scipy.cluster.hierarchy as shc
import scipy.spatial.distance as ssd  # pdist


logger = logging.getLogger()


END_FORMATTING = "\033[0m"
WHITE_BG = "\033[0;30;47m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
RED = "\033[31m"
GREEN = "\033[32m"
MAGENTA = "\033[35m"
BLUE = "\033[34m"
CYAN = "\033[36m"
YELLOW = "\033[93m"
DIM = "\033[2m"


def get_arguments():

    parser = argparse.ArgumentParser(
        prog="compare_snp_prokaion.py", description="Pipeline to compare variants (SNVs) with any non model organism. Specialised in Mycobacterium tuberculosis")

    parser.add_argument('-i', '--input', dest="input_dir", metavar="input_directory",
                        type=str, required=False, help='REQUIRED. Input directory containing all vcf files')

    parser.add_argument('-s', '--sample_list', default=False, required=False,
                        help='File with sample names to analyse instead of all samples')

    parser.add_argument('-d', '--distance', default=50, required=False,
                        help='Minimun distance to cluster groups after comparison')

    parser.add_argument('-c', '--only-compare', dest="only_compare", required=False,
                        default=False, help='Add already calculated snp binary matrix')

    parser.add_argument('-r', '--recalibrate', required=False,
                        type=str, default=False, help='Coverage folder')

    parser.add_argument('-b', '--bam_folder', required=False,
                        type=str, default=False, help='Bam folder')

    parser.add_argument('-w', '--window', required=False, type=int,
                        default=2, help='Number of snps in 10 to discard. Default: 2')

    parser.add_argument('-C', '--complex', required=False,
                        action='store_true', help='Remove complex positions')

    parser.add_argument('-B', '--remove_bed', required=False, type=str,
                        default=False, help='BED file with positions to remove')

    parser.add_argument('-P', '--extract_bed', required=False, type=str,
                        default=False, help='BED file with important positions or genes to annotate')

    parser.add_argument('-R', '--reference', required=False, type=str, default=False,
                        help='Reference fasta file used in original variant calling')

    parser.add_argument("--min_threshold_discard_uncov_sample", required=False, type=float,
                        default=0.6, help="Minimum uncovered genome to discard a sample. Default: 0.6")

    parser.add_argument("--min_threshold_discard_uncov_pos", required=False, type=float,
                        default=0.5, help="Minimum covered position to discard it. Default: 0.5")

    parser.add_argument("--min_threshold_discard_htz_sample", required=False, type=float,
                        default=0.6, help="Minimum heterozygosity to discard a sample. Default: 0.6")

    parser.add_argument("--min_threshold_discard_htz_pos", required=False, type=float,
                        default=0.5, help="Minimum heterozygosity to discard a position. Default: 0.5")

    parser.add_argument("--min_threshold_discard_all_sample", required=False, type=float,
                        default=0.6, help="Minimum inaccuracies to discard a sample. Default: 0.6")

    parser.add_argument("--min_threshold_discard_all_pos", required=False, type=float,
                        default=0.5, help="Minimum inaccuracies to discard a position. Default: 0.5")

    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Name of all the output files, might include path')

    parser.add_argument('-S', '--only_snp', required=False,
                        action='store_true', help='Create the results only with SNPs, removing INDELs')

    arguments = parser.parse_args()

    return arguments


def check_create_dir(path):

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)


def check_file_exists(file_name):
    """
        Check file exist and is not 0 Kb, if not program exit.
    """

    # Retrieve the file info to check if has size > 0
    file_info = os.stat(file_name)

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        logger.info(RED + BOLD + "File: %s not found or empty\n" %
                    file_name + END_FORMATTING)
        sys.exit(1)

    return os.path.isfile(file_name)


def import_tsv_variants(tsv_file, sample, min_total_depth=8, min_alt_dp=8, only_snp=True):

    input_file = os.path.abspath(tsv_file)

    df = pd.read_csv(input_file, sep="\t")
    df = df.drop_duplicates(subset=["POS", "REF", "ALT"], keep="first")

    df = df[((df.TOTAL_DP >= min_total_depth) & (df.ALT_DP >= min_alt_dp))]
    df = df[["REGION", "POS", "REF", "ALT", "ALT_FREQ", "TYPE"]]
    df = df.rename(columns={"ALT_FREQ": sample})

    if only_snp == True:
        df = df[df.TYPE == "snp"]
        df = df.drop(["TYPE"], axis=1)
        return df
    else:
        df = df.drop(["TYPE"], axis=1)
        return df


def import_to_pandas(file_table, header=False, sep='\t'):

    if header == False:
        # Exclude first line, exclusive for vcf
        dataframe = pd.read_csv(file_table, sep=sep, skiprows=[0], header=None)
    else:
        # Use first line as header
        dataframe = pd.read_csv(file_table, sep=sep, header=0)

    return dataframe


def extract_lowfreq(tsv_file, sample, min_total_depth=8, min_alt_dp=8, min_freq_include=0.7, only_snp=True):

    input_file = os.path.abspath(tsv_file)

    df = pd.read_csv(input_file, sep='\t')
    df = df.drop_duplicates(subset=['POS', 'REF', 'ALT'], keep='first')

    df = df[(df.ALT_DP <= min_alt_dp) & (df.ALT_FREQ >= min_freq_include)]

    df = df[['REGION', 'POS', 'REF', 'ALT', 'ALT_FREQ', 'TYPE']]
    df['ALT_FREQ'] = '?'
    df = df.rename(columns={'ALT_FREQ': sample})

    if only_snp == True:
        df = df[df.TYPE == 'snp']
        df = df.drop(['TYPE'], axis=1)
        return df
    else:
        df = df.drop(['TYPE'], axis=1)
        return df


def extract_uncovered(cov_file):

    base_file = os.path.basename(cov_file)
    input_file = os.path.abspath(cov_file)
    sample = base_file.split('.')[0]

    df = pd.read_csv(input_file, sep='\t', header=None)
    df.columns = ['REGION', 'POS', sample]
    df = df[df[sample] == 0]
    df = df.replace(0, '!')
    return df


def ddbb_create_intermediate(variant_dir, coverage_dir, min_freq_discard=0.1, min_alt_dp=8, only_snp=False, samples=False):

    df = pd.DataFrame(columns=["REGION", "POS", "REF", "ALT"])

    # Merge all raw
    for root, _, files in os.walk(variant_dir):
        for name in files:
            if name == "snps.all.ivar.tsv":
                sample = root.split("/")[-1]
                if samples == False:
                    logger.debug("Adding: " + sample)
                    filename = os.path.join(root, name)
                    dfv = import_tsv_variants(
                        filename, sample, min_total_depth=8, min_alt_dp=8, only_snp=only_snp)
                    df = df.merge(dfv, how="outer")
                else:
                    if sample in samples:
                        logger.debug("Adding: " + sample)
                        filename = os.path.join(root, name)
                        dfv = import_tsv_variants(
                            filename, sample, min_total_depth=8, min_alt_dp=8, only_snp=only_snp)
                        if dfv.shape[0] > 0:
                            df = df.merge(dfv, how="outer")
                        else:
                            logger.debug(sample + " HAS NO VARIANTS")
                            sample.remove(sample)

    # Round frequencies
    df = df.round(2)

    # Remove <= 0.1 (parameter in function)
    # IF HANDLE HETEROZYGOUS change this 0 for X or 0.5
    def handle_lowfreq(x): return None if x <= min_freq_discard else x
    df.iloc[:, 4:] = df.iloc[:, 4:].applymap(handle_lowfreq)

    # Drop all NaN rows
    df['AllNaN'] = df.apply(lambda x: x[4:].isnull().values.all(), axis=1)
    df = df[df.AllNaN == False]
    df = df.drop(['AllNaN'], axis=1).reset_index(drop=True)

    # Include poorly covered
    for root, _, files in os.walk(variant_dir):
        for name in files:
            if name == 'snps.all.ivar.tsv':
                filename = os.path.join(root, name)
                sample = root.split('/')[-1]
                if samples == False:
                    logger.debug('Adding lowfreqs: ' + sample)
                    dfl = extract_lowfreq(
                        filename, sample, min_total_depth=8, min_alt_dp=min_alt_dp, only_snp=only_snp)
                    df[sample].update(df[['REGION', 'POS', 'REF', 'ALT']].merge(
                        dfl, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])
                else:
                    if sample in samples:
                        logger.debug('Adding lowfreqs: ' + sample)
                        dfl = extract_lowfreq(
                            filename, sample, min_total_depth=8, min_alt_dp=min_alt_dp, only_snp=only_snp)
                        df[sample].update(df[['REGION', 'POS', 'REF', 'ALT']].merge(
                            dfl, on=['REGION', 'POS', 'REF', 'ALT'], how='left')[sample])

    indel_positions = df[(df['REF'].str.len() > 1) | (
        df['ALT'].str.len() > 1)].POS.tolist()
    indel_len = df[(df['REF'].str.len() > 1) | (
        df['ALT'].str.len() > 1)].REF.tolist()
    indel_len = [len(x) for x in indel_len]

    indel_positions_final = []

    for position, del_len in zip(indel_positions, indel_len):
        indel_positions_final = indel_positions_final + \
            [x for x in range(position - (5 + del_len),
                              position + (5 + del_len))]

    # Include uncovered
    samples_coverage = df.columns.tolist()[4:]
    for root, _, files in os.walk(coverage_dir):
        for name in files:
            if name.endswith('.cov'):
                filename = os.path.join(root, name)
                sample = name.split('.')[0]
                if sample in df.columns[4:]:
                    samples_coverage.remove(sample)
                    logger.debug('Adding uncovered: ' + sample)
                    dfc = extract_uncovered(filename)
                    dfc = dfc[~dfc.POS.isin(indel_positions_final)]
                    df[sample].update(df[['REGION', 'POS']].merge(
                        dfc, on=['REGION', 'POS'], how='left')[sample])

    if len(samples_coverage) > 0:
        logger.info('WARNING: ' + (',').join(samples_coverage) +
                    ' coverage file not found')

    # Asign 0 to rest (Absent)
    df = df.fillna(0)

    # Determine N (will help in poorly covered determination)

    def extract_sample_count(row):
        count_list = [i not in ['!', 0, '0'] for i in row[4:]]
        samples = np.array(df.columns[4:])
        # sample[np.array(count_list)] # Filter array with True/False array

        return (sum(count_list), (',').join(samples[np.array(count_list)]))

    if 'N' in df.columns:
        df = df.drop(['N', 'Samples'], axis=1)
    if 'Position' in df.columns:
        df = df.drop('Position', axis=1)

    df[['N', 'Samples']] = df.apply(
        extract_sample_count, axis=1, result_type='expand')

    df = df[df.N > 0]

    df['Position'] = df.apply(lambda x: ('|').join(
        [str(x['REGION']), str(x['REF']), str(x['POS']), str(x['ALT'])]), axis=1)

    df = df.drop(['REGION', 'REF', 'POS', 'ALT'], axis=1)

    df = df[['Position', 'N', 'Samples'] +
            [col for col in df.columns if col not in ['Position', 'N', 'Samples']]]

    return df


def bed_to_df(bed_file):
    """
    Import bed file separated by tabs into a pandas df
    -Handle header line
    -Handle with and without description (If there is no description adds true or false to annotated df)
    """
    header_lines = 0
    # Handle likely header by checking colums 2 and 3 as numbers
    with open(bed_file, 'r') as f:
        next_line = f.readline().strip()
        line_split = next_line.split(None)  # This split by any blank character
        start = line_split[1]
        end = line_split[2]
        while not start.isdigit() and not end.isdigit():
            header_lines = header_lines + 1
            next_line = f.readline().strip()
            # This split by any blank character
            line_split = next_line.split(None)
            start = line_split[1]
            end = line_split[2]

    if header_lines == 0:
        # delim_whitespace=True
        df = pd.read_csv(bed_file, sep="\t", header=None)
    else:
        df = pd.read_csv(bed_file, sep="\t", skiprows=header_lines,
                         header=None)  # delim_whitespace=True

    df = df.iloc[:, 0:4]
    df.columns = ["#CHROM", "start", "end", "description"]

    return df


def remove_bed_positions(df, bed_file, path_compare):
    bed_df = bed_to_df(bed_file)
    filtered_position = path_compare + ".filter_position.tsv"

    counter = 0

    for _, row in df.iterrows():
        position_number = int(row.Position.split("|")[2])
        if any(start <= position_number <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())):
            logger.info('Position: {} removed found in {}'.format(
                row.Position, bed_file))
            if counter == 0:
                df_filter = df[df.Position == row.Position]
            else:
                df_aux = df[df.Position == row.Position]
                df_filter = pd.concat([df_filter, df_aux], ignore_index=True)
            df = df[df.Position != row.Position]
            counter += 1

    if counter:
        df_filter.to_csv(filtered_position, sep="\t", index=False)

    return df


def extract_bed_positions(df, bed_file, path_compare):
    bed_df = bed_to_df(bed_file)
    filtered_position = path_compare + ".annotated_position.tsv"

    counter = 0

    for _, row in df.iterrows():
        position_number = int(row.Position.split("|")[2])
        if any(start <= position_number <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())):
            logger.info('Position: {} annotated found in {}'.format(
                row.Position, bed_file))
            if counter == 0:
                df_filter = df[df.Position == row.Position]
            else:
                df_aux = df[df.Position == row.Position]
                df_filter = pd.concat([df_filter, df_aux], ignore_index=True)
            counter += 1

    if counter:
        df_filter.to_csv(filtered_position, sep="\t", index=False)


def recheck_variant_rawvcf_intermediate(row, positions, alt_snps, variant_dir, min_cov_low_freq=10):
    """
    CU458896.1	3068036	.	G	A	262.784	.	AB=0.8;ABP=14.7363;AC=1;AF=0.5;AN=2;AO=12;CIGAR=1X;DP=15;DPB=15;DPRA=0;EPP=9.52472;EPPR=3.73412;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=1;NUMALT=1;ODDS=0.121453;PAIRED=0.333333;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=410;QR=108;RO=3;RPL=7;RPP=3.73412;RPPR=9.52472;RPR=5;RUN=1;SAF=4;SAP=5.9056;SAR=8;SRF=1;SRP=3.73412;SRR=2;TYPE=snp	GT:DP:AD:RO:QR:AO:QA:GL	0/1:15:3,12:3:108:12:410:-32.709,0,-5.55972
    """

    sample = row.name
    row_values = row.values
    return_list = [0] * len(positions)

    # Identify vcf
    variant_sample_folder = os.path.join(variant_dir, sample)
    filename = os.path.join(variant_sample_folder, "snps.raw.vcf")

    # Open file and retrieve output

    checked_positions = positions[:]

    with open(filename, 'r') as f:
        for line in f:
            # Parse line/variant in raw VCF
            if line.startswith('#CHROM'):
                headers = line.split('\t')
                # vcf_sample = headers[-1].strip()
            elif not line.startswith('#'):
                line_split = line.split('\t')
                vcf_position = line_split[1]
                params = line_split[-2].split(':')
                value_params = line_split[-1].split(':')
                depth_index = params.index('DP')

                # vcf_reference = line_split[0]
                vcf_alt_base = line_split[4]
                vcf_depth = int(value_params[depth_index])
                # ref_depth_index = params.index('RO')
                alt_depth_index = params.index('AO')

                try:
                    vcf_alt_depth = int(value_params[alt_depth_index])
                except:
                    value_params = [
                        int(x) for x in value_params[alt_depth_index].split(',')]
                    vcf_alt_depth = max(value_params)
                    vcf_alt_base = vcf_alt_base.split(',')[-1]

                vcf_alt_freq = round(vcf_alt_depth/vcf_depth, 2)

                # Process variant information for recalibration
                if vcf_position in positions:
                    if vcf_position in checked_positions:
                        checked_positions.remove(vcf_position)
                    position_index = positions.index(vcf_position)

                    if str(row[position_index]) == '0':
                        # vcf_ref_base = line_split[3]
                        position = positions[position_index]
                        alt_snp = alt_snps[position_index]

                        logger.debug('REC==>ORI:SAMPLE: {}\nPOS: {},  ALT: {}, DP: {}, FREQ: {}\nPOS: {}, ALT: {}, DP: {}, FREQ: {}'.format(
                            sample, vcf_position, vcf_alt_base, vcf_depth, vcf_alt_freq, position, alt_snp, vcf_alt_depth, row[position_index]))

                        try:
                            logger.debug('VALUE:\n{}, AO: {}'.format(
                                (':').join(value_params), alt_depth_index))
                        except:
                            logger.debug('problem in VALUES POS: {}, sample: {}'.format(
                                vcf_position, sample))

                        if vcf_depth <= min_cov_low_freq and vcf_depth > 0:
                            logger.debug('Position: {} LOWDEPTH: {}, DP: {}'.format(
                                vcf_position, vcf_alt_freq, vcf_depth))
                            row[position_index] = '?'
                        else:
                            if vcf_depth >= min_cov_low_freq and alt_snp == vcf_alt_base and vcf_alt_freq >= 0.1:
                                logger.debug('Position: {} RECOVERED from 0 to: {} in sample {}'.format(
                                    vcf_position, vcf_alt_freq, sample))
                                row[position_index] = vcf_alt_freq

                            # Handle INDEL incosistency between alt positions
                            elif (len(alt_snp) > 1) and vcf_alt_freq >= 0.1:
                                logger.debug('Position INDEL: {} RECOVERED from 0 to: {} in sample {}'.format(
                                    vcf_position, vcf_alt_freq, sample))
                                row[position_index] = vcf_alt_freq
                            elif alt_snp in vcf_alt_base:
                                logger.debug('Position COMPLEX: {} RECOVERED from 0 in {} to: {} in {} sample {}'.format(
                                    vcf_position, alt_snp, vcf_alt_freq, vcf_alt_base, sample))
                                row[position_index] = vcf_alt_freq

                elif 'complex' in line:
                    vcf_position = int(vcf_position)
                    positions_complex = [x for x in range(
                        vcf_position - 4, vcf_position + 4)]
                    positions_complex = [str(x) for x in positions_complex]
                    intersection = set(positions_complex).intersection(
                        set(checked_positions))
                    intersection = list(intersection)

                    if len(intersection) > 0:
                        for i in intersection:
                            if i in checked_positions:
                                checked_positions.remove(i)
                            position_index = positions.index(i)
                            if str(row[position_index]) == '0' or str(row[position_index]) == '!':
                                if vcf_depth <= min_cov_low_freq and vcf_depth > 0:
                                    logger.debug('Position: {} LOWDEPTH PRECOMPLEX: {}, DP: {}, sample {}'.format(
                                        i, vcf_alt_freq, vcf_depth, sample))
                                    row[position_index] = '?'
                                else:
                                    logger.debug('INTERSECTION PRECOMPLEX: Position: {}, AF: {}, sample: {}'.format(
                                        i, vcf_alt_freq, sample))
                                    row[position_index] = vcf_alt_freq

    return row


def recalibrate_ddbb_vcf_intermediate(snp_matrix_ddbb_file, variant_dir, min_cov_low_freq=12):
    """
    https://github.com/nalepae/pandarallel
    """

    pandarallel.initialize()

    selected_columns = ['Position', 'N',
                        'Samples', 'CHROM', 'REF', 'POS', 'ALT']

    df = pd.read_csv(snp_matrix_ddbb_file, sep='\t')
    df[['CHROM', 'REF', 'POS', 'ALT']
       ] = df.Position.str.split('|', expand=True)

    df = df[selected_columns +
            [col for col in df.columns if col not in selected_columns]]

    dft = df.transpose()

    samples = [str(x) for x in dft.index if x not in selected_columns]
    POSs = dft.loc['POS', :].tolist()
    # REFs = dft.loc['REF', :].tolist()
    ALTs = dft.loc['ALT', :].tolist()

    # n_samples = len(samples)

    header_df = dft.loc[selected_columns, :]
    samples_df = dft.drop(selected_columns, axis=0)

    samples_df = samples_df.parallel_apply(lambda x: recheck_variant_rawvcf_intermediate(
        x, POSs, ALTs, variant_dir, min_cov_low_freq=min_cov_low_freq), axis=1)

    final_dft = pd.concat([header_df, samples_df], ignore_index=True)

    final_df = final_dft.transpose()
    final_df.columns = selected_columns + samples
    final_df = final_df.drop(['CHROM', 'REF', 'POS', 'ALT'], axis=1)

    def extract_sample_count(row):
        count_list = [i not in ['!', 0, '0'] for i in row[3:]]
        selected_samples = [sample for (sample, true_false) in zip(
            samples, count_list) if true_false]
        return (sum(count_list), (',').join(selected_samples))

    final_df[['N', 'Samples']] = final_df.apply(
        extract_sample_count, axis=1, result_type='expand')

    return final_df


def remove_position_range(df):

    INs = df[df['Position'].str.contains(
        r'\|[ATCG]{3,}\|[0-9]{1,}', regex=True)]

    bed_df = pd.DataFrame(columns=['#CHROM', 'REF', 'start', 'ALT'])
    bed_df[['#CHROM', 'REF', 'start', 'ALT']
           ] = INs['Position'].str.split('|', expand=True)
    bed_df['start'] = bed_df['start'].astype(int)

    # Include only positions in +1 position
    bed_df['start'] = bed_df['start'] + 1
    bed_df['lenREF'] = bed_df['REF'].str.len()

    # Commented to avoid collision with insertions
    bed_df['lenALT'] = bed_df['ALT'].str.len()
    bed_df['end'] = bed_df['start'] + bed_df['lenREF'] - 1

    for _, row in df.iterrows():
        position_number = int(row.Position.split('|')[2])
        if any(start <= position_number <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())) and not re.search('\|[ATCG]{2,}', row.Position):
            logger.debug(
                'Position: {} removed found in INDEL'.format(row.Position))
            df = df[df.Position != row.Position]

    return df


def extract_complex_list(variant_dir, samples=False):

    all_complex = []

    for root, _, files in os.walk(variant_dir):
        for name in files:
            if name == 'snps.all.ivar.tsv':
                sample = root.split('/')[-1]
                filename = os.path.join(root, name)
                if samples == False:
                    logger.debug('Obtaining complex in: ' + sample)
                    df = pd.read_csv(filename, sep='\t')
                    sample_complex = df[~df.OLDVAR.isna()]['POS'].tolist()
                    sample_complex = [int(x) for x in sample_complex]
                    all_complex = all_complex + sample_complex
                else:
                    if sample in samples:
                        logger.debug('Obtaining complex in: ' + sample)
                        df = pd.read_csv(filename, sep='\t')
                        sample_complex = df[~df.OLDVAR.isna()]['POS'].tolist()
                        sample_complex = [int(x) for x in sample_complex]
                        all_complex = all_complex + sample_complex

    return sorted(set(all_complex))


def add_window_distance(vcf_df, window_size=10):
    """
    Add a column indicating the maximum number of SNPs in a windows of 10 or supplied distance
    """

    list_pos = vcf_df.POS.to_list()  # All positions
    set_pos = set(list_pos)  # To set for later comparing
    # Max to iter over positions (independent from reference)
    max_pos = max(vcf_df.POS.to_list())
    all_list = list(range(1, max_pos + 1))  # Create a list to slide one by one

    df_header = 'window_' + str(window_size)
    vcf_df[df_header] = 1  # Create all 1 by default

    # Slide over windows
    for i in range(0, max_pos, 1):
        # This splits the list in windows of determined length
        window_pos = all_list[i:i+window_size]
        set_window_pos = set(window_pos)
        # How many known positions are in every window for later clasification
        num_conglomerate = set_pos & set_window_pos

        if len(num_conglomerate) > 1:
            for i in num_conglomerate:
                # Retrieve index with the known position
                index = vcf_df.index[vcf_df['POS'] == i][0]
                if vcf_df.loc[index, df_header] < len(num_conglomerate):
                    vcf_df.loc[index, df_header] = len(num_conglomerate)


def revised_df(df, out_dir=False, complex_pos=False, min_freq_include=0.7, min_threshold_discard_uncov_sample=0.6, min_threshold_discard_uncov_pos=0.5, min_threshold_discard_htz_sample=0.6, min_threshold_discard_htz_pos=0.5, min_threshold_discard_all_pos=0.5, min_threshold_discard_all_sample=0.6, remove_faulty=True, drop_samples=True, drop_positions=True, windows_size_discard=2):

    if remove_faulty == True:

        uncovered_positions = df.iloc[:, 3:].apply(lambda x:  sum(
            [i in ['!', '?'] for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=1)
        heterozygous_positions = df.iloc[:, 3:].apply(lambda x: sum([(i not in ['!', '?', 0, 1, '0', '1']) and (float(
            i) < min_freq_include) and (float(i) > 0.1) for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=1)

        report_position = pd.DataFrame({'Position': df.Position, 'uncov_fract': uncovered_positions,
                                        'htz_frac': heterozygous_positions, 'faulty_frac': uncovered_positions + heterozygous_positions})
        faulty_positions = report_position['Position'][(report_position.uncov_fract >= min_threshold_discard_uncov_pos) | (
            report_position.htz_frac >= min_threshold_discard_htz_pos) | (report_position.faulty_frac >= min_threshold_discard_all_pos)].tolist()

        uncovered_samples = df.iloc[:, 3:].apply(lambda x: sum(
            [i in ['!'] for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=0)
        heterozygous_samples = df.iloc[:, 3:].apply(lambda x: sum([(i not in ['!', '?', 0, 1, '0', '1']) and (float(
            i) < min_freq_include) and (float(i) > 0.1) for i in x.values])/sum([(i not in [0, '0']) for i in x.values]), axis=0)
        report_samples = pd.DataFrame({'sample': df.iloc[:, 3:].columns, 'uncov_fract': uncovered_samples,
                                       'htz_frac': heterozygous_samples, 'faulty_frac': uncovered_samples + heterozygous_samples})
        faulty_samples = report_samples['sample'][(report_samples.uncov_fract >= min_threshold_discard_uncov_sample) | (
            report_samples.htz_frac >= min_threshold_discard_htz_sample) | (report_samples.faulty_frac >= min_threshold_discard_all_sample)].tolist()

        # Calculate close SNPs/INDELs and remove those with 2 or more mutations in 10bp
        df['POS'] = df.apply(lambda x: x.Position.split('|')[2], axis=1)
        df['POS'] = df['POS'].astype(int)
        df = df.sort_values('POS')
        add_window_distance(df)

        if out_dir:
            out_dir = os.path.abspath(out_dir)
            report_samples_windows = os.path.join(
                out_dir, 'report_windows.tsv')
            report_samples_file = os.path.join(out_dir, 'report_samples.tsv')
            report_faulty_samples_file = os.path.join(
                out_dir, 'faulty_samples.tsv')
            report_positions_file = os.path.join(
                out_dir, 'report_positions.tsv')
            report_faulty_positions_file = os.path.join(
                out_dir, 'faulty_positions.tsv')
            intermediate_cleaned_file = os.path.join(
                out_dir, 'intermediate.highfreq.tsv')

            report_position.to_csv(report_positions_file,
                                   sep='\t', index=False)
            report_samples.to_csv(report_samples_file, sep='\t', index=False)

            with open(report_faulty_samples_file, 'w+') as f:
                f.write(('\n').join(faulty_samples))
            with open(report_faulty_positions_file, 'w+') as f2:
                f2.write(('\n').join(faulty_positions))

            df.to_csv(report_samples_windows, sep='\t')

        clustered_positions = df['POS'][df.window_10 >
                                        windows_size_discard].tolist()

        if drop_positions == True:
            df = df[~df.Position.isin(faulty_positions)]
        if drop_samples == True:
            df = df.drop(faulty_samples, axis=1)

        logger.debug('FAULTY POSITIONS:\n{}\n\nFAULTY SAMPLES:\n{}'.format(
            ("\n").join(faulty_positions), ("\n").join(faulty_samples)))

    if len(clustered_positions) == 0:
        clustered_positions = [0]
    logger.debug('CLUSTERED POSITIONS' + '\n' +
                 (',').join([str(x) for x in clustered_positions]))

    if complex_pos:
        logger.debug('COMPLEX POSITIONS' + "\n" +
                     (',').join([str(x) for x in complex_pos]))
        clustered_positions = list(set(clustered_positions + complex_pos))

    logger.debug('ALL CLOSE POSITIONS' + "\n" +
                 (',').join([str(x) for x in clustered_positions]))

    # Remove close mutations and complex positions
    df = df[~df.POS.isin(clustered_positions)]

    # Remove complex variants in_freq_include
    df['valid'] = df.apply(lambda x: sum(
        [i != '?' and i != '!' and float(i) >= min_freq_include for i in x[3:]]), axis=1)
    df = df[df.valid >= 1]
    df = df.drop('valid', axis=1)
    df = df.drop('window_10', axis=1)
    df = df.drop('POS', axis=1)

    if out_dir != False:
        df.to_csv(intermediate_cleaned_file, sep='\t', index=False)

    df = df.replace('!', 0)
    df = df.replace('?', 1)
    df.iloc[:, 3:] = df.iloc[:, 3:].astype(float)

    # Replace Htz to 0
    # f = lambda x: 1 if x >= min_freq_include else 0 # IF HANDLE HETEROZYGOUS CHANGE THIS 0 for X or 0.5
    # df.iloc[:,3:] = df.iloc[:,3:].applymap(f)

    # Replace Htz with 1
    # IF HANDLE HETEROZYGOUS CHANGE THIS 0 for X or 0.5
    def fn(x): return 1 if x > 0.5 else 0
    df.iloc[:, 3:] = df.iloc[:, 3:].applymap(fn)
    df.N = df.apply(lambda x: sum(x[3:]), axis=1)

    def extract_sample_name(row):
        count_list = [i not in ['!', 0, '0'] for i in row[3:]]
        samples = np.array(df.columns[3:])
        # samples[np.array(count_list)] # Filter array with True/False array
        return ((',').join(samples[np.array(count_list)]))

    df['Samples'] = df.apply(extract_sample_name, axis=1)

    # Remove positions with 0 samples after htz
    df = df[df.N > 0]

    return df


def extract_lowcov(intermediate_df):
    df = pd.read_csv(intermediate_df, sep='\t')

    # Replace Htz with 1
    # IF HANDLE HETEROZYGOUS CHANGE THIS 0 for X or 0.5
    def fn(x):
        if isinstance(x, (int, float)):
            return 1 if x > 0.5 else 0
        else:
            return x

    df.iloc[:, 3:] = df.iloc[:, 3:].applymap(lambda x: fn(float(x)) if isinstance(x, str) and x.replace('.', '', 1).isdigit() else x)

    # Crear un nuevo DataFrame para almacenar las filas filtradas
    df_symbol = pd.DataFrame(columns=df.columns)

    # Iterar sobre las filas del DataFrame
    for index, row in df.iterrows():
        # Filtrar los elementos diferentes de '!' y '?'
        filtered_row = [x for x in row[3:] if x not in ['!', '?']]
        # Verificar si todos los elementos de la fila (excluyendo '?' y '!') son iguales
        if len(set(filtered_row)) == 1:
            # Añadir la fila al nuevo DataFrame
            df_symbol = pd.concat([df_symbol, pd.DataFrame([row])], ignore_index=True)

    df_symbol = df_symbol[df_symbol.apply(lambda row: ('!' in row.values) or ('?' in row.values), axis=1)]

    return df_symbol


def extract_only_snps(revised_df):
    df = pd.read_csv(revised_df, sep='\t')
    # df['ref'] = df['Position'].apply(lambda x: x.split('|')[1])
    # df['alt'] = df['Position'].apply(lambda x: x.split('|')[3])
    # df_snps = df[(df['ref'].str.len() == 1) & (df['alt'].str.len() == 1)]
    df_snps = df[df['Position'].str.split('|', expand=True)[
        [1, 3]].applymap(len).eq(1).all(axis=1)]
    return df_snps


def dendogram_dataframe(dataframe, output_file):

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index

    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T,
                    method='average')  # method='single'

    plt.rcParams['lines.linewidth'] = 8  # Dendrogram line with
    plt.rcParams['xtick.major.size'] = 10  # Only affect to tick (line) size
    plt.rcParams.update({'font.size': 30})  # Increase x tick label size
    # plt.tick_params(labelsize=30)
    plt.figure(figsize=(30, 50))
    plt.ylabel('samples', fontsize=30)
    plt.xlabel('snp distance', fontsize=30)

    shc.dendrogram(Z, labels=labelList, orientation='left', distance_sort='descending',
                   show_leaf_counts=True, color_threshold=10, leaf_font_size=20)

    plt.savefig(output_file, format="png")


def linkage_to_newick(dataframe, output_file):
    """
    Thanks to https://github.com/biocore/scikit-bio/issues/1579
    Input :  Z = linkage matrix, labels = leaf labels
    Output:  Newick formatted tree string
    """

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index

    labelList = dataframe_only_samples.columns.tolist()
    Z = shc.linkage(dataframe_only_samples.T, method='average')

    tree = shc.to_tree(Z, False)

    def buildNewick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            # logger.info("%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick))
            return "%s:%f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0:
                newick = f"):{(parentdist - node.dist)/2}{newick}"
            else:
                newick = ");"
            newick = buildNewick(node.get_left(), newick,
                                 node.dist, leaf_names)
            newick = buildNewick(node.get_right(), ",%s" %
                                 (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            # logger.info(newick)
            return newick

    with open(output_file, 'w') as f:
        f.write(buildNewick(tree, "", tree.dist, labelList))
    return buildNewick(tree, "", tree.dist, labelList)


def matrix_to_rdf(snp_matrix, output_name):

    with open(output_name, 'w+') as fout:
        snp_number = snp_matrix.shape[0]
        first_line = "  ;1.0\n"
        # logger.info(first_line)
        fout.write(first_line)

        snp_list = snp_matrix.Position.tolist()
        snp_list = [x.split('|')[2] for x in snp_list]
        snp_list = " ;".join([str(x) for x in snp_list]) + " ;\n"
        # logger.info(snp_list)
        fout.write(snp_list)

        third_line = ("10;" * snp_number) + "\n"
        # logger.info(third_line)
        fout.write(third_line)

        transposed_snp_matrix = snp_matrix.T

        for index, row in transposed_snp_matrix.iloc[3:, :].iterrows():
            sample_header = ">" + index+";1;;;;;;;\n"
            # logger.info(sample_header)
            fout.write(sample_header)
            snp_row = "".join([str(x) for x in row.tolist()]) + "\n"
            # logger.info(snp_row)
            fout.write(snp_row)

        ref_header = ">REF;1;;;;;;;\n"
        # logger.info(ref_header)
        fout.write(ref_header)
        ref_snp = "0" * snp_number
        # logger.info(ref_snp)
        fout.write(ref_snp)


def matrix_to_common(snp_matrix, output_name):

    max_samples = max(snp_matrix.N.tolist())
    total_samples = len(snp_matrix.columns[3:])

    if max_samples == total_samples:
        with open(output_name, 'w+') as fout:
            common_snps = snp_matrix['Position'][snp_matrix.N == max_samples].astype(
                str).tolist()
            line = "\n".join(common_snps)
            fout.write("Position\n")
            fout.write(line)
    else:
        logger.info("No common SNPs were found")


def pairwise_to_cluster(pw, threshold=0):

    groups = {}
    columns = pw.columns.tolist()
    sorted_df = pw[(pw[columns[0]] != pw[columns[1]]) & (
        pw[columns[2]] <= threshold)].sort_values(by=[columns[2]])

    def rename_dict_clusters(cluster_dict):
        reordered_dict = {}
        for i, k in enumerate(list(cluster_dict)):
            reordered_dict[i] = cluster_dict[k]
        return reordered_dict

    def regroup_clusters(list_keys, groups_dict, both_samples_list):
        # Sum previous clusters
        list_keys.sort()
        new_cluster = sum([groups_dict[key] for key in list_keys], [])
        # Add new cluster
        cluster_asign = list(set(new_cluster + both_samples_list))
        # Remove duped cluster
        first_cluster = list_keys[0]
        groups_dict[first_cluster] = cluster_asign
        rest_cluster = list_keys[1:]
        for key in rest_cluster:
            del groups_dict[key]
        groups_dict = rename_dict_clusters(groups_dict)
        return groups_dict

    for _, row in sorted_df.iterrows():
        group_number = len(groups)
        sample_1 = str(row[0])
        sample_2 = str(row[1])
        both_samples_list = row[0:2].tolist()

        if group_number == 0:
            groups[group_number] = both_samples_list

        all_samples_dict = sum(groups.values(), [])

        if sample_1 in all_samples_dict or sample_2 in all_samples_dict:
            # Extract cluster which have the new samples
            key_with_sample = {key for (key, value) in groups.items() if (
                sample_1 in value or sample_2 in value)}

            cluster_with_sample = list(key_with_sample)
            cluster_with_sample_name = cluster_with_sample[0]
            number_of_shared_clusters = len(key_with_sample)

            if number_of_shared_clusters > 1:
                groups = regroup_clusters(
                    cluster_with_sample, groups, both_samples_list)
            else:
                groups[cluster_with_sample_name] = list(
                    set(groups[cluster_with_sample_name] + both_samples_list))

        else:
            groups[group_number] = both_samples_list

    for _, row in pw[(pw[pw.columns[0]] != pw[pw.columns[1]]) & (pw[pw.columns[2]] > threshold)].iterrows():
        sample_1 = str(row[0])
        sample_2 = str(row[1])
        all_samples_dict = sum(groups.values(), [])

        if sample_1 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_1]

        if sample_2 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_2]

    cluster_df = pd.DataFrame(groups.values(), index=list(groups))

    cluster_df_return = cluster_df.stack().droplevel(
        1).reset_index().rename(columns={'index': 'group', 0: 'id'})

    return cluster_df_return


def calculate_N(row):

    return len(row.samples)


def calculate_mean_distance(row, df):

    if row.N > 1:
        list_sample = row.samples
        list_sample = [str(x) for x in list_sample]
        list_sample = [x.split(".")[0] for x in list_sample]
        dataframe = df.loc[list_sample, list_sample]
        stacked_df = dataframe.stack()
        mean_distance = stacked_df.mean(skipna=True)
        min_distance = stacked_df.min(skipna=True)
        max_distance = stacked_df.max(skipna=True)
        return round(mean_distance, 2), min_distance, max_distance
    else:
        return 'NaN', 'NaN', 'NaN'


def matrix_to_cluster(pairwise_file, matrix_file, distance=0):

    output_dir = ('/').join(pairwise_file.split('/')[0:-1])

    logger.info('Reading Matrix')
    dfdist = pd.read_csv(matrix_file, index_col=0, sep='\t', )
    dfdist.columns = dfdist.columns.astype(str)
    dfdist.index = dfdist.index.astype(str)
    logger.info('Reading Pairwise')
    pairwise = pd.read_csv(pairwise_file, sep="\t", names=[
                           'sample_1', 'sample_2', 'dist'])
    logger.info('Creating Clusters')

    clusters = pairwise_to_cluster(pairwise, threshold=distance)

    cluster_summary = clusters.groupby('group')['id'].apply(
        list).reset_index(name='samples')
    cluster_summary['N'] = cluster_summary.apply(calculate_N, axis=1)
    cluster_summary = cluster_summary.sort_values(by=['N'], ascending=False)

    logger.info('Reseting group number by length')
    sorted_index = cluster_summary.index.to_list()
    sorted_index.sort()
    sorted_index = [x + 1 for x in sorted_index]
    cluster_summary['group'] = sorted_index
    cluster_summary = cluster_summary.sort_values(by=['N'], ascending=False)

    cluster_summary[['mean', 'min', 'max']] = cluster_summary.apply(
        lambda x: calculate_mean_distance(x, dfdist), axis=1, result_type="expand")

    final_cluster = cluster_summary[["group", "samples"]].explode(
        "samples").reset_index(drop=True)
    final_cluster = final_cluster.sort_values(by=['group'], ascending=True)

    final_cluster_file = os.path.join(
        output_dir, "group_table_" + str(distance) + ".tsv")
    cluster_summary_file = os.path.join(
        output_dir, "group_summary_" + str(distance) + ".tsv")

    cluster_summary.to_csv(cluster_summary_file, sep='\t', index=False)
    final_cluster.to_csv(final_cluster_file, sep='\t', index=False)


def snp_distance_matrix(dataframe, output_matrix, output_pairwise):

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(
        dataframe_only_samples.T, metric="hamming")  # dataframe.T means transposed
    snp_distance_df = pd.DataFrame(hamming_distance * len(dataframe_only_samples.index),
                                   index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns)  # Add index
    snp_distance_df = snp_distance_df.astype(int)
    pairwise = snp_distance_df.stack().reset_index(name='distance').rename(
        columns={'level_0': 'sample_1', 'level_1': 'sample_2'})

    snp_distance_df.to_csv(output_matrix, sep='\t', index=True)
    pairwise.to_csv(output_pairwise, sep='\t', header=False, index=False)


def hamming_distance_matrix(dataframe, output_file):

    dataframe_only_samples = dataframe.set_index(dataframe['Position']).drop(
        ['Position', 'N', 'Samples'], axis=1)  # Extract three first colums and use 'Position' as index
    hamming_distance = pairwise_distances(
        dataframe_only_samples.T, metric="hamming")  # dataframe.T means transposed
    hamming_distance_df = pd.DataFrame(
        hamming_distance, index=dataframe_only_samples.columns, columns=dataframe_only_samples.columns)  # Add index

    hamming_distance_df.to_csv(output_file, sep='\t', index=True)


def ddtb_compare(final_database, distance=0):

    database_file = os.path.abspath(final_database)
    check_file_exists(database_file)
    presence_ddbb = import_to_pandas(database_file, header=True)

    output_path = database_file.split('.')[0]

    logger.info('Output path is: ' + output_path)
    logger.info('\n' + BLUE + BOLD + 'Comparing all samples in ' +
                database_file + END_FORMATTING)

    # Calculate SNP distance for all and save file
    logger.info(CYAN + "SNP distance" + END_FORMATTING)
    snp_dist_file = output_path + ".snp.tsv"
    pairwise_file = output_path + ".snp.pairwise.tsv"
    snp_distance_matrix(presence_ddbb, snp_dist_file, pairwise_file)

    # Calculate hamming distance for all and save file
    logger.info(CYAN + "Hamming distance" + END_FORMATTING)
    hmm_dist_file = output_path + ".hamming.tsv"
    hamming_distance_matrix(presence_ddbb, hmm_dist_file)

    """
    # Represent pairwise snp distance for all and save file
    logger.info(CYAN + "Drawing distance" + END_FORMATTING)
    prior_represent = datetime.datetime.now()
    png_dist_file = output_path + ".snp.distance.png"
    # clustermap_dataframe(presence_ddbb, png_dist_file)
    after_represent = datetime.datetime.now()
    logger.info("Done with distance drawing in: %s" %
                (after_represent - prior_represent))
    """

    # Represent dendrogram snp distance for all and save file
    logger.info(CYAN + "Drawing dendrogram" + END_FORMATTING)
    png_dend_file = output_path + ".snp.dendrogram.png"
    dendogram_dataframe(presence_ddbb, png_dend_file)

    # Output a Newick file distance for all and save file
    logger.info(CYAN + "Newick dendrogram" + END_FORMATTING)
    newick_file = output_path + ".nwk"
    linkage_to_newick(presence_ddbb, newick_file)

    # Output a binary snp matrix distance in rdf format
    logger.info(CYAN + "rdf format" + END_FORMATTING)
    rdf_file = output_path + ".rdf"
    matrix_to_rdf(presence_ddbb, rdf_file)

    # Output a list of all common snps in group compared
    logger.info(CYAN + "Common SNPs" + END_FORMATTING)
    common_file = output_path + ".common.txt"
    matrix_to_common(presence_ddbb, common_file)

    # Output files with group/cluster assigned to samples
    logger.info(CYAN + "Assigning clusters" + END_FORMATTING)
    matrix_to_cluster(pairwise_file, snp_dist_file, distance=distance)


if __name__ == '__main__':

    args = get_arguments()

    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split('/')[-1]
    check_create_dir(output_dir)

    ##### LOGGING #####

    # Create log file with date and time

    right_now = str(datetime.datetime.now())
    right_now_full = '_'.join(right_now.split(' '))
    log_filename = group_name + '_' + right_now_full + '.log'
    log_folder = os.path.join(output_dir, 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info('########## COMPARE SNPs ##########')
    logger.info(args)

    group_compare = os.path.join(output_dir, group_name)
    compare_snp_matrix = group_compare + '.tsv'
    input_dir = os.path.abspath(args.input_dir)

    out_variant_dir = os.path.join(input_dir, 'Variants')
    out_stats_dir = os.path.join(input_dir, 'Stats')
    out_stats_coverage_dir = os.path.join(
        out_stats_dir, 'Coverage')  # Subfolder

    if args.sample_list:
        sample_file = os.path.abspath(args.sample_list)
        with open(sample_file, 'r') as f:
            content = f.read()
            sample_list = content.split('\n')
            sample_list = [x.strip() for x in sample_list]
            sample_list = [x for x in sample_list if x != '']

    if args.only_compare == False:

        ##### SNPs COMPARISON using tsv variant files #####

        logger.info('\n\n' + BLUE + BOLD + 'STARTING COMPARISON IN GROUP: ' +
                    group_name + END_FORMATTING + '\n')

        today = str(datetime.date.today())
        check_create_dir(output_dir)
        folder_compare = today + "_" + group_name
        path_compare = os.path.join(output_dir, folder_compare)
        check_create_dir(path_compare)
        full_path_compare = os.path.join(path_compare, group_name)

        compare_snp_matrix_recal = full_path_compare + ".revised.final.tsv"
        compare_snp_matrix_recal_intermediate = full_path_compare + ".revised_intermediate.tsv"
        compare_snp_matrix_recal_mpileup = full_path_compare + \
            ".revised_intermediate_vcf.tsv"
        compare_snp_matrix_INDEL_intermediate = full_path_compare + \
            ".revised_INDEL_intermediate.tsv"
        compare_only_snps = full_path_compare + "_ONLY_SNPs.revised.tsv"

        # Create intermediate

        prior = datetime.datetime.now()

        recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(
            out_variant_dir, out_stats_coverage_dir, min_freq_discard=0.2, min_alt_dp=12, only_snp=False, samples=sample_list)
        # recalibrated_snp_matrix_intermediate.to_csv(
        #     compare_snp_matrix_recal_intermediate, sep='\t', index=False)

        after = datetime.datetime.now()
        print(('\n' + "Done with function ddbb_create_intermediate in: %s" %
               (after - prior) + "\n"))

        # Remove SNPs from BED file (PE/PPE)

        if args.remove_bed:

            prior = datetime.datetime.now()

            recalibrated_snp_matrix_intermediate = remove_bed_positions(
                recalibrated_snp_matrix_intermediate, args.remove_bed, full_path_compare)

            after = datetime.datetime.now()
            print(('\n' + "Done with function remove_bed_positions in: %s" %
                   (after - prior) + "\n"))

        recalibrated_snp_matrix_intermediate.to_csv(
            compare_snp_matrix_recal_intermediate, sep="\t", index=False)

        # Recalibrate intermediate with VCF

        prior = datetime.datetime.now()

        recalibrated_snp_matrix_mpileup = recalibrate_ddbb_vcf_intermediate(
            compare_snp_matrix_recal_intermediate, out_variant_dir, min_cov_low_freq=10)
        recalibrated_snp_matrix_mpileup.to_csv(
            compare_snp_matrix_recal_mpileup, sep='\t', index=False)

        after = datetime.datetime.now()
        print(('\n' + "Done with function recalibrate_ddbb_vcf_intermediate in: %s" %
               (after - prior) + "\n"))

        # Remove SNPs located within INDELs

        prior = datetime.datetime.now()

        compare_snp_matrix_INDEL_intermediate_df = remove_position_range(
            recalibrated_snp_matrix_mpileup)
        compare_snp_matrix_INDEL_intermediate_df.to_csv(
            compare_snp_matrix_INDEL_intermediate, sep='\t', index=False)

        after = datetime.datetime.now()
        print(("Done with function remove_position_range in: %s" %
               (after - prior) + "\n"))

        # Extract all positions marked as complex

        prior = datetime.datetime.now()

        complex_variants = extract_complex_list(out_variant_dir)
        logger.debug('Complex positions in all samples:\n{}'.format(
            (','.join([str(x) for x in complex_variants]))))

        after = datetime.datetime.now()
        print(("Done with function extract_complex_list in: %s" %
               (after - prior) + "\n"))

        # Extract all low coverage o not covered positions

        prior = datetime.datetime.now()

        symbol_file = full_path_compare + '_symbol_lowcov.tsv'
        symbol_lowcov = extract_lowcov(compare_snp_matrix_INDEL_intermediate) # It is made by the INDEL_intermediate.tsv, taking 0 and 1 into account, can also be made with intermediate.highfreq.tsv
        symbol_lowcov.to_csv(symbol_file, sep='\t', index=False)

        after = datetime.datetime.now()
        print(("Done with function extract_lowcov in: %s" %
               (after - prior) + "\n"))

        # Clean all faulty positions and samples for final table

        prior = datetime.datetime.now()

        if args.complex:
            remove_complex_positions = complex_variants
        else:
            remove_complex_positions = False

        recalibrated_revised_INDEL_df = revised_df(compare_snp_matrix_INDEL_intermediate_df, path_compare, complex_pos=remove_complex_positions, min_freq_include=0.7, min_threshold_discard_uncov_sample=args.min_threshold_discard_uncov_sample, min_threshold_discard_uncov_pos=args.min_threshold_discard_uncov_pos, min_threshold_discard_htz_sample=args.min_threshold_discard_htz_sample,
                                                   min_threshold_discard_htz_pos=args.min_threshold_discard_htz_pos, min_threshold_discard_all_pos=args.min_threshold_discard_all_pos, min_threshold_discard_all_sample=args.min_threshold_discard_all_sample, remove_faulty=True, drop_samples=True, drop_positions=True, windows_size_discard=args.window)

        recalibrated_revised_df = recalibrated_revised_INDEL_df[~recalibrated_revised_INDEL_df['Position'].isin(symbol_lowcov['Position'])]

        recalibrated_revised_df.to_csv(
            compare_snp_matrix_recal, sep='\t', index=False)

        if args.only_snp:
            compare_only_snps_df = extract_only_snps(
                compare_snp_matrix_recal)
            compare_only_snps_df.to_csv(
                compare_only_snps, sep="\t", index=False)

        after = datetime.datetime.now()
        print(("Done with function revised_df in: %s" % (after - prior) + "\n"))

        # Matrix to pairwise and nwk

        prior = datetime.datetime.now()

        ddtb_compare(compare_snp_matrix_recal, distance=args.distance)

        if args.only_snp:
            ddtb_compare(compare_only_snps, distance=args.distance)

        after = datetime.datetime.now()
        print(("Done with function ddtb_compare in: %s" % (after - prior) + "\n"))

        # Annotated SNPs from BED file (genes or positions of interest)

        prior = datetime.datetime.now()

        if args.extract_bed:
            annotated_snps_final = extract_bed_positions(
                recalibrated_revised_INDEL_df, args.extract_bed, full_path_compare)

        after = datetime.datetime.now()
        print(("Done with function extract_bed_positions in: %s" %
               (after - prior) + "\n"))

        logger.info('\n' + MAGENTA + BOLD + 'COMPARISON FINISHED IN GROUP: ' +
                    group_name + END_FORMATTING + '\n')

        logger.info("\n" + MAGENTA + BOLD +
                    "##### END OF ONT VARIANT CALLING PIPELINE #####" + "\n" + END_FORMATTING)

    else:

        prior = datetime.datetime.now()

        compare_matrix = os.path.abspath(args.only_compare)
        ddtb_compare(compare_matrix, distance=args.distance)

        after = datetime.datetime.now()
        print(("Done with function ddtb_compare in: %s" % (after - prior) + "\n"))

        logger.info('\n' + MAGENTA + BOLD + 'COMPARISON FINISHED IN GROUP: ' +
                    group_name + END_FORMATTING + '\n')

        logger.info("\n" + MAGENTA + BOLD +
                    "##### END OF ONT VARIANT CALLING PIPELINE #####" + "\n" + END_FORMATTING)
