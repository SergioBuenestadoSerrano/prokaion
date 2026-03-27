#!/usr/bin/env python

from cmd import Cmd
import os
import sys
import re
import subprocess
import shutil
import logging
import datetime
import pandas as pd
import numpy as np
from statistics import mean
import multiprocessing
import concurrent.futures
from pandarallel import pandarallel


logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com), Pedro J. Sola (pedroscampoy@gmail.com)
Version = 0
Created: 24 March 2021

TODO:
    Adapt check_reanalysis
    Check file with multiple arguments
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""

# Colors and formatting

"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
"""

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


### Executing functions ###


def execute_subprocess(cmd, isShell=False, isInfo=False):
    """
    https://crashcourse.housegordon.org/python-subprocess.html
    https://docs.python.org/3/library/subprocess.html 
    Execute and handle errors with subprocess, outputting stderr instead of the subprocess CalledProcessError
    """

    logger.debug("")
    logger.debug(cmd)

    if cmd[0] == "samtools" or cmd[0] == "bwa":
        prog = " ".join(cmd[0:2])
        param = cmd[3:]
    else:
        prog = cmd[0]
        param = cmd[1:]

    try:
        command = subprocess.run(
            cmd, shell=isShell, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if command.returncode == 0:
            logger.debug(
                GREEN + DIM + "Program %s successfully executed" % prog + END_FORMATTING)
        else:
            logger.info(RED + BOLD + "Command %s FAILED\n" % prog + END_FORMATTING + BOLD + "with parameters: " + END_FORMATTING + " ".join(
                param) + "\n" + BOLD + "EXIT-CODE: %d\n" % command.returncode + "ERROR:\n" + END_FORMATTING + command.stderr.decode().strip())

        if isInfo:
            logger.info(command.stdout.decode().strip())
        else:
            logger.debug(command.stdout.decode().strip())

        logger.debug(command.stderr.decode().strip())

    except OSError as e:
        sys.exit(RED + BOLD + "Failed to execute program '%s': %s" % (prog,
                                                                      str(e)) + END_FORMATTING)


### Manipulation of files and paths ###


def check_create_dir(path):
    # exists = os.path.isfile(path)
    # exists = os.path.isdir(path)

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)


def check_file_exists(file_name):
    """
    Check file exist and is not 0Kb, if not program exit.
    """

    # Retrieve the file into to check if has size > 0
    file_info = os.stat(file_name)

    if not os.path.isfile(file_name) or file_info.st_size == 0:
        logger.info(RED + BOLD + "File: %s not found or empty\n" %
                    file_name + END_FORMATTING)
        sys.exit(1)
    return os.path.isfile(file_name)


def check_remove_file(file_name):
    """
    Check file exist and remove it.
    """

    if os.path.exists(file_name):
        os.remove(file_name)


def extract_read_list(input_dir):

    input_dir = os.path.abspath(input_dir)
    all_files = []

    for root, _, files in os.walk(input_dir):
        if root == input_dir:  # This only apply to parent folder, not subdirectories
            for name in files:
                filename = os.path.join(root, name)
                # is_files = re.match(r".*\.f(ast)*[q5](\.gz)*", name)
                is_files = re.match(r".*\.(f(ast)*[q5]|pod5)(\.gz)*", name)
                if is_files:
                    all_files.append(filename)

    all_files = sorted(all_files)

    return all_files


def extract_sample_list(file):

    basename_file = os.path.basename(file)
    basename_file = basename_file.split(".")[0]

    return basename_file


def check_reanalysis(output_dir, samples_to_analyze):

    output_dir = os.path.abspath(output_dir)

    new_samples = []

    variant_dir = os.path.join(output_dir, "Variants")
    compare_dir = os.path.join(output_dir, "Compare")

    previous_files = [variant_dir, compare_dir]

    # Check how many folders exist
    file_exist = sum([os.path.exists(x)
                      for x in previous_files])  # True = 1, False = 0

    # Handle reanalysis: First time; reanalysis or realysis with aditional samples
    if file_exist > 0:  # Already analysed

        previous_samples_list = os.listdir(variant_dir)

        if len(samples_to_analyze) == len(previous_samples_list):
            logger.info(
                MAGENTA + "\nPrevious analysis detected, no new sequences added\n" + END_FORMATTING)
        else:
            new_samples = set(samples_to_analyze) - set(previous_samples_list)
            logger.info(MAGENTA + "\nPrevious analysis detected, " +
                        str(len(new_samples)) + " new sequences added\n" + END_FORMATTING)

    return list(new_samples)


def file_to_list(file_name):

    list_F = []

    file_name_abs = os.path.abspath(file_name)
    with open(file_name_abs, "r") as f:
        for line in f:
            list_F.append(line.strip())

    return list_F


def remove_low_quality(input_dir, output_dir, cov20=70, unmapped_per=25, min_hq_snp=8, type_remove="Uncovered"):

    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))

    output_dir = os.path.abspath(output_dir)
    out_bam_dir = os.path.join(output_dir, "Bam")
    out_consensus_dir = os.path.join(output_dir, "Consensus")
    out_stats_dir = os.path.join(output_dir, "Stats")
    out_stats_bamstats_dir = os.path.join(out_stats_dir, "Bamstats")
    out_stats_coverage_dir = os.path.join(out_stats_dir, "Coverage")

    folder_list = [out_bam_dir, out_consensus_dir, out_stats_dir,
                   out_stats_bamstats_dir, out_stats_coverage_dir]

    uncovered_dir = os.path.join(output_dir, type_remove)  # Uncovered or Mixed
    uncovered_dir_variants = os.path.join(uncovered_dir, "Variants")
    uncovered_bam_dir = os.path.join(uncovered_dir, "Bam")
    uncovered_consensus_dir = os.path.join(uncovered_dir, "Consensus")
    uncovered_stats_dir = os.path.join(uncovered_dir, "Stats")
    uncovered_stats_bam_dir = os.path.join(uncovered_stats_dir, "Bamstats")
    uncovered_stats_cov_dir = os.path.join(uncovered_stats_dir, "Coverage")

    uncovered_folder_list = [uncovered_dir, uncovered_bam_dir, uncovered_consensus_dir,
                             uncovered_stats_dir, uncovered_stats_bam_dir, uncovered_stats_cov_dir]
    [check_create_dir(x) for x in uncovered_folder_list]

    uncovered_samples = []

    for root, _, files in os.walk(output_dir):
        # Any previous file created except for Table for mixed samples
        # and Species for both uncovered and mixed
        if root.endswith("Stats"):
            for name in files:
                filename = os.path.join(root, name)
                # Retreive coverage stats
                if name.endswith("overal.stats.tab"):
                    coverage_stat_file = filename
                    stats_df = pd.read_csv(coverage_stat_file, sep="\t")
                    stats_df = stats_df.fillna(0)

                    # Handle repeated names
                    stats_df["HQ70_SNP"] = stats_df["HQ70_SNP"].astype(str)

                    def f(x):
                        return (x if x.replace(".", "", 1).isdigit() else max(x.strip("()").split(",")))

                    stats_df["HQ70_SNP"] = stats_df.apply(
                        lambda x: f(x.HQ70_SNP), axis=1)
                    stats_df["HQ70_SNP"] = stats_df["HQ70_SNP"].astype(float)

                    # Store samples under any of the parameters indicated
                    uncovered_samples = stats_df["#SAMPLE"][(stats_df["COV>20X"] <= cov20) | (stats_df["UNMAPPED_PROP"] >= unmapped_per) | (
                        stats_df["HQ70_SNP"] < min_hq_snp)].tolist()
                    # print(uncovered_samples)

                    # Create a df with only covered to replace the original
                    covered_df = stats_df[~stats_df["#SAMPLE"].isin(
                        uncovered_samples)]
                    # print(covered_df)
                    covered_df.to_csv(coverage_stat_file,
                                      sep="\t", index=False)

                    # Create a df with uncovered
                    uncovered_df = stats_df[stats_df["#SAMPLE"].isin(
                        uncovered_samples)]
                    uncovered_table_filename = right_now_full + "_uncovered.summary.tab"
                    uncovered_table_file = os.path.join(
                        uncovered_stats_dir, uncovered_table_filename)

                    if len(uncovered_samples) > 0:
                        uncovered_df.to_csv(
                            uncovered_table_file, sep="\t", index=False)

                elif name.endswith(".coverage.summary.tab"):
                    covstats_df = pd.read_csv(filename, sep="\t")
                    final_covstat = filename

    uncovered_samples = [str(x) for x in uncovered_samples]
    # print(uncovered_samples)

    def_covstats_df = covstats_df[~covstats_df["#SAMPLE"].isin(
        uncovered_samples)]
    def_covstats_df.to_csv(final_covstat, sep="\t", index=False)

    logger.debug("Uncovered_samples: ")
    logger.debug(uncovered_samples)

    fastq = extract_read_list(input_dir)
    # print(fastq)

    sample_list_F = []

    for sample in fastq:
        sample = extract_sample_list(sample)
        # sample_name = sample.split('_')[1]
        sample_list_F.append(sample)

    # print(sample_list_F)

    # MOVE Fastq
    if len(uncovered_samples) > 0:
        for uncovered_sample in uncovered_samples:
            try:
                uncovered_index = sample_list_F.index(uncovered_sample)
                # print(uncovered_index)
                destination_fastq = os.path.join(
                    uncovered_dir, fastq[uncovered_index].split("/")[-1])
                # print(destination_fastq)
                logger.debug("Moving FAULTY Fastq {} TO {}".format(
                    fastq[uncovered_index], destination_fastq))
                # print(fastq[uncovered_index])
                shutil.move(fastq[uncovered_index], destination_fastq)
            except:
                logger.info("ERROR: No uncovered detected")

    # Move Variant Folder
    for root, _, files in os.walk(output_dir):
        if root.endswith("Variants") and not "Uncovered" in root:
            # print(root)
            if len(uncovered_samples) > 0:
                for uncovered in uncovered_samples:
                    filename = os.path.join(root, uncovered)
                    # print(filename)
                    destination_file = os.path.join(
                        uncovered_dir_variants, uncovered)
                    # print(destination_file)
                    if not os.path.exists(destination_file):
                        logger.debug("Moving FAULTY folder {} TO {}".format(
                            filename, destination_file))
                        shutil.move(filename, destination_file)
                    else:
                        logger.debug(
                            "{} already exist".format(destination_file))

    # Move Bam, Consensus and Stats files to Uncovered folder
    for outdir, uncovdir in zip(folder_list, uncovered_folder_list[1:]):
        for root, _, files in os.walk(outdir):
            # print(root)
            if root == outdir:
                for name in files:
                    filename = os.path.join(outdir, name)
                    # print(filename)
                    sample = filename.split("/")[-1].split(".")[0]
                    # print(sample)
                    if sample in uncovered_samples:
                        destination_file = os.path.join(uncovdir, name)
                        # print(destination_file)
                        logger.debug("Moving FAULTY file {} TO {}".format(
                            filename, destination_file))
                        shutil.move(filename, destination_file)

    return uncovered_samples


### Processing Reference files ###


def optimal_chunk_size(reference):
    """
    https://www.bioinformatics.nl/emboss-explorer/
        http://www.sacs.ucsf.edu/Documentation/emboss/infoseq.html
    """

    input_reference = os.path.abspath(reference)
    input_folder = os.path.dirname(reference)
    chunks_size = "chunks_size"
    chunks_size_path = os.path.join(input_folder, "chunks_size")

    cmd_infoseq = ["infoseq", "-auto", "-only", "-length", "-noheading",
                   "-odirectory", input_folder, "-outfile", chunks_size, input_reference]
    # print(cmd_infoseq)
    execute_subprocess(cmd_infoseq, isShell=False)

    df = pd.read_csv(chunks_size_path, header=None)
    size = df.loc[df.index[0]]
    return size


def samtools_faidx(reference):
    # samtools faidx reference.fa

    input_reference = os.path.abspath(reference)
    fai_reference = input_reference + ".fai"

    if os.path.isfile(fai_reference):
        logger.info(fai_reference + " already EXISTS")
    else:
        logger.info(GREEN + "Indexing " + fai_reference + END_FORMATTING)
        cmd_faidx = "samtools", "faidx", reference
        execute_subprocess(cmd_faidx, isShell=False)


def create_reference_chunks(reference):

    ref_size = optimal_chunk_size(reference)
    min_freebayes_chunk_size = 1000
    chunk_cpu = multiprocessing.cpu_count() - 2
    num_chunks = round(
        max(min_freebayes_chunk_size, int(ref_size) / chunk_cpu))

    input_reference = os.path.abspath(reference)
    input_folder = os.path.dirname(reference)
    out_reference_file = os.path.join(
        input_folder, "reference." + str(num_chunks) + ".regions")
    fai_reference = input_reference + ".fai"

    if os.path.isfile(out_reference_file):
        logger.info(out_reference_file + " already EXISTS")
    else:
        logger.info(GREEN + "Creating " + out_reference_file + END_FORMATTING)
        cmd_chunks = "fasta_generate_regions.py {} {} > {}".format(
            fai_reference, num_chunks, out_reference_file)
        execute_subprocess(cmd_chunks, isShell=True)

    return out_reference_file


### BAM Variant ###


def extract_indels(input_vcf):

    input_vcf = os.path.abspath(input_vcf)
    vcf_dir = ("/").join(input_vcf.split("/")[0:-1])
    output_indel_vcf = os.path.join(vcf_dir, "snps.indel.vcf")

    with open(output_indel_vcf, "w+") as fout:
        with open(input_vcf, "r") as f:
            for line in f:
                if "TYPE=ins" in line or "TYPE=del" in line:
                    fout.write(line)


def merge_vcf(snp_vcf, indel_vcf):

    snp_vcf = os.path.abspath(snp_vcf)
    indel_vcf = os.path.abspath(indel_vcf)

    vcf_dir = ("/").join(snp_vcf.split("/")[0:-1])
    output_complete_vcf = os.path.join(vcf_dir, "snps.all.vcf")

    with open(output_complete_vcf, "w+") as fout:
        with open(snp_vcf, "r") as f1:
            for line in f1:
                fout.write(line)
        with open(indel_vcf, "r") as f2:
            for line in f2:
                if not line.startswith("#"):
                    fout.write(line)


def create_bamstat(input_bam, output_file, threads=36):

    cmd_bamstat = "samtools flagstat --threads {} {} > {}".format(
        str(threads), input_bam, output_file)
    # print(cmd_bamstat)
    execute_subprocess(cmd_bamstat, isShell=True)


def create_coverage(input_bam, output_file):

    cmd_coverage = "samtools depth -aa {} > {}".format(input_bam, output_file)
    # print(cmd_coverage)
    execute_subprocess(cmd_coverage, isShell=True)


def calculate_cov_stats(file_cov):
    sample = file_cov.split("/")[-1].split(".")[0]
    df = pd.read_csv(file_cov, sep="\t", names=["#CHROM", "POS", "COV"])
    unmapped_pos = len(df.POS[df.COV == 0].tolist())
    pos_0_10 = len(df.POS[(df.COV > 0) & (df.COV <= 10)].tolist())
    pos_high10 = len(df.POS[(df.COV > 10)].tolist())
    pos_high20 = len(df.POS[(df.COV > 20)].tolist())
    pos_high30 = len(df.POS[(df.COV > 30)].tolist())
    pos_high50 = len(df.POS[(df.COV > 50)].tolist())
    pos_high100 = len(df.POS[(df.COV >= 100)].tolist())
    pos_high500 = len(df.POS[(df.COV >= 500)].tolist())
    pos_high1000 = len(df.POS[(df.COV >= 1000)].tolist())
    total_pos = df.shape[0]
    unmapped_prop = "%.2f" % ((unmapped_pos / total_pos) * 100)
    prop_0_10 = "%.2f" % ((pos_0_10 / total_pos) * 100)
    prop_high10 = "%.2f" % ((pos_high10 / total_pos) * 100)
    prop_high20 = "%.2f" % ((pos_high20 / total_pos) * 100)
    prop_high30 = "%.2f" % ((pos_high30 / total_pos) * 100)
    prop_high50 = "%.2f" % ((pos_high50 / total_pos) * 100)
    prop_high100 = "%.2f" % ((pos_high100 / total_pos) * 100)
    prop_high500 = "%.2f" % ((pos_high500 / total_pos) * 100)
    prop_high1000 = "%.2f" % ((pos_high1000 / total_pos) * 100)

    mean_cov = "%.2f" % (df.COV.mean())

    return (sample, mean_cov, unmapped_prop, prop_0_10, prop_high10, prop_high20, prop_high30, prop_high50, prop_high100, prop_high500, prop_high1000)


def obtain_group_cov_stats(directory, group_name):

    directory_path = os.path.abspath(directory)
    samples_to_skip = []
    previous_stat = False

    output_group_name = group_name + ".coverage.summary.tab"
    output_file = os.path.join(directory_path, output_group_name)

    if os.path.exists(output_file):
        previous_stat = True
        df_stat = pd.read_csv(output_file, sep="\t")
        samples_to_skip = df_stat["#SAMPLE"].tolist()
        logger.debug("Skipped samples for coverage calculation:" +
                     (",").join(str(samples_to_skip)))

    columns = ["#SAMPLE", "MEAN_COV", "UNMAPPED_PROP", "COV1-10X",
               "COV>10X", "COV>20X", "COV>30X", "COV>50X", "COV>100X", "COV>500X", "COV>1000X"]

    files_list = []

    for root, _, files in os.walk(directory):
        for name in files:
            if name.endswith(".cov"):
                filename = os.path.join(root, name)
                sample = name.split(".")[0]
                # df[columns] = df.apply(calculate_cov_stats(filename), axis=1, result_type='expand')
                if not sample in samples_to_skip:
                    files_list.append(filename)

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as executor:
        dfs = executor.map(calculate_cov_stats, files_list)
    df = pd.DataFrame(dfs, columns=columns)

    if previous_stat:
        df = pd.concat([df_stat, df], ignore_index=True, sort=True)
        df = df[columns]
        df.to_csv(output_file, sep="\t", index=False)
    else:
        df = df[columns]
        df.to_csv(output_file, sep="\t", index=False)


def extract_snp_count(out_variant_dir, filename_out):

    filename_out = str(filename_out)
    if "." in filename_out:
        filename_out = filename_out.split(".")[0]

    raw_var_folder = os.path.join(out_variant_dir, filename_out)
    filename = os.path.join(raw_var_folder, "snps.all.ivar.tsv")

    if os.path.exists(filename):
        df = pd.read_csv(filename, sep="\t")
        df = df.drop_duplicates(subset=["POS", "REF", "ALT"], keep="first")
        high70_quality_snps = df["POS"][(df.TOTAL_DP >= 30) & 
                                      (df.ALT_FREQ >= 0.7) & 
                                      (df.TYPE == "snp") & 
                                      ~df.OLDVAR.isin(['complex', 'mnp'])].tolist()
        high80_quality_snps = df["POS"][(df.TOTAL_DP >= 30) & 
                                      (df.ALT_FREQ >= 0.8) & 
                                      (df.TYPE == "snp") & 
                                      ~df.OLDVAR.isin(['complex', 'mnp'])].tolist()
        high90_quality_snps = df["POS"][(df.TOTAL_DP >= 30) & 
                                      (df.ALT_FREQ >= 0.9) & 
                                      (df.TYPE == "snp") & 
                                      ~df.OLDVAR.isin(['complex', 'mnp'])].tolist()
        htz_snps = df["POS"][(df.TOTAL_DP >= 30) & 
                             (df.ALT_FREQ < 0.7) & 
                             (df.ALT_FREQ >= 0.3) & 
                             (df.TYPE == "snp") & 
                             ~df.OLDVAR.isin(['complex', 'mnp'])].tolist()
        indels = df["POS"][(df.TOTAL_DP >= 30) & 
                           (df.ALT_FREQ >= 0.7) & 
                           ((df.TYPE == "ins") | (df.TYPE == "del"))].tolist()
        return (len(high70_quality_snps), len(high80_quality_snps), len(high90_quality_snps), len(htz_snps), len(indels))
    else:
        logger.debug("FILE " + filename + " NOT FOUND")
        return None


def extract_mapped_reads(out_stats_bamstats_dir, filename_out):

    filename_out = str(filename_out)
    if "." in filename_out:
        filename_out = filename_out.split(".")[0]
    filename = os.path.join(out_stats_bamstats_dir, filename_out + ".bamstats")

    if os.path.exists(filename):
        reads_mapped = 0
        perc_mapped = 0
        # properly_paired = 0
        # paired_percentage = 0
        with open(filename, "r") as f:
            logger.debug("File bamstat: {}".format(filename))
            for line in f:
                if "mapped" in line and "%" in line:
                    reads_mapped = line.split(" ")[0]
                    perc_mapped = line.split("(")[-1].split("%")[0]
                # elif 'properly paired' in line:
                #     properly_paired = line.split(' ')[0]
                #     paired_percentage = line.split('(')[-1].split('%')[0]

        # logger.debug((',').join(
        #     [reads_mapped, perc_mapped, properly_paired, paired_percentage]))

        # logger.debug(len([x for x in [reads_mapped, perc_mapped,
        #                               properly_paired, paired_percentage] if x != 0]))

        if len([x for x in [reads_mapped, perc_mapped] if x != 0]):
            return int(reads_mapped), float(perc_mapped)
        else:
            return 0, 0
    else:
        logger.info("FILE " + filename + " NOT FOUND")
        return None


# def extract_n_consensus(out_consensus_dir, filename_out):

#     filename_out = str(filename_out)
#     if "." in filename_out:
#         filename_out = filename_out.split(".")[0]

#     filename = os.path.join(out_consensus_dir, filename_out + ".fa")

#     if os.path.exists(filename):
#         with open(filename, "r") as f:
#             content = f.read()
#             content_list = content.split("\n")
#             # sample_fq = content_list[0].strip('>')
#             # In case fasta is in several lines(not by default)
#             sequence = ("").join(content_list[1:]).strip()
#             all_N = re.findall(r"N+", sequence)
#             if all_N:
#                 leading_N = re.findall(r"^N+", sequence)
#                 tailing_N = re.findall(r"N+$", sequence)
#                 length_N = [len(x) for x in all_N]
#                 individual_N = [x for x in length_N if x == 1]
#                 mean_length_N = mean(length_N)
#                 sum_length_N = sum(length_N)
#                 total_perc_N = sum_length_N / len(sequence) * 100
#                 return (len(all_N), len(individual_N), len(leading_N), len(tailing_N), sum_length_N, total_perc_N, mean_length_N)
#             else:
#                 return (0, 0, 0, 0, 0, 0, 0)

#     else:
#         logger.info("FILE " + filename + " NOT FOUND")
#         return None


def obtain_overal_stats(out_stats_dir, output_dir, group):
    pandarallel.initialize()

    samples_to_skip = []
    previous_stat = False

    overal_stat_file = os.path.join(out_stats_dir, group + ".overal.stats.tab")
    # print(overal_stat_file)
    out_variant_dir = os.path.join(output_dir, "Variants")
    out_stats_bamstats_dir = os.path.join(out_stats_dir, "Bamstats")
    out_consensus_dir = os.path.join(output_dir, "Consensus")

    columns = ["#SAMPLE", "MEAN_COV", "UNMAPPED_PROP", "COV1-10X",
               "COV>10X", "COV>20X", "COV>30X", "COV>50X", "COV>100X", "COV>500X", "COV>1000X"]

    if os.path.exists(overal_stat_file):
        previous_stat = True
        df_stat = pd.read_csv(overal_stat_file, sep="\t")
        samples_to_skip = df_stat["#SAMPLE"].tolist()
        logger.debug("Skipped samples for coverage calculation:" +
                     (",").join(str(samples_to_skip)))

    for root, _, files in os.walk(out_stats_dir):
        for name in files:
            if name.endswith("coverage.summary.tab"):
                # print(name)
                filename = os.path.join(root, name)
                # print(filename)
                df = pd.read_csv(filename, sep="\t")
                df = df[~df["#SAMPLE"].isin(samples_to_skip)]
                # print(df)
                if df.shape[0] > 0:
                    df[["HQ70_SNP", "HQ80_SNP", "HQ90_SNP", "HTZ_SNP", "INDELS"]] = df.parallel_apply(lambda x: extract_snp_count(
                        out_variant_dir, x["#SAMPLE"]), axis=1, result_type="expand")
                    df[["reads_mapped", "perc_mapped"]] = df.parallel_apply(lambda x: extract_mapped_reads(
                        out_stats_bamstats_dir, x["#SAMPLE"]), axis=1, result_type="expand")
                    # df[["N_groups", "N_individual", "N_leading", "N_tailing", "N_sum_len", "N_total_perc", "N_mean_len"]] = df.parallel_apply(
                    #     lambda x: extract_n_consensus(out_consensus_dir, x["#SAMPLE"]), axis=1, result_type="expand")

    if previous_stat:
        df = pd.concat([df_stat, df], ignore_index=True, sort=True)
        df = df[columns + [col for col in df.columns if col != "#SAMPLE" and col != "MEAN_COV" and col != "UNMAPPED_PROP" and col !=
                           "COV1-10X" and col != "COV>10X" and col != "COV>20X" and col != "COV>30X" and col != "COV>50X" and col != "COV>100X" and col != "COV>500X" and col != "COV>1000X"]]
        df.to_csv(overal_stat_file, sep="\t", index=False)
    else:
        df = df[columns + [col for col in df.columns if col != "#SAMPLE" and col != "MEAN_COV" and col != "UNMAPPED_PROP" and col !=
                           "COV1-10X" and col != "COV>10X" and col != "COV>20X" and col != "COV>30X" and col != "COV>50X" and col != "COV>100X" and col != "COV>500X" and col != "COV>1000X"]]
        df.to_csv(overal_stat_file, sep="\t", index=False)


### VCF processing ###


def import_VCF42_freebayes_to_tsv(vcf_file, sep="\t"):

    vcf_file = os.path.abspath(vcf_file)
    tsv_file = (".").join(vcf_file.split(".")[:-1]) + ".tsv"

    headers = []
    extra_fields = ["TYPE", "DP", "RO", "AO"]
    with open(tsv_file, "w+") as fout:
        with open(vcf_file, "r") as f:
            next_line = f.readline().strip()
            while next_line.startswith("#"):
                next_line = f.readline().strip()
                if next_line.startswith("#CHROM"):
                    headers = next_line.split("\t")

        headers = headers[:7] + extra_fields + ["OLDVAR"]
        fout.write(("\t").join(headers) + "\n")

        with open(vcf_file, "r") as f:
            for line in f:
                extra_field_list = []
                # and not 'complex' in line and not 'mnp' in line
                if not line.startswith("#"):
                    line_split = line.split(sep)[:8]
                    info = line_split[-1].split(";")
                    for field in extra_fields:
                        extra_field_list.append(
                            [x.split("=")[-1] for x in info if field in x][0])
                    if "OLDVAR" in line:
                        extra_field_list.append(
                            [x.split("=")[-1] for x in info if "OLDVAR" in x][0].split(",")[0])
                    output_line = ("\t").join(
                        line_split[:7] + extra_field_list)
                    fout.write(output_line + "\n")


def import_tsv_freebayes_to_df(tsv_file, sep="\t"):

    tsv_file = os.path.abspath(tsv_file)

    df = pd.read_csv(tsv_file, sep=sep)

    df.rename(columns={"#CHROM": "REGION", "RO": "REF_DP",
                       "DP": "TOTAL_DP", "AO": "ALT_DP", "QUAL": "ALT_QUAL"}, inplace=True)

    df["REF_FREQ"] = df["REF_DP"] / df["TOTAL_DP"]
    df["ALT_FREQ"] = df["ALT_DP"] / df["TOTAL_DP"]

    df = df.sort_values(by=["POS"]).reset_index(drop=True)

    return df[["REGION", "POS", "ID", "REF", "ALT", "ALT_QUAL", "FILTER", "TOTAL_DP", "TYPE", "REF_DP", "ALT_DP", "REF_FREQ", "ALT_FREQ", "OLDVAR"]]


def vcf_to_ivar_tsv(input_vcf, output_tsv):

    input_tsv = (".").join(input_vcf.split(".")[:-1]) + ".tsv"

    import_VCF42_freebayes_to_tsv(input_vcf)

    df = import_tsv_freebayes_to_df(input_tsv)
    df.to_csv(output_tsv, sep="\t", index=False)


### Consensus ###


def ivar_consensus(filename_bam_out, out_consensus_dir, filename_out, min_quality=5, min_frequency_threshold=0.6, min_depth=5, uncovered_character="N"):
    # Find another solution, if we set q>7 an error occur "Segmentation fault", they are trying to fix it.
    """
    ivar consensus
        Usage: samtools mpileup -aa -A -R -d 0 -Q 0 <input.bam> | ivar consensus -p <prefix> 
        Note : samtools mpileup output must be piped into ivar consensus
        Input Options    Description
           -q    Minimum quality score threshold to count base (Default: 20)
           -t    Minimum frequency threshold(0 - 1) to call consensus. (Default: 0)
                 Frequently used thresholds | Description
                 ---------------------------|------------
                                          0 | Majority or most common base
                                        0.2 | Bases that make up atleast 20% of the depth at a position
                                        0.5 | Strict or bases that make up atleast 50% of the depth at a position
                                        0.9 | Strict or bases that make up atleast 90% of the depth at a position
                                          1 | Identical or bases that make up 100% of the depth at a position. Will have highest ambiguities
           -m    Minimum depth to call consensus(Default: 10)
           -k    If '-k' flag is added, regions with depth less than minimum depth will not be added to the consensus sequence. Using '-k' will override any option specified using -n 
           -n    (N/-) Character to print in regions with less than minimum coverage(Default: N)
        Output Options   Description
           -R    Ignore RG tags. Treat all reads in one BAM as one sample.
           -p    (Required) Prefix for the output fasta file and quality file
    """

    prefix = os.path.join(out_consensus_dir, filename_out)

    cmd_consensus = "samtools mpileup -aa -A -R -d 0 -B -Q 0 {} | ivar consensus -p {} -q {} -t {} -m {} -n {}".format(
        filename_bam_out, prefix, min_quality, min_frequency_threshold, min_depth, uncovered_character)
    # print(cmd_consensus)
    execute_subprocess(cmd_consensus, isShell=True)


def replace_consensus_header(input_fasta):
    with open(input_fasta, "r+") as f:
        content = f.read()
        header = content.split("\n")[0].strip(">")
        new_header = header.split("_")[1].strip()
        content = content.replace(header, new_header)
        f.seek(0)
        f.write(content)
        f.truncate()


### Annotation ###


def rename_reference_snpeff(input_vcf, chrom_filename, new_ref_name="Chromosome"):

    with open(input_vcf) as f:
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            next_line = f.readline()
        old_ref = f.readline().split("\t")[0].strip()

    with open(chrom_filename, "w+") as fo:
        with open(input_vcf, "r") as f:
            content = f.read()
            new_content = content.replace(old_ref, new_ref_name)
            fo.write(new_content)


def snpeff_execution(vcf_file, annot_file, database=False):

    with open(vcf_file, "r") as f:
        content_list = f.read().split("\n")
        lines = len(content_list)

    if lines > 1:
        cmd_snpeff = ["snpEff", "-noStats", database, vcf_file]
        with open(annot_file, "w+") as outfile:
            # Calculate coverage and save it in the output file
            subprocess.run(cmd_snpeff, stdout=outfile,
                           stderr=subprocess.PIPE, check=True, universal_newlines=True)
    else:
        with open(annot_file, "w+") as outfile:
            outfile.write("No annotation found")


def import_annot_to_pandas(vcf_file, sep='\t'):
    """
    Order several annoattion by:
    Putative impact: Effects having higher putative impact are first.
    Effect type: Effects assumed to be more deleterious effects first.
    Canonical transcript before non-canonical.
    Marker genomic coordinates (e.g. genes starting before first)
    https://pcingola.github.io/SnpEff/se_inputoutput/
    Parse vcf outputted by snpEFF which adds the ANN field
    Dependences: calculate_ALT_AD
                calculate_true_ALT
    """
    header_lines = 0
    with open(vcf_file) as f:
        first_line = f.readline().strip()
        if first_line == 'No annotation found':
            return pd.read_csv(vcf_file, sep=sep)
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            # logger.info(next_line)
            next_line = f.readline()

    # Use first line as header
    df = pd.read_csv(vcf_file, sep=sep, skiprows=[
                     header_lines], header=header_lines)

    ann_headers = ['Allele',
                   'Annotation',
                   'Annotation_Impact',
                   'Gene_Name',
                   'Gene_ID',
                   'Feature_Type',
                   'Feature_ID',
                   'Transcript_BioType',
                   'Rank',
                   'HGVS.c',
                   'HGVS.p',
                   'cDNA.pos / cDNA.length',
                   'CDS.pos / CDS.length',
                   'AA.pos / AA.length',
                   'ERRORS / WARNINGS / INFO']

    df['Type'] = df.apply(lambda row: 'snp' if len(row['REF']) == 1 and len(row['ALT']) == 1 else 'indel', axis=1)

    anlelle_headers = ['DP', 'REF_DP', 'ALT_DP']
    df[anlelle_headers] = pd.DataFrame([[None, None, None]] * len(df), columns=anlelle_headers)

    def assign_values(row):

        info_parts = row.INFO.split(';') if isinstance(row.INFO, str) else []

        def get_value(prefix):
            for part in info_parts:
                if part.startswith(prefix):
                    return part.split('=')[1]
            return None

        if row['Type'] == 'snp':
            values = [get_value('DP='), get_value('RO='), get_value('AO=')]
        elif row['Type'] in ['indel', 'del', 'ins']:
            values = [get_value('DP='), get_value('RO='), get_value('AO=')]
        else:
            values = [None, None, None]

        return pd.Series(values)

    df[anlelle_headers] = df.apply(assign_values, axis=1)

    for head in anlelle_headers:
        df[head] = df[head].str.split("=").str[-1]

    df['ALT_FREQ'] = pd.to_numeric(df['ALT_DP']) / pd.to_numeric(df['DP'])

    # return df

    df['TMP_ANN_16'] = df['INFO'].apply(
        lambda x: ('|').join(x.split('|')[0:15]))

    df.INFO = df.INFO.str.split("ANN=").str[-1]

    df = df.join(df.pop('INFO')
                   .str.strip(',')
                   .str.split(',', expand=True)
                   .stack()
                   .reset_index(level=1, drop=True)
                   .rename('INFO')).reset_index(drop=True)

    df['TMP_ANN_16'] = df['INFO'].apply(
        lambda x: ('|').join(x.split('|')[0:15]))
    df[ann_headers] = df['TMP_ANN_16'].str.split('|', expand=True)
    df['HGVS.c'] = df['HGVS.c'].str.split(".").str[-1]
    df['HGVS.p'] = df['HGVS.p'].str.split(".").str[-1].replace('', '-')

    df.drop(["INFO", "TMP_ANN_16", "REF_DP", "ALT_DP", "Allele", "ERRORS / WARNINGS / INFO", "cDNA.pos / cDNA.length"], inplace=True, axis=1)

    return df


def annotate_snpeff(input_vcf_file, output_annot_file, database="Mycobacterium_tuberculosis_h37rv"):

    # Execute snpEff
    snpeff_execution(input_vcf_file, output_annot_file, database=database)

    # Clean vcf format and remove info
    annot_df = import_annot_to_pandas(output_annot_file)
    annot_df.to_csv(output_annot_file, sep="\t", index=False)


def bed_to_df(bed_file):
    """
    Import bed file separated by tabs into a pandas df
    -Handle header line
    -Handle with and without description (If there is no description adds true or false to annotated df)
    """

    header_lines = 0

    # Handle likely header by checking colums 2 and 3 as numbers
    with open(bed_file, "r") as f:
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
        # Delim_whitespace=True
        df = pd.read_csv(bed_file, sep="\t", header=None)
    else:
        df = pd.read_csv(bed_file, sep="\t", skiprows=header_lines,
                         header=None)  # Delim_whitespace=True

    df = df.iloc[:, 0:4]
    df.columns = ["#CHROM", "start", "end", "description"]

    return df


def add_bed_info(bed_df, position):
    """
    Identify a position within a range
    credits: https://stackoverflow.com/questions/6053974/python-efficiently-check-if-integer-is-within-many-ranges
    """

    # dict_position = bed_to_dict(bed_file)
    if any(start <= position <= end for (start, end) in zip(bed_df.start.values.tolist(), bed_df.end.values.tolist())):
        description_out = bed_df.description[(
            bed_df.start <= position) & (bed_df.end >= position)].values[0]
        return description_out
    else:
        return None


def annotate_bed(tsv_df, bed_files):
    with open(tsv_df, "r") as f:
        content = f.read().strip()
        if content == "No annotation found":
            return pd.DataFrame(columns=["POS", "REF", "ALT", "INFO"])
        else:
            df = pd.read_csv(tsv_df, sep="\t")
            # Extract file name and use it as header
            variable_list = [x.split("/")[-1].split(".")[0] for x in bed_files]

            for variable_name, bed_file in zip(variable_list, bed_files):
                logger.info("Annotating BED: {}".format(bed_file))
                bed_annot_df = bed_to_df(bed_file)
                df[variable_name] = df["POS"].apply(
                    lambda x: add_bed_info(bed_annot_df, x))

            return df


def import_VCF_to_pandas(vcf_file):

    header_lines = 0

    with open(vcf_file) as f:
        first_line = f.readline().strip()
        next_line = f.readline().strip()
        while next_line.startswith("##"):
            header_lines = header_lines + 1
            # Logger.info(next_line)
            next_line = f.readline()

    if first_line.startswith("##"):
        df = pd.read_csv(vcf_file, sep="\t", skiprows=[
                         header_lines], header=header_lines)

        df["ALT"] = df["ALT"].str.upper()
        df["REF"] = df["REF"].str.upper()
        # Check INFO
        if "INFO" in df.columns:
            return df
        else:
            last_column = df.columns[-1]
            df = df.rename(columns={last_column: "INFO"})
            return df
    else:
        logger.info("This vcf file is not properly formatted")
        sys.exit(1)


def annotate_vcf(tsv_df, vcfs, include_ref=False):
    df = pd.read_csv(tsv_df, sep="\t")

    for vcf in vcfs:
        logger.info("Annotating VCF: {}".format(vcf))
        header = (".").join(vcf.split("/")[-1].split(".")[0:-1])
        dfvcf = import_VCF_to_pandas(vcf)
        if include_ref:
            dfvcf = dfvcf[["POS", "REF", "ALT", "INFO"]]
        else:
            dfvcf = dfvcf[["POS", "ALT", "INFO"]]
        dfvcf = dfvcf.rename(columns={"INFO": header})
        df = df.merge(dfvcf, how="left")

    return df


def user_annotation(tsv_file, output_file, vcf_files=[], bed_files=[]):

    bed_df = annotate_bed(tsv_file, bed_files)
    vcf_df = annotate_vcf(tsv_file, vcf_files)

    df = bed_df.merge(vcf_df)

    df.to_csv(output_file, sep="\t", index=False)


def checkAA(snpEffRow, dfAnnot, gene):
    df = dfAnnot
    df['aaAnnot'] = df['aa'] + ":" + df['annot']
    presence_list = [annot in snpEffRow for annot in dfAnnot.aa]

    if any(":" in annot for annot in dfAnnot.annot):
        for idx, row in dfAnnot.iterrows():
            if ":" in row.annot:
                annot_split = row.annot.split(":")
                if row.aa in snpEffRow and annot_split[0] in gene:
                    return annot_split[-1]
    else:
        annotation_list = np.array(df.aaAnnot.tolist())
        return (',').join(annotation_list[np.array(presence_list)])


def annotate_aa(annot_file, aas):

    df = pd.read_csv(annot_file, sep="\t")
    df = df.drop_duplicates(subset=["POS", "ALT"], keep="first")

    for aa in aas:
        header = (".").join(aa.split("/")[-1].split(".")[0:-1])
        dfaa = pd.read_csv(aa, sep="\t", names=['aa', 'annot'])
        if not header in df.columns:
            print("ANNOTATING AA: {}".format(aa))
            df['HGVS.p'] = df['HGVS.p'].astype(str)
            df[header] = df.apply(lambda x: checkAA(
                x['HGVS.p'], dfaa, x['Gene_Name']), axis=1)
        else:
            print("SKIPPED AA: {}".format(aa))

    return df


def user_annotation_aa(annot_file, output_file, aa_files=[]):
    with open(annot_file, "r") as f:
        content = f.read().strip()

        if content == "No annotation found":
            logger.debug("{} file has NO Annotation".format(annot_file))
            with open(output_file, "w+") as fout:
                fout.write("No annotation found")
        else:
            df = annotate_aa(annot_file, aa_files)
            # Filter SnpEFF output with aa annotations
            # df.drop_duplicates(subset=['HGVS.p'], keep='first', inplace=True)
            df.to_csv(output_file, sep="\t", index=False)


def make_blast(query_fasta, database, sample, output_folder, db_type="nucl", query_type="nucl", evalue=0.0001, threads=36):

    blast_command = "blastn" if query_type == "nucl" else "blastp"
    database_name = database.split("/")[-1].split(".")[0]
    output_database_tmp = os.path.join(
        output_folder, database_name + ".blast.tmp")
    output_blast = os.path.join(
        output_folder, sample + "." + database_name + ".blast")
    blastdb_cmd = ["makeblastdb", "-in", database,
                   "-out", output_database_tmp, "-dbtype", db_type]

    logger.info((",").join(blastdb_cmd))

    execute_subprocess(blastdb_cmd)

    blast_cmd = [blast_command, "-query", query_fasta, "-db", output_database_tmp, "-out", output_blast, "-evalue", str(evalue), "-num_threads", str(
        threads), "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"]

    logger.info((",").join(blast_cmd))

    execute_subprocess(blast_cmd)


def kraken(sample, report, kraken2_db, krona_html, threads=34):

    cmd_kraken = "kraken2 --db {} --use-names --threads {} --report {} --gzip-compressed {} --output -".format(
        kraken2_db, str(threads), report, sample)

    # print(cmd_kraken)
    execute_subprocess(cmd_kraken, isShell=True)

    cmd_krona = "ktImportTaxonomy -m 3 -q 2 -t 5 {} -o {}".format(
        report, krona_html)

    # print(cmd_krona)
    execute_subprocess(cmd_krona, isShell=True)


def mash_screen(sample, out_dir, mash_database, winner=True, threads=16):
    """
    https://mash.readthedocs.io/en/latest/index.html
    https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh #MASH refseq database
    wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
    mash screen -w -p 4 ../refseq.genomes.k21s1000.msh 4_R1.fastq.gz 4_R2.fastq.gz > 4.winner.screen.tab
    identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment
    """

    if not os.path.isfile(mash_database):
        logger.info(RED + BOLD + "Mash database can't be found\n" + END_FORMATTING + "You can download it typing:\n\
            wget https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh")
        sys.exit(1)

    cmd_mash = "mash screen -w -p {} {} {} > {}".format(
        threads, mash_database, sample, out_dir)

    # print(cmd_mash)
    execute_subprocess(cmd_mash, isShell=True)
