#!/usr/bin/env python

import os
import re
import logging
import argparse
import sys
import subprocess
import datetime
import gzip
import shutil

# Local application imports

from misc_prokaion import (
    check_create_dir, extract_read_list, extract_sample_list, execute_subprocess)


logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com)
Version = 0
Created: 24 March 2021

TODO:
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
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


def get_arguments():

    parser = argparse.ArgumentParser(
        prog="guppy_prokaion.py", description="Pipeline to basecalling and barcoding fast5 files from MinION sequencing")

    parser.add_argument("-i", "--input", dest="input_dir", metavar="input_directory",
                        type=str, required=True, help="REQUIRED. Input directory containing all fast5 files")

    parser.add_argument("-o", "--output", type=str, required=True,
                        help="REQUIRED. Output directory to extract all results")

    parser.add_argument('--basecall', type=str, default='dorado', required=True, help="REQUIRED. Program to use in data preprocessing (basecaller and barcoder)")

    parser.add_argument('--model', type=str, default='~/Scripts/git_repos/Dorado/Models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0', required=False, help='The basecaller model to run. Default dna_r10.4.1_e8.2_400bps_sup@v5.2.0 | dna_r10.4_e8.1_sup.cfg')

    parser.add_argument("-s", "--samples", metavar="Samples", type=str, required=False,
                        help="Sample list for conversion from barcode to samples ID")

    parser.add_argument("-b", "--require_barcodes_both_ends", required=False, action="store_true",
                        help="Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification")

    parser.add_argument("--barcode_kit", type=str, required=False, default="EXP-NBD104",
                        help="Kit of barcodes used [EXP-NBD104|SQK-RBK110-96|EXP-NBD196]. Default: EXP-NBD104")

    parser.add_argument('-g', '--gpu', required=False, default=False,
                        action="store_true", help='Specify GPU device: "auto", or "cuda:<device_id>"')

    parser.add_argument("-t", "--threads", type=int, dest="threads", required=False,
                        default=30, help="Threads to use. Default: 30)")

    parser.add_argument("--num_callers", type=int, dest="num_callers",
                        required=False, default=8, help="Number of parallel basecallers. Default: 8")

    parser.add_argument('--chunks', type=int, dest='chunks',
                        required=False, default=1536, help='Maximum chunks per runner. Default: 1536')

    parser.add_argument("--records_per_fastq", type=int, dest="records_per_fastq",
                        required=False, default=0, help="Maximum number of records per fastq")

    parser.add_argument("-q", "--min_base_quality", type=int, dest="min_quality", required=False,
                        default=8, help="Filter on a minimum average read quality score. Default: 8")

    parser.add_argument("-rq", "--min_read_quality", type=int, dest="min_read_quality",
                             required=False, default=8, help="Filter on a minimum average read quality score. Default: 8")

    parser.add_argument("--trim", type=str, dest="trim",
                             required=False, default="adapters", help="Specify what to trim. Options are 'none', 'all', and 'adapters'. The default behaviour is to trim all detected adapters, primers, and barcodes. Choose 'adapters' to just trim adapters. The 'none' choice is equivelent to using --no-trim. Note that this only applies to DNA. RNA adapters are always trimmed. Default: adapters")

    parser.add_argument("--length", type=int, dest="length", required=False,
                        default=500, help="Filter on a minimum read length. Default: 500bp")

    parser.add_argument("--headcrop", type=int, dest="headcrop", required=False,
                        default=20, help="Trim n nucleotides from start of read. Default: 20")

    parser.add_argument("--tailcrop", type=int, dest="tailcrop", required=False,
                        default=20, help="Trim n nucleotides from end of read. Default: 20")

    arguments = parser.parse_args()

    return arguments


def pod5_conversion(input_dir, out_pod5):

    """
    https://github.com/nanoporetech/pod5-file-format
    """

    cmd_conversion = "pod5 convert fast5 {}/*.fast5 --output {} --one-to-one {} --threads {}".format(input_dir, out_pod5, input_dir, str(args.threads))

    print(cmd_conversion)
    execute_subprocess(cmd_conversion, isShell=True)


def basecalling_dorado(input_dir, out_basecalling_dir):

    """
    https://github.com/nanoporetech/dorado
    """

    fastq_basecalled = os.path.join(out_basecalling_dir, 'Dorado_basecalled.fastq')

    cmd_basecalling = "dorado basecaller {} {} --trim {} --emit-fastq > {}".format(str(args.model), input_dir, args.trim, fastq_basecalled)

    print(cmd_basecalling)
    execute_subprocess(cmd_basecalling, isShell=True)


def demux_dorado(fastq_basecalled, out_barcoding_dir, require_barcodes_both_ends=False, barcode_kit="EXP-NBD104"):

    if require_barcodes_both_ends:
        logger.info(
            GREEN + BOLD + "Barcodes are being used at both ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = "--barcode-both-ends"
    else:
        logger.info(
            YELLOW + BOLD + "Barcodes are being used on at least 1 of the ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = ""

    cmd_demux = "dorado demux {} --kit-name {} --output-dir {} --emit-fastq {} --emit-summary".format(fastq_basecalled, barcode_kit, out_barcoding_dir, require_barcodes_both_ends)

    print(cmd_demux)
    execute_subprocess(cmd_demux, isShell=True)


def reorganize_demux(out_barcoding_dir):

    for filename in os.listdir(out_barcoding_dir):
        if filename.endswith('.fastq'):
            fastq_file_path = os.path.join(out_barcoding_dir, filename)
            # with open(fastq_file_path, 'rb') as f_in:
            #     with gzip.open(os.path.join(out_barcoding_dir, filename + '.gz'), 'wb') as f_out:
            #         shutil.copyfileobj(f_in, f_out)
            # os.remove(fastq_file_path)

            if 'barcode' in filename:
                barcode_folder = os.path.join(out_barcoding_dir, filename.split('_')[2].split('.')[0])
                os.makedirs(barcode_folder, exist_ok=True)
                shutil.move(fastq_file_path, barcode_folder)


def basecalling_ion(input_dir, out_basecalling_dir, config='dna_r9.4.1_450bps_fast.cfg', records=0):

    # -i: Path to input fast5 files
    # -s: Path to save fastq files
    # -c: Config file to use > https://community.nanoporetech.com/posts/guppy-v5-0-7-release-note (fast // hac // sup)
    # -x: Specify GPU device: 'auto', or 'cuda:<device_id>'
    # --num_callers: Number of parallel basecallers to Basecaller, if supplied will form part
    # --gpu_runners_per_device: Number of runners per GPU device.
    # --cpu_threads_per_caller: Number of CPU worker threads per basecaller
    # --chunks_per_runner: Maximum chunks per runner
    # --compress_fastq: Compress fastq output files with gzip
    # --records_per_fastq: Maximum number of records per fastq file, 0 means use a single file (per worker, per run id)

    if args.gpu != False:
        logger.info(
            GREEN + 'Basecalling executing on GPU device' + END_FORMATTING)
        gpu_device = "auto"
        cmd = ['guppy_basecaller', '-i', input_dir, '-s', out_basecalling_dir, '-c',
               config, '-x', gpu_device, '--records_per_fastq', str(records), '--compress_fastq']
    else:
        logger.info(
            YELLOW + 'Basecalling executing on CPU device' + END_FORMATTING)
        cmd = ['guppy_basecaller', '-i', input_dir, '-s', out_basecalling_dir, '-c', config, '--num_callers', str(args.num_callers), '--chunks_per_runner', str(
            args.chunks), '--cpu_threads_per_caller', str(args.threads), '--records_per_fastq', str(records), '--compress_fastq']

    print(cmd)
    execute_subprocess(cmd, isShell=False)


def barcoding_ion(out_basecalling_dir, out_barcoding_dir, require_barcodes_both_ends=False, barcode_kit="EXP-NBD104", threads=30):

    # -i: Path to input files
    # -r: Search for input file recursively
    # -s: Path to save files
    # -t: Number of worker threads
    # --num_barcoding_threads: Number of worker threads to use for barcoding.
    # -x: Specify GPU device to accelerate barcode detection: 'auto', or 'cuda:<device_id>'.

    # --fastq_out: Output Fastq files
    # --compress_fastq: Compress fastq output files with gzip
    # --barcode_kits: Space separated list of barcoding kit(s) or expansion kit(s) to detect against. Must be in double quotes
    # --require_barcodes_both_ends: Reads will only be classified if there is a barcode above the min_score at both ends of the read
    # --records_per_fastq: Maximum number of records per fastq file, 0 means use a single file (per worker, per run id)
    # --allow_inferior_barcodes: Reads will still be classified even if both the barcodes at the front and rear (if applicable) were not the best scoring barcodes above the min_score.

    # --detect_barcodes: Detect barcode sequences at the front and rear of the read.
    # --detect_adapter: Detect adapter sequences at the front and rear of the read.
    # --detect_primer: Detect primer sequences at the front and rear of the read.
    # --enable_trim_barcodes: Enable trimming of barcodes from the sequences in the output files. By default is false, barcodes will not be trimmed.
    # --trim_adapters: Trim the adapters from the sequences in the output files.
    # --trim_primers: Trim the primers from the sequences in the output files.

    # --min_score_barcode_front: Minimum score to consider a front barcode to be a valid barcode alignment (Default: 60).
    # --min_score_barcode_rear: Minimum score to consider a rear barcode to be a valid alignment (and min_score_front will then be used for the front only when this is set).

    # if args.gpu != False:
    #     gpu_device = "-x auto "
    # else:
    #     gpu_device = ""

    if require_barcodes_both_ends:
        logger.info(
            GREEN + BOLD + "Barcodes are being used at both ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = "--require_barcodes_both_ends"
    else:
        logger.info(
            YELLOW + BOLD + "Barcodes are being used on at least 1 of the ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = ""

    cmd = ["guppy_barcoder", "-i", out_basecalling_dir, "-s", out_barcoding_dir, "-r", require_barcodes_both_ends,
           "--barcode_kits", barcode_kit, "-t", str(threads), '--detect_barcodes', '--min_score_barcode_front', str(85), '--min_score_barcode_rear', str(85), '--enable_trim_barcodes', '--detect_mid_strand_barcodes', '--min_score_barcode_mid', str(60), '--detect_primer', '--trim_primers', '--detect_adapter', '--trim_adapters', "--fastq_out", "--compress_fastq"]

    print(cmd)
    execute_subprocess(cmd, isShell=False)


def rename_files(output_samples):

    with open(output_samples, "w") as bc_output:
        for fastq_file in sum_files:
            with open(fastq_file, "r") as f_in:
                shutil.copyfileobj(f_in, bc_output)


def ONT_QC_filtering(output_samples, filtered_samples):

    # -c: Write on standard output, keep the original files unchanged
    # -q: Filter on a minimum average read quality score
    # --minlength: Filter on a minimum read length
    # --headcrop: Trim n nucleotides from start of read
    # --tailcrop: Trim n nucleotides from end of read

    cmd_filtering = "chopper -i {} -q {} --minlength {} --headcrop {} --tailcrop {} --threads {} | gzip > {}".format(output_samples, str(args.min_read_quality), str(args.length), str(args.headcrop), str(args.tailcrop), str(args.threads), filtered_samples)

    # print(cmd_filtering)
    execute_subprocess(cmd_filtering, isShell=True)

def ONT_quality(filtered_sample, out_qc, threads=30):

    # --fastq_rich: Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy with additional information concerning with channel and time
    # --N50: Show the N50 mark in the read length histogram

    cmd_QC = ["NanoPlot", "--fastq", filtered_sample, "--N50", "-o", out_qc, "-t", str(threads)]

    # print(cmd_QC)
    execute_subprocess(cmd_QC, isShell=False)


if __name__ == "__main__":

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split("/")[-1]
    check_create_dir(output_dir)

    # Logging
    # Create log file with date and time

    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = group_name + "_" + right_now_full + ".log"
    log_folder = os.path.join(output_dir, "Logs")
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s:%(message)s")

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    logger.info(
        "\n" + BLUE + "############### START PROCESSING FAST5 FILES ###############" + END_FORMATTING + "\n")
    logger.info(args)


    # Declare folders created in pipeline and key files

    out_basecalling_dir = os.path.join(output_dir, "Basecalling")
    check_create_dir(out_basecalling_dir)
    basecalling_summary = os.path.join(
        out_basecalling_dir, "sequencing_summary.txt")
    out_barcoding_dir = os.path.join(output_dir, "Barcoding")
    check_create_dir(out_barcoding_dir)
    barcoding_summary = os.path.join(
        out_barcoding_dir, "barcoding_summary.txt")
    out_samples_dir = os.path.join(output_dir, "Samples_Fastq")
    check_create_dir(out_samples_dir)
    out_samples_filtered_dir = os.path.join(out_samples_dir, "Filtered_Fastq")
    check_create_dir(out_samples_filtered_dir)
    fastq_basecalled = os.path.join(out_basecalling_dir, 'Dorado_basecalled.fastq')
    # out_correction_dir = os.path.join(out_samples_dir, 'Corrected')
    # check_create_dir(out_correction_dir)
    out_qc_dir = os.path.join(output_dir, "Quality")
    check_create_dir(out_qc_dir)

    # Obtain all fast5 files from folder

    fast5 = extract_read_list(args.input_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fast5:
        sample = extract_sample_list(sample)
        sample_list.append(sample)

    logger.info("\n" + CYAN + "{} Samples will be analysed: {}".format(
        len(sample_list), ",".join(sample_list)) + END_FORMATTING)


    ############### START PIPELINE ###############

    # Basecalling #

    ### Dorado

    if args.basecall == 'dorado':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + "STARTING BASECALLING WITH DORADO" + END_FORMATTING + "\n")

        if os.path.exists(fastq_basecalled):

            logger.info("\n" + YELLOW + BOLD +
                        "Ommiting BASECALLING" + END_FORMATTING + "\n")

        else:

            if any(file.endswith('.fast5') for file in os.listdir(args.input_dir)):

                out_pod5_dir = os.path.join(input_dir, "Pod5")
                check_create_dir(out_pod5_dir)

                pod5_conversion(input_dir, out_pod5_dir)

                basecalling_dorado(out_pod5_dir, out_basecalling_dir)

            elif any(file.endswith('.pod5') for file in os.listdir(args.input_dir)):

                basecalling_dorado(input_dir, out_basecalling_dir)

            else:
                sys.exit(f"Error: The files in {input_dir} are not in .fast5 or pod5 format")

            after = datetime.datetime.now()
            print(("Done with function basecalling_dorado in: %s" % (after - prior) + "\n"))


    ### Guppy

    elif args.basecall == 'guppy':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + "STARTING BASECALLING WITH GUPPY" + END_FORMATTING + "\n")

        if any(file.endswith('.fast5') for file in os.listdir(args.input_dir)):

            if os.path.isfile(basecalling_summary):
                logger.info("\n" + YELLOW + BOLD + "Ommiting BASECALLING" + END_FORMATTING + "\n")

            else:
                basecalling_ion(input_dir, out_basecalling_dir, config=args.model, records=args.records_per_fastq)

            for root, _, files in os.walk(out_basecalling_dir):
                for name in files:
                    if name.startswith('guppy_basecaller_log'):
                        log_file = os.path.join(out_basecalling_dir, name)
                        os.remove(log_file)

        else:
            sys.exit(f"Error: The files in {input_dir} are not in .fast5 format")

        after = datetime.datetime.now()
        print(("Done with function basecalling_ion in: %s" % (after - prior) + "\n"))


    # Barcoding #

    ## Dorado

    if args.basecall == 'dorado':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + BOLD + "STARTING BARCODING" + END_FORMATTING)

        if any(f.endswith('sequencing_summary.txt') for f in os.listdir(out_barcoding_dir)):
            logger.info("\n" + YELLOW + BOLD +
                        "Ommiting BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
        else:
            logger.info("\n" + GREEN +
                        "STARTING BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
            demux_dorado(fastq_basecalled, out_barcoding_dir, barcode_kit=args.barcode_kit, require_barcodes_both_ends=args.require_barcodes_both_ends)
            reorganize_demux(out_barcoding_dir)

        after = datetime.datetime.now()
        print(("Done with function demux_dorado & reorganize_demux in: %s" % (after - prior) + "\n"))


    ## Guppy

    elif args.basecall == 'guppy':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + BOLD + "STARTING BARCODING" + END_FORMATTING)

        if os.path.isfile(barcoding_summary):
            logger.info("\n" + YELLOW + BOLD +
                        "Ommiting BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
        else:
            logger.info("\n" + GREEN +
                        "STARTING BARCODING/DEMULTIPLEX" + END_FORMATTING + "\n")
            barcoding_ion(out_basecalling_dir, out_barcoding_dir, barcode_kit=args.barcode_kit,
                        threads=args.threads, require_barcodes_both_ends=args.require_barcodes_both_ends)

        after = datetime.datetime.now()
        print(("Done with function barcoding_ion in: %s" % (after - prior) + "\n"))


    # Read Filtering

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD +
                "STARTING SAMPLE FILTERING" + END_FORMATTING)

    if args.samples == None:
        logger.info(
            '\n' + GREEN + 'Filtering samples' + END_FORMATTING)
        for root, _, files in os.walk(out_barcoding_dir):
            for subdirectory in _:
                if subdirectory.startswith('barcode'):
                    barcode_dir = os.path.join(root, subdirectory)
                    for root2, _, files2 in os.walk(barcode_dir):
                        if len(files2) > 1:
                            barcode_path = root2
                            sample = barcode_path.split('/')[-1]
                            # print(sample)
                            output_samples = os.path.join(
                                out_samples_dir, sample + '.fastq')
                            # print(output_samples)
                            filtered_samples = os.path.join(
                                out_samples_filtered_dir, sample + '.fastq.gz')
                            # print(filtered_samples)

                            logger.info('\n' + BLUE + BOLD +
                                        sample + END_FORMATTING)

                            sum_files = []
                            for name in files2:
                                filename = os.path.join(barcode_path, name)
                                # print(filename)
                                sum_files.append(filename)
                            logger.info(MAGENTA + BOLD + "Processing {} files in {}".format(
                                len(sum_files), sample) + END_FORMATTING)

                            if os.path.isfile(filtered_samples):
                                logger.info(
                                    YELLOW + sample + ' sample already renamed and filtered' + END_FORMATTING)
                            else:
                                logger.info(
                                    GREEN + 'Renaming sample ' + sample + END_FORMATTING)
                                rename_files(output_samples)

                                logger.info(
                                    GREEN + 'Filtering sample ' + sample + END_FORMATTING)
                                ONT_QC_filtering(
                                    output_samples, filtered_samples)
                        else:
                            None

    else:
        logger.info(
            '\n' + GREEN + 'Filtering & Renaming' + END_FORMATTING)
        with open(args.samples, 'r') as f:
            for line in f:
                barcode, sample = line.split('\t')
                sample = sample.strip()

                barcode_path = None
                for root, dirs, _ in os.walk(out_barcoding_dir):
                    for d in dirs:
                        if d == barcode:
                            barcode_path = os.path.join(root, d)
                            break
                    if barcode_path:
                        break

                if barcode_path is None:
                    logger.warning(YELLOW + 'Barcode folder not found for: ' + barcode + END_FORMATTING)
                    continue

                output_samples = os.path.join(out_samples_dir, sample + '.fastq')
                filtered_samples = os.path.join(out_samples_filtered_dir, sample + '.fastq.gz')

                logger.info('\n' + BLUE + BOLD + sample + END_FORMATTING)

                sum_files = []
                for root, _, files in os.walk(barcode_path):
                    for name in files:
                        sum_files.append(os.path.join(root, name))

                logger.info(MAGENTA + BOLD + "Processing {} files in {}".format(
                    len(sum_files), sample) + END_FORMATTING)

                if os.path.isfile(filtered_samples):
                    logger.info(YELLOW + sample + ' sample already renamed and filtered' + END_FORMATTING)
                else:
                    logger.info(GREEN + 'Renaming sample ' + sample + END_FORMATTING)
                    rename_files(output_samples)
                    logger.info(GREEN + 'Filtering sample ' + sample + END_FORMATTING)
                    ONT_QC_filtering(output_samples, filtered_samples)

    after = datetime.datetime.now()
    print(('\n' + "Done with function rename_files & ONT_QC_filtering in: %s" %
           (after - prior) + "\n"))

    # Quality Check

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD +
                "QUALITY CHECK IN RAW" + END_FORMATTING + '\n')

    for root, _, files in os.walk(out_samples_filtered_dir):
        for name in files:
            filtered_sample = os.path.join(root, name)
            # print(filtered_sample)
            out_qc = os.path.join(
                out_qc_dir, os.path.basename(filtered_sample.split(".")[0]))
            check_create_dir(out_qc)
            # print(out_qc)
            report = [x for x in os.listdir(out_qc) if "NanoPlot-report" in x]
            report_file = os.path.join(out_qc, "".join(report))
            # print(report)

            if os.path.isfile(report_file):
                logger.info(YELLOW + report_file +
                            " EXIST\nOmmiting QC for sample " + name + END_FORMATTING)
            else:
                logger.info(GREEN + "Checking quality in sample " +
                            name + END_FORMATTING)
                ONT_quality(filtered_sample, out_qc, threads=args.threads)

    after = datetime.datetime.now()
    print(("\n" + "Done with function ONT_quality in: %s" % (after - prior) + "\n"))

    logger.info("\n" + MAGENTA + BOLD +
                "#####END OF ONT DATA PROCESSING PIPELINE #####" + END_FORMATTING + "\n")
