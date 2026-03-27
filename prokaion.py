#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging
import shutil


# Third party imports
import argparse
import subprocess
import datetime
import gzip
import multiprocessing
import pandas as pd


# Local application imports

from misc_prokaion import (check_create_dir, check_file_exists, extract_read_list, extract_sample_list, execute_subprocess, check_reanalysis, file_to_list, samtools_faidx, create_reference_chunks, extract_indels, merge_vcf, vcf_to_ivar_tsv, create_bamstat,
                           create_coverage, obtain_group_cov_stats, obtain_overal_stats, ivar_consensus, replace_consensus_header, remove_low_quality, rename_reference_snpeff, annotate_snpeff, user_annotation, user_annotation_aa, make_blast, kraken, mash_screen)

from compare_snp_prokaion import (ddbb_create_intermediate, recalibrate_ddbb_vcf_intermediate,
                                  remove_position_range, extract_complex_list, revised_df, ddtb_compare, remove_bed_positions, extract_only_snps, extract_bed_positions, extract_lowcov)


"""
=============================================================
HEADER
=============================================================
Institution: IiSGM
Author: Sergio Buenestado-Serrano (sergio.buenestado@gmail.com) & Pedro J. Sola (pedroscampoy@gmail.com)
Version = 0
Created: 22 November 2021

TODO:
    Check program is installed (dependencies)
================================================================
END_OF_HEADER
================================================================
"""


# COLORS AND AND FORMATTING

"""
http://ozzmaker.com/add-colour-to-text-in-python/
The above ANSI escape code will set the text colour to bright green. The format is;
\033[  Escape code, this is always the same
1 = Style, 1 for normal.
32 = Text colour, 32 for bright green.
40m = Background colour, 40 is for black.
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE = '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

logger = logging.getLogger()


# ARGUMENTS

def get_arguments():

    parser = argparse.ArgumentParser(
        prog='prokaion.py', description='Pipeline to call variants (SNVs) with any non model organism. Specialised in Mycobacterium tuberculosis (MTB)')

    input_group = parser.add_argument_group('Input', 'Input parameters')

    input_group.add_argument('-i', '--input', dest='input_dir', metavar='input_directory',
                             type=str, required=True, help='REQUIRED. Input directory containing all fast5 files')

    input_group.add_argument('-o', '--output', type=str, required=True,
                             help='REQUIRED. Output directory to extract all results')

    input_group.add_argument('-t', '--threads', type=int, dest='threads',
                             required=False, default=30, help='Threads to use. Default: 30)')

    basecall_group = parser.add_argument_group('Basecall', 'Basecall parameters')

    basecall_group.add_argument('-s', '--samples', metavar='Samples', type=str,
                             required=False, help='Sample list for conversion from barcode to samples ID')

    basecall_group.add_argument('--basecall', type=str, default='dorado', required=True, help="REQUIRED. Program to use in data preprocessing (basecaller and barcoder)")

    basecall_group.add_argument('--model', type=str, default='~/Scripts/git_repos/Dorado/Models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0', required=False, help='The basecaller model to run. Default dna_r10.4.1_e8.2_400bps_hac@v4.3.0 | dna_r10.4_e8.1_hac.cfg')

    basecall_group.add_argument('-b', '--require_barcodes_both_ends', required=False, action='store_true',
                             help='Require barcodes at both ends. By default it only requires the barcode at one end for the sequences identification')

    basecall_group.add_argument('--barcode_kit', type=str, required=False, default='SQK-RBK110-96',
                             help='Kit of barcodes used [SQK-RBK110-96|EXP-NBD196|SQK-NBD112-24]. Default: SQK-RBK110-96')

    basecall_group.add_argument('-g', '--gpu', dest='gpu', required=False, default=False,
                             action="store_true", help='Specify GPU device: "auto", or "cuda:<device_id>"')

    basecall_group.add_argument('--num_callers', type=int, dest='num_callers',
                             required=False, default=8, help="Number of parallel basecallers. Default: 8")

    basecall_group.add_argument('--chunks', type=int, dest='chunks', required=False,
                             default=1536, help='Maximum chunks per runner. Default: 1536')

    basecall_group.add_argument('--records_per_fastq', type=int, dest='records_per_fastq',
                             required=False, default=0, help='Maximum number of records per fastq')

    basecall_group.add_argument("-rq", "--min_read_quality", type=int, dest="min_read_quality",
                             required=False, default=15, help="Filter on a minimum average read quality score. Default: 15")

    basecall_group.add_argument("--trim", type=str, dest="trim",
                             required=False, default="adapters", help="Specify what to trim. Options are 'none', 'all', and 'adapters'. The default behaviour is to trim all detected adapters, primers, and barcodes. Choose 'adapters' to just trim adapters. The 'none' choice is equivelent to using --no-trim. Note that this only applies to DNA. RNA adapters are always trimmed. Default: adapters")

    basecall_group.add_argument("--headcrop", type=int, dest="headcrop", required=False,
                             default=20, help="Trim n nucleotides from start of read. Default: 20")

    basecall_group.add_argument("--tailcrop", type=int, dest="tailcrop", required=False,
                             default=20, help="Trim n nucleotides from end of read. Default: 20")

    basecall_group.add_argument("--length", type=int, dest="length", required=False,
                             default=150, help="Filter on a minimum read length. Default: 150bp")

    varcal_group = parser.add_argument_group('Varcal', 'Varcal parameters')

    varcal_group.add_argument("-sample", "--sample", metavar="sample",
                              type=str, required=False, help="Sample to identify further files")

    varcal_group.add_argument("-L", "--sample_list", type=str, required=False,
                              help="Sample names to analyse only in the file supplied")

    varcal_group.add_argument("-bayes", "--bayes", required=False, action="store_true",
                              help="Variant Calling is done with freebayes-parallel")

    varcal_group.add_argument('-r', '--reference', metavar="reference",
                              type=str, required=True, help='REQUIRED. File to map against')

    varcal_group.add_argument("--ploidy", type=int, dest="ploidy", required=False,
                              default=1, help="Sets the default ploidy for the analysis")

    varcal_group.add_argument("-amplicon", "--amplicon", required=False, action="store_true",
                              help="Mapping performed using ngmlr")

    species_group = parser.add_argument_group(
        "Species determination", "Species databases")

    species_group.add_argument("--kraken2", dest="kraken2_db", type=str,
                               default=False, required=False, help="Kraken2 database")

    species_group.add_argument("--mash_db", dest="mash_db", type=str, required=False,
                               default=False, help="MASH NCBI annotation containing bacterial database")

    variant_group = parser.add_argument_group(
        "Variant Calling", "Variant Calling parameters")

    variant_group.add_argument("-f", "--min_allele_frequency", type=int, dest="min_allele_frequency", required=False,
                               default=0.2, help="Minimum fraction of observations supporting an alternate allele. Default: 0.2")

    variant_group.add_argument("-q", "--min_base_quality", type=int, dest="min_quality", required=False,
                               default=12, help="Exclude alleles from analysis below threshold. Default: 12")

    variant_group.add_argument("-m", "--min_mapping_quality", type=int, dest="min_mapping", required=False,
                               default=60, help="Exclude alignments from analysis below threshold. Default: 60")

    variant_group.add_argument("-freq", "--min_frequency", type=int, dest="min_frequency", required=False,
                               default=0.7, help="Minimum fraction of observations to call a base. Default: 0.7")

    variant_group.add_argument("-d", "--min_depth", type=int, dest="min_depth",
                               required=False, default=8, help="Minimum depth to call a base. Default: 8")

    quality_group = parser.add_argument_group(
        'Quality parameters', "Parameters for diferent Quality conditions")

    quality_group.add_argument('-c', '--coverage20', type=int, default=70, required=False,
                               help='Minimum percentage of coverage at 20x to classify as uncovered. Default: 70')

    quality_group.add_argument('-u', '--unmapped', type=int, default=25, required=False,
                               help='Minimum percentage of unmapped reads to classify as uncovered. Default: 25')

    quality_group.add_argument('-n', '--min_snp', type=int, required=False,
                               default=75, help='SNP number to pass quality threshold. Default: 75 HQ SNP')

    annot_group = parser.add_argument_group(
        'Annotation', 'Parameters for variant annotation')

    annot_group.add_argument('-B', '--annot_bed', type=str, default=[],
                             required=False, action='append', help='BED file to annotate')

    annot_group.add_argument('-V', '--annot_vcf', type=str, default=[],
                             required=False, action='append', help='VCF file to annotate')

    annot_group.add_argument("-F", "--annot_fasta", type=str, default=[],
                             required=False, action="append", help="FASTA file to annotate")

    annot_group.add_argument('-A', '--annot_aa', type=str, default=[],
                             required=False, action='append', help='Aminoacid file to annotate')

    annot_group.add_argument('-R', '--remove_bed', type=str, default=False,
                             required=False, help='BED file with positions to remove')

    annot_group.add_argument('-P', '--extract_bed', required=False, type=str,
                             default=False, help='BED file with important positions or genes to annotate')

    annot_group.add_argument('--snpeff_database', type=str, required=False,
                             default=False, help='snpEFF annotation database')

    compare_group = parser.add_argument_group(
        'Compare', 'parameters for compare_snp')

    compare_group.add_argument('-S', '--only_snp', required=False,
                               action='store_true', help='Create the results only with SNPs, removing INDELs')

    compare_group.add_argument("-w", "--window", required=False, type=int,
                               default=2, help="Number of SNPs in 10bp to discard. Default: 2")

    compare_group.add_argument("--min_threshold_discard_uncov_sample", required=False, type=float,
                               default=0.6, help="Minimum uncovered genome to discard a sample. Default: 0.6")

    compare_group.add_argument("--min_threshold_discard_uncov_pos", required=False,
                               type=float, default=0.55, help="Minimum covered position to discard it. Default: 0.55")

    compare_group.add_argument("--min_threshold_discard_htz_sample", required=False, type=float,
                               default=0.6, help="Minimum heterozygosity to discard a sample. Default: 0.6")

    compare_group.add_argument("--min_threshold_discard_htz_pos", required=False, type=float,
                               default=0.5, help="Minimum heterozygosity to discard a position. Default: 0.5")

    compare_group.add_argument("--min_threshold_discard_all_sample", required=False, type=float,
                               default=0.6, help="Minimum inaccuracies to discard a sample. Default: 0.6")

    compare_group.add_argument("--min_threshold_discard_all_pos", required=False, type=float,
                               default=0.55, help="Minimum inaccuracies to discard a position. Default: 0.55")

    compare_group.add_argument('-D', '--distance', type=int, default=50, required=False,
                               help='Minimun distance to cluster groups after comparison')

    arguments = parser.parse_args()

    return arguments


### Functions from Guppy_prokaion.py ###

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
    #     gpu_device = '-x auto '
    # else:
    #     gpu_device = ''

    if require_barcodes_both_ends:
        logger.info(
            GREEN + BOLD + "Barcodes are being used at both ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = "--require_barcodes_both_ends"
    else:
        logger.info(
            YELLOW + BOLD + "Barcodes are being used on at least 1 of the ends" + END_FORMATTING + "\n")
        require_barcodes_both_ends = ""

    cmd = ["guppy_barcoder", "-i", out_basecalling_dir, "-s", out_barcoding_dir, "-r", require_barcodes_both_ends,
           "--barcode_kits", barcode_kit, "-t", str(threads), '--detect_barcodes', '--min_score_barcode_front', str(85), '--min_score_barcode_rear', str(85), '--enable_trim_barcodes', '--detect_mid_strand_barcodes', '--min_score_barcode_mid', str(90), '--detect_primer', '--trim_primers', '--detect_adapter', '--trim_adapters']

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


def ONT_quality(output_samples, out_qc, threads=30):

    # --fastq_rich: Data is in one or more fastq file(s) generated by albacore, MinKNOW or guppy with additional information concerning with channel and time
    # --N50: Show the N50 mark in the read length histogram

    cmd_QC = ["NanoPlot", "--fastq", output_samples,
              "--N50", "-o", out_qc, "-t", str(threads)]

    # print(cmd_QC)
    execute_subprocess(cmd_QC, isShell=False)


### Functions from Varcal_prokaion.py ###

def minimap2_mapping(HQ_filename, filename_bam_out, reference):
    """
    https://github.com/lh3/minimap2
        # Oxford Nanopore genomic reads
        minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam
    http://www.htslib.org/doc/samtools.html
    """

    # -a: Output in the SAM format
    # -x: Preset (always applied before other options; see minimap2.1 for details) []
    #    - map-pb/map-ont - PacBio CLR/Nanopore vs reference mapping
    #    - map-hifi - PacBio HiFi reads vs reference mapping
    #    - ava-pb/ava-ont - PacBio/Nanopore read overlap
    #    - asm5/asm10/asm20 - asm-to-ref mapping, for ~0.1/1/5% sequence divergence
    #    - splice/splice:hq - long-read/Pacbio-CCS spliced alignment
    #    - sr - genomic short-read mapping
    # -t: Number of threads

    # -b: Output BAM
    # -S: Ignored (input format is auto-detected)
    # -F: Only include reads with none of the FLAGS in INT present
    # --threads: Number of additional threads to use

    cmd_minimap2 = "minimap2 -ax map-ont {} {} | samtools view -bS -F 4 - | samtools sort -o {}".format(
        reference, HQ_filename, filename_bam_out)
    # print(cmd_minimap2)
    execute_subprocess(cmd_minimap2, isShell=True)

    cmd_indexing = "samtools", "index", filename_bam_out
    # print(cmd_indexing)
    execute_subprocess(cmd_indexing, isShell=False)


def ngmlr_mapping(HQ_filename, filename_bam_out, reference, threads=30):
    """
    https://github.com/philres/ngmlr
    """

    # -r <file>, --reference <file>
    # -q <file>, --query <file>
    # -o <string>, --output <string>
    # -t <int>, --threads <int>
    # -x <pacbio, ont>, --presets <pacbio, ont>
    # -i <0-1>, --min-identity <0-1> - Alignments with an identity lower than this threshold will be discarded [0.65]

    cmd_ngmlr = "ngmlr -t {} -r {} -q {} -x ont | samtools view -bS -F 4 - | samtools sort -o {}".format(
        threads, reference, HQ_filename, filename_bam_out)
    # print(cmd_ngmlr)
    execute_subprocess(cmd_ngmlr, isShell=True)

    cmd_indexing = "samtools", "index", filename_bam_out
    # print(cmd_indexing)
    execute_subprocess(cmd_indexing, isShell=False)


def freebayes_variant(reference, filename_bam_out, output_vcf, threads=36, frequency=0.1, ploidy=1, base_qual=7, map_qual=60):
    """
    https://github.com/freebayes/freebayes
        # Freebayes-parallel
        freebayes-parallel <(fasta_generate_regions.py {fai_reference} {chunks}) {threads} {args} > {output}
    """

    # --region_file: Genome partitioning according to computational threads.
    # --haplotype_length: Allow haplotype calls with contiguous embedded matches of up to this length
    # --use-best-n-alleles: Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores
    # --min-alternate-count: Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position
    # --min-alternate-fraction: Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position
    # -p: Sets the default ploidy for the analysis to N
    # --min-coverage:
    # -q: Exclude alleles from analysis if their supporting base quality is less than Q
    # -m: Exclude alignments from analysis if they have a mapping quality less than Q
    # --strict-vcf: Generate strict VCF format (FORMAT/GQ will be an int)

    region_file = create_reference_chunks(reference)

    cmd_bayes = "freebayes-parallel {} {} -f {} --haplotype-length 0 --use-best-n-alleles 1 --min-alternate-count 0 --min-alternate-fraction {} -p {} --min-coverage 1 -q {} -m {} --strict-vcf {} > {}".format(
        region_file, str(threads), reference, str(frequency), str(ploidy), str(base_qual), str(map_qual), filename_bam_out, output_vcf)
    print(cmd_bayes)
    execute_subprocess(cmd_bayes, isShell=True)


def freebayes(reference, filename_bam_out, output_vcf, frequency=0.1, ploidy=1, base_qual=7, map_qual=60):

    cmd_bayes = "freebayes -f {} --haplotype-length 0 --use-best-n-alleles 1 --min-alternate-count 0 --min-alternate-fraction {} -p {} --min-coverage 1 -q {} -m {} --strict-vcf {} > {}".format(
        reference, str(frequency), str(ploidy), str(base_qual), str(map_qual), filename_bam_out, output_vcf)
    # print(cmd_bayes)
    execute_subprocess(cmd_bayes, isShell=True)


def bcftool_filter(output_raw_vcf, output_vcf):
    """
    https://samtools.github.io/bcftools/bcftools.html
        # bcftools view: View, subset and filter VCF files by position and filtering expression
        # bcftools annotate: Add or remove annotations
    """

    # --include: include sites for which Expression is true
    # --remove: list of annotations to remove

    # cmd_bcf_view = "bcftools view --include 'QUAL >= 10 && FMT/DP >= 1 && (FMT/AO)/(FMT/DP) >= 0.1' {} -o {}".format(output_raw_vcf, output_vcf)
    # # print(cmd_bcf_view)
    # execute_subprocess(cmd_bcf_view, isShell=True)

    # cmd_bcf_annot = 'bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL {} -o {}'.format(output_vcf, output_vcf_filt)
    # # print(cmd_bcf_annot)
    # execute_subprocess(cmd_bcf_annot, isShell=True)

    cmd_bcf = "bcftools view --include 'QUAL >= 10 && FMT/DP >= 1 && (FMT/AO)/(FMT/DP) >= 0.1' {} | bcftools annotate --remove ^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL -o {}".format(
        output_raw_vcf, output_vcf)
    # print(cmd_bcf)
    execute_subprocess(cmd_bcf, isShell=True)


def snippy_sub(output_vcf_filt, output_vcf_sub):
    """
    https://github.com/tseemann/snippy/tree/master/bin
        # snippy-vcf_extract_subs: Convert MNP,COMPLEX into SNP and ignore INS,DEL
    """

    cmd_snippy_subs = "snippy-vcf_extract_subs {} > {}".format(
        output_vcf_filt, output_vcf_sub)
    # print(cmd_snippy_subs)
    execute_subprocess(cmd_snippy_subs, isShell=True)


######################################################################
########################### START PIPELINE ###########################
######################################################################


if __name__ == "__main__":

    args = get_arguments()

    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output)
    group_name = output_dir.split("/")[-1]
    check_create_dir(output_dir)

    reference = os.path.abspath(args.reference)

    # Declare folders created in pipeline and key files

    out_basecalling_dir = os.path.join(input_dir, "Basecalling")
    check_create_dir(out_basecalling_dir)
    basecalling_summary = os.path.join(
        out_basecalling_dir, "sequencing_summary.txt")
    out_barcoding_dir = os.path.join(input_dir, "Barcoding")
    check_create_dir(out_barcoding_dir)
    barcoding_summary = os.path.join(
        out_barcoding_dir, "barcoding_summary.txt")
    out_samples_dir = os.path.join(input_dir, "Samples_Fastq")
    check_create_dir(out_samples_dir)
    fastq_basecalled = os.path.join(out_basecalling_dir, 'Dorado_basecalled.fastq')
    # out_samples_filtered_dir = os.path.join(out_samples_dir, "Filtered_Fastq")
    # check_create_dir(out_samples_filtered_dir)
    # out_correction_dir = os.path.join(out_samples_dir, 'Corrected')
    # check_create_dir(out_correction_dir)

    out_qc_dir = os.path.join(output_dir, "Quality")
    check_create_dir(out_qc_dir)

    out_species_dir = os.path.join(output_dir, "Species")
    check_create_dir(out_species_dir)

    out_bam_dir = os.path.join(output_dir, "Bam")
    check_create_dir(out_bam_dir)

    out_variant_dir = os.path.join(output_dir, "Variants")
    check_create_dir(out_variant_dir)

    out_consensus_dir = os.path.join(output_dir, "Consensus")
    check_create_dir(out_consensus_dir)

    out_stats_dir = os.path.join(output_dir, "Stats")
    check_create_dir(out_stats_dir)
    out_stats_bamstats_dir = os.path.join(
        out_stats_dir, "Bamstats")  # subfolder
    check_create_dir(out_stats_bamstats_dir)
    out_stats_coverage_dir = os.path.join(
        out_stats_dir, "Coverage")  # subfolder
    check_create_dir(out_stats_coverage_dir)

    out_annot_dir = os.path.join(output_dir, "Annotation")
    check_create_dir(out_annot_dir)
    out_annot_snpeff_dir = os.path.join(out_annot_dir, "snpeff")  # subfolder
    check_create_dir(out_annot_snpeff_dir)
    out_annot_user_dir = os.path.join(out_annot_dir, "user")  # subfolder
    check_create_dir(out_annot_user_dir)
    out_annot_user_aa_dir = os.path.join(out_annot_dir, "user_aa")  # subfolder
    check_create_dir(out_annot_user_aa_dir)
    out_annot_blast_dir = os.path.join(out_annot_dir, "blast")  # subfolder
    check_create_dir(out_annot_blast_dir)

    out_compare_dir = os.path.join(output_dir, "Compare")
    check_create_dir(out_compare_dir)

    samtools_faidx(args.reference)

    create_reference_chunks(args.reference)

    # Logging
    # Create log file with date and time

    today = str(datetime.date.today())
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

    # Dorado

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

    # Guppy

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


    # Barcoding

    ## Dorado

    if args.basecall == 'dorado':

        prior = datetime.datetime.now()

        logger.info("\n" + GREEN + BOLD + "STARTING BARCODING" + END_FORMATTING)

        if any(f.endswith('unclassified.fastq') for f in os.listdir(out_barcoding_dir)):
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
                                output_dir, sample + '.fastq.gz')
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
                # print(barcode,sample)
                barcode_path = os.path.join(out_barcoding_dir, barcode)
                # print(barcode_path)
                output_samples = os.path.join(
                    out_samples_dir, sample.strip() + '.fastq')
                # print(output_samples)
                filtered_samples = os.path.join(
                    output_dir, sample.strip() + '.fastq.gz')

                logger.info('\n' + BLUE + BOLD + sample + END_FORMATTING)

                sum_files = []
                for root, _, files in os.walk(barcode_path):
                    for name in files:
                        filename = os.path.join(barcode_path, name)
                        # print(filename)
                        sum_files.append(filename)
                    logger.info(MAGENTA + BOLD + "Processing {} files in {}".format(
                        len(sum_files), sample) + END_FORMATTING)

                    if os.path.isfile(filtered_samples):
                        logger.info(
                            YELLOW + sample + ' sample already renamed and filtered' + END_FORMATTING)
                    else:
                        logger.info(GREEN + 'Renaming sample ' +
                                    sample + END_FORMATTING)
                        rename_files(output_samples)

                        logger.info(GREEN + 'Filtering sample ' +
                                    sample + END_FORMATTING)
                        ONT_QC_filtering(output_samples, filtered_samples)

    after = datetime.datetime.now()
    print(('\n' + "Done with function rename_files & ONT_QC_filtering in: %s" %
           (after - prior) + "\n"))

    # Quality Check

    prior = datetime.datetime.now()

    logger.info("\n" + GREEN + BOLD +
                "QUALITY CHECK IN RAW" + END_FORMATTING + '\n')

    for root, _, files in os.walk(output_dir):
        for name in files:
            if name.endswith('.fastq.gz'):
                # raw_sample = os.path.join(root, name) ## (output_samples)
                filtered_samples = os.path.join(root, name) ## (output_dir)
                # print(raw_sample)
                out_qc = os.path.join(
                    out_qc_dir, os.path.basename(filtered_samples.split(".")[0]))
                check_create_dir(out_qc)
                # print(out_qc)
                report = [x for x in os.listdir(
                    out_qc) if "NanoPlot-report" in x]
                report_file = os.path.join(out_qc, "".join(report))
                # print(report)

                if os.path.isfile(report_file):
                    logger.info(YELLOW + report_file +
                                " EXIST\nOmmiting QC for sample " + name + END_FORMATTING)
                else:
                    logger.info(GREEN + "Checking quality in sample " +
                                name + END_FORMATTING)
                    ONT_quality(filtered_samples, out_qc, threads=args.threads)

    after = datetime.datetime.now()
    print(("\n" + "Done with function ONT_quality in: %s" % (after - prior) + "\n"))

    logger.info("\n" + MAGENTA + BOLD +
                "#####END OF ONT DATA PROCESSING PIPELINE #####" + END_FORMATTING + "\n")

    logger.info(
        "\n" + BLUE + "############### START VARIANT CALLING ###############" + END_FORMATTING + "\n")

    logger.info(args)

    # Obtain all fastq files from folder

    fastq = extract_read_list(output_dir)

    # Check how many files will be analysed

    sample_list = []

    for sample in fastq:
        sample = extract_sample_list(sample)
        # sample = sample.split('_')[1]
        sample_list.append(sample)

    # logger.info('\n' + CYAN + '{} Samples will be analysed: {}'.format(
    #     len(sample_list), ', '.join(sample_list)) + END_FORMATTING)

    # Check if there are samples to filter out

    sample_list_F = []
    if args.sample_list == None:
        logger.info("\n" + "No samples to filter" + "\n")
        for sample in fastq:
            sample = extract_sample_list(sample)
            sample_list_F.append(sample)
    else:
        logger.info("Samples will be filtered")
        sample_list_F = file_to_list(args.sample_list)

    new_samples = check_reanalysis(args.output, sample_list_F)

    logger.info(CYAN + "\n%d samples will be analysed: %s" %
                (len(sample_list_F), ",".join(sample_list_F)) + END_FORMATTING + '\n')

    logger.info(CYAN + "\n%d NEW samples will be analysed: %s" %
                (len(new_samples), ",".join(new_samples)) + END_FORMATTING + '\n')

    new_sample_number = 0

    for sample in fastq:
        # Extract sample name
        sample = extract_sample_list(sample)
        args.sample = sample

        if sample in sample_list_F:
            sample_number = str(sample_list_F.index(sample) + 1)
            sample_total = str(len(sample_list_F))

            if sample in new_samples:
                new_sample_number = str(int(new_sample_number) + 1)
                new_sample_total = str(len(new_samples))
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample + " (" + sample_number + "/" +
                            sample_total + ")" + " (" + new_sample_number + "/" + new_sample_total + ")" + END_FORMATTING)
            else:
                logger.info("\n" + WHITE_BG + "STARTING SAMPLE: " + sample +
                            " (" + sample_number + "/" + sample_total + ")" + END_FORMATTING)

            # if not os.path.isfile(output_final_vcf):
            HQ_filename = os.path.join(output_dir, sample + ".fastq.gz")
            # print(HQ_filename)
            # filename_out = sample.split('.')[0].split('_')[1]
            filename_out = sample
            # print(filename_out)

            logger.info("\n" + GREEN + BOLD +
                        "STARTING ANALYSIS FOR SAMPLE " + filename_out + END_FORMATTING + '\n')

            ##### SPECIES DETERMINATION #####

            prior = datetime.datetime.now()

            # Species determination with kraken2 and its standard database and visualization with ktImportTaxonomy from kronatools kit

            sample_species_dir = os.path.join(out_species_dir, sample)
            # print(sample_species_dir)
            check_create_dir(sample_species_dir)
            report = os.path.join(sample_species_dir, sample)
            krona_html = os.path.join(report + ".html")
            mash_output = os.path.join(report + ".screen.tab")

            if args.kraken2_db != False:
                if os.path.isfile(krona_html):
                    logger.info(
                        YELLOW + krona_html + " EXIST\nOmmiting species determination with Kraken2 for " + sample + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Species determination with Kraken2 for sample " + sample + END_FORMATTING)
                    kraken(HQ_filename, report, args.kraken2_db,
                           krona_html, threads=args.threads)
            else:
                logger.info(
                    YELLOW + BOLD + "No Kraken database suplied, skipping specie assignation in group " + group_name + END_FORMATTING)

            # Species determination with mash and its bacterial database

            if args.mash_db != False:
                if os.path.isfile(mash_output):
                    logger.info(
                        YELLOW + mash_output + " EXIST\nOmmiting species determination with Mash screen for " + sample + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Species determination with Mash for sample " + sample + END_FORMATTING)

                    mash_screen(HQ_filename, mash_output,
                                args.mash_db, winner=True, threads=args.threads)

                    # Name the columns of the mash output and sort them in descending order by identity
                    output_sort_species = pd.read_csv(mash_output, sep='\t', header=None, names=[
                                                      'Identity', 'Share-hashes', 'Median-multiplicity', 'p-value', 'ID accession', 'Organism']).sort_values(by=['Identity'], ascending=False)
                    output_sort_species.to_csv(
                        mash_output, sep='\t', index=None)
            else:
                logger.info(
                    YELLOW + BOLD + "No MASH database suplied, skipping specie assignation in group " + group_name + END_FORMATTING)

            after = datetime.datetime.now()
            print(("Done with function kraken & mash_screen in: %s" %
                   (after - prior) + "\n"))

            ##### MAPPING #####

            # Mapping with minimap2, sorting Bam and indexing it (also can be made with bwa index & bwa mem -x ont2d)

            prior = datetime.datetime.now()

            filename_bam_out = os.path.join(
                out_bam_dir, filename_out + ".sort.bam")
            filename_bai_out = os.path.join(
                out_bam_dir, filename_out + ".sort.bam.bai")
            # print(filename_bam_out)

            if os.path.isfile(filename_bai_out):
                logger.info(YELLOW + filename_bam_out +
                            " EXIST\nOmmiting mapping for " + filename_out + END_FORMATTING)
            else:
                logger.info(GREEN + "Mapping sample " +
                            filename_out + END_FORMATTING)

                if args.amplicon:
                    ngmlr_mapping(HQ_filename, filename_bam_out,
                                  reference, threads=args.threads)
                else:
                    minimap2_mapping(HQ_filename, filename_bam_out,
                                     reference=args.reference)

            after = datetime.datetime.now()
            print(f"Done with function {'ngmlr_mapping' if args.amplicon else 'minimap2_mapping'} in: %s" % (
                after - prior) + "\n")

            ##### VARIANT CALLING #####

            # Variant calling with freebayes-parallel (also can be made with nanopolish, we should use nanopolish index & nanopolish variants)

            if filename_bam_out.endswith("bam"):
                # print(filename_bam_out)
                sample_variant_dir = os.path.join(
                    out_variant_dir, filename_out)
                # print(sample_variant_dir)
                check_create_dir(sample_variant_dir)
                output_raw_vcf = os.path.join(
                    sample_variant_dir, "snps.raw.vcf")
                # print(output_raw_vcf)

                prior = datetime.datetime.now()

                if os.path.isfile(output_raw_vcf):
                    logger.info(
                        YELLOW + output_raw_vcf + " EXIST\nOmmiting Variant Calling for sample " + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Starting Variant Calling for sample " + filename_out + END_FORMATTING)
                    if args.amplicon:
                        freebayes(args.reference, filename_bam_out, output_raw_vcf,
                                      frequency=args.min_allele_frequency, ploidy=args.ploidy, base_qual=args.min_quality, map_qual=args.min_mapping)
                    else:
                        freebayes_variant(args.reference, filename_bam_out, output_raw_vcf, threads=args.threads,
                                        frequency=args.min_allele_frequency, ploidy=args.ploidy, base_qual=args.min_quality, map_qual=args.min_mapping)

                after = datetime.datetime.now()
                print(("Done with function freebayes_variant in: %s" %
                       (after - prior) + "\n"))

                # Filtering the raw variant calling by quality, depth and frequency with bcftools. Also extracting complex variations and MNP with snippy-vcf_extract_subs

                output_vcf = os.path.join(sample_variant_dir, "snps.vcf")
                # print(output_vcf)
                output_vcf_filt = os.path.join(
                    sample_variant_dir, "snps.filt.vcf")
                # print(output_vcf_filt)
                output_vcf_sub = os.path.join(
                    sample_variant_dir, "snps.subs.vcf")
                # print(output_vcf_sub)

                prior = datetime.datetime.now()

                if os.path.isfile(output_vcf_sub):
                    logger.info(
                        YELLOW + output_vcf_sub + " EXIST\nOmmiting Variant Calling filter in " + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Variant Calling filtering in sample " + filename_out + END_FORMATTING)
                    bcftool_filter(output_raw_vcf, output_vcf)
                    snippy_sub(output_vcf, output_vcf_sub)

                after = datetime.datetime.now()
                print(("Done with function bcftools_filter & snippy_sub in: %s" % (
                    after - prior) + "\n"))

                # Variant format combination, extracting INDELs and combining with subs.vcf (SNPs, MNPs and complex)

                out_variant_indel_sample = os.path.join(
                    sample_variant_dir, "snps.indel.vcf")
                # print(out_variant_indel_sample)
                out_variant_all_sample = os.path.join(
                    sample_variant_dir, "snps.all.vcf")
                # print(out_variant_all_sample)
                chrom_filename = os.path.join(
                    sample_variant_dir, "snps.all.chromosome.vcf")
                # print(chrom_filename)

                prior = datetime.datetime.now()

                if os.path.isfile(out_variant_indel_sample):
                    logger.info(YELLOW + out_variant_indel_sample +
                                " EXIST\nOmmiting INDEL filtering in " + filename_out + END_FORMATTING)
                else:
                    logger.info(GREEN + "Filtering INDELs in " +
                                filename_out + END_FORMATTING)
                    extract_indels(output_vcf)

                if os.path.isfile(out_variant_all_sample):
                    logger.info(YELLOW + out_variant_all_sample +
                                "EXIST\nOmmiting VCF combination for sample " + filename_out + END_FORMATTING)
                else:
                    logger.info(GREEN + "Combining VCF in " +
                                filename_out + END_FORMATTING)
                    merge_vcf(output_vcf_sub, out_variant_indel_sample)

                after = datetime.datetime.now()
                print(("Done with function extract_indels & merge_vcf in: %s" %
                       (after - prior) + "\n"))

                # Variant format adaptation

                out_variant_tsv_file = os.path.join(
                    sample_variant_dir, "snps.all.ivar.tsv")
                # print(out_variant_tsv_file)

                prior = datetime.datetime.now()

                if os.path.isfile(out_variant_tsv_file):
                    logger.info(YELLOW + out_variant_tsv_file +
                                " EXIST\nOmmiting format adaptation for " + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Adapting variants format in sample " + filename_out + END_FORMATTING)
                    vcf_to_ivar_tsv(out_variant_all_sample,
                                    out_variant_tsv_file)

                after = datetime.datetime.now()
                print(("Done with function vcf_to_ivar_tsv in: %s" %
                       (after - prior) + "\n"))

                ##### CONSENSUS #####

                # Building consensus fasta file with samtools mpileup and ivar consensus

                out_consensus_file = os.path.join(
                    out_consensus_dir, filename_out + ".fa")
                # print(out_consensus_file)

                prior = datetime.datetime.now()

                if os.path.isfile(out_consensus_file):
                    logger.info(YELLOW + out_consensus_file +
                                " EXIST\nOmmiting Consensus for " + filename_out + END_FORMATTING)
                else:
                    logger.info(
                        GREEN + "Creating Consensus in sample " + filename_out + END_FORMATTING)

                    # Find another solution, if we set q>7 an error occur "Segmentation fault", they are trying to fix it.
                    ivar_consensus(filename_bam_out, out_consensus_dir, filename_out, min_quality=5,
                                   min_frequency_threshold=0.7, min_depth=12, uncovered_character="N")
                    replace_consensus_header(out_consensus_file)

                for root, _, files in os.walk(out_consensus_dir):
                    for name in files:
                        if name.endswith('qual.txt'):
                            qual_consensus = os.path.join(
                                out_consensus_dir, name)
                            os.remove(qual_consensus)

                after = datetime.datetime.now()
                print(("Done with function ivar_consensus & replace_consensus_header in: %s" % (
                    after - prior) + "\n"))

        ##### CREATE STATS AND QUALITY FILTERS #####

        # Create Bamstats

        out_bamstats_name = filename_out + ".bamstats"
        out_bamstats_file = os.path.join(
            out_stats_bamstats_dir, out_bamstats_name)
        # print(out_bamstats_file)
        # print(filename_bam_out)

        prior = datetime.datetime.now()

        if os.path.isfile(out_bamstats_file):
            logger.info(YELLOW + out_bamstats_file +
                        "EXIST\nOmmiting Bamstats for " + filename_out + END_FORMATTING)
        else:
            logger.info(GREEN + "Creating Bamstats in sample " +
                        filename_out + END_FORMATTING)
            create_bamstat(filename_bam_out, out_bamstats_file,
                           threads=args.threads)

        after = datetime.datetime.now()
        print(("Done with function create_bamstat in: %s" %
               (after - prior) + "\n"))

        # Create Coverage

        out_coverage_name = filename_out + ".cov"
        out_coverage_file = os.path.join(
            out_stats_coverage_dir, out_coverage_name)
        # print(out_coverage_file)

        prior = datetime.datetime.now()

        if os.path.isfile(out_coverage_file):
            logger.info(YELLOW + out_coverage_file +
                        " EXIST\nOmmiting Coverage for " + filename_out + END_FORMATTING)
        else:
            logger.info(GREEN + "Creating Coverage in sample " +
                        filename_out + END_FORMATTING)
            create_coverage(filename_bam_out, out_coverage_file)

        after = datetime.datetime.now()
        print(("Done with function create_coverage in: %s" %
               (after - prior) + "\n"))

    # Coverage Output summary

    prior = datetime.datetime.now()

    logger.info(GREEN + BOLD + "Creating summary report for coverage results in group " +
                group_name + END_FORMATTING)
    obtain_group_cov_stats(out_stats_dir, group_name)

    # Reads and Variants output summary

    logger.info(GREEN + BOLD + "Creating overal summary report in group " +
                group_name + END_FORMATTING)
    obtain_overal_stats(out_stats_dir, output_dir, group_name)

    after = datetime.datetime.now()
    print(("Done with function obtain_group_cov_stats & obtain_overal_stats in: %s" % (
        after - prior) + "\n"))

    # Remove Uncovered

    prior = datetime.datetime.now()

    logger.info(GREEN + "Removing low quality samples in group " +
                group_name + END_FORMATTING)

    uncovered_samples = remove_low_quality(
        output_dir, output_dir, cov20=args.coverage20, unmapped_per=args.unmapped, min_hq_snp=args.min_snp, type_remove="Uncovered")

    if len(uncovered_samples) > 1:
        logger.info(RED + BOLD + "Uncovered samples: " +
                    (",").join(uncovered_samples) + END_FORMATTING)
    else:
        logger.info(GREEN + "NO uncovered samples found" + END_FORMATTING)

    after = datetime.datetime.now()
    print(("Done with function remove_low_quality in: %s" % (after - prior) + "\n"))

    ##### ANNOTATION #####

    logger.info("\n\n" + BLUE + BOLD + "STARTING ANNOTATION IN GROUP: " +
                group_name + END_FORMATTING + "\n")

    # Annotation with SnpEFF

    prior = datetime.datetime.now()

    if args.snpeff_database != False:
        for root, _, files in os.walk(out_variant_dir):
            for name in files:
                # print(name)
                if name == "snps.all.vcf":
                    sample = root.split("/")[-1]
                    # print(sample)
                    filename = os.path.join(root, name)
                    # print(filename)
                    chrom_filename = os.path.join(
                        root, "snps.all.chromosome.vcf")
                    out_annot_file = os.path.join(
                        out_annot_snpeff_dir, sample + ".annot")
                    # print(out_annot_file)

                    if os.path.isfile(out_annot_file):
                        logger.info(YELLOW + DIM + out_annot_file +
                                    " EXIST\nOmmiting SnpEFF Annotation for sample " + sample + END_FORMATTING)
                    else:
                        logger.info(
                            GREEN + "Annotating sample with SnpEFF: " + sample + END_FORMATTING)
                        rename_reference_snpeff(filename, chrom_filename)
                        annotate_snpeff(
                            chrom_filename, out_annot_file, database=args.snpeff_database)

    else:
        logger.info(YELLOW + BOLD + "No SnpEFF database suplied, skipping annotation in group " +
                    group_name + END_FORMATTING)

    after = datetime.datetime.now()
    print(("Done with function rename_reference_snpeff & annotate_snpeff in: %s" % (
        after - prior) + "\n"))

    # Annotation for user defined (bed & vcf annot)

    prior = datetime.datetime.now()

    if not args.annot_bed and not args.annot_vcf:
        logger.info(
            YELLOW + BOLD + "Ommiting User Annotation, no BED or VCF files supplied" + END_FORMATTING)
    else:
        for root, _, files in os.walk(out_variant_dir):
            for name in files:
                # print(name)
                if name == "snps.all.ivar.tsv":
                    sample = root.split("/")[-1]
                    # print(sample)
                    logger.info(
                        "User bed/vcf annotation in sample {}".format(sample))
                    filename = os.path.join(root, name)
                    # print(filename)
                    out_annot_file = os.path.join(
                        out_annot_user_dir, sample + ".tsv")
                    user_annotation(
                        filename, out_annot_file, vcf_files=args.annot_vcf, bed_files=args.annot_bed)

    after = datetime.datetime.now()
    print(("Done with function user_annotation in: %s" % (after - prior) + "\n"))

    # Annotation for user aa defined (aminoacid annot)

    prior = datetime.datetime.now()

    if not args.annot_aa:
        logger.info(
            YELLOW + BOLD + "Ommiting User aa annotation, no AA files supplied" + END_FORMATTING)
    else:
        for root, _, files in os.walk(out_annot_snpeff_dir):
            if root == out_annot_snpeff_dir:
                for name in files:
                    if name.endswith(".annot"):
                        sample = name.split(".")[0]
                        logger.info(
                            "User aa annotation in sample {}".format(sample))
                        filename = os.path.join(root, name)
                        out_annot_aa_file = os.path.join(
                            out_annot_user_aa_dir, sample + ".tsv")

                        if os.path.isfile(out_annot_aa_file):
                            user_annotation_aa(
                                out_annot_aa_file, out_annot_aa_file, aa_files=args.annot_aa)
                        else:
                            user_annotation_aa(
                                filename, out_annot_aa_file, aa_files=args.annot_aa)

    after = datetime.datetime.now()
    print(("Done with function user_annotation_aa in: %s" % (after - prior) + "\n"))

    # Annotation with blast (fasta annot)

    prior = datetime.datetime.now()

    if not args.annot_fasta:
        logger.info(
            YELLOW + BOLD + "Ommiting User FASTA annotation, no FASTA files supplied" + END_FORMATTING)
    else:
        for root, _, files in os.walk(out_consensus_dir):
            for name in files:
                if name.endswith(".fa"):
                    filename = os.path.join(root, name)
                    sample = root.split("/")[-1]
                    logger.info(
                        "User FASTA annotation in sample {}".format(sample))

                    for db in args.annot_fasta:
                        make_blast(filename, db, sample, out_annot_blast_dir, db_type="nucl",
                                   query_type="nucl", evalue=0.0001, threads=args.threads)

    after = datetime.datetime.now()
    print(("Done with function make_blast in: %s" % (after - prior) + "\n"))

    ##### COMPARISON #####

    # SNPs comparison using tsv variant files

    logger.info("\n\n" + GREEN + BOLD + "STARTING COMPARISON IN GROUP: " +
                group_name + END_FORMATTING + "\n")

    folder_compare = today + "_" + group_name
    path_compare = os.path.join(out_compare_dir, folder_compare)
    check_create_dir(path_compare)
    full_path_compare = os.path.join(path_compare, group_name)

    compare_snp_matrix_recal = full_path_compare + ".revised.final.tsv"
    compare_snp_matrix_recal_intermediate = (
        full_path_compare + ".revised_intermediate.tsv")
    compare_snp_matrix_recal_mpileup = (
        full_path_compare + ".revised_intermediate_vcf.tsv")
    compare_snp_matrix_INDEL_intermediate = (
        full_path_compare + ".revised_INDEL_intermediate.tsv")
    compare_only_snps = full_path_compare + "_ONLY_SNPs.revised.tsv"

    # Create intermediate

    prior = datetime.datetime.now()

    recalibrated_snp_matrix_intermediate = ddbb_create_intermediate(
        out_variant_dir, out_stats_coverage_dir, min_freq_discard=args.min_allele_frequency, min_alt_dp=args.min_depth, only_snp=False)
    # recalibrated_snp_matrix_intermediate.to_csv(compare_snp_matrix_recal_intermediate, sep="\t", index=False)

    after = datetime.datetime.now()
    print(("Done with function ddbb_create_intermediate in: %s" %
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
        compare_snp_matrix_recal_mpileup, sep="\t", index=False)

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

    # Extract all low coverage o not covered positions

    prior = datetime.datetime.now()

    symbol_file = full_path_compare + '_symbol_lowcov.tsv'
    symbol_lowcov = extract_lowcov(compare_snp_matrix_INDEL_intermediate) # It is made by the INDEL_intermediate.tsv, taking 0 and 1 into account, can also be made with intermediate.highfreq.tsv
    symbol_lowcov.to_csv(symbol_file, sep='\t', index=False)

    after = datetime.datetime.now()
    print(("Done with function extract_lowcov in: %s" %
               (after - prior) + "\n"))

    # Extract all positions marked as complex

    prior = datetime.datetime.now()

    complex_variants = extract_complex_list(out_variant_dir)
    logger.debug('Complex positions in all samples:\n{}'.format(
        (','.join([str(x) for x in complex_variants]))))

    after = datetime.datetime.now()
    print(("Done with function extract_complex_list in: %s" %
           (after - prior) + "\n"))

    # Clean all faulty positions and samples for final table

    prior = datetime.datetime.now()

    recalibrated_revised_INDEL_df = revised_df(compare_snp_matrix_INDEL_intermediate_df, path_compare, complex_pos=complex_variants, min_freq_include=0.7, min_threshold_discard_uncov_sample=args.min_threshold_discard_uncov_sample, min_threshold_discard_uncov_pos=args.min_threshold_discard_uncov_pos, min_threshold_discard_htz_sample=args.min_threshold_discard_htz_sample,
                                               min_threshold_discard_htz_pos=args.min_threshold_discard_htz_pos, min_threshold_discard_all_pos=args.min_threshold_discard_all_pos, min_threshold_discard_all_sample=args.min_threshold_discard_all_sample, remove_faulty=True, drop_samples=True, drop_positions=True, windows_size_discard=args.window)

    recalibrated_revised_df = recalibrated_revised_INDEL_df[~recalibrated_revised_INDEL_df['Position'].isin(symbol_lowcov['Position'])]

    recalibrated_revised_df.to_csv(
            compare_snp_matrix_recal, sep='\t', index=False)

    if args.only_snp:
        compare_only_snps_df = extract_only_snps(
            compare_snp_matrix_recal)
        compare_only_snps_df.to_csv(compare_only_snps, sep="\t", index=False)

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
