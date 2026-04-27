# Prokaion: Whole-genome sequencing analysis pipeline for prokaryotes (specialized in _Mycobacterium tuberculosis_)

`Prokaion` is a modular, in-house pipeline designed for the analysis of nanopore WGS data from any prokaryotic microorganism (specialized in _Mycobacterium tuberculosis_). It performs basecalling and barcoding with guppy or Dorado, quality control, variant calling, taxonomic assignment, variant annotation and comparative genomics in a reproducible and structured way.

---

## 📥 Clone the Repository

To get started, clone the `prokaion` repository from GitHub:

```bash
git clone https://github.com/MG-IiSGM/prokaION
cd prokaion
```

---

## 🔧 Installation

This pipeline is built to run in a Conda environment. Use the provided `prokaion.yml` file to install all dependencies:

```bash
conda env create -f prokaion.yml
conda activate prokaion
```

> Channels used: `bioconda`, `conda-forge`, `anaconda`, `r`, `defaults`, `etetoolkit`.

---

## 🚀 Usage

Run the pipeline from the command line using the main script:

```bash
python prokaion.py -i input_folder -o output_folder --basecall dorado --model last_model --barcode_kit kit_used -s Summary_file -S -g -rq 20 -T 100 -r reference.fa --snpeff_database Mycobacterium_tuberculosis_h37rv --mash_database bacteria_mash.msh --kraken2 Kraken2_database/ -R repeats_phage_coord.bed -V lin_tbdb_MTBanc.vcf -V Resistances_OMS.vcf -A Resistance_OMS.aa
```

### 🔹 Parameters

| Parameter           | Description                                                                                                                    |
| ------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `-i`                | Input folder with paired-end FASTQ files (`*_1.fastq.gz`, `*_2.fastq.gz`)                                                      |
| `-o`                | Output folder                                                                                                                  |
| `--basecall`        | Basecalled and barcoded program (`Guppy for R9`, `Dorado for R10`)                                                             |
| `--model`           | Model used, dorado download --models-directory to download all models and used the last version of dna_r10.4.1_e8.2_400bps_sup |
| `--barcode_kit`     | Kit of barcodes used [SQK-RBK110-96 / EXP-NBD196 / SQK-NBD112-24]                                                              |
| `-s`                | Sample list for conversion from barcode to samples ID                                                                          |
| `-g`                | Use GPU device                                                                                                                 |
| `-rq`               | Filter on a minimum average read quality score.                                                                                |
| `-T`                | Number of threads (optional)                                                                                                   |
| `-r`                | Reference genome in FASTA format                                                                                               |
| `--snpeff_database` | SnpEff database name (optional, for variant annotation)                                                                        |
| `--mash_database`   | Mash database for taxonomic classification (optional)                                                                          |
| `--kraken2`         | Kraken2 database path (optional)                                                                                               |
| `-R`                | BED file with coordinates to remove (e.g. repeats, phages) (optional)                                                          |
| `-V`                | VCF file(s) with custom variant annotations (lineage/resistance markers)                                                       |
| `-A`                | File with codon-level resistance annotations (requires SnpEff output)                                                          |

Additional parameters can be applied. These can be checked in the script help `prokaion.py --help`

---

### 📄 Example of Optional Input Files

The following optional inputs can be provided to enrich the variant annotation and filtering process:

#### 🔹 `-R` or `-B` : BED file with regions to exclude (`-R`) (e.g., repeats or phage insertions) or annotate any position within a region (`-B`)

This file defines genomic coordinates (e.g. in BED format) that will be excluded from downstream analysis, such as repetitive or poorly mappable regions, or annotated in the _M. tuberculosis_ genome.

**Example:**

| Chromosome | Start  | End    | INFO |
| ---------- | ------ | ------ | ---- |
| MTB_anc    | 33582  | 33794  | +    |
| MTB_anc    | 80185  | 80373  | +    |
| MTB_anc    | 103710 | 104663 | -    |

---

#### 🔹 `-V` : VCF file with user-defined variant annotations (e.g., lineage or resistance markers)

This file should contain variants of interest to annotate VCF outputs during the analysis, such as lineage-defining mutations.

**Example:**

| CHROM      | POS  | ID  | REF | ALT | QUAL | FILTER | FORMAT | INFO           |
| ---------- | ---- | --- | --- | --- | ---- | ------ | ------ | -------------- |
| Chromosome | 1131 | .   | C   | A   | 20   | PASS   | .      | lineage4.2.2.1 |
| Chromosome | 4206 | .   | C   | T   | 20   | PASS   | .      | lineage4.1.3   |
| Chromosome | 9944 | .   | A   | C   | 20   | PASS   | .      | lineage4.2.2   |

---

#### 🔹 `-A` : Codon-level resistance annotations

This file allows you to specify amino acid substitutions linked to resistance phenotypes. These will be matched against variant calls (e.g. via SnpEff output).

**Example:**

| Substitution | Gene:Drug       |
| ------------ | --------------- |
| Asp435Phe    | rpoB:Rifampicin |
| Asp435Tyr    | rpoB:Rifampicin |
| Asp435Val    | rpoB:Rifampicin |

> These inputs are optional, but if `-A` is specified, SnpEff annotation must be run (`--snpeff_database`) to enable codon matching.

---

## 📂 Output Structure

The pipeline produces a structured directory with the following subfolders:

- `Quality/` – FastQC reports for raw data.
- `Stats/`
  - `Coverage/` – per-sample `.cov` files for genome-wide depth.
  - `Bamstats/` – alignment stats in `.bamstats` format.
  - `coverage.summary.tab` – summary of coverage thresholds.
  - `overall.stats.tab` – includes high-quality SNPs, indels, heterozygosity and high-quality sequences.
- `Variants/` – one folder per sample with BAM, BAI, raw and filtered VCFs.
- `Species/` – taxonomic assignment results (if `--mash_database` or `--kraken2` used).
- `Uncovered/` – samples failing the 70% >20X coverage threshold (configurable).
- `Annotation/`
  - `snpeff/` – SnpEff-annotated VCFs.
  - `user/` – VCFs annotated with custom variant lists from `-V`.
  - `user_aa/` – codon-level resistance annotations from `-A`.
- `Compare/` – inter-sample SNV comparisons, phylogeny (`.nwk`), dendrograms and intermediate data.

> For ambiguous or low-confidence regions (e.g. heterozygous SNPs), manual review with a genome browser is recommended.

---

## ✅ Quality Control

- **FastQC** is used to assess raw data quality.
- **Filtering and trimming** are done via Snippy and custom wrappers.
- **Coverage and alignment statistics** are generated with `samtools` and integrated into summary reports.

---

## 📦 Reproducibility and Portability

We strongly recommend using containers for reproducibility. If you are using a shared cluster or HPC environment, you can build a Docker or Singularity image based on the `prokaion.yml` environment file.

---

## 🧪 Example (Minimal)

```bash
python prokaion.py -i . -o . -T 16 -r H37Rv.fasta
```

---

## 📬 Contact

For questions or contributions, please contact the authors or open an issue in the GitHub repository.

**Author:** [Sergio Buenestado-Serrano](mailto:sergio.buenestado@gmail.com)  
**Email:** <sergio.buenestado@gmail.com>

---

## 📄 License

This software is distributed under the **GNU GENERAL PUBLIC LICENSE**.

Although the project is distributed under an open license, **collaboration is highly welcome**. Contributions, feedback, and suggestions to improve Prokaion are encouraged.
