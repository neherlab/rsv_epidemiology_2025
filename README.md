# rsv_epidemiology_2025
A collection of code and data used for **"Molecular Epidemiology of Respiratory Syncytial Virus in Switzerland from 2020 to 2024"** publication.

---

## `src/genome_assembly`
This directory contains input data and metadata used for genome assembly and downstream analysis of RSV samples.

### `data/` Folder Structure

```
data/
├── annotation.tsv
├── primers/
├── ref_genomes/
└── samples/
```

---

#### `annotation.tsv`

Sample metadata file.  
**Columns:**
- `NGS_ID`: Full sample ID (matches FASTQ filenames)
- `Virus Type`: Either `HRSVA` or `HRSVB`
- `amplicon_types`: Comma-separated amplicon types (`500bp, 1000bp, 2000bp` in the original analysis)
- `Pools`: `Pool 1, Pool 2` in the original analysis, different for separate even/odd amplicon pool experiments

---

#### `primers/`

CSV files specifying primer schemes used for different virus types and amplicon sizes.

Naming: `rsv_{a|b}_primers_{500bp|1000bp|2000bp|revised}{|_even|_odd}.csv`

- RSV A or B;
- Amplicon sizes; revised for pooled together from the original analysis;
- empty for combined odd & even amplicons, different otherwise.

---

#### `ref_genomes/`

Reference genome FASTA files:
- `rsva_ref.fasta`
- `rsvb_ref.fasta`

*(Additional BWA index files may be created here, can be cleaned with `cleanup`)*

---

#### `samples/`

One folder per sample (named by `NGS_ID` in `annotation.tsv`). Each contains:

```
<NGS_ID>/
├── <NGS_ID>_R1.fastq.gz
└── <NGS_ID>_R2.fastq.gz
```

---

### Configuration (`config.yaml`)

The `config.yaml` file specifies various parameters for running the pipeline. Below is a breakdown of the options in the current configuration:

| Parameter         | Description                                                                 | Example                   |
|-------------------|-----------------------------------------------------------------------------|---------------------------|
| `base_dir`        | The base directory where the project is located.                            | `"."` (current directory)  |
| `res_dir_name`    | The name of the directory where results will be stored.                     | `"results"`                |
| `ann_filename`    | The name of the annotation file in `src/genome_assembly/data`. This file containis metadata for the samples. | `"annotation.tsv"`         |
| `read_length`     | TrimGalore min length.                     | `80`                       |
| `min_cov`         | The minimum coverage threshold for variants.                               | `100`                      |
| `freq_threshold`  | The frequency threshold: bases with the main allele of lower values are treated as undefined                     | `0.9`                      |
| `allow_cleanup`   | A flag that determines if intermediate files should be cleaned up after execution. Set to `True` to enable cleanup. | `True`                     |

### Example Commands

Here are some example commands run from `src/genome_assembly`:

#### 1. **Run the genome assembly pipeline**
```bash
snakemake --configfile config.yaml --cores 1 --use-conda -p
```

#### 2. **Cleanup intermediate files**
```bash
snakemake cleanup --cores 1
```

---
