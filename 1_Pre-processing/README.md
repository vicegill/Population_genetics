# Pre-processing pipeline

We pre-process our samples using the snakemake-based [greenpipe](https://github.com/moiexpositoalonsolab/grenepipe/wiki). 

To run the whole pre-processing pipeline just type:

`
conda activate grenepipe
snakemake --conda-frontend mamba --directory /dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2022_SwallowWGS/Results/1_Pre-processing_v3/ --profile profiles/slurm --conda-prefix ~/conda-envs`

Here we describe in detail the parmeters used for each step of the pre-processing pipeline.

### Pipeline overview

0. [Set-up](#0-set-up)
1. [Trimming](#1-Trimming)
2. [Mapping](#2-Mapping)
3. [Dedup](#3-Dedup)
4. [Indel realignment](#4-Indel-realignment)
5. [Calling](#5-Calling)
6. [Filtering](#6-Filtering)

---

## 0. Set-up

For the pipeline to run we install miniconda3 locally and we install and activate the greepipe environment following recommendations listed [here](https://github.com/moiexpositoalonsolab/grenepipe/wiki/Setup-and-Usage). 

Subsequently, we modify the config.yaml file that contains all configurations for the data and the tools to use. The modify version of this file used in this project can be found in [here](https://github.com/pnunezv/shallow_WGS_project/tree/main/1_Pre-processing/grenepipe-0.10.0/config.yaml) and on the results folder on the lrz server: `/dss/dsslegfs01/pr53da/pr53da-dss-0034/projects/2022_SwallowWGS/1_Pre-processing_v3/config.yaml`

### Input data

#### Sample information

The configuration (config.yaml) expects data: samples to point to a tab-separated table that lists all sample names and the absolute paths to their fastq files, our is stored on: `/dss/dsshome1/04/ra96her/shallow_project/Data/sample_data/sample_files.tsv` which looks like this: 

```
sample          unit    platform        fq1     fq2
18CH024 1       ILLUMINA        /dss/dsslegfs01/pr53da/pr53da-dss-0034/rawdata/WGS/Illumina/CH_18_R_18CH024_BL_ADL_M__SRR21484715__R1__14f60d1f84bc56cbc59242044bff626d.fq.gz     /dss/dsslegfs01/pr53da/pr53da-dss-0034/rawdata/WGS/Illumina/CH_18_R_18CH024_BL_ADL_M__SRR21484715__R2__c24f10115b1830c812bf872ff5b72a10.fq.gz
18CH024 2       ILLUMINA        /dss/dsslegfs01/pr53da/pr53da-dss-0034/rawdata/WGS/Illumina/CH_18_R_18CH024_BL_ADL_M__SRR21484714__R1__f164da55fdd6293f0eb60cdb930d5667.fq.gz     /dss/dsslegfs01/pr53da/pr53da-dss-0034/rawdata/WGS/Illumina/CH_18_R_18CH024_BL_ADL_M__SRR21484714__R2__18f0be171285e6efcc5460c6b604afa7.fq.gz
18CH081 1       ILLUMINA        /dss/dsslegfs01/pr53da/pr53da-dss-0034/rawdata/WGS/Illumina/CH_18_R_18CH081_BL_ADL_M__SRR21484703__R1__cf9b0fd7879be342be3e8d5acdf59bc8.fq.gz     /dss/dsslegfs01/pr53da/pr53da-dss-0034/rawdata/WGS/Illumina/CH_18_R_18CH081_BL_ADL_M__SRR21484703__R2__f5e0413bde0d1770835e74c62d9d60d0.fq.gz
```

#### Reference genome

On the other hand, the config file also expects a reference genome in our case we are using 
*Hirundo rustica* (Barn swallow) reference genome stored in: `/dss/dsslegfs01/pr53da/pr53da-dss-0034/assemblies/Hirundo.rustica/genome/v2/Hirundo.rustica_genome_bHirRus1_PacBio.Illumina.Hi-C.DLS_v2.fasta`. This reference genome is available to download in [NCBI](https://www.ncbi.nlm.nih.gov/data-hub/genome/GCF_015227805.1/).

##### Assembly stats

|                       |            |
| --------------------- | ---------- |
| Genome size           | 1.1 Gb     |
| Number of chromosomes | 39         |
| Number of scaffolds   | 617        |
| Scaffold N50          | 76.2 Mb    |
| Scaffold L50          | 5          |
| Number of contigs     | 1,719      |
| Contig N50            | 2.8 Mb     |
| Contig L50            | 100        |
| GC percent            | 42.5       |
| Assembly level        | Chromosome |

### Cluster set-up

Since our pipeline is running on a cluster environment we adapt a slurm profile present in the greenpipe docs. In particular we modify the [cluster_config.yaml](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/profiles/slurm/cluster_config.yaml). This file contains per-rule customization, for example, execution times, memory limits, partitions and user names to use for the submission, etc. 

**NOTE:** We modify [slumr-status.py](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/profiles/slurm/slumr-status.py) to adpat it to the bioHPC cluster requiremets. In essence, we specify cluster and user name for `job status check -M biohpc_gen -u ra96her` 

## 1. Trimming

Removal of adapter sequences in a process called read trimming. To perform this removal we used [Cutadapt software](https://cutadapt.readthedocs.io/en/stable/). Rule used to run this step can be found [here](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/rules/trimming-cutadapt.smk)

We used the following parameters: 

| Parameter | Description                                                                                                                                        |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------- |
| -a        | Regular 3’ adapter - Foward                                                                                                                        |
| -g        | Regular 5’ adapter - Foward                                                                                                                        |
| -A        | Regular 3’ adapter - Reverse                                                                                                                       |
| -G        | Regular 5’ adapter -  Reverse                                                                                                                      |
| -q  10    | Assume that quality values in the FASTQ file are encoded as ascii(quality + N). This needs to be set to 64 for some very old Illumina FASTQ files. |

> Cutadapt removes adapter sequences from high-throughput sequencing reads.
> Martin M.
> EMBnet journal. 2011.
> doi:10.14806/ej.17.1.200

## 2. Mapping

Read mapping is the process to align the reads on a reference genomes. To map the reads we used [bwa mem](https://bio-bwa.sourceforge.net/bwa.shtml). Rule used to run this step can be found [here](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/rules/mapping-bwa-mem.smk). Moreover, this rule use a snakemake template wrapper thar can be found [here](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/rules/mapping-bwa-mem.smk)

We used the following parameters: 

| Parameter    | Description                                                      |
| ------------ | ---------------------------------------------------------------- |
| -M           | Mark shorter split hits as secondary (for Picard compatibility). |
| -R "'@RG\\t" | Complete read group header line                                  |

> Fast and accurate short read alignment with Burrows-Wheeler transform.
> Li H, Durbin R.
> Bioinformatics. 2009.
> [doi:10.1093/bioinformatics/btp324](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=false)

Aditionally, after mapping we sort our mapped reads using [samtools](http://www.htslib.org/doc/samtools-sort.html). The following parameters are used:

| Parameter | Description                                                                                            |
| --------- | ------------------------------------------------------------------------------------------------------ |
| -m 4G     | Approximately the maximum required memory per thread                                                   |
| -l 9      | Set the desired compression level for the final output file. 9 (best compression but slowest to write) |
| -@ 20     | Set number of sorting and compression threads                                                          |

> The Sequence Alignment/Map format and SAMtools.
> Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, et al.
> Bioinformatics. 2009.
> [doi:10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688?login=false)

## 3. Dedup

We remove PCR duplicates using [picard tools](https://gatk.broadinstitute.org/hc/en-us/articles/360057439771-MarkDuplicates-Picard). This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Rule used to run this step can be found [here](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/rules/duplicates-picard.smk). Moreover, this rule use a snakemake template wrapper thar can be found [here](https://github.com/snakemake/snakemake-wrappers/blob/master/bio/picard/markduplicates/wrapper.py)

We used the following parameters: 

| Parameter                               | Description                                                                                                               |
| --------------------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| REMOVE_DUPLICATES=true                  | If true do not write duplicates to the output file instead of writing them with appropriate flags set.                    |
| MAX_RECORDS_IN_RAM=150000               | When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. |
| MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 | Maximum number of file handles to keep open when spilling read ends to disk.                                              |

> Picard toolkit.
> Broad Institute; 2018.
> GitHub repository, online: [http://broadinstitute.github.io/picard/](https://broadinstitute.github.io/picard/)

## 4. Indel realingnment

Reads near detected indels are realigned to remove alignment artifacts. This step is apparently integrate within gatk4 processing steps, there is no need to add it manually, [check](https://github.com/moiexpositoalonsolab/grenepipe/issues/2).

> The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data.
> McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, et al.
> Genome Res. 2010.
> [doi:10.1101/GR.107524.110](https://genome.cshlp.org/content/20/9/1297)

## 5. Calling

Reads near detected indels are realigned to remove alignment artifacts. To perform this removal we used [GATK](https://gatk.broadinstitute.org/hc/en-us). Rule used to run this step can be found [here](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/rules/calling-haplotypecaller.smk)
Moreover, this rule use a snakemake template wrapper thar can be found [here](https://github.com/snakemake/snakemake-wrappers/blob/master/bio/gatk/haplotypecaller/wrapper.py)

We used the following parameters: 

| Parameter            | Description                                                           |
| -------------------- | --------------------------------------------------------------------- |
| --intervals          | One or more genomic intervals over which to operate                   |
| -mbq 15              | Minimum base quality required to consider a base for calling          |
| -ERC GVCF            | Mode for emitting reference confidence scores. On wrapper             |
| -hets 0.0037         | Heterozygosity value used to compute prior likelihoods for any locus. |
| -indelHeterozygosity | Heterozygosity for indel calling. There is no infoo                   |

*NOTE*: To obtain the heterozygosity values we ran the following comands using awk:

```
vcf=/dss/dsslegfs01/pr53da/pr53da-dss-0034/vcf/hirundo_rustica+smithii.allsites.final.auto.snps.miss02.maf05.ingroup.vcf.gz

#For SNPs
awk -F\0\/0: '!/^ *#/ {total += NF-1; count++} END { print total/count }' 
-0.00374433 #This makes no sense, so we:
```

*BE CAREFUL HERE! We obtain this estimates only for autosomes!!! We may need to run a separate run for chr W & Z.*
*ALSO THERE IS NO INDEL INFO *

The expected heterozygosity value used to compute prior probability that a locus is non-reference. The default priors are for provided for humans: het = 1e-3 which means that the probability of N samples being hom-ref at a site is: 1 - sum_i_2N (het / i). That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population of organisms would differ from each other (one being A and the other B) at a rate of 1 in 100 bp. 

*NOTE2*: Intervals are given as an input. But I don't get from where they are obtaining them. Also, if we were to calculate the indel realignement this we intervals would be obtain from there. 

The next step in the pipeline is Combine the GVCF, this will be done using the same rule as the haplotype caller. But the wrapper used is [this](https://github.com/snakemake/snakemake-wrappers/blob/master/bio/gatk/combinegvcfs/wrapper.py). Fot this step we decided not to add extra parameters. 

Afterwards, we proceed to genotype the vcf. For this step the same rule as before is used, and the wrapper can be found [here](https://github.com/snakemake/snakemake-wrappers/tree/master/bio/gatk/genotypegvcfs) 

| Parameter                 | Description                                                                      |
| ------------------------- | -------------------------------------------------------------------------------- |
| -hets 0.015               | Heterozygosity value used to compute prior likelihoods for any locus.            |
| -indelHeterozygosity 0.01 | Heterozygosity for indel calling. .                                              |
| -allSites                 | Genotype all sites                                                               |
| -stand_call_conf 15       | The minimum phred-scaled confidence threshold at which variants should be called |

> The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data.
> McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, et al.
> Genome Res. 2010.
> [doi:10.1101/GR.107524.110](https://genome.cshlp.org/content/20/9/1297)

> The Variant Call Format and VCFtools.
> Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, et al. 
> Bioinformatics. 2011.
> [doi:10.1093/bioinformatics/btr330](https://doi.org/10.1093/bioinformatics/btr330)

## 6. Filtering

As a last step of this pipeline we need to filter our vcf files. To perform this filtering we used gatk. The rule used to run this step can be found [here](https://github.com/pnunezv/shallow_WGS_project/blob/main/1_Pre-processing/grenepipe-0.10.0/rules/filtering.smk). 

The first step in the filtering pipeline is to selct variants. the wrapper for this can be found [here](https://github.com/snakemake/snakemake-wrappers/blob/master/bio/gatk/selectvariants/wrapper.py). No extra parameters were added. Next it comes the actual filtering, wrapper used in this step can be found [here](https://github.com/snakemake/snakemake-wrappers/blob/master/bio/gatk/variantfiltration/wrapper.py). For this we chose to do a hard filtering following [gatk recommendations](https://gatkforums.broadinstitute.org/gatk/discussion/6925/understanding-and-adapting-the-generic-hard-filtering-recommendations) .  We used the next filter expression: 

SNPs

```
QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0
```

Indels

```
QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0
```

Where:

| Filter                | Description                                                                                                                                                                                 |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| QD < 2.0              | quality by depth (QUAL/DP),                                                                                                                                                                 |
| FS > 60               | fisher strand. Phred-scaled probability that there is strand bias at the site. gatk recommended hard filter.                                                                                |
| MQ < 40               | root mean square mapping quality over all the reads at the site.                                                                                                                            |
| MQRankSum < -12.5     | MappingQualityRankSumTest. A value close to zero is best and indicates little difference between the mapping qualities.                                                                     |
| ReadPosRankSum < -8.0 | the u-based z-approximation from the Rank Sum Test for site position within reads. It compares whether the positions of the reference and alternate alleles are different within the reads. |
| SOR > 3.0             | StrandOddsRatio. Estimate strand bias. Most variants have an SOR value less than 3.                                                                                                         |
| QUAL < 25.0           | quality                                                                                                                                                                                     |

*NOTE*: We did not add a mask setting. 

> The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data.
> McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, et al.
> Genome Res. 2010.
> [doi:10.1101/GR.107524.110](https://genome.cshlp.org/content/20/9/1297)

The resulting VCF file from the grenepipe pipeline was further filtered using vcftools (v0.1.16) (Danecek et al., 2011). Indels were removed (-remove-indels) and only biallelic SNPs with a minimum quality of 20 and depth from 5 to 35 were kept (--min-alleles 2 --max-alleles 2 -minQ 20  --min-meanDP 5 --max-meanDP 35). Moreover, we imposed a threshold of 0.05 for minor allele frequency and removed  95% of missing data (--max-missing 0.95 --maf 0.05). Filtered VCF yields a total of 17,248,742 sites. 

Script: ``` 1_filter_vcftools.sh ```

Additionally, we check the quality of the data with vcftools stats, using ```0_checkvcf_vcftools.sh```
We plot this stats using the Rmarkdown script: ```rustica.allsites/rustica.allsites.Rmd```