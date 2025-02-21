# Rapid-cloning-of-Yr6

This repository shows a modified pipeline inspired by [MutChromSeq](https://github.com/steuernb/MutChromSeq/) which was used for rapid cloning of strip rust resistant locus *Yr6* based on RNA-seq data of EMS-mutagenized mutants.

## Data preperation

This pipeline requires RNA-seq data of multiple independent susceptiable mutants. Any commonly used aligner is acceptable. We opted to employ [STAR](https://alexdobin.github.io/STAR/) to map the reads of mutants to CDC Stanley wheat genome assembly (https://www.ebi.ac.uk/ena/browser/view/GCA_902810665). In order to exclude nonspecific alignments as much as possible, while also mitigating the impact of SNPs present in the reference sequence, we used the parameter `--outFilterMismatchNmax 2`, which removes alignments with more than 2 substitutions. The potential PCR-duplicated reads were removed using [sambamba](https://github.com/biod/sambamba).  

## Variant calling

Variant calling was perform with [minipileup](https://github.com/lh3/minipileup) and the result was sorted and compressed by [bcftools](https://samtools.github.io/bcftools/).

```bash
minipileup -f stanley.fa -vc <path/to/samples.bam> | \
		bcftools sort -Oz -o SNP.sort.vcf.gz
	bcftools index SNP.sort.vcf.gz
```
**Note :** We did not utilize variant calling tools such as GATK, in order to minimize the introduction of false negatives due to hard filtering. The false positives resulting from this can be overcome through mutual validation among independent mutants.

## Vcf file chopping

Inheriting the core concept of methods such as MutChromSeq, when mutations from multiple independent mutants converge on the same gene, transcript, or genomic segment, it indicates an association/cosegregation of the gene with the mutant phenotype (i.e. stripe rust susceptibility). However, the probability of non-associated mutation sites enriching on the same ginomic segment increases with the length of the segment. We developed the python script, `vcfChopperMultipler.v2.py`, to "chop" the chromosomes in the vcf file into overlapping segments for subsequent analysis.

**Note:** `pysam` is required by this script.

```bash
python vcfChopperMultipler.v2.py \
	-i SNP.sort.vcf.gz -j 40 \
    > SNP.chop.vcf
```

Here's the parameters of the script.

| Parameter | Description |
| --- | --- |
| -h | Print help information |
| -w | Work path, default is "./" |
| -i | Input vcf file, must be bgzipped and with a suffix of "vcf.gz" |
| -l | Length of interval, default is 50000. |
| -p | Length of overlap, default is 5000. |
| -t | Offset of the first interval, default is 0. |
| -f | File of a list of chromosomes to process parallelly. Choromesomes not mentioned will be deprecated. |
| -j | Number of blocks to be processed parallely. Default is 1 (with out parallely processing). |

## MutCandidator

We developed the python script, `MutCandidator.v3.py`, to filter and compare the SNP sites among mutants and finally identify the ginomic segment carrying the candidate resistant gene.

```bash
python MutCandidator.v3.py -i SNP.chop.vcf -n 10 -a 0.3 -s
```

Here's the parameters of the script.

| Parameter | Description |
| --- | --- |
| -h | Print help information |
| -w | Work path, default is "./" |
| -i | Input vcf file, usually the chopped vcf file. |
| -n | Minimum number of mutants to report a contig. Default is 2. |
| -c | Minimum depth for a site to be regarded. Default is 10. |
| -a | Maximum reference allele frequency of a SNP to be reported. Default is 0.01. |
| -z | Number of mutants that are allowed to have SNP at a same site. Default is 2. |
| -m | Maximum SNP number per contig per mutant able to be retained after filter. Default is 0 (unlimited). |
| -s | Strict mode, with only C->T or G->A type permitted. |

The script will first filter the SNPs based on parameters `c`, `a`, `z`, and `m`, and then outputs a filtered vcf file named *SNP.chop.vcf.a0.3c10z2.filter.vcf* in which the excluded sites are denoted as `.`. This process may take tens of minutes. Note that the filename represents the parameters given. If this filtered vcf file is absent (e.g. you changed parameter `a`), the script will regenerate it; otherwise, the script will read the existing file and identify the candidate segments based on the parameter `n`, which will be very fast. The candidates will be reported in a file named *SNP.chop.vcf.n10a0.3c10z2.s.report.txt*

In our study, the reports says:

```
Candidate contig: LR865758.1_707445000
10 mutants
	TW4-4.dedup.bam:
		Pos: 12639	G->A	depth: 15
		Pos: 38707	G->A	depth: 705

	TW4-2.dedup.bam:
		Pos: 23997	C->T	depth: 208
		Pos: 38522	C->T	depth: 278

	TW4-14A.dedup.bam:
		Pos: 24503	G->A/G	depth: 117/1

	TW4-8.dedup.bam:
		Pos: 24667	C->T	depth: 112

	TW4-3.dedup.bam:
		Pos: 25157	C->T	depth: 193

	TW4-1.dedup.bam:
		Pos: 40206	C->T	depth: 324

	TW4-5.dedup.bam:
		Pos: 40914	C->T	depth: 685

	TW4-10A.dedup.bam:
		Pos: 40993	G->A/G	depth: 168/1

	TW4-19A.dedup.bam:
		Pos: 41093	C->T	depth: 353

	T2W4-6-LGM26650_L1.dedup.bam:
		Pos: 41946	G->A	depth: 373
```

The reports tells you the candidate segment and the information of SNP sites. The name of the segment `LR865758.1_707445000` indicates that it is a segment starting from the 707,445,001st base on chromosome LR865758.1 (chr7B), and the positions of each SNP locus are calculated relative to the corresponding segment. That is, the position of a SNP site on the chromosome should be $\text{Pos} + 707445000$.

