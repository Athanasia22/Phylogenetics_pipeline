# Phylogenetics_pipeline
My data is from zooplankton!
Data analysis is conducted for genes: COI, 18S, ITS


1. [Project Overview](#project-overview)
2. [Fetch Sequences](#fetch-sequences)
3. [Multiple Sequence Alignment](#multiple-sequence-alignment)
4. [Preparation for Quality Control (QC)](#preparation-for-quality-control-(QC))
5. [QC and filtering](#QC-and-filtering)
6. [Trimming with trimal](#trimming-with-trimal)
7. [Check effects of trimming](#check-effects-of-trimming)
8. [Tree construction](#tree-construction)
9. [Tree visualization](#tree-visaulization)
10. [Supplemetary data creation](#supplementary-data-creation)

## 1. Project Overview
This project was built to analyze the copepod species of the genus <i>Pseudodiaptomus</i> but can be used for any sequences.


## 2. Fetch sequences 

Before using the script download the outgroup manually:

<b><i>Add an outgroup for root!</b></i> Can be from the same family but different genus or closely related family. BlastN your data, select to show 5000 sequences and pick one from the list.

> **Use the script**: `API_NCBI.ipynb`
<p>This script uses API requests to fetch sequences and metadata from NCBI</p>

1. Adds your fetched sequences, your own data and the outgroup.
2. Includes all names of the gene you're working on
3. If complete genomes exist, adds only the part you're interested in

Usual databases: NCBI, BOLD (for COI). <b>Best:</b> download data from the whole genus

><i><b>Note</b>: your local files and outgroup must include the name of the gene and organism to be appended. New folders will be created! </i>


## 3. Multiple sequence alignment
There are many tools: MAFFT, MUSCLE, Clustal Omega, T-Coffee etc. 

1. Perform multiple sequence alignment (MSA) in the raw data
2. Cut the sequences at the boundaries of your own data if you have any
3. Optional: perform MSA again

## 4. Preparation for Quality Control (QC)

If you have sequences of ITS, you will likely have the <b>whole ITS1-5.8S-ITS2 operon</b>. Before filtering and extra trimming it's best that you have identified the 5.8S area so that it's easier to find it again when more columns are removed.

The 5.8S gene is bout 160-170bp. Whith this in mind:

1. Open the alignment in a viewer like Jalview
2. Find a sequence among them that is annotated in NCBI and copy the 5.8S gene area <b>without gaps</b>
3. Search for this known sequence
4. Adjust the edges: if you find conserved nucleotides among your sequences include them BUT remember the overall length

## 5. QC and filtering

><b>Use the script </b>: QC_haplotype_phylogenetics.R

<b>This script:</b>

1. Calculates the percentage of gaps and ambiguous bases (N)
2. Removes sequences have more than 45% gaps
3. Identifies sequences that are 100% identical and collapses them in two steps: Before trimming and after
4. Fixes the reading frame for coding genes (dividable to 3, no stop codons)
5. Creates partiotion file for coding genes and for operons
6. Creates fasta for trimming, final fasta and "receipts"

## 6. Trimming with trimal

After the initial haplotype collapsing, perform trimming with trimal using these steps in our HPC:
<i> you can change the names in []. Remember to <b>remove the []</b></i>

```bash
module load  gcc/14.2.0 miniconda3/24.7.1
source $CONDA_PROFILE/conda.sh
conda create [trimal]
conda activate [trimal]
conda install bioconda::trimal
trimal -in [alignment_for_trimal].fasta -out [alignment_trimmed].fasta -automated1
```
<i>! Remember to be in the folder in which you have your alignment_for_trimal.fasta </i>

<b>After this continue with the QC_haplotype_phylogenetics.R pipeline </b>

## 7. Check effects of trimming

To check the effects of trimming, variable site, GC contect, gaps etc. we can use AMAS through our HPC.
<i> you can change the names in []. Remember to <b>remove the []</b></i>

```bash
module load  gcc/14.2.0 miniconda3/24.7.1
source $CONDA_PROFILE/conda.sh
conda create [amas]
conda activate [amas]
conda install bioconda::amas
```

Create a summary of the three fasta files

```bash
conda run -n amas AMAS.py summary -f fasta -d dna -i [] -o []
```
The reason why we use AMAS.py and not amas is because sometimes it's downloaded in a different name, to check use these in your selected conda environment:

```bash
echo $PATH
grep -i amas
```
Use the answer, in this case: AMAS.py

## 8. Tree construction

For Maximum Likelihood (ML) trees, the most famous and fast algorith is the IQTREE. We use that in the HPC.
In the Aristotelis HPC of Aristotle Univeristy of Thessaloniki it's already installed.

```bash
module load gcc/14 iq-tree/2.3.2
```
If you need to download it:
<i> you can change the names in []. Remember to <b>remove the []</b></i>

```bash
module load  gcc/14.2.0 miniconda3/24.7.1
source $CONDA_PROFILE/conda.sh
conda create [iqtree]
conda activate [iqtree]
conda install bioconda::iqtree
```

<h3>IQTREE commands</h3>

<p><b>-m MFP </b>: Model Finder Pro, determins best model for tree construction</p>
<p><b>-m MFP+MERGE </b>: Model Finder Pro and greedy algorith MERGE, finds best model for each partition and megres them if statistically better</p>
<p><b>-m TESTONLY </b>: makes statistics without tree construction</p>
<p><b>-B or -bb </b>: Ultrafast Bootstrap (UFboot) </p>
<p><b>-alrt </b>: SH-aLRT (Shimodaira-Hasegawa approximate Likelihood Ratio Test)</p>
<p><b>-bnni </b>: optimize UFBoot trees by nearest neighbor interchange (NNI) on corresponding bootstrap alignments</p>
<p><b>-nt AUTO </b>: automatically determines the optimal number of cores </p>
<p><b>-st  </b>: sequence type</p>
<p><b>-s  </b>: sequences file</p>
<p><b>-p </b>: partition file</p>
<p><b>-pre </b>: name of the files created</p>


1. Run the -m TESTONLY first to check for sequences that fail the chi2 test. If they are not sequences of interest you can remove them.

```bash
iqtree2 -s [alignment].fasta -m TESTONLY
```
for the coding genes with partition specifically:
<b> Change the translation table accordingly, NT2AA5 is for invertebrate mitochondrial </b>

```bash
iqtree2 -s [alignment].fasta -m TESTONLY -st NT2AA5
```

<b>For the trees without -m TESTONLY</b>

2. Run ML trees first without partition (partition is used for coding genes or for the ITS "operon").
<i>For publications you could use bootstrap of 10,000, but for an initial tree to be faster you can just use 1,000</i>

```bash
iqtree2 -s [alignment].fasta -m MFP -B 1000 -alrt 2000 -nt AUTO -pre [no_partition] -bnni
```

3. Run the partition model

```bash
iqtree2 -s [alignment].fasta -p [alignment_partition].nex -m MFP+MERGE -B 1000 -bnni -alrt 2000 -nt AUTO -pre [partition]
```

4. Check the AIC, BIC etc. statistics from the .iqtree file to decide on the tree (partition or not)
5. Download the .treefile

### Tree visualization

><b>Use the script</b>: visualization_metadata_phylogenetics.R

This script:
1. Uses the .treefile and metadata file created from the QC_haplotype_phylogenetics.R script
2. Roots the tree using the Haplotype name
3. Integrates the metadata to create a title: <b>Accession number</b> <i> Species name</i> | Country [n=(number of collapsed haplotypes)]
4. Depicts Bootstrap and SH-aLRT (Shimodaira-Hasegawa approximate Likelihood Ratio Test) into numbers or dots with a legend

### Supplementary data creation

Use the _representative_map.tsv and copy the Representative accession number and Collapsed haplotypes.
MOdify accordingly (names of columns, comma instead of semicolon etc)