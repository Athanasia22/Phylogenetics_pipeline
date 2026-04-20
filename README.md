# Phylogenetics_pipeline
My data is from zooplankton!

Data analysis is conducted for genes: COI, 18S, ITS

1. [Project Overview](#project-overview)
2. [Fetch Sequences](#fetch-sequences)

## Project Overview
This project analyzes copepod species of the genus <i>Pseudodiaptomus</i>


## Fetch sequences 
Usual databases: NCBI, BOLD (for COI)

<b>Best:</b> download data from the whole genus, because especially for COI, a phylogenetic tree for species only is generally not good.

<b><i>Add an outgroup for root!</b></i> Can be from the same family but different genus or closely related family. BlastN your data, select to show 5000 sequences and pick one from the list.



<h2>2. Bulild your alignment</h2>
Add your main sequences, your own data and the outgroup.
You can do this manually in a fasta file (easiest), or in an app like MEGA11, Unipro UGENE, JalView
or with coding.

1. Make sure to include all names of the gene you're working on:
   
   <p> a. (Pseudodiaptomus[Organism]) AND (COI OR COX1 OR "cytochrome oxidase subunit I" OR "cytochrome c oxidase subunit 1")</p>
   <p> b. (Pseudodiaptomus[Organism]) AND (18S OR "small subunit ribosomal RNA" OR SSU) </p>
   <p> c. (Pseudodiaptomus[Organism]) AND (ITS OR "internal transcribed spacer" OR "ITS1" OR "ITS2") </p>
2. Download from NCBI: Send to--> Complete Record--> File... FASTA
3. If you have any complete genomes and are a fraction of all of your sequences you can manually download for them the coding sequences for them and keep only the gene you like. Then, edit the headers and insert them into the fasta with the rest of the sequences 


<h2>3. Multiple sequence alignment</h2>
There are many tools out there: MAFFT, MUSCLE, Clustal Omega, T-Coffee etc. Most used ones are MAFFT and ClustalW.


In JalView or other fasta editors we can check the occupancy of the sequences and trim the start and end for positions for example under a threshold of occupancy. Or use trimmming algorithms, such as TrimAl. For 18S sometimes the sequences differ greatly in length depending on the sequencing machine. For 18S since it's non coding we can more freely trim the sequences.

The first steps to clearing up the data are these:


1. Peform MSA in an editor app e.g. Jalview with e.g. MAFFT
2. Cut all the sequences at the length of your own sequence (If you don't have any reference sequence, cut where most sequences start/end).
Using a custom code:
3. Remove sequences that don't cover at least 45% of your sequence
4. Collapse them into haplotypes

10. The final fasta and "receipts" will get created

11. To check the effects of trimming, variable site, GC contect, gaps etc.:

We can use AMAS and download through BioConda to our HPC.

```bash
module load  gcc/14.2.0 miniconda3/24.7.1
source $CONDA_PROFILE/conda.sh
conda create []
conda activate []
conda install bioconda::amas
```

Create a summary of the three fasta files

```bash
conda run -n amas AMAS.py summary -f fasta -d dna -i []  -o []
```
The reason why we use AMAS.py and not amas is because sometimes it's downloaded in a different name, to check use these in your selected conda environment:

```bash
echo $PATH
grep -i amas
```
Use the answer, in this case: AMAS.py


Create the partition nexus file </h1>


<h1> Going to tree-construction: </h1>

1. Run ML trees first without partition (partition is used for coding genes only)
   
2. Check if some sequences fail the chi2 test and remove them (If we removed sequences with only a few bases and many gaps it's likely that not many sequences will fail the chi2 test. But since it depends on the composition of the sequences we inserted it heavily depends on which sequences we are comparing each time. So, if some sequences fail in one test might not fail in another.)
3. Run again the non-partition if neccesary and then run the partition
4. Check the AIC, BIC etc. statistics to decide on the tree (partition or not)




<h3>Coding genes</h3>
For codon partitioning we need to have exact length of sequences for all sequences that is dividable to 3.

<h2>Make nexus file</h2>
If we use IQ TREE without codon partition we don't need a nexus file. If we do, or use other programms like Mr Bayes, we do need a nexus file. But we do need a nexus file for codon partition.