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
11. [Bibliography](#bibliography)

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
>you can change the names in []. Remember to <b>remove the []</b>

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

## 9. Tree visualization

><b>Use the script</b>: visualization_metadata_phylogenetics.R

This script:
1. Uses the .treefile and metadata file created from the QC_haplotype_phylogenetics.R script
2. Roots the tree using the Haplotype name
3. Integrates the metadata to create a title: <b>Accession number</b> <i> Species name</i> | Country [n=(number of collapsed haplotypes)]
4. Depicts Bootstrap and SH-aLRT (Shimodaira-Hasegawa approximate Likelihood Ratio Test) into numbers or dots with a legend

## 10. Supplementary data creation

Use the _representative_map.tsv and copy the Representative accession number and Collapsed haplotypes.
MOdify accordingly (names of columns, comma instead of semicolon etc)


## 11. Bibliography

<details>
<summary><b>Click to expand Software Citations</b></summary>

R Core Team (2025). _R: A Language and Environment for Statistical
  Computing_. R Foundation for Statistical Computing, Vienna, Austria.
  <https://www.R-project.org/>.

Wickham H, Hester J, Bryan J (2026). _readr: Read Rectangular Text
  Data_. doi:10.32614/CRAN.package.readr
  <https://doi.org/10.32614/CRAN.package.readr>, R package version
  2.2.0, <https://CRAN.R-project.org/package=readr>.

Paradis E, Schliep K (2019). “ape 5.0: an environment for modern
  phylogenetics and evolutionary analyses in R.” _Bioinformatics_,
  *35*, 526-528. doi:10.1093/bioinformatics/bty633
  <https://doi.org/10.1093/bioinformatics/bty633>.

Paradis E (2010). “pegas: an R package for population genetics with
  an integrated-modular approach.” _Bioinformatics_, *26*, 419-420.
  doi:10.1093/bioinformatics/btp696
  <https://doi.org/10.1093/bioinformatics/btp696>.


Guangchuang Yu. (2022). Data Integration, Manipulation and
  Visualization of Phylogenetic Trees (1st edition). Chapman and
  Hall/CRC. <https://doi:10.1201/9781003279242>

Shuangbin Xu, Lin Li, Xiao Luo, Meijun Chen, Wenli Tang, Li Zhan, Zehan Dai, Tommy T. Lam, Yi Guan, Guangchuang Yu.      Ggtree: A serialized data object for    visualization of a phylogenetic tree and annotation data iMeta 2022, 4(1):e56. <https:doi:10.1002/imt2.56>

Guangchuang Yu. Using ggtree to visualize data on tree-like
structures. _Current Protocols in Bioinformatics_, 2020, 69:e96. <https:doi:10.1002/cpbi.96>

Guangchuang Yu, Tommy Tsan-Yuk Lam, Huachen Zhu, Yi Guan. Two methods for mapping and visualizing associated data on phylogeny using ggtree. _Molecular Biology and Evolution 2018_, 35(2):3041-3043. <https:doi:10.1093/molbev/msy194>

Guangchuang Yu, David Smith, Huachen Zhu, Yi Guan, Tommy Tsan-YukLam. ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. _Methods in Ecology and Evolution 2017_, 8(1):28-36. <https:doi:10.1111/2041-210X.12628>

H. Wickham. ggplot2: Elegant Graphics for Data Analysis.
  _Springer-Verlag New York_, 2016.

Wickham H, François R, Henry L, Müller K, Vaughan D (2026). _dplyr: A
Grammar of Data Manipulation_. doi:10.32614/CRAN.package.dplyr <https://doi.org/10.32614/CRAN.package.dplyr>, R package version 1.2.0, <https://CRAN.R-project.org/package=dplyr>.

Wickham H (2025). _stringr: Simple, Consistent Wrappers for Common
  String Operations_. doi:10.32614/CRAN.package.stringr
  <https://doi.org/10.32614/CRAN.package.stringr>, R package version 1.6.0, <https://CRAN.R-project.org/package=stringr>.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Abbassi, L., Gebremedhin, M. T., & de Hoon, M. J. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423. <https://doi.org/10.1093/bioinformatics/btp163>

Eric W Sayers, Evan E Bolton, J Rodney Brister, Kathi Canese, Jessica Chan, Donald C Comeau, Ryan Connor, Kathryn Funk, Chris Kelly, Sunghwan Kim, Tom Madej, Aron Marchler-Bauer, Christopher Lanczycki, Stacy Lathrop, Zhiyong Lu, Francoise Thibaud-Nissen, Terence Murphy, Lon Phan, Yuri Skripchenko, Tony Tse, Jiyao Wang, Rebecca Williams, Barton W Trawick, Kim D Pruitt, Stephen T Sherry, Database resources of the national center for biotechnology information, _Nucleic Acids Research_, Volume 50, Issue D1, 7 January 2022, Pages D20–D26, <https://doi.org/10.1093/nar/gkab1112>


The pandas development team. (2026). pandas-dev/pandas: Pandas (v3.0.2). Zenodo. <https://doi.org/10.5281/zenodo.19340003>

Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. CreateSpace, Scotts Valley, CA.
<https://www.python.org/>

Capella-Gütiérrez, S., Silla-Martínez, J. M. & Gabaldón, T. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. _Bioinformatics_ 25, 1972–1973 (2009). 

Minh, B. Q. et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic _Inference in the Genomic Era. Mol. Biol. Evol. 37_, 1530–1534 (2020).

Kalyaanamoorthy, S., Minh, B. Q., Wong, T. K. F., von Haeseler, A. & Jermiin, L. S. ModelFinder: Fast model selection for accurate phylogenetic estimates. _Nat. Methods 14_, 587–589 (2017).

Chernomor, O., von Haeseler, A. & Minh, B. Q. Terrace Aware Data Structure for Phylogenomic Inference from Supermatrices. _Syst. Biol._ 65, 997–1008 (2016).

Hoang, D. T., Chernomor, O., von Haeseler, A., Minh, B. Q. & Vinh, L. S. UFBoot2: Improving the Ultrafast Bootstrap Approximation. _Mol. Biol. Evol. 35_, 518–522 (2018).

</details>