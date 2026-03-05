# Phylogenetics_pipeline
My data is from zooplankton!

Data analysis is conducted from genes: COI, 18S, ITS

<h2>1. Choose your data.</h2>
Usual databases: NCBI

<b>Moto:</b> go big or go home

<b>Best:</b> download data from the whole genus, because especially for COI, a phylogenetic tree for species only is generally not good.


<b><i>Add an outgroup for root!</b></i>

<h2>2. Bulild your alignment</h2>
Add your main sequences, your own data and the outgroup.
You can do this manually in a fasta file, or in an app like MEGA11, Unipro UGENE, JalView
or with coding

<h2>3. Multiple sequence alignment</h2>
There are many tools out there: MAFFT, MUSCLE, Clustal Omega, T-Coffee etc

<h2>4. Trim your alignment</h2>
Sometimes the start and end of the alignment can be messy.

In JalView or other fasta editors we can check the occupancy of the sequences and trim the start and end for positions for example under 80% occupancy

<h3>Coding genes</h3>
For codon partitioning we need to have exact length of sequences for all sequences that is dividable to 3.


