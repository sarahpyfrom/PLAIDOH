# PLAIDOH LncRNA Annotation Software

Supplementary information and a description of **PLAIDOH** output files can be found in Pyfrom et al, 2018 (PLAIDOH: A Novel Approach for Functional Prediction of Long Non-Coding RNAs Identifies Cancer-Specific LncRNA Activities). We suggest you take a look at Table 1!

## Installation

After cloning or downloading **PLAIDOH** from github, you should have the following directory structure: 

```
|-- Default_Files
    |-- Biomartquery_hg19.txt.zip
    |-- ChIP_Files
	|-- <lots of tissue filenames>
    |-- EnhancerAtlas_Enhancers.bed
    |-- K562_and_MCF7_POL2RAChIApet.txt.zip
    |-- Nuclear_Fraction_GM12878_hg19.txt
    |-- RBPs
	|-- <lots of RBP filenames>
    |-- SEA_SuperEnhancers.bed
    |-- hESCTADboundariesCombined.bed
|-- Example_Input_File.txt
|-- PLAIDOH.pl
|-- README.md
```

### Dependencies
PLAIDOH has few dependencies, but does require a few common bioinformatics tools/languages on your PATH. Most people will likely already have them installed:
 - [bedtools2](http://bedtools.readthedocs.io/en/latest/index.html)
 - [perl5](https://www.perl.org/)
 - [R](https://www.r-project.org/)

## Quick Start

To test your dependencies and create some necessary default files, PLAIDOH should be run on the included test file:

```Bash
perl PLAIDOH.pl --cellENH A549 --cellSE A549 --cellCHIP A549 --RBP NONO Example_Input_File.txt
```

This process should is only necessary once. It should generate two files:

```
Default_Files/RNABindingProtein_eCLIP_h19.txt
Default_Files/ENCODE_ChIPseq_pValues.txt
```

And will unzip: 
```
K562_and_MCF7_POL2RAChIApet.txt.zip
Biomartquery_hg19.txt.zip
```

To run **PLAIDOH** on your own dataset, you should use the following command with the appropriate options to specify your own datasets: 

```Bash
perl PLAIDOH.pl [options] <Input tab-delimited bed file with columns: #Chr Start Stop Name Type(lncRNA, antisense_RNA or protein_coding) Sample1_Expression Sample2_Expression etc >
                
                Please Note: The first line MUST start with a #.
		Input expression file and all user-generated bed
                files should be sorted using the following command if you wish
                to use any default PLAIDOH input files:
                sort -k1,1 -k2,2n in.bed > in.sorted.bed
                
        OPTIONS
        
        -s filename for Super-enhancer file
           --cellSE Optional selection to pick a specific cell line
                variable from the SE list
           
       -e filename for Enhancer data
           --cellENH Optional selection to pick a specific cell line
                variable from the Enhancer list
       
       -f filename for Nuclear Fraction file
       
       -b filename for Biomart Query file
        
        -c filename for Chip-seq data
           --cellCHIP Optional selection to pick a specific cell line
                variable from the ChipSeq list
        
        -p filename for CHIA-pet data
        
        -r filename for RBP file
           --RBP Optional selection to pick a specific cell line
                variable from the RBP list
```
