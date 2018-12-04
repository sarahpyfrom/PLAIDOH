# PLAIDOH LncRNA Annotation Software
`Manuscript submitted. Code made public for reviewer access, but in Beta-testing until manuscript review and publication. Please feel free to get in touch with requests or comments: sarahpyfrom@gmail.com`

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

This process is only necessary once. It should create two Default files:

```
Default_Files/RNABindingProtein_eCLIP_h19.txt
Default_Files/ENCODE_ChIPseq_pValues.txt
```

And will unzip: 
```
K562_and_MCF7_POL2RAChIApet.txt.zip
Biomartquery_hg19.txt.zip
```

The test run should also create all of PLAIDOH's standard output files for the Example input file:

`PLAIDOH_OUTPUT_Example_Input_File.txt` This is the primary PLAIDOH output file described in detail in the table below.
 
`PLAIDOH_RUN_LOG.txt` Each time PLAIDOH.pl is run, this file is updated with a new entry, recording the input files and options selected by the user and including "COMPLETED SUCCESSFULLY" with a date and time at which the program finished running if PLAIDOH successfully completed its run and all calculations. 

`lncs_Example_Input_File.txt` contains all entries from the input file with "lncRNA" or "antisense_RNA" in the "Type" column.

`protein_coding_Example_Input_File.txt` contains all entries from the input file with "protein_coding" in the "Type" column.

`misc_Example_Input_File.txt` contains all entries from the input file with anything OTHER THAN "protein_coding", "lncRNA" or "antisense_RNA" in the "Type" column.

`User_Selected_Enhancers.txt`, `User_Selected_Super_Enhancers.txt`, and `User_Selected_RBPs.txt`, are the subset of each datatype that the user's options selected from the chosen Enhancer, Super enhancer, and RBP files. 

`lncRNAs_bound_by_RBPs.txt` is output from a bedtools intersect between lnc_Input_Filename.txt and any RBPs selected by the user.

`lncs_and_proteins_intersect_Example_Input_File.txt` is output from a bedtools intersect between lncs_Example_Input_File.txt and protein_coding_Example_Input_File.txt.

`lncRNAs_in_SuperEnhancers.txt` is output from a bedtools intersect between lncs_Example_Input_File.txt and any User_Selected_Super_Enhancers.txt.


PLAIDOH also outputs three R-generated plots of your data: 

`RankedCisRegulatoryScores.jpg` is a graph plotting every Cis-Regulatory score for every LncRNA/Coding gene Pair (LCP) identified by PLAIDOH, ranked from lowest to highest. A linear regression line is plotted in red and can be used to determine the best cut-off score for a likely cis-regulatory lncRNA. 

`RankedEnhancerScores.jpg`  is a graph plotting every Enhnacer score for every LCP identified by PLAIDOH, ranked from lowest to highest. A linear regression line is plotted in red and can be used to determine the best cut-off score for a likely enhancer-associated lncRNA. 

`CisRegulatoryScorebyEnhancerScore.jpg` is a graph plotting the Cis-regulatory and Enhancer score for each LCP.
 

To run **PLAIDOH** on your own dataset, you should use the following command with the appropriate options to specify your own datasets: 

```
    USAGE

        perl PLAIDOH.pl [options] <Input Expression File> 
                
    REQUIRED INPUT FILE
        
        Tab-delimited file with columns: #Chr  Start   Stop   Name   Type(lncRNA, antisense_RNA or protein_coding) Sample1_Expression Sample2_Expression etc
        
        The input file should include all lncRNA and protein-coding transcripts of interest.
        
        Please see Example_Input_File.txt for an Example of a properly-formatted Input file. 
         
                Please Note: Input expression file and all user-generated bed
                files should be sorted using the following command if you wish
                to use any default PLAIDOH input files:
                sort -k1,1 -k2,2n in.txt > in_sorted.txt
                
    OPTIONS
        
       -s filename for Super-enhancer file.
            Default file: SEA_SuperEnhancers.bed
        
        --cellSE Optional selection to pick a specific cell line
               from the super enhancer list
           
       -e filename for Enhancer data.
            Default file: EnhancerAtlas_Enhancers.bed
           
        --cellENH Optional selection to pick a specific cell line
               from the Enhancer list
       
       -f finename for Nuclear Fraction file
            Default file: Nuclear_Fraction_GM12878_hg19.txt
       
       -b filename for Biomart Query file (includes GO terms and other annotatin options)
            Default file: Biomartquery_hg19.txt
        
       -c filename for Chip-seq data
            Default file: ENCODE_ChIPseq_pValues.txt
            
        --cellCHIP Optional selection to pick a specific cell line
                from the ChipSeq list
        
       -p filename for CHIA-pet data
             Default file: K562_and_MCF7_POL2RAChIApet.txt   
        
       -r filename for RBP file
            Default file: RNABindingProtein_eCLIP_h19.txt
        
        --RBP Optional selection to pick a specific RBP
                variable from the RBP list

                
        Detailed descriptions for all input files and their sources can be found in the Methods
        secition of Pyfrom, Luo and Payton, 2018 BMC Genomics.
```

## PLAIDOH Output File

PLAIDOH outputs a several supporting files and a single final file with the name "Output_$input_filename". It has 42 tab-delimited columns, which are described in the following table. Each line contains information about a single lncRNA and one coding gene within 400kb of the lncRNA:

Column Number | Column Heading | Description
----------------|-------------|-----------
1 | LINE_NUMBER | Unique numerical designation for each lncRNA and protein pair (LCP)
2 | CHRlnc | Chromosome of lncRNA
3 | STARTlnc | Downstream chromosomal position of lncRNA
4 | STOPlnc | Upstreamstream chromosomal position of lncRNA
5 | NAMElnc | Name of lncRNA
6 | TYPElnc | "lncRNA"
7 | NUMBERlnc | Unique numerical designation for each lncRNA; different for each isoform of a given lncRNA included in input
8 | CHRcoding | Chromosome of protein
9 | STARTcoding | Downstream chromosomal position of protein
10 | STOPprot | Upstreamstream chromosomal position of protein
11 | NAMEcoding | Name of protein
12 | TYPEcoding | "protein_coding"
13 | NUMBERcoding | Unique numerical designation for each protein; different for each isoform of a given protein included in input
14 | LNC_CODING_OVERLAP | "Overlap" if current lncRNA and protein overlap by 1+ base pairs. "Distal" if there is no overlap.
15 | DIST_L_to_C | Shortest distance between lncRNA and coding-gene as determined by START and STOP positions from input.
16 | LNC_AVERAGE | Average expression of lncRNA across all input samples.
17 | PROTEIN_AVERAGE | Average expression of coding-gene across all input samples.
18 | PEARSON | Pearson correlation between lncRNA and protein expression as calculated by R cor.test package.
19 | PEARSONp | p-value of above pearson correlation.
20 | SPEARMAN | Spearman correlation between lncRNA and protein expression as calculated by R cor.test package.
21 | SPEARMANp | p-value of above spearman correlation.
22 | CHIA | "0" if no evidence of ChIA-PET interaction between lncRNA and coding-gene. Any >0 value represents the score of interaction (between 100 and 1000 if default input is used).
23 | TAD | 1 if both lncRNA and protein are both at least partly within the same TAD, 0 if entirety of lncRNA or protein is in a separate TAD.
24 | INTERGENIC | 1 if lncRNA does not overlap with ANY "protein_coding" gene in the input file. 0 if lncRNA overlaps with ANY "protein_coding" gene in the input file.
25 | CHIPSEQ | Full list of ChIP-seq peaks and associated scores that overlap any part of the lncRNA. 
26 | H3K4ME3 | Highest -log10 p-value of any overlapping H3K4ME3 ChIP-seq peak. 0 if no overlap.
27 | H3K27AC | Highest -log10 p-value of any overlapping H3K27AC ChIP-seq peak. 0 if no overlap.
28 | H3K4ME1 | Highest -log10 p-value of any overlapping H3K4ME1 ChIP-seq peak. 0 if no overlap.
29 | RNA Binding Proteins | List of all RNA Binding proteins and their binding score (200 or 1000) bound to lncRNA.
30 | NEAR_ENH | Shortest distance between lncRNA and closest enhancer as determined by START and STOP positions from input files.
31 | NEAR_SE | Shortest distance between lncRNA and closest super enhancer as determined by START and STOP positions from input files.
32 | GO | Gene Ontology for Protein. Matched by name of protein. Empty if exact match for protein not found in GO input file.
33 | L_strand | "1" if lncRNA on + strand. "-1" if lncRNA on negative strand. Empty if exact match for protein not found in GO input file.
34 | C_strand | "1" if protein on + strand. "-1" if protein on negative strand. Empty if exact match for protein not found in GO input file.
35 | TSStoTSSdist | Distance from TSS of lncRNA to TSS of protein. Only available if both lncRNA and protein names can be found in GO input file.
36 | LNC_NUC_FRACTION | Fraction of lncRNA transcript found in nucleus. If default files are used, this refers to the fraction of reads found in the nucleus vs cytoplasm of GM12878 cell line from ENCODE.
37 | CODING_NUC_FRACTION | Fraction of protein transcript found in nucleus. If default files are used, this refers to the fraction of reads found in the nucleus vs cytoplasm of GM12878 cell line from ENCODE.
38 | SENSE_CATEGORY | "Distal" = The lncRNA and coding gene are non-overlapping AND (separated by at least 2kb OR they are on the same strand). "Antisense_Proximal" = The lncRNA and coding gene's TSS' are less than 2kb apart on opposing strands, but the transcripts do not overlap. "Antisense_Overlap" = The lncRNA and coding gene's TSS' are less than 2kb apart on opposing strands, and the transcripts do overlap by one bp or greater. " "Overlap" = The lncRNA and coding gene overlap but their TSS' are greater than 2kb apart OR they are on the same strand.
39 | ADJUSTED_SPEARMAN | Adjusted Spearman p-value of lncRNA/protein correlation to correct for multiple tests. 
40 | SCORE | Cis-regulatory score calculated by PLAIDOH
41 | ENHANCER_SCORE | Enhancer Score calculated by PLAIDOH
42 | FRACTION_SCORE | For each possible lncRNA/RBP combination a “Fraction Score” was calculated. If there is no evidence that the lncRNA and RBP interact based on the ENCODE eCLIP data, the pair is given a score of 0. If an RBP binds a given lncRNA and the sub-cellular localization of both the lncRNA and corresponding RBP are both determined to be Nuclear or both are Cytoplasmic, the pair is given a score of 2 or -2, respectively. If the lncRNA is primarily Nuclear and the RBP does not have a Nuclear localization by Immunofluorescence (or vice versa), the overlapping region is given a score of 1 or -1. If both the RBP and lncRNA are considered to be present in both the nucleus and cytoplasm (as described above) the pair is given a score of 3. 
