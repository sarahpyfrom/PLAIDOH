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

This process is only necessary once. It should generate two files:

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

```
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
           --RBP Optional selection to pick a specific RBP
                variable from the RBP list
```

###PLAIDOH Output File

PLAIDOH outputs a several supporting files and a single final file with the name "Master_$input_filename". It has 42 tab-delimited columns, which are described in the following table. Each line contains information about a single lncRNA and one coding gene within 400kb of the lncRNA:

Column Number	Column Heading	Description
1	LINE_NUMBER	Unique numerical designation for each lncRNA and protein pair (LCP)
2	CHRlnc	Chromosome of lncRNA
3	STARTlnc	Downstream chromosomal position of lncRNA
4	STOPlnc	Upstreamstream chromosomal position of lncRNA
5	NAMElnc	Name of lncRNA
6	TYPElnc	"lncRNA"
7	NUMBERlnc	Unique numerical designation for each lncRNA; different for each isoform of a given lncRNA included in input
8	CHRcoding	Chromosome of protein
9	STARTcoding	Downstream chromosomal position of protein
10	STOPprot	Upstreamstream chromosomal position of protein
11	NAMEcoding	Name of protein
12	TYPEcoding	"protein_coding"
13	NUMBERcoding	Unique numerical designation for each protein; different for each isoform of a given protein included in input
14	LNC_CODING_OVERLAP	"Overlap" if current lncRNA and protein overlap by 1+ base pairs. "Distal" if there is no overlap.
15	DIST_L_to_C	Shortest distance between lncRNA and coding-gene as determined by START and STOP positions from input.
16	LNC_AVERAGE	Average expression of lncRNA across all input samples.
17	PROTEIN_AVERAGE	Average expression of coding-gene across all input samples.
18	PEARSON	Pearson correlation between lncRNA and protein expression as calculated by R cor.test package.
19	PEARSONp	p-value of above pearson correlation.
20	SPEARMAN	Spearman correlation between lncRNA and protein expression as calculated by R cor.test package.
21	SPEARMANp	p-value of above spearman correlation.
22	CHIA	"0" if no evidence of ChIA-PET interaction between lncRNA and coding-gene. Any >0 value represents the score of interaction (between 100 and 1000 if default input is used).
23	TAD	1 if both lncRNA and protein are both at least partly within the same TAD, 0 if entirety of lncRNA or protein is in a separate TAD.
24	INTERGENIC	1 if lncRNA does not overlap with ANY "protein_coding" gene in the input file. 0 if lncRNA overlaps with ANY "protein_coding" gene in the input file.
25	CHIPSEQ	Full list of ChIP-seq peaks and associated scores that overlap any part of the lncRNA. 
26	H3K4ME3	Highest -log10 p-value of any overlapping H3K4ME3 ChIP-seq peak. 0 if no overlap.
27	H3K27AC	Highest -log10 p-value of any overlapping H3K27AC ChIP-seq peak. 0 if no overlap.
28	H3K4ME1	Highest -log10 p-value of any overlapping H3K4ME1 ChIP-seq peak. 0 if no overlap.
29	RNA Binding Proteins	List of all RNA Binding proteins and their binding score (200 or 1000) bound to lncRNA.
30	NEAR_ENH	Shortest distance between lncRNA and closest enhancer as determined by START and STOP positions from input files.
31	NEAR_SE	Shortest distance between lncRNA and closest super enhancer as determined by START and STOP positions from input files.
32	GO	Gene Ontology for Protein. Matched by name of protein. Empty if exact match for protein not found in GO input file.
33	L_strand	"1" if lncRNA on + strand. "-1" if lncRNA on negative strand. Empty if exact match for protein not found in GO input file.
34	C_strand	"1" if protein on + strand. "-1" if protein on negative strand. Empty if exact match for protein not found in GO input file.
35	TSStoTSSdist	Distance from TSS of lncRNA to TSS of protein. Only available if both lncRNA and protein names can be found in GO input file.
36	LNC_NUC_FRACTION	Fraction of lncRNA transcript found in nucleus. If default files are used, this refers to the fraction of reads found in the nucleus vs cytoplasm of GM12878 cell line from ENCODE.
37	CODING_NUC_FRACTION	Fraction of protein transcript found in nucleus. If default files are used, this refers to the fraction of reads found in the nucleus vs cytoplasm of GM12878 cell line from ENCODE.
38	SENSE_CATEGORY	"Distal" = The lncRNA and coding gene are non-overlapping AND (separated by at least 2kb OR they are on the same strand). "Antisense_Proximal" = The lncRNA and coding gene's TSS' are less than 2kb apart on opposing strands, but the transcripts do not overlap. "Antisense_Overlap" = The lncRNA and coding gene's TSS' are less than 2kb apart on opposing strands, and the transcripts do overlap by one bp or greater. " "Overlap" = The lncRNA and coding gene overlap but their TSS' are greater than 2kb apart OR they are on the same strand.
39	ADJUSTED_SPEARMAN	Adjusted Spearman p-value of lncRNA/protein correlation to correct for multiple tests. 
40	SCORE	Cis-regulatory score calculated by CREPE
41	ENHANCER_SCORE	Enhancer Score calculated by CREPE
42	FRACTION_SCORE	For each possible lncRNA/RBP combination a “Fraction Score” was calculated. If there is no evidence that the lncRNA and RBP interact based on the ENCODE eCLIP data, the pair is given a score of 0. If an RBP binds a given lncRNA and the sub-cellular localization of both the lncRNA and corresponding RBP are both determined to be Nuclear or both are Cytoplasmic, the pair is given a score of 2 or -2, respectively. If the lncRNA is primarily Nuclear and the RBP does not have a Nuclear localization by Immunofluorescence (or vice versa), the overlapping region is given a score of 1 or -1. If both the RBP and lncRNA are considered to be present in both the nucleus and cytoplasm (as described above) the pair is given a score of 3. 
