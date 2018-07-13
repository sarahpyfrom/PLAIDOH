README for PLAIDOH lncRNA Annotation Software

After downloading from github, PLAIDOH is all ready to go! <--- Currently a lie! Must manually unzip the two zip files. All else should be golden. 
Supplementary information and a description of PLAIDOH output files can be found in Pyfrom et al, 2018 (PLAIDOH: A Novel Approach for Functional Prediction of Long Non-Coding RNAs Identifies Cancer-Specific LncRNA Activities). We suggest you take a look at Table 1!

0. Downloaded Files

When you Download PLAIDOH from github, it should come with the following directories/files: 

Default_Files
	Biomartquery_hg19.txt.zip
	ChIP_Files
		<lots of tissue filenames>
	EnhancerAtlas_Enhancers.bed
	K562_and_MCF7_POL2RAChIApet.txt.zip
	Nuclear_Fraction_GM12878_hg19.txt
	RBPs
		<lots of RBP filenames>
	SEA_SuperEnhancers.bed
	hESCTADboundariesCombined.bed
Example_Input_File.txt
PLAIDOH.pl
README.txt



1. Quick Start

To test your dependencies and implementation, I recommend you run the tester file the first time you run PLAIDOH using the following command:


perl PLAIDOH.pl --cellENH A549 --cellSE A549 --cellCHIP A549 --RBP NONO Example_Input_File.txt



The first time you run PLAIDOH, it may take some additional time as it checks for all the necessary default files. This process should only happen once and, upon generating those files, it will run the chosen command/query. It should generate two files:

Default_Files/RNABindingProtein_eCLIP_h19.txt
Default_Files/ENCODE_ChIPseq_pValues.txt

and will unzip: 
K562_and_MCF7_POL2RAChIApet.txt.zip
Biomartquery_hg19.txt.zip

To run PLAIDOH on your own dataset, you should use the following statement: 

perl PLAIDOH_JULY042018.pl [options] <Input tab-delimited bed file with columns: #Chr Start Stop Name Type(lncRNA, antisense_RNA or protein_coding) Sample1_Expression Sample2_Expression etc >
                
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
       
       -f finename for Nuclear Fraction file
       
       -b filename for Biomart Query file
        
        -c filename for Chip-seq data
           --cellCHIP Optional selection to pick a specific cell line
                variable from the ChipSeq list
        
        -p filename for CHIA-pet data
        
        -r filename for RBP file
           --RBP Optional selection to pick a specific cell line
                variable from the RBP list
