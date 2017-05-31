# prepABEL
A simple script combining BASH line commands to grab data from dosage VCF files coming from Michigan Imputation Server and prepare input files to run survival analyses with ProbABEL, a software to do survival analysis with dosage data. prepABEL prepares the following ProbABEL input files:
- Genomic Predictor file (.gp), 
- INFO file (.info),
- MAP (.map; optional for ProbABEL)

prepABEL chunks dose VCF files according to the number of rows and transpose the SNP x individuals data to an individual x SNP matrix file (Fig. 1), in the format required by ProbABEL.

prepABEL is programmed to be run on supercomputer resources, but can be easily adapted to run in a stand-alone context. This can be done modifying the number of lines to be processed (see "number_lines" variable; default value: 50000). An example script to run prepABEL for autosomal chromosomes (1 to 22) is also provided here.

prepABEL may also be modified to chunk the output files into smallers files if lack of memory is a limit to run ProbABEL with large genomic predictor files. An example script to perfom the chunking is provided too.

prepABEL is being presented in the ATS 2017 Conference in Washington, 19-24 May, 2017. It has been developed by our staff at the Genomics Divion of the Institute of Energy and Renewable Energy (ITER), Tenerife, Canary Islands, SPAIN.

prepABEL will be available for downloading in June 1st, 2017. Please, feel free to contact me at "jlorsal" @ "gmail.com".

![prepABEL layout](https://github.com/genomicsITER/prepABEL/blob/master/prepABEL_layout.png)
Figure 1. PrepABEL is able to transpose vcf data into a genomic predictor file (a), and capture selected columns to prepare an INFO SNP (b) and MAP (c) files to conduct survival analysis.

