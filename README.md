# prepABEL
A simple script combining BASH line commands to grab data from dosage VCF files coming from Michigan Imputation Server and prepare input files to run survival analyses with ProbABEL, a software to do survival analysis with dosage data. prepABEL prepares the following ProbABEL input files:
- Genomic Predictor file (.gp), 
- INFO file (.info),
- MAP (.map; optional for ProbABEL)

prepABEL chunks dose VCF files according to the number of rows and transpose the SNP x individuals data to an individual x SNP matrix file (Fig. 1), in the format required by ProbABEL.

prepABEL is programmed to be run on supercomputer resources, but can be easily adapted to run in a stand-alone context. This can be done modifying the number of lines to be processed (see "number_lines" variable; default value: 50000). An example script to run prepABEL for autosomal chromosomes (1 to 22) is also provided here.

prepABEL may also be modified to chunk the output files into smallers files if lack of memory is a limit to run ProbABEL with large genomic predictor files. An example script to perfom the chunking is provided too.

prepABEL has been presented in the ATS 2017 Conference in Washington, 19-24 May, 2017. It has been developed by our staff at the Genomics Divion of the Institute of Energy and Renewable Energy (ITER), Tenerife, Canary Islands, SPAIN.

prepABEL code is available below. Please, feel free to contact me at "jlorsal" @ "gmail.com" for any further information. 

<p align="center">
  <img src="https://github.com/genomicsITER/prepABEL/blob/master/prepABEL_layout.png" width="auto"/>
</p>

<b>Figure 1</b>. PrepABEL is able to transpose vcf data into a genomic predictor file (a), and capture selected columns to prepare an INFO SNP (b) and MAP (c) files to conduct survival analysis.


# Citation

Please, cite prepABEL as follows:

# PrepABEL: A Bioinformatic Tool for Preprocessing Imputed Genotype Dosage Data for Survival Analyses

J. Lorenzo-Salazar, S.F. Ma, J. M. Oldham, I. Noth, and C. Flores

American Thoracic Society International Conference, A6796

May 19-24, 2017, Washington DC, USA


# Code
```Bash
#!/bin/bash
# --------------------------------------------------------------------------------------------------------
# The following code is to be executed in a HPC infraestructure working with SLURM
# Delete code between dashed lines if a local environment will be used instead
# Test this code with a chromosome to check that everything runs smoothly
# Read the first element of the command line array to grab the chromosome number from prepABEL_launcher.sh
chromosome=chr$1  # $1 grabs the first element of shell command line array
#SBATCH -J prepABEL
#SBATCH -p batch
#SBATCH -N 1
#SBATCH --tasks-per-node=16
#SBATCH -t 1-00:00:00
#SBATCH --mem=60000
#SBATCH -o prepABEL_%j.out
#SBATCH -e prepABEL_%j.err
#SBATCH -D /path-to-folder/slurm/
##########################################################
#Send job to a queue: sbatch prepABEL.sh
#########################################################
## Source modules in the actual user profile
source /etc/profile.d/profile.modules.sh
## Load ProABEL module
module load probAbel/0.5
# --------------------------------------------------------------------------------------------------------

######################################################################################################################### #
# Preparation of files required by probABEL to run Survival Analysis on dosage VCF files from Michigan Imputation Server  #
######################################################################################################################### #
# Description of pipeline:
# 1. Prepare the SNP info file for probABEL
# 2. Inspect imputed VCF file and chunk into smaller files (decide whether you lately will want chunked VCFs 
#    with or without headers)
# 3. Extract the list of individuals
# 4. Build up a map file with ids, positions, ref and alt alleles for probABEL
# 5. Iterate over the chunked VCF small files, and collect dosage data
# 6. Paste all intermediate files to get a final Genomic Predictor file for probABEL
# 7. Remove unnecessary files
# 8. Auxiliary file operations: replace missing data with 'NA' and replace tabulators with spaces, if needed 
#    (uncomment lines)

# Unless you want to process everything within the same batch, write down the name of the chromosome imputed VCF file
# Do not include file extensions
#chromosome="chr1" # <--- UPDATE this variable with each chromosome under analysis if running in local mode

echo "============================================================================"
echo "=== 				                      prepABEL 				                       ==="
echo "===	 			                           v1	 			                           ==="
echo "===      A program to prepare files for ProbABEL survival analyses       ==="
echo "===  with dosage VCF imputation files from Michigan Imputation Server    ==="
echo "===         See more at: https://github.com/genomicsITER/prepABEL        ==="
echo "===                           Genomics Division                          ==="
echo "===       Institute of Technology and Renewable Energy (ITER) 2017       ==="
echo "============================================================================"
echo ""
echo ""

# INPUTS
root1=/path-to-base-dir/
input=$root1/$chromosome
file=$chromosome
gp_filename=$chromosome
phenotypes=$root1/phenotypes

# OUTPUTS
root2=/path-to-base-dir/
output=$root2/$chromosome
coxph=$root2/coxph

echo "================================ Path to files ================================"
echo "Analyzing chromosome $chromosome"
echo "Root path: $root1"
echo "Output files saved to $output"
echo ""


# ################### #
# 1. Prepare SNP file #
# ################### #

echo "============================= Preparing SNP file ============================="
echo "Base filename is... " $file

cut -f1-7 $input/$file.info > $output/probabel.$file.info

echo ""


# ##################################################### #
# 2. Inspect VCF file and chunking in smaller VCF files #
# ##################################################### #

# To process dose.VCF from Michigan Imputation Server, you must deal with very large files, 
# problaby containing thousands of millions of data
# A better strategy will be to chunk the imputed VCF file into smaller files
# After collecting the necessary data, all chunked prepared files will be pasted altogether

echo "================================ Inspection of VCF files ================================"
echo "Inspecting the size of the imputed VCF..."

# Count the lines in file
# a=$(wc -l $file | awk '{print $1}')
# echo "The number of lines in the VCF file is: " $a

# A conservative number of lines will be 75.000 for my laptop, considering 1000 individuals x 87500 lines = 
# 87.500.000 variant values per file!
# <<<<<< THIS NUMBER IS VERY IMPORTANT TO AVOID MEMORY LIMITATIONS >>>>>>
number_lines=50000

# Grab the header
echo "Grabbing the header of $file.dose.vcf"
head -n 1000 $input/$file.dose.vcf | grep "^#" > $output/header

# Count the number of lines in header
number_of_lines_in_header=$(wc -l $output/header | awk '{print $1}')
row_to_start_grabbing=$(( $number_of_lines_in_header + 1))

# Grab the non header lines
echo "Grabbing the non header lines of $file.dose.vcf"
echo "Warning: this will take a time depending on your VCF size"

grep -v "^#" $input/$file.dose.vcf > $output/variants

# Results for the abstract: take a 20 GB file and remove the header
# Command: time { grep -v "^#" $file.dose.vcf > variants2; }
# real	6m38.662s
# user	0m2.560s
# sys	0m25.852s

# Results for a paper/abstract: take a 20 GB file and remove the header
# A file for a big chromosome could be about 100 GB, which represents about 30 "real" minutes
# Command: time { tail -n +$row_to_start_grabbing chr21.dose.vcf > variants; }
# real	6m52.265s
# user	0m2.232s
# sys	0m34.892s

echo "Created variants file with no header"

# Split into chunks with $c lines
echo "Splitting the VCF into chunks of smaller size..."
split -d -l $number_lines -a 3 $output/variants $output/variants_

# Reattach the header to each and clean up
echo "Reattaching the VCF header to each chunk VCF and clean up..."


# ################################## #
# 3. Extract the list of individuals #
# ################################## #

# Extract the individuals from the header
echo "Extracting individuals from the VCF file as a row-vector..."

# Take the last row of header, as determined by "number_of_lines_in_header"
awk -v hd=$number_of_lines_in_header 'NR==hd {for (i=10 ; i<=NF ; i++) {printf("%s\t", $i)} print ""}' $output/header > $output/individuals

# Transpose the row-vector of individuals as a column-vector, and insert an index as column 1 and the string "MLDOSE" 
# as column 3
# Reference for transposing a file: http://stackoverflow.com/questions/1729824/transpose-a-file-in-bash

echo "Transposing the row-vector of individuals into a column-vector, inserting indexes and strings..."

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>x { x = NF }
END {    
    for(j=1; j<=x; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print j"->"str" MLDOSE" # Prepare a matrix with columns 1, 2, and 3 following probABEL specifications
    }
}' $output/individuals > $output/t_individuals

echo ""


# ################################# #
# 4. Prepare a MAP file for probABEL #
# ################################# #

echo "================================ Preparing MAP file ================================"

# Use 'variants' file (it does not have a header) and select columns to build up a map file

echo "rs position REF(0) ALT(1)" > $output/guide
awk '{print $3,$2,$4,$5}' $output/variants > $output/probabel.$file.map_
cat $output/guide $output/probabel.$file.map_ > $output/probabel.$file.map
rm -f $output/guide $output/probabel.$file.map_

echo ""


# ###################################################################################### #
# 5. Iterate over the number of chunked VCF files and prepare the GENOMIC PREDICTOR file #
# ###################################################################################### #

# For each chunked VCF small file

for file in $output/variants_*; do
echo ""
echo "*************************** File in Process ***************************************"
echo "Working file is $file"

# Save chunked VCF filename into a temporal variable
f=$file

# Read the dosage VCF file, line by line, and extract the dosage field information for each SNP 
# at each individual. Write the data into a SNP by individual matrix:
echo "Reading dosage VCF data..."

# If chunked VCF have header, add "NR>10"
# awk -F: -v f=5 'NR>10{{for(i=f ; i<=NF ; i+=2) printf("%s\t",$i)} print""}' $file > $file.snps_x_individuals
# If chunked VCF do not have headers, start to grab lines from the very first line
# If necessary, you may replace tabs by spaces in the output
awk -F: -v k=5 '{{for(i=k ; i<=NF ; i+=2) printf("%s\t",$i)} print""}' $file > $file.snps_x_individuals

# Transpose the SNP by individual matrix into a individual by SNP matrix ('individuals_x_snps' file)
echo "Transposing dosage VCF data..."

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>x { x = NF }
END {    
    for(j=1; j<=x; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' $file.snps_x_individuals > $file.individuals_x_snps

# Remove intermediate files
rm -f $f
rm -f $file.snps_x_individuals

echo ""

# END of loop for chunked-files
# ------------------------------------------------------------------------------------------------------------------ #
# Warning: before "done" there must be a ";" to close the for-loop
echo "" ;
done


# ################################################################################# #
# 6. Paste all intermediate files to get a final Genomic Predictor file for probABEL #
# ################################################################################# #

# Paste the matrix of index-individual-"MLDOSE" and dosage matrix to prepare the probABEL 'genomic_predictor' file.
echo "Pasting individuals and dose matrices into a genomic prediction file..."
paste $output/t_individuals $output/variants_*.individuals_x_snps > $output/probabel.$gp_filename.gp

echo ""


# ################################################################################# #
# 7. Paste all intermediate files to get a final Genomic Predictor file for probABEL #
# ################################################################################# #

echo "================================ Removal of temporal files ================================"

# Remove intermediate files
echo "Removing intermediate files..."
echo ""
rm -f $output/individuals $output/t_individuals $output/header $output/variants_*.individuals_x_snps $output/$file.individuals_x_snps 

# Delete the big files
rm -f $output/variants


# ######################################################################################## #
# 8. Auxiliar file operations: replace missing data with "NA" and replace tabs with spaces #
# ######################################################################################## #

echo "================================ Replace tabs with spaces. Ending of prepABEL ================================"

# Replace missing data with 'NA' within phenotype file. If necessary, uncomment next line
#awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "NA" }; 1' $phenotypes/phenotypes_revised.csv > $phenotypes/phenotypes_revised_noNA.csv

# Replace tabs with spaces, if needed
echo "Replacing tabs by spaces in ProABEL prepared files..."

awk '{$1=$1}1' OFS=" " $output/probabel.$gp_filename.info > $output/probabel.$gp_filename.spaces.info

awk '{$1=$1}1' OFS=" " $output/probabel.$gp_filename.gp > $output/probabel.$gp_filename.spaces.gp

awk '{$1=$1}1' OFS=" " $output/probabel.$gp_filename.map > $output/probabel.$gp_filename.spaces.map


# Discard files with tabs and keep files with spaces as delimiters
rm -f $output/probabel.$gp_filename.info $output/probabel.$gp_filename.gp $output/probabel.$gp_filename.map

# Rename ProbABEL files
mv $output/probabel.$gp_filename.spaces.info $output/probabel.$gp_filename.info
mv $output/probabel.$gp_filename.spaces.gp $output/probabel.$gp_filename.gp
mv $output/probabel.$gp_filename.spaces.map $output/probabel.$gp_filename.map

echo "================================ End of Script ================================"


# ########################################### #
# 9. Perform Survival Analysis using ProbABEL #
# ########################################### #

# Run Cox Proportional Relative Hazard model using variants dosage data, SNP information, MAP information, 
# and phenotypes, using all covariates declared (see ProbABEL manual at http://www.genabel.org/manuals/ProbABEL) 
pacoxph -p $phenotypes/probabel.phenotypes \
-d $output/probabel.$gp_filename.gp \
-i $output/probabel.$gp_filename.info \
-m $output/probabel.$gp_filename.map \
--allcov \
-o $coxph/probabel.$gp_filename.coxph

# End of script
'''
