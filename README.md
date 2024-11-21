![alt text](https://github.com/Singh-Lab/Dyscovr/blob/main/Dyscovr_Logo.png "Logo")
# Dyscovr (dyscovr.princeton.edu)
#### Mona Singh Lab
#### Lewis-Sigler Institute for Integrative Genomics, Princeton University

Uncover links between mutated driver genes and dysregulation of putative target genes across the genome, within and across 19 cancer types.

# Publication Information
Please cite as: Geraghty, S.E., Boyer, J.A., Fazel-Zarandi, M. Arzouni, N., Ryseck, R., McBride, M.J., Parsons, L.R., Rabinowitz, J.D., and Singh M. (2024). Integrative Computational Framework, Dyscovr, Links Mutated Driver Genes to Expression Dysregulation Across 19 Cancer Types. bioRxiv. https://doi.org/10.1101/2024.11.20.624509

# Source Files
* Input data files for the TCGA were downloaded from the GDC Data Portal Repository: https://portal.gdc.cancer.gov/repository. See Suppl. Table 6 for information about file downloads.
* Input data files for METABRIC (https://doi.org/10.1038/nature10983) were downloaded from cBioPortal (https://www.cbioportal.org/study/summary?id=brca_metabric).
* After files are downloaded, they are run through the various preprocessing pipelines found in the data_processing folder. These steps vary by data type, but generally include filtering and normalization procedures.

# Preprocessing
1. After initial preprocessing, files are processed according to the premodel_processing.R file. This file merges patient- and sample-level clinical data frames and ensures that all data files are subsetted to only those patient samples that have all data types of interest (e.g. mutation, CNA, expression, methylation, and clinical data).
2. The functions in generate_protein_input_tables.R are run in order to determine which driver genes are mutated at a sufficient frequency (e.g. 5%) and to create a table(s) with those drivers for use in the following step. To create your own driver gene input file, simply create a 2-column file with "swissprot_ids" and "ensg_ids", with each row containing the Uniprot or Swissprot and Ensembl IDs for each driver gene of interest.
3. Create a file with a list of target genes of interest (either all putative genes, or a subset; the smaller the list of genes, the faster the software will run). This file should take the format of a 2-column data from with column names "ensg" and "swissprot", with each row containing the Uniprot/Swissprot and Ensembl IDs for each target gene of interest.

# Running Dyscovr
## Install Dependencies
Install dependencies using install.packages() and BiocManager::install(). See STAR Methods Software and Algorithms table for list of dependencies.
## Clone Dyscovr Git Repository 
Clone the repository into a location with ample storage for input and results files.
```
git clone https://github.com/scamilli97/Dyscovr.git
```
## Create Directory Hierarchy for Input and Output Files
In repository, create directories ```input_files```, ```output_files```, and ```output_visualizations``` (e.g. ```mkdir input_files```).

### Input Files
1. Within ```input_files```, create a separate directory for each cancer type of interest, e.g. ```BRCA```, ```PanCancer```, etc. Within each of these cancer type directories, create a separate directory to store each data type, including ```Clinical```, ```CNA```, ```Expression```, ```Methylation```, ```Mutation``` and ```Sample``` (they should have these names and are case-sensitive). Store your preprocessed data files in these locations. Also within ```input_files/CANCER_TYPE```, create a directory called ```target_lists```, where you store the target list file you created in Preprocessing part 3. Store your file of input drivers directly in the ```input_files``` directory.
2. Within input_files, create another directory called ```lm_input_tables```. This is where Dyscovr will store the model input tables it creates from your preprocessed data. Within this directory, create another named by the run you wish to perform (e.g. nonsynonymous, silent, etc.). You will be able to specify this directory name later as an argument.
3. If running per-cancer, e.g. across all cancer types individually, within the PanCancer directory create a separate folder for each cancer type (e.g. ```ACC```, ```BLCA```, etc.). Each of these folders should have its own ```lm_input_tables``` directory containing at least one run folder (see step 2).

### Output Files
Within ```output_files```, create a matching directory with cancer type name, e.g. ```PanCancer```, containing a directory with your run name, e.g. ```nonsynonymous```. This is where your overall output files from Dyscovr will be written. Again, if running per-cancer, you need a nested directory as such: ```PanCancer/ACC/nonsynonymous/```, with a separate directory for each cancer type name you are investigating.

### Visualization (Optional)
Dyscovr defaults to create histograms of the p-value distribution and coefficient distribution for a given run. Within ```output_visualizations```, create a matching directory with cancer type name, e.g. ```PanCancer```, containing a directory with your run name, e.g. ```nonsynonymous```.

## Create Dyscovr Input Tables
Create input tables for the Dyscovr framework. This is performed in a separate step in order to shorten runtime when multiple variations of Dyscovr (e.g. with varying arguments) are run on the same input data. This step relies on helper functions found in dyscovr_helper_functions.R and creates two types of input files, stored separately to conserve memory: A. Features independent of the target gene: a table with features that are independent of the given target gene (including clinical features and driver-specific features). Only one of these files is created per cancer type. B. Features dependent on the target gene: specifically, target gene expression, mutation status, CNA status, and methylation status. One of these files is created per target gene in consideration. 
### Adjust Arguments
Dyscovr is run using SLURM job allocation software on a HPC. View the ```create_dyscovr_input_tables.sh``` file to adjust arguments as desired. Detailed descriptions of each argument with defaults are found at the top of ```create_dyscovr_input_tables.R```. Default memory allocation is 64GB, and allocating less may cause an out-of-memory error. CPUs-per-task is set to 8 for purposes of multithreading and speed improvements.
### Run via the Command Line
Once the batch file is adjusted appropriately, ```cd``` to the main Dyscovr directory and call ```sbatch create_dyscovr_input_tables.sh``` from the command line.

## Run the Dyscovr Model
The dyscovr.R file runs an individual regression per target gene, and for each target gene, binds together the two types of input files (the target-independent feature file and the target-specific feature file, see 'Input Files') to use as input to regression. This file outputs merged data tables, consisting of regression output from all target genes together. This included uncorrected output (no ID conversion, no multiple hypothesis testing correction, and all feature results included) as well as corrected output for only driver mutation (_MUT) or CNA (_CNA) features, with driver-target pairings ordered by q-value. Variables excluded during multicollinearity checks are output to a separate file.
### Adjust Arguments
Dyscovr is run using SLURM job allocation software on a HPC. View the ```run_dyscovr.sh``` file to adjust arguments as desired. Detailed descriptions of each argument with defaults are found at the top of ```dyscovr.R```. Default memory allocation is 64GB, and allocating less may cause an out-of-memory error. CPUs-per-task is set to 8 for purposes of multithreading and speed improvements.
### Run via the Command Line
Once the batch file is adjusted appropriately, ```cd``` to the main Dyscovr directory and call ```sbatch run_dyscovr.sh``` from the command line. Outfiles will be found in the directories you created in 'Create Directory Hierarchies...' and specified in the arguments of ```run_dyscovr.sh```.

NOTE: When running on all genes in the human genome as targets, run on multiple, smaller input files for memory purposes (allgene_targets_pt1.csv, allgene_targets_pt2.csv, etc.). Post-processing is required to bind these results files back together and perform multiple hypothesis testing correction on the full set of tests; see ```recombine_multipart_output.R```.
