# GWAS
This repository is dedicated to host the in-house scripts developed and used for GWAS.


# Abbreviation
| ABVR     | Detail                                  |
| -------- | --------------------------------------- |
| js       | Joint calling                           |
| gs       | Genes selected only                     |
| vep      | Variant effect predictor output         |
| qc       | Quality control                         |
| mscvo    | Most severe consequence variant only    |
| OUD      | Opioid Use Disorder                     |

# Let's install singularity to enable use of tools available as sigularity-image 
## We used Ubuntu 24.04 (LTS)
First of all I added the singularity repository using
```
wget -O- http://neuro.debian.net/lists/noble.us-tn.libre | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkps://keyserver.ubuntu.com 0xA5D32F012649A5A9
```
follow the instructions
```
sudo apt update
```
after completion use
```
sudo apt install singularity
```
It installed singularity version 4.1.1, you can check your's using
singularity version


# Required Tools
* PLINK
* VEP (using an .sif)
I followed the instructions provided at link below
https://useast.ensembl.org/info/docs/tools/vep/script/vep_download.html#singularity
to get the vep.sif, install and then setup the cache for homo_sapiens reference genome assembly GRCh38
We used $PROJECTDIR/tools for vep_data instead of $HOME provided in the instrctions. Where PROJECTDIR is the 
main-folder of this project.

* samtools (bcftools) v1.20

# Large-files from groups website
Large files including variant calling file (.vcf) and raw STRING databse PPIs file and completely normalized PPIs file, and also singularity images of the vep and samtools are made available. These files can not be uploaded to GitHub rather made available from groups webpage.
Download large file from http://compbio.clemson.edu/media/download/large-files.zip
Unzip and copy the folders inside unzipped folder to the folders with same name in OUDgwas folder

# Workflow
![Workflow](./image/oud-gwas-workflow.png)

# Case-control association for variants with coding consequence on reported 12 genes.
* The 12 genes that are already reported in literature to harbor variants associated to OUD are listed in `data/reported-12-genes.txt` and corresponding variants filtered from the all variants are saved in folder `only-reported-12-gene` and all related association analysis results are contained in it.

Covered steps: 1.1 to 1.4.1.

To run this analysis use:
```
bash only-reported-12-genes-mscvo-analysis.bash
```
We did not found any variant with significance from this analaysis.

# Case-control association for all variants on reported 12 genes (including intronic).
Covers steps &amp; 1.5.1 shown in workflow and all results are stored in `only-reported-12-gene`.

To run this analysis use:
```
bash only-reported-12-genes-all-varinats-analysis.bash
```
Based on raw p-value < 0.05 we found some varinats, though upon multiple hypothesis testing none of the variants reached significance. Which is due to very small sample size used here. We considered variants only based on raw p-value in this step.
The final varinats list after association alalysis is also kept in `association-results` folder.

# Variant effect prediction of all variants, filtering of *mscvo* variants, and case-control association
All the related files will be created in the folder `genome-wide-association`.
Covers steps: 2.1 to 2.6
Final association results is also kept in `association-results` folder for easy access and comparison to what one gets after running the script.

To run this analysis use:
```
bash genome-wide-mscvo-analysis.bash
```

running this script completely may take several hours on normal WorkStations.

# Co-occurance analysis
Covers steps: 2.7 &amp; 2.7.1
The co-occurance analysis is performed using the python code provided in the notebook
install all the dependencies in a dedicated conda-environment before running the code.

Use jupyter-lab to see and run the code contained in the notebook `GWAS-py-v01.ipynb`
Final association results is also kept in `association-results` folder.


# PPIs network analysis
Covers steps 3.1 to 3.4
Preprocessing of the PPIs in dbSTRING, normalization and maping gene_ids to gene_names is also done.
Use jupyter-lab to see and run the code contained in the notebook `GWAS-py-v01.ipynb`

# Network of Genes-of-Interest PPIs construction and analysis.
The PPIs network is constructed and all the relevant data and Cytoscape session is saved in folder `PPI-net`.

# Enrichment query and results
Function enrichment query and results are proved in the folder `enrichemnt-results`.