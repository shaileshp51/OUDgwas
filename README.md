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

# Let's install singularity to enable use of tools available as sigularity-image 
## We used Ubuntu 24.04 (LTS)
First of all I added the singularity repository using
wget -O- http://neuro.debian.net/lists/noble.us-tn.libre | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkps://keyserver.ubuntu.com 0xA5D32F012649A5A9

follow the instructions
sudo apt update

after completion use
sudo apt install singularity

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

* samtools (bcftools)

# Workflow
```
mermaid
graph TB
A1[all genomic variants] -->|step 1.1: filter<br>in reported genes| B(step 1.2: variant effect prediction)
    B --> C{step 1.3: Coding consequence}
    C -->|Yes| D[step 1.4: case-control association] 
    C -->|No| E[step 1.5: case-control association]
    D --> F{step 1.4.1<br>check <br>significance<br>raw p-value < 0.05} -->|YES|H(None)
    E --> G{step 1.5.1<br>check <br>significance<br>raw p-value < 0.01} -->|YES|I(13 intronic variants)
A1 -->B1(step 2.1: variant effect prediction)
B1 -->C1{step 2.2: Coding consequence} --> |YES| D1[step 2.3: case-control association] 
D1 -->|step 2.4: filter <br>alternate-allele <br>in cases only| F1{step 2.5: check <br>significance<br>raw p-value < 0.05} -->|YES|H3{step 2.6: raw p-value < 0.005} --> |YES|H1(varinats<br>on<br>11 genes)
F1 --> |YES|G1(step 2.7: Co-occurance) --> |"step 2.7.1: filter<br>alter-alleles-count<br>sum<br>in cases only &ge; 10"|I1(multi-variant genes)
H1 --> J1(Genes-of-Interest)
I1 --> J1
J1 --> F3

A3(PPIs-STRING-hs) --> B3(step 3.1: preprocess) 
B3 --> C3(step 3.2: normalize<br>combined_score)
C3 --> D3(step 3.3: plot distribution<br>combined_score<br>normalized_score<br>analyze)
D3 --> E3("step 3.4: filter<br>high-confidence PPIs:<br>combined_scode &gt; 800<br> or <br>normalized_score &gt; 0.8")
E3 --> F3(setp 3.5: select GOIs for subnetwork<br>with Cytoscape)
F3 --> G3(step 3.6: extend subnet: <br>also include <br>immideate neighbors of GOIs<br> with Cytoscape)
G3 --> I3(step 3.7: further extend subnet:<br> also include nodes with<br> two or more neighbors in selected nodes<br>with Cytoscape)
I3 --> J3(step 3.8: cluster subnet into <br> connected components <br> with Cytoscape)
```
J3 --> K3(setp 3.9: select genes <br>in largest compnent:<br>say extended-GOIs<br> with Cytoscape)
K3 --> L3(step 3.10: gProfile: functional enrichment analysis)
L3 --> M3(step 3.11: report hightly enriched<br> terms relevant to OUD) 
```
