#!/bin/bash


shopt -s expand_aliases
source ~/.bashrc
source ~/.bash_aliases
#source activate /opt/skp/anaconda3/envs/py311_genomics


# How to query SNP in dabase
# ((14[Chromosome]) AND 79404467[Base Position]) 

GENE_EXTENSION=5000

# Change the directory to the root-dir of the project
PROJDIR=`pwd` # for MyWorkStation

export APPTAINER_BINDPATH="$PROJDIR"



bcftools="singularity exec $PROJDIR/tools/samtools-1.20.sif bcftools"

SUBSET="only-reported-12-genes"
#<<'STEP1'
# Create a bed file listing only 12-genes that are reported in literature to have association with OUD
# with each genes bounday is extended by $GENE_EXTENSION bases
awk -v "gext=$GENE_EXTENSION" '{printf "%s\t%d\t%d\n", $2, $3-gext, $4+gext}' $PROJDIR/data/reported-12-genes.txt > $PROJDIR/$SUBSET/reported-12-genes.bed
# extract the variants related to the reoprted-12-genes 
$bcftools view $PROJDIR/data/combined.jc.vcf.gz -q 0.05:minor -R $PROJDIR/$SUBSET/reported-12-genes.bed -o $PROJDIR/$SUBSET/combined.jc.gs.vcf
#STEP1

#<<'STEP2'
# bzip and index the .vcf containing the genes selected
$bcftools view $PROJDIR/$SUBSET/combined.jc.gs.vcf -q 0.05:minor -Oz -o $PROJDIR/$SUBSET/combined.jc.gs.vcf.gz
$bcftools index $PROJDIR/$SUBSET/combined.jc.gs.vcf.gz
#STEP2

vcfprefix="combined.jc.gs"
# <<'STEP5'
# step 5. Predict the effects of variants using Ensemble Variant Effect Predictor
singularity exec $PROJDIR/tools/vep.sif vep --fork 16 --dir $PROJDIR/tools/vep_data \
    --max_af -i $PROJDIR/$SUBSET/${vcfprefix}.vcf --format vcf \
    --force_overwrite -o $PROJDIR/$SUBSET/${vcfprefix}.vep.vcf --vcf --cache
# STEP5


# Lets filter only most-severe-consequence-variants, see consequence-list to filter in filter-variants.txt
conseq="_mscvo" #"_mscvo"

#<<'STEP6'
$bcftools view --header $PROJDIR/$SUBSET/${vcfprefix}.vep.vcf > $PROJDIR/$SUBSET/header.txt
egrep -f $PROJDIR/data/filter-most-severe-consequence-variants.txt $PROJDIR/$SUBSET/${vcfprefix}.vep.vcf > $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}.vcf

awk '{if($1=="#[1]CHROM"||$1=="chrX"||$1=="chrY"||$1=="chr1"||$1=="chr2"||$1=="chr3"||$1=="chr4"||$1=="chr5"||$1=="chr6"||$1=="chr7"||$1=="chr8"||$1=="chr9"||$1=="chr10"||$1=="chr11"||$1=="chr12"||$1=="chr13"||$1=="chr14"||$1=="chr15"||$1=="chr16"||$1=="chr17"||$1=="chr18"||$1=="chr19"||$1=="chr20"||$1=="chr21"||$1=="chr22"){print $0}}' $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}.vcf > $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}1.vcf
cat $PROJDIR/$SUBSET/header.txt $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}1.vcf > $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}2.vcf

$bcftools reheader $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}2.vcf --header $PROJDIR/$SUBSET/header.txt > $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}3.vcf

# Lets check the correctness of the filtered vcf file
$bcftools view $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}3.vcf -Oz -o $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}3.vcf.gz

mv $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}3.vcf $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}.vcf
mv $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}3.vcf.gz $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}.vcf.gz

rm $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}2.vcf

$bcftools  +split-vep $PROJDIR/$SUBSET/${vcfprefix}.vep${conseq}.vcf -H -s worst -f '%CHROM\t%POS\t%ID[\t%Consequence][\t%IMPACT][\t%Gene][\t%GT]\n'  > $PROJDIR/$SUBSET/annot${conseq}.txt
awk '{if($1=="#[1]CHROM"||$1=="chrX"||$1=="chrY"||$1=="chr1"||$1=="chr2"||$1=="chr3"||$1=="chr4"||$1=="chr5"||$1=="chr6"||$1=="chr7"||$1=="chr8"||$1=="chr9"||$1=="chr10"||$1=="chr11"||$1=="chr12"||$1=="chr13"||$1=="chr14"||$1=="chr15"||$1=="chr16"||$1=="chr17"||$1=="chr18"||$1=="chr19"||$1=="chr20"||$1=="chr21"||$1=="chr22"){print $0}}' $PROJDIR/$SUBSET/annot${conseq}.txt > $PROJDIR/$SUBSET/annot${conseq}2.txt
mv $PROJDIR/$SUBSET/annot${conseq}2.txt $PROJDIR/$SUBSET/annot_${vcfprefix}${conseq}.txt
#STEP6

#<<'COMMENT'
# Quality Control criteria
params_GENO=0.1 # Missingness per SNP
params_MIND=0.1 # Missingness per indidividual
params_MAF=0.05 # Minor allele frequency
params_HWE=0.0000001 # Hardy-Weinberg threshold
input_FILT_VCF="${vcfprefix}.vep${conseq}.vcf"
input_PHENO_FILT="sample-phenotypes-plink.txt"
params_BFILE="${vcfprefix}.vep${conseq}.qc"
# 
$PROJDIR/tools/plink_linux_x86_64_20231211/plink \
    --keep-allele-order \
    --allow-no-sex \
    --nonfounders \
    --make-bed \
    --double-id \
    --geno ${params_GENO} \
    --pheno $PROJDIR/data/${input_PHENO_FILT} \
    --maf ${params_MAF} \
    --vcf $PROJDIR/$SUBSET/${input_FILT_VCF} \
    --out $PROJDIR/$SUBSET/${params_BFILE} \
    --hwe ${params_HWE} 

$PROJDIR/tools/plink_linux_x86_64_20231211/plink --bfile $PROJDIR/$SUBSET/${params_BFILE} --recode --out $PROJDIR/$SUBSET/${params_BFILE}


# Run association with plink
#plink case-control 
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --assoc \
    --ci 0.95 --out $PROJDIR/$SUBSET/${params_BFILE}_assoc
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --assoc \
    counts --ci 0.95 --out $PROJDIR/$SUBSET/${params_BFILE}_assoc-count                               
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --assoc \
    --adjust --ci 0.95 --out $PROJDIR/$SUBSET/${params_BFILE}_assoc-adj                                
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --logistic \
    --ci 0.95 --out $PROJDIR/$SUBSET/${params_BFILE}_assoc-logistic


#COMMENT
