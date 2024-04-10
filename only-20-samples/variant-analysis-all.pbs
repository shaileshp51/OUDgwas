#!/bin/bash
#PBS -N opi-gen
#PBS -l select=1:ncpus=2:mem=64gb:mpiprocs=1,walltime=72:00:00
# module purge
# module add biocontainers
# module add bcftools/1.17

#cd ${PBS_O_WORKDIR}
source ~/.bashrc
source ~/.bash_aliases
source activate /opt/skp/anaconda3/envs/py311_genomics

GENE_EXTENSION=5000
PROJROOT="/project/ealexov/opioid" # for  Palmetto
PROJROOT="/home/skp/ClemsonWork/Projects/running/Opioid" # for MyWorkStation
PROJDIR="${PROJROOT}/spanday" # for Palmetto
PROJDIR="/home/skp/ClemsonWork/Projects/running/Opioid" # for MyWorkStation
SUBSET="only-20-samples"
export APPTAINER_BINDPATH="$PROJROOT"


<<'STEP1'
# Filter the variants for only first 20 samples, as the last 
# four samples were added to minimum sample count for sequenceing

bcftools view -S selected-sample.txt $PROJDIR/combined.jc.vcf.gz > $PROJDIR/$SUBSET/combined.jc.vcf

# compress the selected-samples vcf 
bcftools view $PROJDIR/$SUBSET/combined.jc.vcf -Oz -o $PROJDIR/$SUBSET/combined.jc.vcf.gz
bcftools index $PROJDIR/$SUBSET/combined.jc.vcf.gz
STEP1

#<<'STEP2'
awk -v "gext=$GENE_EXTENSION" '{printf "%s\t%d\t%d\n", $2, $3-gext, $4+gext}' genes-selected.txt  > genes-selected.bed
# Step 3. extract the variants related to the genes-selected (Genes of Interest)
bcftools view $PROJDIR/$SUBSET/combined.jc.vcf.gz -q 0.03:minor -R $PROJDIR/$SUBSET/genes-selected.bed -o $PROJDIR/$SUBSET/combined.jc.gs.vcf
#STEP2

#<<'STEP4'
#step 4. bzip and index the .vcf containing the genes selected
bcftools view $PROJDIR/$SUBSET/combined.jc.gs.vcf -q 0.03:minor -Oz -o $PROJDIR/$SUBSET/combined.jc.gs.vcf.gz
bcftools index $PROJDIR/$SUBSET/combined.jc.gs.vcf.gz
#STEP4

#<<'STEP5'
# step 5. Predict the effects of variants using Ensemble Variant Effect Predictor
singularity exec $PROJDIR/tools/vep.sif vep --fork 16 --dir $PROJDIR/tools/vep_data \
    --max_af -i $PROJDIR/$SUBSET/combined.jc.gs.vcf --format vcf \
    --force_overwrite -o $PROJDIR/$SUBSET/combined.jc.gs.vep.vcf --vcf --cache
#STEP5

# Lets filter only most-severe-consequence-variants, see consequence-list to filter in filter-variants.txt
conseq="_mscvo" 
bcftools view --header $PROJDIR/$SUBSET/combined.jc.gs.vep.vcf > $PROJDIR/$SUBSET/combined.jc.gs.vep${conseq}.vcf
egrep -f $PROJDIR/$SUBSET/filter-variants.txt $PROJDIR/$SUBSET/combined.jc.gs.vep.vcf >> $PROJDIR/$SUBSET/combined.jc.gs.vep${conseq}.vcf
# Lets check the correctness of the filtered vcf file
bcftools view $PROJDIR/$SUBSET/combined.jc.gs.vep${conseq}.vcf -Oz $PROJDIR/$SUBSET/combined.jc.gs.vep${conseq}.vcf.gz
#<<'STEP6'
# step 6. convert .vcf to PLINK format
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --noweb --const-fid 0 --maf 0.03 --keep-allele-order \
    --vcf $PROJDIR/$SUBSET/combined.jc.gs.vep.vcf \
    --recode --out $PROJDIR/$SUBSET/combined.jc.gs.vep.plink
#STEP6

# From Vijay
# Quality Control criteria
params_GENO=0.1 # Missingness per SNP
params_MIND=0.1 # Missingness per indidividual
params_MAF=0.03 # Minor allele frequency
params_HWE=0.0000001 # Hardy-Weinberg threshold
input_FILT_VCF="combined.jc.gs.vep${conseq}.vcf"
input_PHENO_FILT="sample-phenotypes-plink.txt"
params_BFILE="combined.jc.gs.vep${conseq}.qc"
# 
$PROJDIR/tools/plink_linux_x86_64_20231211/plink \
    --keep-allele-order \
    --allow-no-sex \
    --nonfounders \
    --make-bed \
    --double-id \
    --geno ${params_GENO} \
    --pheno ${input_PHENO_FILT} \
    --maf ${params_MAF} \
    --vcf $PROJDIR/$SUBSET/${input_FILT_VCF} \
    --out $PROJDIR/$SUBSET/${params_BFILE}

# Options excluded for now
#    --mind ${params_MIND} \
#    --hwe ${params_HWE} \

$PROJDIR/tools/plink_linux_x86_64_20231211/plink --bfile ${params_BFILE} --recode --out ${params_BFILE}


# Run association with plink
#plink case-control search in manual
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --assoc \
    --out $PROJDIR/$SUBSET/${params_BFILE}_assoc
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --assoc \
    --ci 0.95 --out $PROJDIR/$SUBSET/${params_BFILE}_assoc-ci.95
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --assoc \
    --counts --out $PROJDIR/$SUBSET/${params_BFILE}_assoc-count                               
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --assoc \
    --adjust --out $PROJDIR/$SUBSET/${params_BFILE}_assoc-adj                                
$PROJDIR/tools/plink_linux_x86_64_20231211/plink --file $PROJDIR/$SUBSET/${params_BFILE} --allow-no-sex --logistic \
    --out $PROJDIR/$SUBSET/${params_BFILE}_assoc-logistic


# Run association using GEMMA
$PROJDIR/tools/gemma-0.98.5-linux-static-AMD64 \
    -bfile $PROJDIR/$SUBSET/${params_BFILE} \
    -gk 1 \
    -o RelMat

$PROJDIR/tools/gemma-0.98.5-linux-static-AMD64 \
    -bfile $PROJDIR/$SUBSET/${params_BFILE} \
    -k output/RelMat.cXX.txt \
    -lmm 2 \
    -o GWASresults.lmm

mv $PROJDIR/$SUBSET/output $PROJDIR/$SUBSET/output${conseq}
