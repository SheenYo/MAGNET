#!/bin/bash 

LocSource=$(pwd) 

source /$LocSource/ConfigFiles/Tools.config
source /$LocSource/ConfigFiles/Thresholds.config

cd OUTPUT_DIR/Stage3_GWAS

echo -e "
#######################################
# Stage 3: GWAS has started now       #
#######################################
"

#cd .. /Stage3_GWAS
chmod 770 ../Stage2_GenoImpute/Data_SNPfile*.raw

export pheno
export phenofile
export covars
export fixed

python <<END
import glob
import os
import optparse

parser = optparse.OptionParser()
parser.add_option("-f", "--factor", action="store", dest="fact", help="Factor X")
(options, args) = parser.parse_args()

inputfiles = glob.glob("../Stage2_GenoImpute/Data_SNPfile*.raw")

print len(inputfiles)

pheno2=os.environ["pheno"]
phenofile=os.environ["phenofile"]
covars=os.environ["covars"]
fixed=os.environ["fixed"]

i = 1
for inputfile in inputfiles:
    jobname = "job_"+str(i)+"_factor.sh"
    f = open(jobname,"w")
    f.write("R --no-save --args ../Stage2_GenoImpute/Data_SNPfile"+str(i)+".raw "+str(phenofile)+" "+str(pheno2)+" "+str(covars)+" "+str(fixed)+" ResultsMerged"+str(pheno2)+"_"+str(i)+".txt<  ../../Scripts/Code5_Regression.R")
    f.close()


    i = i+1
END

chmod 770 job_*.sh

njobs=$(ls -l ../Stage2_GenoImpute/Data_SNPfile*.raw|wc -l)


for j in $(seq 12 12 $njobs); do
    lower=$(expr "$j" - 11)
    upper=$j    
    
    for i in $(seq $lower $upper); do 
    ./job_"$i"_factor.sh > prg"$i".out  2>&1 &
	sleep 1
    done;
    wait
done;

wait
# final run

remainder=$(expr $njobs % 12)
lower=$(expr $njobs - $remainder)
for i in $(seq $lower $njobs); do
    ./job_"$i"_factor.sh > prg"$i".out  2>&1 &
    sleep 1
done;
wait

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e "
#######################################################################################
# 				Post regression 				      #
#######################################################################################
"

cat ResultsResults*> tmp.txt

#header=$(head -1 tmp.txt)

sed 's/\"//g' tmp.txt>tmp2.txt

sed '/snp/d' tmp2.txt>tmp3.txt

sed 's/X//g' tmp3.txt>tmp4.txt

echo -e "SNP\tbeta\tse\tt\tP\tadj.pvalues" | cat - tmp4.txt > tmp5.txt

awk {'print $1 "\t" $5'} tmp5.txt>SNP_Pvalue$pheno

mv tmp5.txt RegressionOutput$pheno.txt

ImputedData=../Stage2_GenoImpute/Merged_FinalQC_SNPs_Data

# To find which SNPs belong to which Chromosme, we merge the Regression output file with the final bim file
# Merged_FinalQC_SNPs_DE.bim and statistical summary of Regression summary

sed 's/:/./' $ImputedData.bim>File2 # Since SNPs which are represented as chr:bp are write as chr.bp so we substitute ":" with "."

sed 's/ /\t/g' File2>File3


WD=$(pwd)

awk 'FNR==NR{a[$1]=$0;next}{if ($2 in a){p=$2;$2="";print a[p] "\t" $0} else {print p=$2 "\t" "NA" "\t" "NA" "\t" "NA" "\t" "NA" "\t" "NA" "\t" "NA" "\t" $1 "\t" $3 "\t" $4 "\t" $5 "\t" $6}}' RegressionOutput$pheno.txt File3>File4$pheno

echo -e "SNP\tbeta\tse\tt\tP\tadj.pvalues\tCHR\t\CHR2\tGP\tBP\tA1\tA2" | cat - File4$pheno > Complete_SNP_info2$pheno

sed 's/ /\t/g' Complete_SNP_info2$pheno>Complete_SNP_info$pheno

echo -e "SNP\tbeta\tse\tt\tP\tadj.pvalues\tCHR\tCHR2\tGP\tBP\tA1\tA2" | cat - Complete_SNP_info$pheno>Complete_SNP_info_2$pheno

sed 's/ /\t/g' Complete_SNP_info2$pheno>Complete_SNP_info_$pheno


ColManhattan="#c51b8a"
pheno=$pheno
R --no-save --slave --args $pheno $ColManhattan $snpsOfInterest < ../../../Scripts/Manhattanplot.R 



######################################################################################
# 14) Gene enrichment analysis(MAGMA)						     #
######################################################################################


#For MAGMA annotation file is needed which can be generated using the reference genotype file

#For a windown size of 5
pheno=$pheno
#1-Generating Annotation file
magma --annotate window=$windowSize \
      --snp-loc $MagmaSNPloc \
      --gene-loc $MagmaGeneloc \
      --out $pheno.annot
#Time required elapsed: 00:00:57

#2-Gene set enrichment was performed in bataches chromosomewise and submitted to server


python <<END
import glob
import os
import optparse

parser = optparse.OptionParser()
parser.add_option("-f", "--factor", action="store", dest="fact", help="Factor X")
(options, args) = parser.parse_args()

inputfiles = glob.glob("../Stage2_GenoImpute/Merged_FinalQC_SNPs_Data*.fam")

print len(inputfiles)



for i in range(1,23):
    jobname ="../../Scripts/MAGMAChr"+str(i)+".sh"	
    f = open(jobname,"w")
    f.write("magma --bfile $MagmaRef \
        --pval SNP_Pvalue$pheno N=$MagmaN \
        --batch "+str(i)+" chr \
        --gene-settings adap-permp \
        --gene-annot $pheno.annot.genes.annot \
        --out Magma_$pheno")
    f.close()

with open("MAGMAChr"+str(i)+".sh", "a") as myfile:
    myfile.write("\nexit")
    myfile.close()
END


chmod 755 ../../Scripts/MAGMA*.sh


for((i=1; i<=22;i++)); do

echo "start MAGMA Chromosome $i"

     ./../../Scripts/MAGMAChr$i.sh &
     sleep 1
done;

echo "all processes started"

wait

echo "all processes ended"
    
rm MAGMAChr*.sh

chmod 770 *
echo  "3-Merge All the batches"

magma --merge Magma_$pheno --out Magma_$pheno

#4-Make a tab separated file out of the wide space file of MagmaOut.genes.out

sed -e 's/  */\t/g'  Magma_$pheno.genes.out > Magma_$pheno.genes_tabseparated.txt

awk -v Pval=$MagmaPERMP '$10<=Pval {print $0}' Magma_$pheno.genes_tabseparated.txt > tmp

head -1 Magma_$pheno.genes_tabseparated.txt |sed 's/  */\t/g'| cat - tmp > $pheno.SignificantGenes.txt ## attach header

Nsig=($(wc -l $pheno.SignificantGenes.txt))

echo -e  "overal $Nsig genes passed permutation p-value treshold of $MagmaPERMP\n"
echo -e  "File saved to $pheno.Significant.txt"

cat $pheno.SignificantGenes.txt



