#!/bin/bash 

LocSource=$(pwd)

source /$LocSource/ConfigFiles/Tools.config
source /$LocSource/ConfigFiles/Thresholds.config

cd OUTPUT_DIR/Stage2_GenoImpute

#Check if the user performed MAGNET Stage1 or has provided a quality checked file

if [ $PerformStage1=="No" ]
then
	echo -e "User provided quality checked file will be used \n"
	QCFile=$QualityCheckedFile
else
	QCFile=../Stage1_GenoQC/FinalQC_Study #variable assigned to the final QC file from Stage1_GenoQC.sh
fi


echo -e "
#####################################################
# 1)Check if only autosomal chromosomes are present #
#####################################################
"

echo -e "Only autosomal chromosomes will be selected from the final quality check file provided from Stage1_GenoQC.sh \n "

if [ -e $QCFile.bim ]
then
	awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' $QCFile.bim > snp.txt
else
	echo -e "Please check the existence of $QCFile \n"

sleep 2
fi

echo -e '
######################################################
'

$PLINK --bfile $QCFile --extract snp.txt --make-bed --out QC_autosomal

echo -e '
######################################################
'

$PLINK --bfile $QCFile --list-duplicate-vars ids-only suppress-first #List the duplicated SNPs 

echo -e '
######################################################
'

$PLINK --bfile $QCFile --exclude plink.dupvar --make-bed --out QC_autosomal_for_StrandCheck #Delete the duplicated SNPs

echo -e '
######################################################
'

awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' QC_autosomal_for_StrandCheck.bim > snp.txt

$PLINK   --bfile QC_autosomal_for_StrandCheck --extract snp.txt --make-bed --out QC_autosomal


echo -e "
#########################################################
# 2)Perform liftOver, if genome build is other than hg19#
#########################################################
"

if [[ $LIFTOVER=="No" ]]
then
	echo -e "Your data is already in hg19 format \n"
	
else

	echo -e							".........Updating the map positions of plink files....... \n"
	awk {'print $4 "\t" $3'} $hg19SNPs >  map.file
	$PLINK --bfile QC_autosomal --update-map map.file --make-bed --out QC_automsomal_Updated

	echo -e							 ".........Map positions are updated.........\n"
	wait
	echo -e							 "........Find all the SNPs whose bp position is 0........ \n"

	awk '$4 == 0 {print}' QC_automsomal_Updated.bim > not_annotated.bim


	echo ".....$(cat not_annotated.bim|wc -l) SNPs had unannotated positions"


	echo -e "Preparing file for liftOver \n"    

	#Fill the middle column with 0 as we do not have the starting position of the snp and need to put "chr" before chromsome number, for preparing a bed file for liftover

	awk {'print $1 "\t" $4 "\t" $2'} QC_automsomal_Updated.bim |awk '$1 = $1 FS "0"'|awk '$1="chr"$1'> liftOverFile_Data.bed


	echo -e	"Make sure you use the correct build to liftover, please select your required chain file to update the genome build from the config file, by default hg18 is lifted over to hg19 \n"

	$liftOver liftOverFile_Data.bed $ChainToChoose Data_lifted_afterhg.bed Data_unlifted_afterhg.bed
	wait
	echo -e "$(cat Data_lifted_afterhg.bed|wc -l) SNPs have been lifted \n"

	SizeLifted=($(cat Data_lifted_afterhg.bed | wc -l))

	if [ $SizeLifted -gt 0 ]
	then
		mv Data_lifted_afterhg.bed Data_liftedoutput_map #Rename the final lifted file
		cat Data_unlifted_afterhg.bed>Data_unliftedoutput_map #Remove the unlifted SNPs
		sed '/#/d' Data_unlifted_afterhg.bed>Data_unliftedoutput_map1 #Remove comments from the liftOver output of unlifted SNPs
		awk {'print $4'} Data_unliftedoutput_map1>Data_unliftedoutput_map3 #get rs names of unlifted snps; includes snps already in hg19 position 
		awk {'print $4 "\t" $3'} Data_liftedoutput_map>Final_map
		$PLINK --bfile QC_automsomal_Updated --update-map Final_map --make-bed --out Data_liftedSNPs
	else
		echo "all SNPs mapped nothing to lift over or SNPs were already in hg19 annotation" 
		mv QC_automsomal_Updated.bim  QC_automsomal.bim
		mv QC_automsomal_Updated.bed  QC_automsomal.bed
		mv QC_automsomal_Updated.fam  QC_automsomal.fam
		echo "ready for QC"
	fi


echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e "
###########################################################
# 2a) Remove all the SNPs which have 0 bp position  	    #
###########################################################
"

	awk 'int($4)==0' Data_liftedSNPs.bim|wc -l 

	awk 'int($4)==0' Data_liftedSNPs.bim|awk '{print $2}'> Data_noSNP_coordinates.txt

	read count x <<< $(wc -l Data_noSNP_coordinates.txt)

	if [ $count -gt 0 ]
	then
		$PLINK --bfile Data_liftedSNPs --exclude Data_noSNP_coordinates.txt --make-bed --out Data_liftedSNPs_1
	else
		mv Data_liftedSNPs.fam Data_liftedSNPs_1.fam
		mv Data_liftedSNPs.bim Data_liftedSNPs_1.bim
		mv Data_liftedSNPs.bed Data_liftedSNPs_1.bed
	fi



echo -e 				"...........................$count SNPs were not mapped and were excluded from the analysis\n..........................."


echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e "
#######################################################
# 2b) Remove all SNPs which are not annotated in hg19 #
#######################################################
"

	echo -e "Recheck missing hg19 annotations \n"

	if [[ -e Data_liftedSNPs_1.bim ]]
	then
		awk {'print $2 "\t" $4'} Data_liftedSNPs_1.bim>Data_snps
		echo -e "DataSNPs written \n"
	else
		echo -e "Please check the existence of Data_liftedSNPs_1.bim files \n"
	fi

echo -e '
######################################################
'

	if  [[ -e map.file ]]
	then
    		mv map.file hg19_snps
    	echo -e "hg19file exists \n"
	else
    		awk {'print $4 "\t" $3'} $hg19SNPs> hg19_snps
    		echo -e "hg19SNps written \n"
	fi


echo -e '
######################################################
'

	if [ -e Data_snps_sorted ]
	then
    		echo "Data sorted file exists \n"
	else
    		sort Data_snps>Data_snps_sorted
    	echo -e "Data SNPs sorted \n"
	fi


echo -e '
######################################################
'


	if [ -e hg19_snps_sorted ]
	then
    		echo "hg19 sorted file exists \n"
	else
    		sort hg19_snps>hg19_snps_sorted
		echo -e "hg19 SNPs sorted \n"
	fi


echo -e '
######################################################
'


	comm -23 Data_snps_sorted hg19_snps_sorted>SNPs_notInHG19 # Find bp positions of SNPs in study data that are different from hg19 bp position of SNPs
 

	echo "$(wc -l SNPs_notInHG19) SNPs have a different bp position than hg19, so they will be removed from the main file \n"


	read count x <<< $(wc -l SNPs_notInHG19) ## remove SNPS not HG19 annotated and not lifted 

	if [ $count -gt 0 ]
	then
		$PLINK --bfile Data_liftedSNPs_1 --exclude SNPs_notInHG19 --make-bed --out Data_lifted_hg19only
		echo "$count SNPs have been removed from the main file \n"
	else
		cp  Data_liftedSNPs_1.fam QC_autosomal.fam
		cp  Data_liftedSNPs_1.bim QC_autosomal.bim
		cp  Data_liftedSNPs_1.bed QC_autosomal.bed
		echo -e "no SNPs have been removed all fine \n"
    
	fi



echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

	echo -e "Excluding SNPs with duplicated positions \n"

	sort -k3n Data_lifted_hg19only.bim |  uniq  -f2 -D | cut -f2 > Duplicated_SNPs.txt


	$PLINK --bfile Data_lifted_hg19only --exclude Duplicated_SNPs.txt --make-bed --out QC_autosomal

fi

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e "
##########################################################
# 3) Extract affected individuals if provided		   #
##########################################################
"

if [ $ExtractAffected=="Yes" ]
then
	max=$(awk '{print $6}' QC_autosomal.fam | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print total/count, max, min}'| awk -F ' ' '{print $2}')
	min=$(awk '{print $6}' QC_autosomal.fam | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print total/count, max, min}'| awk -F ' ' '{print $3}')
	if [ ! "$min" -eq "$max" ] #If Plink fam file contains affected indiviudal statuses
	then
		threshold=2
		echo -e "Select only the affected individuals"
		awk -v threshold="$threshold" '$6 == threshold' QC_autosomal.fam|awk '{print $1 "\t" $2}'> Affected_Data_liftedSNPs.txt
		$PLINK --bfile QC_autosomal --keep Affected_Data_liftedSNPs.txt --make-bed --out QC_affected

	elif [ -e $AffectedInds ] #If there is already a list of affected individuals occuring
	then
		$PLINK --bfile QC_autosomal --keep $AffectedInds --make-bed --out QC_affected

	elif [ "$min" -eq "$max" ] && [ ! -e $AffectedInds ]
	then
		echo -e "Please provide a file with affected individuals only, with two columns consisting of FID and IID respectively otherwise all individuals will be selected for performing further analysis"
	    	cp  QC_autosomal.fam QC_affected.fam
	    	cp  QC_autosomal.bim QC_affected.bim
	    	cp  QC_autosomal.bed QC_affected.bed
	fi

	else
	echo -e "All individuals are selected for further analysis \n"
fi


echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"


echo -e "
######################################################
# 4) Stand Check 				       #
######################################################
"

if [ -e QC_affected.bim ]
then
	for ((i=1;i<=22;i++))
	do
	$PLINK --bfile QC_affected --chr "$i" --make-bed --out Affected_Data_liftedSNPs"$i"
	done
fi




echo -e '
########### Performing Strand check ###################
'

START=$(date +%s)
for ((i=1;i<=22;i++))
do
$SHAPEIT -check \
        -B Affected_Data_liftedSNPs"$i" \
        -M $MapFile/genetic_map_chr"$i"_combined_b37.txt \
        --input-ref $ShapeitRefHaps/ALL_1000G_phase1integrated_v3_chr"$i"_impute.hap.gz \
        $ShapeitRefLegend/ALL_1000G_phase1integrated_v3_chr"$i"_impute.legend \
        $ShapeitRefSample/ALL_1000G_phase1integrated_v3.sample \
        --output-log gwas.alignments"$i"
done
END=$(date +%s) #~3 mins


echo -e '
######################################################
'


chmod 770 *

for ((i=1;i<=22;i++))
do
	if [[ -e gwas.alignments"$i".snp.strand ]]
	then
	threshold="Strand"
	awk -v threshold="$threshold" '$2 == threshold' gwas.alignments"$i".snp.strand|awk '{print $4}'| sort | uniq>SNP_to_flip_Chr"$i"
	else
	echo -e "Affected_Data_liftedSNPs"$i" did not had enough variants to proceed \n"
	fi
done

wait

echo -e '
######################################################
'


echo -e "Flip the SNPs which have Strand problems \n"

for ((i=1;i<=22;i++))
do
	if [[ -e SNP_to_flip_Chr"$i" ]]
	then
	$PLINK --bfile Affected_Data_liftedSNPs"$i" --flip SNP_to_flip_Chr"$i" --make-bed --out Set1_Affected_Data_liftedSNPs"$i" &
	else
	mv Affected_Data_liftedSNPs"$i".bed Set1_Affected_Data_liftedSNPs"$i".bed
	mv Affected_Data_liftedSNPs"$i".bim Set1_Affected_Data_liftedSNPs"$i".bim
	mv Affected_Data_liftedSNPs"$i".fam Set1_Affected_Data_liftedSNPs"$i".fam
	fi
done
wait


echo -e '
######################################################
'

for ((i=1;i<=22;i++))
do
	if [[ ! -e Set1_Affected_Data_liftedSNPs"$i".bim ]]
	then
	echo "Affected_Data_liftedSNPs"$i" did not had enough variants to proceed \n"
	fi
done

echo -e '
######################################################
'

echo -e "Check if the shapeit reports "Check" Shapeit that have reported "Missing" \n"


for ((i=1;i<=22;i++))
do
threshold="Missing"
awk -v threshold="$threshold" '$2 == threshold' gwas.alignments"$i".snp.strand|awk '{print $4}'| sort | uniq>SNP_to_Exclude_Chr"$i" &
done
wait


echo -e '
######################################################
'

echo -e "Exclude those SNPs which have strand issues as Missing \n"

for ((i=1;i<=22;i++))
do
	if [[ -e SNP_to_Exclude_Chr"$i" ]]
	then
	$PLINK --bfile Set1_Affected_Data_liftedSNPs"$i" --exclude SNP_to_Exclude_Chr"$i" --make-bed --out Set2_Affected_Data_liftedSNPs"$i" &
	else
	mv Set1_Affected_Data_liftedSNPs"$i".bed Set2_Affected_Data_liftedSNPs"$i".bed
	mv Set1_Affected_Data_liftedSNPs"$i".bim Set2_Affected_Data_liftedSNPs"$i".bim	
	mv Set1_Affected_Data_liftedSNPs"$i".fam Set2_Affected_Data_liftedSNPs"$i".fam
	fi
done
wait


echo -e '
######################################################
'



echo -e "Run the final strand check to see if there are still any complains"

for ((i=1;i<=22;i++))
do
echo "
#######################################################
#######################chr$i###########################
#######################################################

"
	if [[ -e Set2_Affected_Data_liftedSNPs"$i".bim ]]
	then
	$SHAPEIT -check \
        -B Set2_Affected_Data_liftedSNPs"$i" \
        -M $MapFile/genetic_map_chr"$i"_combined_b37.txt \
        --input-ref $ShapeitRefHaps/1000GP_Phase3/1000GP_Phase3_chr"$i".hap.gz \
	            $ShapeitRefHaps/1000GP_Phase3/1000GP_Phase3_chr"$i".legend.gz \
	            $ShapeitRefHaps/1000GP_Phase3/1000GP_Phase3.sample \
        --output-log gwas.alignments"$i" &
	else
	echo -e "Please check the existence of Set2_Affected_Data_liftedSNPs"$i" before proceeding \n"
	sleep 1
	fi
done
wait




echo -e '
######################################################
'

echo " 
###################################################### 
###### Exclude all the problematic snps###############
######################################################"

for ((i=1;i<=22;i++))
do
	if [[ -e gwas.alignments"$i".snp.strand.exclude ]]
	then
	echo "chr$i"
	$SHAPEIT -check \
        -B Set2_Affected_Data_liftedSNPs"$i" \
        -M $MapFile/genetic_map_chr"$i"_combined_b37.txt \
	--noped \
        --input-ref $ShapeitRefHaps/1000GP_Phase3/1000GP_Phase3_chr"$i".hap.gz \
	            $ShapeitRefHaps/1000GP_Phase3/1000GP_Phase3_chr"$i".legend.gz \
	            $ShapeitRefHaps/1000GP_Phase3/1000GP_Phase3.sample \
        --exclude-snp gwas.alignments"$i".snp.strand.exclude &
	fi
done
wait



echo -e "..................................................................Final QC before submitting to pre-phasing................................................................."

for ((i=1;i<=22;i++))
do
	$PLINK --bfile Set2_Affected_Data_liftedSNPs"$i" --geno $GENO --maf $MAF --hwe $HWE --make-bed --out tmp

echo -e '
######################################################
'

	$PLINK --bfile tmp --mind $MIND --make-bed --out   All_Affected_DataSNPs"$i"


done
wait


echo -e "
#####################################################
#5) Prephasing					      #
#####################################################
"



cd ../../Scripts

echo -e ".....................................................................Preparing jobs for Pre-phasing \n................................................................."

for ((i=1;i<=22;i++))
do
	if [[ -e ../OUTPUT_DIR/Stage2_GenoImpute/All_Affected_DataSNPs"$i".bim ]]
	then

python << END
import glob
import os
import optparse

parser = optparse.OptionParser()
parser.add_option("-f", "--factor", action="store", dest="fact", help="Factor X")
(options, args) = parser.parse_args()

inputfiles = glob.glob("../OUTPUT_DIR/Stage2_GenoImpute/All_Affected_DataSNPs*.fam")

print len(inputfiles)

import sys
sys.path.insert(0, '../ConfigFiles')

i = 1
for inputfile in inputfiles: 
    jobname = "shapeitrest"+str(i)+".sh"	
    f = open(jobname,"w")
    f.write("$SHAPEIT -B All_Affected_DataSNPs"+str(i)+" -M $MapFile/genetic_map_chr"+str(i)+"_combined_b37.txt --duohm --exclude-snp gwas.alignments"+str(i)+".snp.strand.exclude -O gwas_Data.chr"+str(i)+".phased -T 8")
    f.close()
    i = i+1

END
	else
	echo -e "Please check the existence of All_Affected_DataSNPs"$i" file \n"
	fi
done

echo -e "------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"


echo -e ".....................................................................Starting Pre-phasing \n................................................................."


chmod a+x *.sh
#Create the file for running at the server

cd ../OUTPUT_DIR/Stage2_GenoImpute



cd OUTPUT_DIR/Stage2_GenoImpute


chmod 770 *

for i in $(seq 22); do
./../../Scripts/shapeitrest"$i".sh &
done

wait

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"


echo -e ".....................................................................Pre-phasing Ended \n................................................................."


echo -e "
#####################################################
# 6) Convert all the outputs to VCF format	      #
#####################################################
"


chmod 770 *

for ((i=1;i<=22;i++))
do
	$SHAPEIT -convert \
        --input-haps gwas_Data.chr"$i".phased \
        --output-vcf gwas_Data.chr"$i".vcf &
done
wait



echo -e "Done converting to vcf \n"

cd OUTPUT_DIR/Stage2_GenoImpute


echo -e ".....................................................................Preparing jobs for submitting to Minimac for Imputation \n................................................................."

echo -e "
######################################################
# 7) Run Minimac.py				       #	
######################################################
"


cd ../../Scripts

python minimac.py

wait

chmod 770 *.sh

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"


echo -e ".....................................................................Starting Minimac for Imputation \n................................................................."

cd ../OUTPUT_DIR/Stage2_GenoImpute



for ((i=1;i<=22;i++))
do
    ./../../Scripts/phase2_"$i".sh &
    sleep 1
done
wait

echo "done Phase 2"

echo -e '
######################################################
'



cd OUTPUT_DIR/Stage2_GenoImpute

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"


echo -e "
######################################################
#8)Find SNPs which passed QC threshold of imputation#
######################################################
"

#Change the spacing

for((i=1;i<=22;i++))
do
	if [[ -e Gwas.Chr"$i"_Study.Imputed.Output.info ]]
	then
	sed -i $'s/\t/ /g' Gwas.Chr"$i"_Study.Imputed.Output.info &
	fi
done
wait

echo -e '
######################################################
'

echo -e "Filter the SNPs based on quality threshold \n"

for((i=1;i<=22;i++))
do
	if [[ -e Gwas.Chr"$i"_Study.Imputed.Output.info ]]
	then
	threshold=$thresholdImp
	awk -v threshold="$threshold" '$7 >= threshold' Gwas.Chr"$i"_Study.Imputed.Output.info| awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}' > tmp"$i" &
	else
	echo -e "WARNING!!!! No info file exists for Chromosome "$i" , please check if there was any error reported before proceeding forward \n"
	#sleep 2	
	fi
done
wait

echo -e '
######################################################
'

echo -e 'Updating SNP names. VCF files contains SNP identifiers as ids combined with chr:snpid:ph.pos:bp:allele1_allele2'

i=22
Nr=$(cat tmp"$i"| awk -v FS='\t' '{print $1}'| head -2 | tail -n 1|grep -o ':'|wc -l)


if [[ $Nr>1 ]]
then
	for((i=1;i<=22;i++))
	do
		if [[ -e tmp"$i" ]]
		then
		awk '{print $1}' tmp"$i" | cut -f1,2 -d':' > SNPnames$i
		paste tmp"$i" SNPnames"$i"|awk '{print $14 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13}'>testset"$i"
		awk '{if ($1 ~ /rs/) {a=sub(/\:.*$/,"",$1); print $0;} else {print $0}}' testset"$i">testset2"$i"
		grep -E "^[0-9]|^rs"  testset2"$i" > Chr"$i"_StudyImputed.txt &
		else
		echo -e "WARNING!!!! No tmp file exists for Chromosome "$i" , please check if there was any error reported before proceeding forward \n"
		fi
	done
else
	for((i=1;i<=22;i++))
	do
	grep -E "^[0-9]|^rs"  tmp"$i" > Chr"$i"_StudyImputed.txt &
	done
fi
wait



echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e '
######################################################
'

echo -e "Convert the vcf files to PLINK files \n"

for((i=1;i<=22;i++))
do
	if [[ -e Gwas.Chr"$i"_Study.Imputed.Output.dose.vcf.gz ]]
	then
	$PLINK --vcf Gwas.Chr"$i"_Study.Imputed.Output.dose.vcf.gz --chr "$i" --make-bed --out Chr"$i"_Data_Imputed_plinkables &
	else
	echo -e "WARNING!!!! No vcf file exists for Chromosome "$i" , please check if there was any error reported before proceeding forward \n"
	fi
done
wait




#VCF files contains SNP identifiers as ids combined with chr:snpid:ph.pos:bp:allele1_allele2, this needs to be separated by the following commands

Nr=$(cat Chr"$i"_Data_Imputed_plinkables.bim| awk -v FS='\t' '{print $2}'| head -2 | tail -n 1|grep -o ':'|wc -l)

if [[ $Nr>1 ]]
then
#Removing everything after second colon (works)
	for((i=1;i<=22;i++))
	do
	awk '{print $2}' Chr"$i"_Data_Imputed_plinkables.bim | cut -f1,2 -d':' > SNPnames$i
	done


# combining files
	paste Chr"$i"_Data_Imputed_plinkables.bim SNPnames"$i"|awk '{print $1 "\t" $7 "\t" $3 "\t" $4 "\t" $5 "\t" $6}'>testset"$i"

	awk '{if ($2 ~ /rs/) {a=sub(/\:.*$/,"",$2); print $0;} else {print $0}}' testset"$i">Chr"$i"_Data_Imputed_plinkables2.bim
	mv Chr"$i"_Data_Imputed_plinkables.bed Chr"$i"_Data_Imputed_plinkables2.bed
	mv Chr"$i"_Data_Imputed_plinkables.fam Chr"$i"_Data_Imputed_plinkables2.fam
else
	mv Chr"$i"_Data_Imputed_plinkables.bed Chr"$i"_Data_Imputed_plinkables2.bed
	mv Chr"$i"_Data_Imputed_plinkables.fam Chr"$i"_Data_Imputed_plinkables2.fam
	mv Chr"$i"_Data_Imputed_plinkables.bim Chr"$i"_Data_Imputed_plinkables2.bim
fi



chmod 770*

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e '
######################################################
'

echo -e "Select only biallelic and the ones passing threshold using PLINK \n"

for((i=1;i<=22;i++))
do
	if [[ -e Chr"$i"_StudyImputed.txt ]]
	then
	$PLINK --bfile Chr"$i"_Data_Imputed_plinkables2 --biallelic-only strict --extract Chr"$i"_StudyImputed.txt --make-bed --out Chr"$i"_Data_Imputed_biallelic_Filtered &
	else
	echo -e "WARNING!!!! No Imputed file exists for Chromosome "$i" , please check if there was any error reported before proceeding forward \n"
	fi
done
wait


echo -e '
######################################################
'

echo -e "
#######################################################
# 9)Convert all the output SNPs "bp:chr" to snps names#
#######################################################
"
#We used the latest build of the time when MAGNET was written i.e dbSNP 151 ( 09-Dec-2018)


echo -e "Extracting the rsid names of the filtered files... \n"

for((i=1;i<=22;i++))
do
	if [[ -e Chr"$i"_Data_Imputed_biallelic_Filtered.bim ]]
	then
	awk 'FNR==NR{a[$4]=$0;next}{if ($2 in a){p=$2;$2="";print a[p],$0} else {print p="NA" " " $4 " " $2 " " $2 " " $1 " " $3 " " $4 " " $5 " " $6}}'  $AnnotationFile/Comp_chr"$i"_wHead.txt Chr"$i"_Data_Imputed_biallelic_Filtered.bim>Annotated_Final_Imputed_Data_Chr"$i"SNPs.bim 
	else
	echo "Chr"$i"_Data_Imputed_biallelic_Filtered.bim does not exist"
	fi
done
wait

echo -e '
######################################################
'

#Arrange the file as in bim format
for ((i=1;i<=22;i++)); 
do
	if [[ -e Annotated_Final_Imputed_Data_Chr"$i"SNPs.bim ]]
	then
	awk '{print $5 " " $3 " " $6 " " $2 " " $8 " " $9}' Annotated_Final_Imputed_Data_Chr"$i"SNPs.bim>tmp"$i"SNPs.bim &
	else
	echo -e "Please check the existence of Annotated_Final_Imputed_Data_Chr"$i"SNPs.bim before proceeding \n"
	sleep 2
	fi
done
wait

echo -e '
######################################################
'


for ((i=1;i<=22;i++)); 
do
	if [[ -e tmp"$i"SNPs.bim ]]
	then
	mv tmp"$i"SNPs.bim Chr"$i"_Data_Imputed_biallelic_Filtered.bim &
	else
	echo -e "Please check the existence of tmp"$i"SNPs.bim before proceeding \n"
	sleep 2
	fi
done
wait

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e "-------------------------------------------------------------------------------Perform Final QC of imputed SNPs----------------------------------------------------------------------------- \n"



echo -e '
##########################################################
#10) Perform QC of imputed SNPs			   #
##########################################################
'

for ((i=1;i<=22;i++)); 
do
	if [[ -e Chr"$i"_Data_Imputed_biallelic_Filtered.bim ]]
	then
	$PLINK --bfile Chr"$i"_Data_Imputed_biallelic_Filtered  --geno $GENO  --hwe $HWE --maf $MAF --make-bed --out tmp2"$i"
	else
	echo -e "Please check the existence of Chr"$i"_Data_Imputed_biallelic_Filtered.bim before proceeding \n"
	sleep 2
	fi
done
wait

echo -e '
######################################################
'

for ((i=1;i<=22;i++)); 
do
	if [[ -e tmp2"$i".bim ]]
	then
	$PLINK --bfile tmp2"$i" --mind $MIND --make-bed --out Final_Imputed"$i" >> ../../LOG/Stage1_GenoQC/plinkoutput_QC1_Chr"$i".log &
	else
	echo -e "Please check the existence of tmp2"$i".bim before proceeding \n"
	sleep 2
	fi
done
wait

echo -e '
######################################################
'


for ((i=1;i<=22;i++)); 
do
rm Annotated_Final_Imputed_Data_Chr"$i"SNPs*
done


echo -e '
######################################################
'

touch listForMergeSNPs_plink.txt

ls -1v Final_Imputed*.bed>a
ls -1v Final_Imputed*.bim>b
ls -1v Final_Imputed*.fam>c

paste a b c>listForMergeSNPs_plink.txt


echo -e "-----------------------------------------------------Merge All the Chromosome SNPs--------------------------------------- \n"

echo -e '
#########################################################
# 11) Merge All the Chromosome SNPs			  #
#########################################################
'

if [[ -e listForMergeSNPs_plink.txt ]]
then
PlinkFinal=$(echo *  | head -n1 listForMergeSNPs_plink.txt |cut -f1|cut -d. -f1)
else
echo -e "Please check the existence of listForMergeSNPs_plink.txt before proceeding \n"
fi

echo -e '
######################################################
'

sed '1d' listForMergeSNPs_plink.txt>SNPsToMerge.txt

$PLINK --bfile $PlinkFinal --merge-list SNPsToMerge.txt --out Merged_Imputed_FinalQC_SNPs_Data

echo -e "If merging not successful exclude any 3+ alleles \n"

if [[ -e Merged_Imputed_FinalQC_SNPs_Data.missnp ]]
then
	for ((i=1;i<=22;i++)); 
	do
	$PLINK --bfile Final_Imputed"$i" --exclude Merged_Imputed_FinalQC_SNPs_Data.missnp --make-bed --out tmp_chr"$i" &
	done
	wait

echo -e '
######################################################
'

#Rename the tmp files to Final Imputed Files
	for ((i=1;i<=22;i++))
	do
		if [[ -e tmp_chr"$i".bim ]]
		then
		mv tmp_chr"$i".bim Final_Imputed"$i".bim
		mv tmp_chr"$i".bed Final_Imputed"$i".bed
		mv tmp_chr"$i".fam Final_Imputed"$i".fam
		fi
	done

echo -e '
######################################################
'

#Save the list of SNPs with 3+ snps

	mv Merged_Imputed_FinalQC_SNPs_Data.missnp SNPs_with_3plusalleles

#Merge them again 

	$PLINK --bfile $PlinkFinal --merge-list SNPsToMerge.txt --out Merged_FinalQC_SNPs_Data


else
	mv Merged_Imputed_FinalQC_SNPs_Data.bim Merged_FinalQC_SNPs_Data.bim
	mv Merged_Imputed_FinalQC_SNPs_Data.bed Merged_FinalQC_SNPs_Data.bed
	mv Merged_Imputed_FinalQC_SNPs_Data.fam Merged_FinalQC_SNPs_Data.fam

fi

echo -e '
######################################################
'

wc -l Merged_FinalQC_SNPs_Data.bim 

if [[ -e Merged_FinalQC_SNPs_Data.bim ]]
then
	echo -e "--------------For further analysis we need to split snp file into chunks of 5000 and report the minor allele frequency of the SNPS--------------\n"
else
	echo "Imputation did not run successfully please check back the log files"
fi

echo -e '
######################################################
'
echo -e '
######################################################
# 12) Split into chunks of 5000,		       #
# Find the minor allele frequency of the SNPS	       #
######################################################
'


awk '{print $2}' Merged_FinalQC_SNPs_Data.bim >AllSNPs.txt


NSNPs=($(wc -l AllSNPs.txt)) #Total number of SNPs

DivFiles=$[$NSNPs/$chunkSize] #Files after dividing into chunks of 5000

if [[($(expr $NSNPs % $chunkSize) -gt 1)]]; then #If there is a remainder existing after dividing it into chunks of 5000
	NFiles=$(expr $DivFiles + 1) #If exists then add 1 to the number of main files
	else Nfiles=$DivFiles #Else let the number as it is
fi


digi=${#NFiles}

x=$NFiles
split -a $digi -l 5000 -d --additional-suffix=.txt AllSNPs.txt splitFiles



chmod 775 splitFiles*.txt

Files=($(ls splitFiles*.txt))

for (( i=0; i<=$x; i++ ))
do
echo ${Files[i]}
	if [[ -e ${Files[i]} ]]
	then
	$PLINK --bfile Merged_FinalQC_SNPs_Data --extract ${Files[i]} --make-bed --out tmp"$i"
	fi
wait
done

for (( i=0; i<=$x; i++ ))
do
	if [[ -e tmp"$i".bim ]]
	then
	echo ${Files[i]}
	$PLINK --bfile tmp"$i" --recodeA --out Data_SNPfile"$i" 
	rm tmp"$i"
	fi
done


echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e "----------------------------------------------------------------Stage 2 Imputation completed successfully-------------------------------------------------------------------------- \n"































