#!/bin/bash

MAGNET=$(pwd)

cd $MAGNET

LocSource=$(pwd)

echo -e "
######################################
# Sourcing the configuration files   #
######################################
"


if [ -e /$LocSource/ConfigFiles/Tools.config ]
then
	source /$LocSource/ConfigFiles/Tools.config
else
	echo -e "Please check the path of configuration files, the program exists \n"
exit
fi

if [ -e /$LocSource/ConfigFiles/Thresholds.config ]
then
	source /$LocSource/ConfigFiles/Thresholds.config
else
	echo -e "Please check the path of configuration files, the program exists \n"
exit
fi

echo -e "
#######################################
# Stage 1: Quality Control	      #
#######################################
"

#### Read the arguments passed

if [[ ! -e $installationDir ]] #Check if the installation directory is existing
then
    	mkdir -p $installationDir/OUTPUT_DIR/Stage1_GenoQC #If not create directory
	chmod 770 *
else
	echo -e "...Installation directory already exists, please make sure you are not overwriting any previous results... \n"
fi
wait


sleep 1

echo -e "...Checking if Data files exist \n"


Alleles=(A G C T) #Save alleles in ACGT format
AllelesRec=(1 2 3 4) #Save alleles in 1234 format


if [ -e $SamplesToQC.fam ]
then
	AllelesInData=$(cut -f5 $SamplesToQC.bim | sort | uniq) #Save the alleles format provided in the data
else
	if [ -e $SamplesToQC.ped ]
	then
	AllelesInData=$(cut -d ' ' -f7 $SamplesToQC.ped | sort | uniq)
	fi
fi

read -a arr2 <<< $AllelesInData



#Create an array
declare -a AllelesInData2

AllelesInData2=($AllelesInData)


if [[ $a == [a-z] ]]
then
	echo "...Alleles will be converted to upper case"
fi


A=${arr2[@]}; # Save the array with alleles in study data in variable A
B=${Alleles[@]}; # Save the array with alleles in ACGT format in variable B
C=${AllelesRec[@]} # Save the array with alleles in 1234 format in variable C



if [[ $A == $B || $C ]]
then #Check if the study alleles are ACGT/1234
	echo -e "...Alleles are correctly coded... \n"
else
	echo -e "...Alles are not correctly coded the code will exit, please provide ACGT coding and run the program again... \n"
exit
fi

echo -e '
######################################################
'
if [[ $A == [a-z] ]]
then
  touch AllelesFile.txt
  awk '{ print $2 "\t" $5 "\t" $6 "\t" toupper ($5) "\t" toupper ($6) }'  $SamplesToQC.bim > AllelesFile.txt
  $PLINK --bfile $SamplesToQC --update-alleles AllelesFile.txt --make-bed --out $SamplesToQC
fi

sleep 2

echo -e '
######################################################
'

if [[ $A == $C ]]
then
	$PLINK --bfile $SamplesToQC --make-bed --alleleACGT --out $SamplesToQC
fi
wait

if [ -e $SamplesToQC.fam ]
then
	echo -e ".....$SamplesToQC.fam exists \n"
	convert="false"
else
	if [ -e $SamplesToQC.ped ]
	then
	echo -e"..... $SamplesToQC.ped exists \n"
	convert="true"
	else
	echo -e ".....$SamplesToQC.fam or $SamplesToQC.ped  does not exist; program exits \n"
	exit
	fi
fi
wait


if [ $convert ==  "true" ]
then
	if [ -e $SamplesToQC.map ]
	then
	$PLINK --file $SamplesToQC --make-bed --out $installationDir/OUTPUT_DIR/Stage1_GenoQC/QC1.1
	else
	echo -e "$SamplesToQC.map file does not exist. Program exits \n"
	exit
  fi
wait

else
	if [ -e $SamplesToQC.bim ]
	then
	echo -e ".....$SamplesToQC.bim exists \n"
	else
	echo -e ".....$SamplesToQC.bim does not exist; program exits \n"
	exit
	fi
wait

	if [ -e $SamplesToQC.bed ]
	then
	echo -e ".....$SamplesToQC.bed exists \n"
	cp $SamplesToQC.bed $installationDir/OUTPUT_DIR/Stage1_GenoQC/QC1.1.bed
	cp $SamplesToQC.bim $installationDir/OUTPUT_DIR/Stage1_GenoQC/QC1.1.bim
	cp $SamplesToQC.fam $installationDir/OUTPUT_DIR/Stage1_GenoQC/QC1.1.fam
	else
	echo -e ".....$SamplesToQC.bed does not exist; program exits \n"
	exit
  fi
fi
wait

echo -e '
######################################################
'


cd $installationDir/OUTPUT_DIR/Stage1_GenoQC

read ind filename <<< $(wc -l QC1.1.fam)
read snps filename <<< $(wc -l QC1.1.bim)

echo -e "...Further calculations will be done in output directory \n"


echo -e "...Your input file has" $ind  "samples and" $snps  "SNPs \n"

sleep 2


echo -e '
#######################################################################################
# Based on QC thresholds provided in config file QC will be performed of genotype data#
#######################################################################################
'
echo -e '
######################################################
#	Checking missingness, HWE and MAF    	             #
######################################################
'

echo -e "...If no QC options provided default filtering criteria will be applied: \n Callrate>95% \n Genotyping-rate>95% \n MAF>0.02 \n HWE>10e-8 \n"
$PLINK --bfile QC1.1  --make-bed --out QC1.1

echo -e '
######################################################
'

if [ -e QC1.1.bed ]
then
	$PLINK --bfile QC1.1  --geno $GENO --hwe $HWE --maf $MAF --make-bed --out tmp
	sleep 1
else
	echo -e "No QC1.1 file exits, please check the file existence. Program exits \n"
	exit
fi
wait

sleep 2

echo -e '
######################################################
'

if [ -e tmp.bed ]
then
	echo -e	"...Checking rate of missingness per individual... \n"
	$PLINK --bfile tmp --mind $MIND --make-bed --out QC1
	sleep 1
else
	echo -e "No tmp file exits, please check the file existence. Program exits \n"
	exit
fi

echo -e '
######################################################
'


rm tmp* #remove extra files
sleep 1

if [ -e QC1.irem ]
then
	read  irem filename <<< $(wc -l QC1.irem)
	sleep 1
else
	irem=0
fi

read SNPs filename <<< $(wc -l QC1.bim)

Diff=$(expr $snps - $SNPs)

echo -e "....." $irem  "samples were removed due to genotyping freq < 95% \n
....." $Diff  "SNPS were removed due to call rates  < 95%; MAF <0.02 or significant (10e-8) deviation from HWE \n"

sleep 2

echo -e '
######################################################
'


echo -e "...Report the missingness per SNP and individual wise...\n"

if [ -e QC1.bed ]
then
	$PLINK --bfile QC1 --missing --out QC1_report
	sleep 1
else
	echo -e "No QC1 file exits, please check the file existence. Program exits \n"
	exit
fi
sleep 1
wait

echo -e '
######################################################
'

echo -e "...Report the hardy weinberg distribution...\n"

if [ -e QC1.bed ]
then
	$PLINK --bfile QC1 --hardy --out PreQC_hardy
	sleep 1
else
	echo -e "No QC1 file exits, please check the file existence. Program exits \n"
	exit
fi
sleep 2
wait


echo -e '
######################################################
'


echo -e "...Report the PreQC Allele Frequency... \n"

if [ -e QC1.bed ]
then
	$PLINK --bfile QC1 --freq --out PreQC_AlleleFreq
	sleep 1
else
	echo -e "No QC1 file exits, please check the file existence. Program exits \n"
	exit
fi
sleep 2
wait

echo -e '
######################################################
'

echo -e "...Report the Inbreeding.. \n"

if [ -e QC1.bed ]
then
	$PLINK --bfile QC1 --het --out PreQC_Inbreeding
	sleep 1
else
	echo -e "No QC1 file exits, please check the file existence. Program exits \n"
	exit
fi
sleep 2
wait

echo -e '
######################################################
'



echo -e '
######################################################
#		Performing Sex check                  #
######################################################
'

echo -e '...checking for genderdiscrepancies... \n'


if [ -e QC1.bed ]
then
	$PLINK --bfile QC1 --check-sex --out QC2_Sexcheck
	sleep 2
else
	echo -e "No QC1 file exits, please check the file existence. Program exits \n"
	exit
fi
sleep 1

grep "PROBLEM" QC2_Sexcheck.sexcheck | awk '$3 == 0' > QC2_NAGender.txt

grep "PROBLEM" QC2_Sexcheck.sexcheck | awk '$3 != 0' > QC2_BadGender.txt

read NAgend filename <<< $(wc -l QC2_NAGender.txt)
read Badgend filename <<< $(wc -l QC2_BadGender.txt)

echo -e "...$NAgend samples had no gender annotated...\n
	 ...$Badgend samples had gender indiscrepancies... \n"


sleep 2

echo -e '
######################################################
'


#If there are samples with no gender annotated

if [ $NAgend -gt 0 ]
then
	echo -e "
	>>>> Your file has" $NAgend "unknown samples \n

	>>>> Your data consists of samples without any gender annotated, if you wish to correct for them, please rerun the analysis after correcting them else analysis will continue with them \n"
	#rename QC1 QC2 QC1.*;
	mv QC1.bim QC2.bim
	mv QC1.bed QC2.bed
	mv QC1.fam QC2.fam
	echo -e ".....samples with unknown gender will not be considered further \n"
	NAupdate=1
else
	NAupdate=0
	echo -e '.....No samples with unknown gender present in Dataset \n'
	#rename QC1 QC2  QC1.*;
	mv QC1.bim QC2.bim
	mv QC1.bed QC2.bed
	mv QC1.fam QC2.fam
fi

sleep 2

echo -e '
######################################################
'


#If there are samples with gender discrepancies

if [ $Badgend -gt 0 ]
then
	echo -e '...the following samples had gender discrepancies \n'
	cat QC2_BadGender.txt
	sleep 2
	echo -e "
	>>>> Your data has" $Badgend "samples with gender inconsistencies, if you wish to correct for them, please rerun the analysis after correcting them \n

	>>>> As per default option, samples with gender inconsistencies will be removed \n"

	awk {'print($1 "\t"  $2)'} QC2_BadGender.txt > QC2_BadGender_removed.txt;

	echo -e ".....Samples with genderincongruencies being removed \n";

	$PLINK --bfile QC2 --remove QC2_BadGender_removed.txt --make-bed --out QC3
	sleep 2
	if [ -e QC3.bim ]
	then
		echo -e "........$Badgend Samples with gender incongruencies were succesfully  removed \n";
	else
		echo -e ".....update of gender incongruencies failed pelase check plinkoutput_QC3.log \n";
	fi

else
	echo -e ".....no gender discrepancies detected \n"
	#rename QC2 QC3 QC2.*
	mv QC2.bim QC3.bim
	mv QC2.bed QC3.bed
	mv QC2.fam QC3.fam
fi
sleep 2

echo -e '
######################################################
'

echo -e "...Final Sexcheck... \n"

if [ -e QC3.bim ]
then
	$PLINK --bfile QC3 --check-sex --out QC4_Sexcheck
	sleep 2
else
	echo -e "Please check the existence of QC3 files, else program exits \n"
	exit
fi

grep "PROBLEM" QC4_Sexcheck.sexcheck | awk '$3 == 0' > QC4_NAGender.txt
grep "PROBLEM" QC4_Sexcheck.sexcheck | awk '$3 != 0' > QC4_BadGender.txt
read prob0 x <<< $(wc -l QC4_NAGender.txt)
read probNA x <<< $(wc -l QC4_NAGender.txt)
sleep 2

echo -e '
######################################################
'

echo -e "...Report Missing Gender Handling: \n"

if [ $NAupdate == 1 ]
then
	echo -e ".....$prob0 samples without annotated gender will not be considered further \n";
fi
sleep 1

if [ $NAupdate == 0 ]
then
	echo -e ".....No Sample/s without annotated gender were present in original dataset \n";
fi

echo -e "...Gender issues in current files,consider reviewing these files \n"
sleep 2

echo -e '
######################################################
'


if [ $prob0 -gt 0 ]
then
  echo -e ".....there are still $prob0 gender incongruencies \n"

  echo -e ".....samples written to QC4_Badgender.txt; please check manually; program exits here \n";
  exit
else
	echo -e ".....No gender incongruencies remained after cleaning files \n";
fi
sleep 2

echo -e '
######################################################
'




echo -e '
######################################################
#	Checking for inbreeding and contamination    #
######################################################

'


if [ -e QC3.bim ]
then
	$PLINK --bfile QC3 --het --out Breeding_output
	sleep 2
else
	echo -e "Please check the existence of QC3 files, else program exits \n"
	exit
fi

echo -e '
######################################################
'

awk '$6>0.2'  Breeding_output.het | awk '$6!="F"'  > QC4_inbredIDs.txt
awk '$6<-0.15'  Breeding_output.het > QC4_contaminatedIDs.txt


read cont x  <<< $(wc -l QC4_contaminatedIDs.txt)
read inbrd x <<< $(wc -l QC4_inbredIDs.txt)

echo -e "
...there are $cont potentially contaminated samples with high heterozygosity scores and
...there are $inbrd samples potentially inbred with low heterozygosity scores \n"
sleep 2

echo -e '
######################################################
'



if [ $cont -gt 0 ]
then
	echo -e "...Contaminated samples will be excluded \n"

	echo -e ".....$cont contaminated samples being removed \n"

	awk {'print $1 "\t" $2'} QC4_contaminatedIDs.txt > Contaminated_samples_removed.txt

	$PLINK --bfile QC3 --remove Contaminated_samples_removed --make-bed --out QC4 >plinkoutput_QC4_rm_cont.txt
	sleep 2

	echo -e ".....$samples succesfully removed, samples written to Contaminated_samples_removed.txt \n";
else
	echo -e ".....no sample/s needed to be removed due to contamination \n"
	#rename QC3 QC4 QC3.*;
	mv QC3.bim QC4.bim
	mv QC3.bed QC4.bed
	mv QC3.fam QC4.fam
	echo -e "....QC4 passed \n"
fi
sleep 2

echo -e '
######################################################
'



if [ $inbrd -gt 0 ]
then
	echo -e '>>>The following individuals have high inbred scores \n'
	cat QC4_inbredIDs.txt

	echo -e ".....$inbrd samples will be removed \n"
	awk {'print $1 " " $2'} QC4_inbredIDs.txt > Inbred_samples_removed.txt

	$PLINK --bfile QC4 --remove Inbred_samples_removed.txt --make-bed --out QC5 > ../../LOG/Stage1_GenoQC/plinkoutput_QC5_rmInbred.log
	sleep 2

	if [ -e QC5.bim ]
	then
      	echo -e "...$inbrd samples removed, Samples written to Inbred_samples_removed.txt \n";
	else
	echo -e "...there were some issues with removing inbred samples. Please check plinkoutput_QC5_rmInbred.log \n";
	exit
	fi
else
	echo -e "...no samples needed to be removed due to inbreeding \n"
	#rename QC4 QC5 QC4.*
	mv QC4.bim QC5.bim
	mv QC4.bed QC5.bed
	mv QC4.fam QC5.fam
	echo -e "...QC5 passed \n"
fi
sleep 2

echo -e '
######################################################
'



echo -e "
#######################################################
#	      Checking for Mendel Errors              #
#######################################################
"

#parameter "Set-me missing" zero out specific Mendelian inconsistencies to be used in conjuction with me 1 1 so that no particular snp/individual is removed but to zero out specific genotypes

echo -e "...excluding families with more then 1% mendel errors and SNPs with more then 10% mendel errors \n"

if [ -e QC5.bim ]
then
	$PLINK --bfile QC5 --me $MEFam $MESNP --set-me-missing --make-bed --out QC6
	sleep 2


echo -e '
######################################################
'

	$PLINK --bfile QC6 --mendel --out MendelsummaryQC6
	sleep 2
else
	echo -e "Please check the existence of QC5 files, else program will exit \n"
	exit
fi

echo -e '
######################################################
'


echo -e "...checking again for callrate>95% genotyping-rate>95%  MAF>0.02 hwe>10e-8 \n"

if [ -e QC6.bim ]
then
	$PLINK --bfile QC6  --geno $GENO --hwe $HWE --maf $MAF --make-bed --out tmp
	sleep 2

echo -e '
######################################################
'

	$PLINK --bfile tmp --mind $MIND --make-bed --out  QC7
	sleep 2

	read SNPs filename <<< $(wc -l QC7.bim)
	snps=$(cat QC5.bim | wc -l)
	Diff=$(expr $snps - $SNPs)

	Indiv=$(cat QC7.fam | wc -l)
	Allind=$(cat QC5.fam | wc -l)

	irem=$(expr $Allind - $Indiv)
else
	echo -e "Please check the existence of QC6 files, program will exit \n"
fi

echo -e '
######################################################
'


echo -e "....." $irem  "sample/s were removed due to Mendel error rates \n"
echo -e "....." $Diff  "SNPS were removed due to call rates  < 95%; MAF <0.02 or significant (10e-8) deviation from HWE \n"
sleep 1


# clean all files

#rename QC7 QC8 QC7.*

mv QC7.bim QC8.bim
mv QC7.bed QC8.bed
mv QC7.fam QC8.fam

read ind filename <<< $(wc -l  QC8.fam)
read snps filename <<< $(wc -l  QC8.bim)

echo -e "
...Your Final file has" $ind  "samples and" $snps  "SNPs that passed QC
...Files are stored as $PLINK binary files QC_final.* \n"

wait


echo -e "
#########################################################
#	  Calculating Identity by descent               #
#########################################################
"


if [ -e QC8.bim ]
then
	$PLINK --bfile QC8  --genome --out QC8_Relcheck
  sleep 2
else
	echo -e "Check the existence of QC8 files. Program exits \n"
	exit
fi

echo -e '
######################################################
'
echo -e "...Relatedness scores across all individuals are being calculated... \n"
wait

if [ -e  QC8_Relcheck.genome ]
then
    echo -e "...Relcheck is succesfully finished.... \n"
else
    echo -e "...there were some major issues with your data. For details check QC8_Relcheck.log... \n"
fi
sleep 2

echo -e "...IBD scores successfully calculated... \n"

awk '$9 >0.8 && $10 >0.8'  QC8_Relcheck.genome > QC8_Relcheck.duplicates

echo -e '
######################################################
'


read dups x <<< $(wc -l QC8_Relcheck.duplicates)

echo -e "...there are/is $dups duplicated sample/s. Samples written to QC8_Relcheck.duplicates \n..."

sleep 2
## Decisions on how to deal with families, user should provide additional files to correct the errors if any


echo ">>> One duplicated individual will be deleted"

awk {'print $1 " " $2'} QC8_Relcheck.duplicates > QC8_duplicatesDeleted.txt;


if [ ! -z $DuplicatedIds ]
	then
	Duplicates=$DuplicatedIds
	else
	Duplicates=QC8_duplicatesDeleted.txt
	fi
$PLINK --bfile QC8 --remove $Duplicates --make-bed --out QC9 > ../../LOG/Stage1_GenoQC/plinkoutput_QC8_rmDupls.log
sleep 2

echo -e '
######################################################
'


if [ -e QC9.bed ]
then
	echo -e "...$(cat $Duplicates| wc -l) Sample/s was/were deleted. Sample/s written to QC8_duplicatesDeleted.txt \n"
else
	echo -e "...there were some issues with deleting duplicated samples, for more details please check plinkoutput_QC8.1_rmDupls.log \n"
fi

echo -e '
######################################################
'

mv QC8.bim QC9.bim
mv QC8.bed QC9.bed
mv QC8.fam QC9.fam


sleep 2

echo -e "...Checking for unrelated Parent offspring duos... \n"

$PLINK --bfile QC9 --genome rel-check --out QC9_Relcheck > plinkoutput_QC9_PO_Relcheck
sleep 2

echo -e '
######################################################
'


awk '($5=="PO"  &&  ($7>0.1 ||  $8<0.8 || $9>0.1))'  QC9_Relcheck.genome > Not_related_POs.txt

awk '{print $1}' Not_related_POs.txt | sort -u > temp

grep -f temp QC9_Relcheck.genome > Not_related_fams.txt

echo -e '
######################################################
'

NnotPO=$(cat temp | wc -l);

echo -e "... there were $NnotPO not related parent offspring duos. Families are written to Not_related_fams.txt... \n"
wait

echo -e ">>> Deleting families with discordant relation status..... \n"

grep -f temp QC9.fam | awk '{print $1 " " $2}' > RmIndividuals_QC9.txt


if [ -e RmIndividuals_QC9.txt ]
then
	$PLINK --bfile QC9 --remove RmIndividuals_QC9.txt --make-bed --out QC10 > ../../LOG/Stage1_GenoQC/plinkoutput_QC9.1_rmNotPO.log
else
	mv QC9.bim QC10.bim
	mv QC9.bed QC10.bed
	mv QC9.fam QC10.fam
	mv QC9.log QC10.log
fi

echo -e '
######################################################
'

if [ -e QC10.bed ]
then
	echo -e ".....$NnotPO Samples were deleted. Samples written to RmIndividuals_QC9.txt \n"
else
	echo -e ".....there were some issues with deleting not related families. For more details please check plinkoutput_QC9.1_rmNotPO.log \n"
	echo -e "... program is terminated \n"
	exit
fi
wait
sleep 2


echo -e '
######################################################
'


echo -e "... calculating significance across families relations \n"

if [ -e QC10.bed ]
then
	$PLINK --bfile QC10 --genome --out QC10_Relcheck > plinkoutput_QC10_Relcheckround2.txt
	sleep 2
else
	echo -e "Please check the existence of QC10 files. Program exits \n"
	exit
fi
wait

echo -e '
######################################################
'


if [ -e  QC10_Relcheck.genome ]
then
    echo -e "....relcheck is succesfully finished.... \n"
else
    echo -e "...there were some major issues with your data. For details check plinkoutput_QC10_relcheck.log... \n"
fi
wait

echo -e '
######################################################
'


awk '($1!=$3 && $10>0.1)' QC10_Relcheck.genome > CrossfamRelations.txt
sleep 2
NCroFamsR=$(cat CrossfamRelations.txt | wc -l)

cat  CrossfamRelations.txt

echo -e "...There are $NCroFamsR 3rd degree or closer Relations across families. Families written to CrossfamRelations.txt. \n"

echo -e '
######################################################
'

echo -e ">>> You can check the output and rerun the analysis otherwise the program by default deleted one of the two families <<< \n"

awk '{print $1}' CrossfamRelations.txt | sort -u > temp

fams=$(cat temp | wc -l)

grep -f temp QC10.fam | awk '{print $1 " " $2}' > RmIndividuals_QC10.txt

if [ -e QC10.bed ]
then
	$PLINK --bfile QC10 --remove RmIndividuals_QC10.txt --make-bed --out QC11 > ../../LOG/Stage1_GenoQC/plinkoutput_QC11.1_rmRelFam.log
else
	mv QC10.bim QC11.bim
	mv QC10.bed QC11.bed
	mv QC10.fam QC11.fam
fi

if [ -e QC11.bed ]
then
	echo -e ".....$fams Families were deleted. Samples written to RmIndividuals_QC8.txt..... \n"
	cat temp
	else
	echo -e ".....there were some issues with deleting related families. For more details please check plinkoutput_QC8.1_rmNotPO.log"
fi


sleep 2

echo -e "
#########################################################
#	 	 Performing Final QC                    #
#########################################################
"

echo -e ".......Perform final QC before proceeding to next stage of the program........  \n"



echo -e ".......Checking again for callrate>95% genotyping-rate>95%  MAF>0.02 hwe>10e-8...... \n"


if [ -e QC11.bed ]
then
	$PLINK --bfile QC11  --geno 0.05 --hwe 10e-8 --maf 0.02 --make-bed --out tmp12
	sleep 2

echo -e '
######################################################
'


	$PLINK --bfile tmp12  --mind 0.05 --make-bed --out QC12 > ../../LOG/Stage1_GenoQC/plinkoutput_QC12.log
	sleep 2
else
	echo -e "Please check the existence of QC11 files. Program exits \n"
	exit
fi

echo -e '
######################################################
'

if [ -e QC12.irem ]
then
	read  irem filename <<< $(wc -l QC12.irem)
else
	irem=0
fi

echo -e '
######################################################
'

read SNPs filename <<< $(wc -l QC12.bim)

Diff=$(expr $snps - $SNPs)

echo -e "....." $irem  "samples were removed to to genotyping freq < 95%..... \n

	....." $Diff  "SNPS were removed due to call rates  < 95%; MAF <0.02 or significant (10e-8) deviation from HWE.... \n"

sleep 2

echo -e '
######################################################
'



read ind filename <<< $(wc -l QC12.fam)
read snps filename <<< $(wc -l QC12.bim)

echo -e "...The file after primary quality check has" $ind  "samples and" $snps  "SNPs... \n and is named as QC12 \n"

echo -e '
######################################################
'

if [ -e QC12.bed ]
then
	$PLINK --bfile QC12 --missing --out PostQC_report
echo -e '
######################################################
'
	$PLINK --bfile QC12 --hardy --out PostQC_Hardyreport
echo -e '
######################################################
'
	$PLINK --bfile QC12 --freq --out PostQC_AlleleFreq
echo -e '
######################################################
'
	$PLINK --bfile QC12 --het --out PostQC_Inbreeding
echo -e '
######################################################
  '
	$PLINK --bfile QC12 --check-sex --out Final_Sexcheck
echo -e '
######################################################
'

	QCFile=QC12
else
	echo -e "Please check the existence of QC12 files. Program exits \n"
	exit
fi
sleep 2

echo -e "
#######################################################
#							#
# 	     MDS plots and IBD analysis		#
#                                     		#
#######################################################
"

echo -e ">>> Provide Path to HAP Map files (bim fam bed coded) for ethnicity check <<<< \n"

hapmap=$Hapmapfile

if [ ! -e $hapmap.fam ] || [ ! -e $hapmap.bim ] || [ ! -e $hapmap.bed ]
then
	echo -e "...$hapmap files do not exist, the program will exit. Please download the hapmap files from the link provided in the user manual and rerun the analysis again... \n"
	exit
else
	echo -e "...$hapmap files found... \n"
fi

sleep 2

echo -e '
######################################################
'



datafileQC=QC12 #Study data

datafileHP=$hapmap


awk '{print $2}' $datafileQC.bim | sort -u > temp1
awk '{print $2}' $datafileHP.bim | sort -u  > temp2

comm -12 temp1 temp2 >overlapping.snps.txt
wait

rm temp*

echo -e '
######################################################
'

$PLINK --bfile $datafileQC --extract overlapping.snps.txt --alleleACGT --make-bed --out localQC #Making sure the allele coding is ACGT


sleep 2

echo -e '
######################################################
'

$PLINK --bfile $datafileHP --extract overlapping.snps.txt --alleleACGT --make-bed --out localHM  #Making sure the allele coding is ACGT

sleep 2

echo -e '
######################################################
'

#Print list of IID and FID for generating MDS plot


awk '{print $1 "\t" $2}' $datafileQC.fam > SampleInds.txt


awk '{print $2 " "  $1 }' localHM.bim > HM3CHR #Hapmap chromosomes
awk '{print $2 " " $4}' localHM.bim > HM3pos #Hapmap SNP base-pair positions
awk '{print $2 " " $3}' localHM.bim > HM3CM #Hapmap centi-morgan

echo -e "...Updating map, chr and cm of the study data... \n"

sleep 2

echo -e '
######################################################
'

$PLINK --bfile localQC  --update-map HM3pos --make-bed  --out localQCp #updated SNP base-pair position
wait
sleep 2

echo -e '
######################################################
'


echo -e "...Updated SNP positions of the study data based on Hapmap data...  \n"

$PLINK --bfile localQCp  --update-map HM3CHR --update-chr --make-bed  --out localQCpc #updated SNP base-pair position and chromosome of the study data
wait
sleep 2


echo -e '
######################################################
'

echo -e "...Updated CHR of the study data based on Hapmap data... \n"


$PLINK --bfile localQCpc  --update-map HM3CM --update-cm --make-bed  --out localQCpcc #updated SNP base-pair position,chromosome and cm of the study data
sleep 2

echo -e "...Updating CM of the study data based on Hapmap data... \n"

echo -e '
######################################################
'

if [ -e localQCpcc-merge.missnp ]
then
    missnps=$(cat localQCpcc-merge.missnp | wc -l)
    if [ $missnps -gt 0 ]
    then
	$PLINK --bfile localQCpc --exclude localQCpcc-merge.missnp --make-bed --out tmp  #exclude SNPs which fail to merge

echo -e '
######################################################
'

	$PLINK --bfile tmp  --update-map HM3CM --update-cm --make-bed  --out localQCpcc
    fi
fi
wait
sleep 2

rm tmp*

echo -e '
######################################################
'

datafile=localQCpcc

awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep -v ambig > localsnplist.txt

echo -e "...Check if there is any inconsistency between the coding of the alleles, note that only alleles with ACGT coding are accepted...\n"


echo -e '
######################################################
'

$PLINK --bfile localHM --extract localsnplist.txt --make-bed --out external   #Extract the SNPs from Hapmap file

sleep 2
echo -e '
######################################################
'

$PLINK --bfile localQCpcc --extract localsnplist.txt --make-bed --out internal  #Extract the SNPs from study data whose bp positions,chromosome and cm positions is updated.

sleep 2
echo -e '
######################################################
'

$PLINK --bfile internal --bmerge external.bed external.bim external.fam --make-bed --out MDSfile
sleep 2
echo -e '
######################################################
'


if [ -e MDSfile-merge.missnp ] #If there are SNPs with merging issues such as strand issues than flip those SNPs
then
	missnps=$(cat MDSfile-merge.missnp | wc -l)
	if [ $missnps -gt 0 ]
    	then
	$PLINK --bfile internal --flip  $LocSource/OUTPUT_DIR/Stage1_GenoQC/MDSfile-merge.missnp --make-bed --out internalf
echo -e '
######################################################
'

	$PLINK --bfile internalf --bmerge external.bed external.bim external.fam --make-bed --out MDSfile2 #Again try to merge the Hamap and study data file with the flipped SNPs
echo -e '
######################################################
  '

		if [[ -e MDSfile2-merge.missnp ]]
		then
		$PLINK --bfile internalf --exclude MDSfile2-merge.missnp --make-bed --out internalf
echo -e '
######################################################
  '

		$PLINK --bfile internalf --bmerge external.bed external.bim external.fam --make-bed --out MDSfile
echo -e '
######################################################
 '
		fi
	fi

fi
sleep 2

$PLINK --bfile MDSfile --cluster --mind .05 --mds-plot 4 --out HM3mds #Create an MDS
sleep 2


echo -e "...Generating all Pre and Post QC plots including MDS... \n"


$R --slave --no-save --args QC1_report.lmiss QC1_report.imiss PreQC_hardy.hwe PreQC_AlleleFreq.frq QC2_Sexcheck.sexcheck PreQC_Inbreeding.het PostQC_report.lmiss PostQC_report.imiss PostQC_Hardyreport.hwe PostQC_AlleleFreq.frq Final_Sexcheck.sexcheck PostQC_Inbreeding.het SampleInds.txt $Siteinfo HM3mds.mds < ../../Scripts/GenerateQC_plots.R

echo -e "...Stage 1 of Genotype QC is successfully completed... \n"

mv QC12.fam FinalQC_Study.fam
mv QC12.bed FinalQC_Study.bed
mv QC12.bim FinalQC_Study.bim

#Remove all the uncessary files

rm external*
rm localsnplist.txt
rm internal*
rm localQCp*
rm HM3CM*
rm HM3pos*
rm HM3CHR*
rm localHM*
rm overlapping.snps.txt
rm QC*
rm Mendel*

echo -e '
##################################################################### Start of Stage2 #####################################################################################
'

LocSource=$MAGNET

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




echo -e '
##################################################################### Start of Stage 3 #####################################################################################
'



LocSource=$MAGNET

source /$LocSource/ConfigFiles/Tools.config
source /$LocSource/ConfigFiles/Thresholds.config

cd OUTPUT_DIR/Stage3_GWAS

#Check if the user performed MAGNET Stage2 or has provided a quality checked file

if [ $PerformStage2=="Yes" ]
then
	echo -e "Stage2_Imputation provided raw files will be used \n"
	echo -e "
	#######################################
	# Stage 3: GWAS has started now       #
	#######################################
	"

	chmod 770 ../Stage2_GenoImpute/Data_SNPfile*.raw

	export pheno
	export phenofile
	export covars
	export fixed
	export Phenotype_Status

python <<END
import glob
import os
import optparse

parser = optparse.OptionParser()
parser.add_option("-f", "--factor", action="store", dest="fact", help="Factor X")
(options, args) = parser.parse_args()

inputfiles = glob.glob("../Stage2_GenoImpute/Data_SNPfile*.raw")

print (len(inputfiles)) #For python 3
#print  len(inputfiles) #For python 2

pheno2=os.environ["pheno"]
phenofile=os.environ["phenofile"]
covars=os.environ["covars"]
fixed=os.environ["fixed"]
Phenotype_Status=os.environ["Phenotype_Status"]

i = 1
for inputfile in inputfiles:
    jobname = "job_"+str(i)+"_factor.sh"
    f = open(jobname,"w")
    f.write("R --no-save --args ../Stage2_GenoImpute/Data_SNPfile"+str(i)+".raw "+str(phenofile)+" "+str(pheno2)+" "+str(covars)+" "+str(fixed)+" ResultsMerged"+str(pheno2)+"_"+str(i)+".txt "+str(Phenotype_Status)+" <  ../../Scripts/Code5_Regression.R")
    f.close()


    i = i+1
END

	chmod 770 job_*.sh

	njobs=$(ls -l ../Stage2_GenoImpute/Data_SNPfile*.raw|wc -l)


	for j in $(seq 12 12 $njobs)
	do
	lower=$(expr "$j" - 11)
	upper=$j
		for i in $(seq $lower $upper);do 
		./job_"$i"_factor.sh > prg"$i".out  2>&1 &
		sleep 1
		done
	wait
	done

	# The remaining jobs

	remainder=$(expr $njobs % 12)
	lower=$(expr $njobs - $remainder)
	for i in $(seq $lower $njobs)
	do
	./job_"$i"_factor.sh > prg"$i".out  2>&1 &
	sleep 1
	done	
	wait

echo -e "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- \n"

echo -e "
#######################################################################################
# 1)Post regression 				      				         #
#######################################################################################
"

	cat ResultsResults*> tmp.txt


	sed 's/\"//g' tmp.txt>tmp2.txt

	sed '/snp/d' tmp2.txt>tmp3.txt

	sed 's/X//g' tmp3.txt>tmp4.txt

	echo -e "SNP\tbeta\tse\tt\tP\tadj.pvalues" | cat - tmp4.txt > tmp5.txt

	awk {'print $1 "\t" $5'} tmp5.txt>SNP_Pvalue$pheno

	mv tmp5.txt RegressionOutput$pheno.txt



	if [ $ImputeData == "Yes" ]
	then
	ImputedData=../Stage2_GenoImpute/Merged_FinalQC_SNPs_Data
	else
	ImputedData=$SamplesToQC
	fi

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


echo -e "
######################################################################################
#2) Gene enrichment analysis(MAGMA)							#
######################################################################################
"

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

#inputfiles = glob.glob("../Stage2_GenoImpute/Merged_FinalQC_SNPs_Data*.fam")

#print len(inputfiles)



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


	for((i=1; i<=22;i++))
	do
	echo "start MAGMA Chromosome $i"
	./../../Scripts/MAGMAChr$i.sh &
	sleep 1
	done

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

else

#If only running Stage3_GWAS
######################################################################################
# 11) Split into chunks of 5000 and Find the minor allele frequency of the SNPS	#
######################################################################################

	cd ../Stage2_GenoImpute

	awk '{print $2}' $SamplesToQC.bim >AllSNPs.txt


	NSNPs=($(wc -l AllSNPs.txt)) #Total number of SNPs

	DivFiles=$[$NSNPs/$chunkSize] #Files after dividing into chunks of 5000

	if [[($(expr $NSNPs % $chunkSize) -gt 1)]]; then #If there is a remainder existing after dividing it into chunks of 5000
	NFiles=$(expr $DivFiles + 1) #If exists then add 1 to the number of main files
	else
	Nfiles=$DivFiles #Else let the number as it is
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
		$PLINK --bfile $SamplesToQC --extract ${Files[i]} --make-bed --out tmp"$i"
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

fi


echo -e '
##################################################################### Start of Stage 4 #####################################################################################
'


cd $MAGNET

LocSource=$MAGNET

source /$LocSource/ConfigFiles/Tools.config
source /$LocSource/ConfigFiles/Thresholds.config


OUTPUTDIR=$MAGNET/OUTPUT_DIR/Stage4_Enrichment/

if [[ ! -e $OUTPUTDIR ]] #Check if the installation directory is existing
then
    	mkdir $MAGNET/OUTPUT_DIR/
	cd $MAGNET/OUTPUT_DIR
	mkdir Stage4_Enrichment #If not create directory
	chmod 770 *
fi
wait


cd $MAGNET/OUTPUT_DIR/Stage4_Enrichment/


if [[ ! -e Output ]]
then
	echo "Creating folder for plots..."
	mkdir Output
	cd Output
	chmod 770 *
fi
wait


cd $MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output
cd ..


if [ $PerformStage3=="No" ]
then
	echo -e " The user provided Gene list would be used...\n"	
	if [ ! -z $GenelistProvided ]
	then
	Genelist=$GenelistProvided
	else
	echo -e "Gene list from Stage 3 will be used"
	fi
else
	Genelist=$GenelistFromStage3
fi


echo -e "Check if Genelist is in the right format"

awk '{print $1}' $Genelist> Genes


re='^[0-9]+$'


if [[ $Genes=$re ]]
	then
	echo "Genes are in correct format."
	else
	echo "Please check the format of Genes"
fi



awk '{print $1}' $KangUnivers> EG
awk '{print $1}' $Genelist>GL
cat EG GL| sort | uniq>$MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output/UniversalGenelist.txt #Kang Genes and MAGMA genes

UniversalGenes=$MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output/UniversalGenelist.txt

rm EG
rm GL


cd $MAGNET

if [[ ! -e tmpRlib ]]
	then
	mkdir -p tmpRlib
	chmod 775 tmpRlib
	else
	echo -e "R library already exists"
fi

R_LIBS="./tmpRlib"
export R_LIBS
R --no-save --slave --args $MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output/ $MAGNET/OUTPUT_DIR/Stage4_Enrichment/OutputDownstream $pheno $SourceData $Genelist $UniversalGenes< $MAGNET/Scripts/Downstream.R









echo -e '################################################End of MAGNET analysis######################################## \n'












