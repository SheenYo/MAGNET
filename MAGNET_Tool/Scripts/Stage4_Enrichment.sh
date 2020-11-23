
MAGNET=$(pwd)

cd $MAGNET

LocSource=$(pwd)

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
