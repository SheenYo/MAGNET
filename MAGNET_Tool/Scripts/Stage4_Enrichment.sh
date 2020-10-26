 
MAGNET=$(pwd) 

cd $MAGNET

LocSource=$(pwd) 

source /$LocSource/ConfigFiles/Tools.config
source /$LocSource/ConfigFiles/Thresholds.config

cd $MAGNET/OUTPUT_DIR/Stage4_Enrichment/


OUTPUTDIR=$MAGNET/OUTPUT_DIR/Stage4_Enrichment/

if [[ ! -e $OUTPUTDIR ]] #Check if the installation directory is existing
	then
    	mkdir $MAGNET/OUTPUT_DIR/
	cd $MAGNET/OUTPUT_DIR
	mkdir Stage4_Enrichment #If not create directory
	chmod 770 *
fi
wait
 

if [[ ! -e Output ]]
	then
	echo "Creating folders for performing GOelite"
	mkdir Output
	cd Output
	mkdir goelite_input
	mkdir goelite_denom
	mkdir goelite_output
	chmod 770 *
	cd $MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output
	else
	cd $MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output
fi
wait



GO="goelite"

if [[ ! -e $GO_input ]]
	then
	mkdir $GO"_input"
	else
	echo "GO input Folder exists"
fi

if [[ ! -e $GO_denom ]]
	then
	mkdir $GO"_denom"
	else
	echo "GO denom Folder exists"
fi


if [[ ! -e $GO_output ]]
	then
	mkdir $GO_output
	else
	echo "GO output Folder exists"
fi

cd ..

echo -e "Please check the existence of Gene list, if you are only running Stage 4, please update the path for Gene list, else by default MAGNET will look into MAGNET/OUTPUT_DIR/Stage3_GWAS/"


echo -e "Check if Genelist is in the right format"

awk '{print $1}' $Genelist> Genes


re='^[0-9]+$'


if [[ $Genes=$re ]]
	then 
	echo "Genes are in correct format."
	else 
	echo "Please check the format of Genes"
fi



awk '{print $1 "\t" "L"}' $Genelist> tmp



echo -e "SourceIdentifier\tSourceCode" | cat - tmp > Output/goelite_input/GOElite.input.txt ## add header for GOElite
#Universal Genes =Kang + genes


awk '{print $1}' $KangUnivers> EG
awk '{print $1}' $Genelist>GL
cat EG GL| sort | uniq>Output/UniversalGenelist.txt #Kang Genes and MAGMA genes


awk '{print $1 "\t" "L"}' Output/UniversalGenelist.txt > tmp2
echo -e "SourceIdentifier\tSourceCode" | cat - tmp2 > Output/goelite_denom/GOElite.universe.txt ## add header for GOElite



##echo "run GOElite"


OutputF=$(pwd)"/Output/goelite_output/"
InputF=$(pwd)"/Output/goelite_input/"
DenomF=$(pwd)"/Output/goelite_denom/"

python $GOelite --update Official --species $GOeliteSpecies --version EnsMart62Plus
python $GOelite --species $GOeliteSpecies --input $InputF --denom $DenomF --output $OutputF 



cd $MAGNET


R_LIBS="./tmpRlib"

if [[ -e $R_LIBS ]]
	then 
	echo "R temproray library exists."
	else 
	echo "R temproray library does not exists, so creating one..."
	mkdir tmpRlib
fi

export R_LIBS
R --no-save --slave --args $MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output/goelite_output $MAGNET/OUTPUT_DIR/Stage4_Enrichment/OutputDownstream $pheno $SourceData $Genelist< $MAGNET/Scripts/Downstream.R

