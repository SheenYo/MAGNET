 
MAGNET=$(pwd) 

cd $MAGNET

LocSource=$(pwd) 

source /$LocSource/ConfigFiles/Tools.config
source /$LocSource/ConfigFiles/Thresholds.config

cd $MAGNET/OUTPUT_DIR/Stage4_Enrichment/

mkdir -p Output
cd Output

mkdir goelite_input
mkdir goelite_denom
mkdir goelite_output

cd ..
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

$PYTHON $GOelite --update Official --species $GOeliteSpecies --version EnsMart62Plus
$PYTHON $GOelite --species $GOeliteSpecies --input $InputF --denom $DenomF --output $OutputF 



cd $MAGNET


mkdir -p tmpRlib 
chmod 775 tmpRlib 
R_LIBS="./tmpRlib"
export R_LIBS
R --no-save --slave --args $MAGNET/OUTPUT_DIR/Stage4_Enrichment/Output/goelite_output $MAGNET/OUTPUT_DIR/Stage4_Enrichment/OutputDownstream $pheno $SourceData $Genelist< $MAGNET/Scripts/Downstream.R


