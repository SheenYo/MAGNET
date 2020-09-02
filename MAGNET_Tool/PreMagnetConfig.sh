#!/bin/bash


#Echo sourcing config and threshold files


chmod -R +xr *
CurrentDir=$(pwd) 

MAGNET=$(pwd) 

chmod -R +xr $MAGNET*


source /$MAGNET/ConfigFiles/Tools.config
source /$MAGNET/ConfigFiles/Thresholds.config

pwd

echo -e '

							     	######################################################################
							      			    	       WELCOME TO MAGNET \n                        

								

											 ||||||||||||      ||||||||||||
											 ||||||||||||      ||||||||||||
											 ||||||||||||      ||||||||||||
 											 |          |      |          |
 											 |          |      |          |
  											 |          |      |          |
											 |          |      |          |
 											 |          ||||||||          |
											 |                            |
 											 |                            |
  											  |                          |
   											   |                        |
   											    |                      |
     											     ||||||||||||||||||||||
				       													     
							        							     		     
					                         		MApping the Genetics of NEuropsychological Trait \n	     
							                     		      v1.0 (2018-8-16) \n			     
			            				                                                                           
																	     
							        #######################################################################


#################################################################
# CHECK IF ALL REQUIRED PROGRAMMES/TOOLS/PACKAGES ARE INSTALLED #
#################################################################
'
sleep 1

echo -e '
#########################################
#	Performing Basic Checks		#
#########################################
'

wget -q --spider http://google.com

if [ $? -eq 0 ]; then
    echo "...Internet is working... \n"
else
    touch "ConnectionError.txt"
    echo "...Check your internet connection, and then try again... The programme will exit now\n">>ConnectionError.txt
    
exit
fi


if [[ $OSTYPE != "linux-gnu" ]]
then
	echo -e "...As some of the tools are only linux based, please run MAGNET on an Linux based computer or contact us... \n"
else
	echo -e "...OS type checked... \n"
fi
wait

sleep 1

if (uname -m|grep "64"); then
	echo -e "...Machine type checked... \n"
else
	echo -e "...Please run MAGNET On a 64 bit Linux machine... \n"
fi
wait 

sleep 1

echo -e '
#########################
# Creating directories	#
#########################
'
if [ ! -e $MAGNET/MAIN_DIR ]
then
	echo -e "...Creating MAIN directory... \n"
	mkdir $MAGNET/MAIN_DIR
	touch $MAGNET/MAIN_DIR/"Warnings.txt"
else
	echo -e "...MAIN directory exists... \n"
fi

sleep 1 


if [ ! -e $MAGNET/INPUT_DIR ]
then
	echo -e "...Creating INPUT directory... \n"
	mkdir $MAGNET/INPUT_DIR
else
	echo -e "...INPUT directory exists... \n"
fi

sleep 1
wait 


if [ ! -e $MAGNET/OUTPUT_DIR ]
then
	echo -e "...Creating OUTPUT directory... \n"
	mkdir $MAGNET/OUTPUT_DIR
	cd $MAGNET/OUTPUT_DIR
	mkdir Stage1_GenoQC
	mkdir Stage2_GenoImpute
	mkdir Stage3_GWAS
	mkdir Stage4_Enrichment
else
	echo -e "...OUTPUT directory exists... \n"
fi

#cd ..

sleep 1 
wait 

if [ ! -e $MAGNET/LOG ]
then
	echo -e "...Creating LOG directory... \n"
	mkdir $MAGNET/LOG
	cd $MAGNET/LOG
	mkdir Stage1_GenoQC
	mkdir Stage2_GenoImpute
	mkdir Stage3_GWAS
	mkdir Stage4_Enrichment
else
	echo -e "...LOG directory exists... \n"
fi

wait 
sleep 1 

echo -e '
#########################################
#  Performing Tools availibility Checks	#
#########################################
'



tools_installed=""
tools_installed_ind=1


echo -e '
##########
#GO-elite#
##########
'

echo -e "...Checking for the Program GOelite, if doesn't exist then it will be downloaded... \n"
	if hash python GOelite 2>/dev/null
	then	
	echo -e "...GO-elite is installed... \n"
	GOelite=$(which GOelite)
	echo GOelite=$GOelite >> $MAGNET/ConfigFiles/Tools.config
	else
	echo -e "...GOelite does not exist... \n Downloading GOelite..."
	#Create main GOelite directory
	mkdir $MAGNET/MAIN_DIR/GOelite
	chmod +xr *
	#Change to GOelite main directory
	cd $MAGNET/MAIN_DIR/GOelite
	#Download GOelite
	#wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/go-elite/GO-Elite_v.1.2.5-Ubuntu_10.04.zip --no-check-certificate
	wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/go-elite/GO-Elite_v.1.2.5-Py.zip --no-check-certificate
	unzip GO-Elite_v.1.2.5-Py.zip
	chmod +xr *

	chmod +xr $MAGNET/MAIN_DIR/GOelite/GO-Elite_v.1.2.5-Py/GO_Elite.py
	#$MAGNET/MAIN_DIR/GOelite/GO-Elite_v.1.2.5-Py/GO_Elite.py > $MAGNET/LOG/GO_Elite.log 2>&1 
	testp=$MAGNET/MAIN_DIR/GOelite/GO-Elite_v.1.2.5-Py/GO_Elite.py | grep 'GO_Elite'
		if [[ ! $testp = 0 ]]
		then
		echo -e 			"\n ...GOelite successfully installed... \n"
		echo GOelite="$MAGNET/MAIN_DIR/GOelite/GO-Elite_v.1.2.5-Py/GO_Elite.py" >> $MAGNET/ConfigFiles/Tools.config
		else
		echo -e "...Error during GO_Elite installation... \n"
		tools_installed_ind=0
		tools_installed=$tools_installed $testp
		fi
	
	fi	

sleep 1
wait 

echo -e '
##########
#liftOver#
##########
'
echo -e "Checking for the Program liftOver, if doesn't exist then it will be downloaded \n"

	if hash liftOver 2>/dev/null
	then 
	echo -e 				"\n ...liftOver is installed... \n"
	liftOver=$(which liftOver)
	echo liftOver=$liftOver >> $MAGNET/ConfigFiles/Tools.config
	else
	echo -e "...You require liftOver but it's not installed, or the path is not correct. \n Downloading liftOver..."
	#Create main liftOver directory
	mkdir $MAGNET/MAIN_DIR/liftOver
	#Change to liftOver main directory
	cd $MAGNET/MAIN_DIR/liftOver
	#Download liftOver
	wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver --no-check-certificate
	#The $MAGNET/MAIN_DIR/liftOver/linux.x86_64/liftOver
	# when you enter the path in the info file, you must replace the value of
	# $MAGNET with the actual path, e.g.
	chmod +xr $MAGNET/MAIN_DIR/liftOver/liftOver
	$MAGNET/MAIN_DIR/liftOver/liftOver > $MAGNET/LOG/liftover.log 2>&1 
	testp='$MAGNET/MAIN_DIR/liftOver/liftOver | grep  'liftOver''
		if [[ ! $testp = 0 ]]
		then
		echo -e 			"\n ...liftOver successfully installed... \n"
		echo "liftOver=$MAGNET/MAIN_DIR/liftOver/liftOver" >> $MAGNET/ConfigFiles/Tools.config
		else
		echo -e "...Error during liftOVer installation... \n"
		tools_installed_ind=0
		tools_installed=$tools_installed $testp
		fi
	
	fi

echo " "
echo " "
echo " "
echo " "

sleep 1 
wait 
echo '
#######
#MAGMA#
#######
'

echo -e "...Checking for the Program MAGMA, if doesn't exist then it will be downloaded... \n"
	if hash magma 2>/dev/null
	then
	echo -e "MAGMA is installed \n"
	magma=$(which magma)
	echo magma=$magma >> $MAGNET/ConfigFiles/Tools.config
	else
	echo -e "You require Magma but it's not installed, or the path is not correct. \n Downloading Magma..."
	#Create main MAGMA directory
	mkdir $MAGNET/MAIN_DIR/magma
	#Change to magma main directory
	cd $MAGNET/MAIN_DIR/magma
	#Download magma
	wget https://ctg.cncr.nl/software/MAGMA/prog/magma_v1.07b.zip --no-check-certificate
	unzip magma_v1.07b.zip
	chmod +xr $MAGNET/MAIN_DIR/magma/magma_v1.07b/magma
	$MAGNET/MAIN_DIR/magma/magma_v1.07b/magma > $MAGNET/LOG/magma.log 2>&1 
	testp='$MAGNET/MAIN_DIR/magma/magma_v1.07b/magma | grep 'magma''
		if [[ ! $testp = 0 ]]
		then
		echo -e 			"\n ...Magma successfully installed... \n"
		echo "magma=$MAGNET/MAIN_DIR/magma/magma_v1.07b/magma" >> $MAGNET/ConfigFiles/Tools.config
		else
		echo "...Error during magma installation..."
		tools_installed_ind=0
		tools_installed=$tools_installed $testp
		fi
	
fi		
echo " "
echo " "
echo " "
echo " "

sleep 1 
wait 

echo '
##########
# Plink  #
##########
'
echo -e "...Checking for the Program plink v1.09 exists in your PATH variable.If this version of program doesn't exist then it will be downloaded... \n"
	if hash plink1.9 2>/dev/null
	then 
	echo -e "...plink exists... \n"
	PLINK=$(which plink1.9)
	echo PLINK=$PLINK >> $MAGNET/ConfigFiles/Tools.config
	else
	echo -e "...You require Plink but it's not installed, or the path is not correct. \n Downloading plink..."
	#Create main plink directory
	mkdir $MAGNET/MAIN_DIR/PLINK
	#Change to plink main directory
	cd $MAGNET/MAIN_DIR/PLINK
	#Download the plink package
	wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20181202.zip --no-check-certificate
	unzip plink_linux_x86_64_20181202.zip
	chmod +xr $MAGNET/MAIN_DIR/PLINK/plink
	$MAGNET/MAIN_DIR/PLINK/plink >> $MAGNET/LOG/PLINK.log 2>&1 
	testp='$MAGNET/MAIN_DIR/PLINK/plink | grep  'Purcell''
		if [[ ! $testp = 0 ]]
		then
		echo -e 			"\n ...plink successfully installed... \n"
		echo "PLINK=$MAGNET/MAIN_DIR/PLINK/plink" >> $MAGNET/ConfigFiles/Tools.config
		else
		echo "...Error during plink installation... \n"
		tools_installed_ind=0
		tools_installed=$tools_installed $testp
		fi
	
fi
		
echo " "
echo " "
echo " "
echo " "

sleep 1
wait 

echo '
##########
#Minimac3#
##########
'
	echo -e "...Downloading minimac..."
	#Create main minimac directory
	mkdir $MAGNET/MAIN_DIR/minimac
	#Change to minimac directory
	cd $MAGNET/MAIN_DIR/minimac
	#Download minimac tool package
	wget ftp://share.sph.umich.edu/minimac3/Minimac3Executable.tar.gz --no-check-certificate
	#Untar the downloaded package
	tar -xzvf Minimac3Executable.tar.gz
	#change directory 
	cd Minimac3Executable/
	chmod +xr *
	./Minimac3 > $MAGNET/LOG/minimac.log 2>&1
	#Try minimac
	itest='./Minimac3 | grep 'Minimac3''
		if [[ ! $itest = 0 ]]
		then
		echo -e 			"\n ...Minimac successfully installed..."
		echo "minimac=$MAGNET/minimac/Minimac3Executable/bin/Minimac3" >> $MAGNET/ConfigFiles/Tools.config
		else
		echo "Error during the minimac installation \n"
		tools_installed_ind=0
		tools_installed=$tools_installed $itest
		fi
	



echo " "
echo " "
echo " "
echo " "

sleep 1 
wait 

echo '
#########
#   R   #
#########
'
echo -e "Checking for R if doesn't exist then it will be downloaded \n"

	if hash R 2>/dev/null
	then
	echo "R is installed"
	R=$(which R)
	echo R=$R >> $MAGNET/ConfigFiles/Tools.config
	else	
	echo -e "You require R but it's not installed, or the path is not correct. \n Please download R >=3.5. and then rerun the analysis">> $MAGNET/MAIN_DIR/Warnings.txt
	fi
	
sleep 1 
wait 

echo '
#########
#SHAPEIT#
#########
'
echo -e "...Checking for the Program SHAPEIT v2.727 exists in your PATH variable.If this version of program doesn't exist then it will be downloaded \n"
	if hash shapeit 2>/dev/null
	then 
	echo -e "...shapeit is installed... \n"
	SHAPEIT=$(which shapeit)
	echo SHAPEIT=$SHAPEIT >> $MAGNET/ConfigFiles/Tools.config
	else
	echo -e "...SHAPEIT2 version v2.727 does not exist \n So Downloading shapeit..."
	mkdir $MAGNET/MAIN_DIR/SHAPEIT
	#Change to SHAPEIT directory
	cd $MAGNET/MAIN_DIR/SHAPEIT
	#Download SHAPEIT tool package
	wget  https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.12.linux.tar.gz --no-check-certificate
	#Untar the downloaded package
	tar -zxvf shapeit.*.tar.gz
	#Try SHAPEIT
	chmod +xr $MAGNET/MAIN_DIR/SHAPEIT/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit
	$MAGNET/MAIN_DIR/SHAPEIT/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit >$MAGNET/LOG/SHAPEIT.log 2>&1 
	shap=$MAGNET/MAIN_DIR/SHAPEIT/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit
		if [[ ! $shap = 0 ]]
		then
			echo -e "\n ...Shapeit successfully installed... \n"
			echo "SHAPEIT=$MAGNET/MAIN_DIR/SHAPEIT/shapeit.v2.904.2.6.32-696.18.7.el6.x86_64/bin/shapeit" >> $MAGNET/ConfigFiles/Tools.config
		else			
			echo -e "Unable to install shapeit \n"
			tools_installed_ind=0
			tools_installed=$tools_installed $shap

		fi
	
	fi

echo " "
echo " "
echo " "
echo " "

echo $tools_installed


if [ $tools_installed_ind = 0 ]
then
	echo "Following tools not installed properly.Check the LOG directory in the install directory for more information">> $MAGNET/MAIN_DIR/Warnings.txt
	echo $tools_installed
else
	echo "...............................................All required tools are installed............................................."
fi



echo -e '
##############################################
#  CHECK IF ALL REQUIRED FILES ARE EXISTING  #												
##############################################
'

echo -e '
########################
#      SampleToQC      #
########################
'

echo -e 'SampleToQC Files \n'

if [ -e $SamplesToQC.fam ]
then
    echo ".....$SamplesToQC.fam exists"
    convert="false"
else
    if [ -e $SamplesToQC.ped ]
    then
	echo "..... $SamplesToQC.ped exists"
	convert="true"
    else
	echo ".....$SamplesToQC.fam or $SamplesToQC.ped  does not exist; program exits">> $MAGNET/MAIN_DIR/Warnings.txt
	exit;
    fi
fi
wait
sleep 1

if [ $convert ==  "true" ]
then
    if [ -e $SamplesToQC.map ]
       then
	  #plink --file $SamplesToQC --make-bed --out QC1.1
	echo "..... $SamplesToQC.map exists"
    else
	echo "$SamplesToQC.map file does not exist. program exits">> $MAGNET/MAIN_DIR/Warnings.txt
	exit;
    fi
else
    if [ -e $SamplesToQC.bim ]
    then
	echo -e ".....$SamplesToQC.bim exists"
    else
	echo -e".....$SamplesToQC.bim does not exist; program exits">> $MAGNET/MAIN_DIR/Warnings.txt
	exit;
    fi

    if [ -e $SamplesToQC.bed ]
    then
	echo -e ".....$SamplesToQC.bed exists"
    else
	echo -e".....$SamplesToQC.bed does not exist; program exits">> $MAGNET/MAIN_DIR/Warnings.txt
	exit;
    fi
    #cp $SamplesToQC.bim QC1.1.bim
    #cp $SamplesToQC.fam QC1.1.fam
    #cp $SamplesToQC.bed QC1.1.bed
fi
wait
sleep 1

echo -e '
######################
#Annotation File hg19#
######################									
'

if [ -e $hg19SNPs ]
	then
	echo -e "............Annotation File hg19 exists \n"
	else
	echo -e "............Annotation File hg19 does not exists \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait
sleep 1 

echo -e '
###################################
#Checking Chain files For liftOver#
###################################
'


if [ -e $Chain16To19 ]
	then
	echo -e "............Chain file exists for hg16Tohg19 \n"
	else
	echo -e "............Chain file exists for hg16Tohg19 does not exists, please refer to manual on how to download it \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait

if [ -e $Chain17To19 ]
	then
	echo -e "............Chain file exists for hg17Tohg19 \n"
	else
	echo -e "............Chain file exists for hg17Tohg19 does not exists, please refer to manual on how to download it \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait

if [ -e $Chain18To19 ]
	then
	echo -e "............Chain file exists for hg18Tohg19 \n............"
	else
	echo -e "............Chain file exists for hg18Tohg19 does not exists, please refer to manual on how to download it \n............">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait
sleep 1 

echo -e '
###################################
#       Map files For Shapeit     #
###################################
'


if [ -e $MapFile ]
	then
	echo -e "............Map file for Shapeit do exists............ \n"
	else
	echo -e "............Map file does not exists for Shapeit, please refer to manual on how to download it \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait

if [ -e $ShapeitRefHaps ]
	then
	echo -e "............Haps file for Shapeit do exists............ \n"
	else
	echo -e "............Haps file does not exists for Shapeit, please refer to manual on how to download it............ \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait

if [ -e $ShapeitRefLegend ]
	then
	echo -e "............Legend file for Shapeit do exists............ \n"
	else
	echo -e "............Legend file does not exists for Shapeit, please refer to manual on how to download it............ \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait

if [ -e $ShapeitRefSample ]
	then
	echo -e "............Sample file for Shapeit do exists............ \n"
	else
	echo -e "............Sample file does not exists for Shapeit, please refer to manual on how to download it............ \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait

if [ -e $AnnotationFile ]
	then
	echo -e "............Annotation file do exists............ \n"
	else
	echo -e "............Annotation file does not exists, please refer to manual on how to download it............ \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait
sleep 1 

echo -e '
###################################
#	  Magma Files 		  #
###################################
'

if [ -e $MagmaSNPloc ]
	then
	echo -e "............Magma SNP Loc file exists............ \n"
	else
	echo -e "............Magma SNP Loc file does not exist, please refer to manual on how to download it............ \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait

if [ -e $MagmaGeneloc ]
	then
	echo -e "............Magma Gene Loc file exists............ \n"
	else
	echo -e "............Magma Gene Loc file does not exists, please refer to manual on how to download it............ \n">> $MAGNET/MAIN_DIR/Warnings.txt
fi
wait
sleep 1

echo -e '
###########
#   PERL  #
###########
'

PERL=$(which perl)


echo PERL=$PERL >> $MAGNET/ConfigFiles/Tools.config

if [ !  -f $PERL ]
    then
    echo -e "............no PERL found $PERL............"
    exit 1	
    else
    echo -e "............PERL found............ \n"
fi
wait

if [ !  -x $PERL ]
    then
    echo -e "............PERL found but it is not executable $PERL............">> $MAGNET/MAIN_DIR/Warnings.txt
    exit 1
fi
wait
sleep 1 

echo -e '
#############
#  PYTHON   #
#############
'

PYTHON=$(which python)

echo "PYTHON=$PYTHON">> $MAGNET/ConfigFiles/Tools.config


if [ !  -f $PYTHON ]
    then
    echo -e "............no PYTHON found $PYTHON............ \n"
    exit 1	
    else 
    echo -e "............Python found, dont panic :-D \n"

    # Check for Python major version. The MAGNET pipeline expects Python2. I fixed some function calls to allow it to work with
    #    Python3, but this is not fully tested, so we warn and suggest using Python2.
    PYTHON3_STRING=$(python --version | fgrep 'Python 3')
    if [ -n "$PYTHON3_STRING" ]; then 
        echo "WARNING: your python interpreter '$PYTHON' is Python v3, which is untested with MAGNET. If you experience problems, set PYTHON to a python2 interpreter in 'Tools.config'."
    fi
fi

echo -e '
##########
#   SH   #
##########
'

sh=$(which bash)

echo "sh=$sh">> $MAGNET/ConfigFiles/Tools.config


if [ !  -f $SH ]
    then
    echo -e "............no SH found $SH............"
    exit 1	
    else
    echo -e "............Shell found............"
fi
wait

if [ !  -x $SH ]
    then
    echo -e "............SH found but it is not executable $SH............">> $MAGNET/MAIN_DIR/Warnings.txt
    exit 1
fi
wait
sleep 1 

echo '
###################################################
#	  Check the job management system	  #
###################################################
'
cd $MAGNET
touch Magnet.sh

if [[ $ENV == "SLURM" ]]
then
	echo -e "#!/bin/bash \n
#SBATCH --ntasks=$tasks_slurm \n
#SBATCH --nodes=$nodes_slurm  \n
#SBATCH --mem-per-cpu=$memory_slurm \n
#SBATCH --time=$time_slurm \n
#SBATCH --partition=$partition_slurm \n
#SBATCH --output=Magnet.out \n
#SBATCh --hint=compute_bound \n
#SBATCH --mail-type=$mail_slurm \n" >> $MAGNET/Magnet.sh
elif [[ $ENV == "PBS" ]]
then
	echo -e " #!/bin/bash \n
#PBS -l procs=$tasks_pbs \n
#PBS -l nodes=$nodes_pbs  \n
#PBS -l pmem=$memory_pbs \n
#PBS -l walltime=$time_pbs \n
#PBS -l q=$partition_pbs \n
#PBS -l o=$outputname \n
#PBS -l m=$mail_pbs \n" >> $MAGNET/Magnet.sh
else
	if [[ $ENV == "Local" ]]
	then
	echo -e "#Running Script locally \n" >> $MAGNET/Magnet.sh
	echo -e "The pipeline requires a server (SLURM or PBS) to perform Stage2 and Stage3, please let us know if you have none of them available \n">> $MAGNET/MAIN_DIR/Warnings.txt
	fi
fi

echo -e "Based on the parameter set in the Config file, the requested job management will be selected, by default "$ENV" will be set"
sleep 1 

echo -e "------------------------------------------------------------------------------------------------------------------------------- \n "

echo -e "Based on the parameter set in the Config file, the requested code will be selected to run, by default "Stage1_GenoQC.sh" will run"

if [[ $Script == "CompleteScript.sh" ]]
then 
	echo -e "
./Scripts/CompleteScript.sh \n" >> $MAGNET/Magnet.sh
elif [[ $Script == "Stage1_GenoQC.sh" ]]
then 
	echo -e "
./Scripts/Stage1_GenoQC.sh \n" >> $MAGNET/Magnet.sh
elif [[ $Script == "Stage2_Imputation.sh" ]]
then
	echo -e "
./Scripts/Stage2_Imputation.sh \n" >> $MAGNET/Magnet.sh
elif [[ $Script == "Stage3_GWAS.sh" ]]
then 
	echo -e "
./Scripts/Stage3_GWAS.sh \n" >> $MAGNET/Magnet.sh
elif [[ $Script == "Stage4_Enrichment.sh" ]]
then 
	echo -e "
./Scripts/Stage4_Enrichment.sh \n" >> $MAGNET/Magnet.sh
fi



echo -e '
#################################
# 	Running MAGNET now      #
#################################
'
chmod -R +xr $MAGNET/*

cd $MAGNET

sbatch Magnet.sh












