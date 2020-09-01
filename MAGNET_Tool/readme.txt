##Welcome to MAGNET##

#For detailed instructions please read the user manual of MAGNET.

#Please make sure you have R installed >=3.5, gunzip unzip utility for linux and then run "PreMagnetConfig.sh" before running MAGNET

#Please download the reference files for MAGNET. Commands to download the reference data needed for each stage are provided below as shell script, please make sure you install them in the #“RefData” folder within MAGNET: 

#Stage 1 QC

#Hapmap plink data
	cd RefData

	wget http://zzz.bwh.harvard.edu/plink/dist/hapmap_r23a.zip
	unzip -n hapmap_r23a.zip 

#and site information

	wget http://zzz.bwh.harvard.edu/plink/dist/hapmap.pop

	echo -e "FID IID Population" | cat - hapmap.pop > Hapmap_siteinfo.txt

#Stage 2 Imputation

# 1) liftOver Chain files
#Depending on the genome build you can download one from the following chain files, by default #hg18tohg19 is set, in case you have a different genome build to be lifted please change the #ChainToChoose #variable in the MAGNET/ConfigFiles/Tools.config e.g if you want to liftover from #hg38 to hg19:

# ChainToChoose=$Chain38To19

	cd RefData

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg16/liftOver/hg16ToHg19.over.chain.gz
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg17/liftOver/hg17ToHg19.over.chain.gz
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz

#2) Hg19SNPs: Hg19 SNPs file containing SNPs rsids, starting and ending bp position, chromosome #number. The user can individually download the files and merge them using the following linux #code: 
       	
	cd RefData

	for ((i =1;i<=22;i++))
	do
	wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/BED/bed_chr_"$i".bed.gz
	done
	chmod 770 *
	gunzip -n bed_chr_*.bed.gz
	cat bed_chr_*.bed>SNPs_all.bed
	rm bed_chr_*

#3) Genetic map file

	cd RefData

	mkdir genetic_map_b37
	cd genetic_map_b37
	for ((i =1;i<=22;i++))
	do
	wget 	http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/genetic_map_chr"$i"_combined_b37.txt
	done
	cd ..

	
#4) Shapeit reference files

	cd RefData
	mkdir ALL_1000G_phase1integrated_v3_impute
	cd ALL_1000G_phase1integrated_v3_impute
	wget https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz
	tar -xvzf 1000GP_Phase3.tgz
	
	cd 1000GP_Phase3

	for ((i =1;i<=22;i++))
	do
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr"$i".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
	done
	cd ../..
	
#5) Annotation file
	

	cd RefData

	mkdir AnnotationFiles_SNP142
	cd AnnotationFiles_SNP142
	wget 	http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp142.txt.gz
	gunzip snp142.txt.gz
	awk '{print $2 " " $4 " " $5}' snp142.txt> Annotation_file.txt # Extract chr, bp and snp id
	
	for((i=1;i<=22;i++))
 	do
 	grep "chr"$i" """ Annotation_file.txt>Annotation_chr"$i".txt
	wait
	awk -F, '{NF=4}1' OFS="\t" Annotation_chr"$i".txt > Annotation_chr"$i"_extraC.txt
	wait
	sed -e 's/chr//g'  Annotation_chr"$i"_extraC.txt|awk '{print $1 " " $2 " " $3 " " $4 " " 	$1":"$2}'>Comp_chr"$i".txt
	wait
	echo $'chr bp_pos snp_name chr:bp' | cat - Comp_chr"$i".txt>Comp_chr"$i"_wHead.txt
	done

	cd ..

#Stage 3 GWAS
	
#MAGMA reference file

	cd RefData
	wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip
	unzip -n g1000_eur.zip

	rm -r g1000_eur.zip 

#MAGMA gene location file
	cd RefData
	wget https://ctg.cncr.nl/software/MAGMA/aux_files/NCBI37.3.zip
	unzip -n NCBI37.3.zip

