import glob
import os
from __future__ import print_function

inputfiles = glob.glob("../OUTPUT_DIR/Stage2_GenoImpute/All_Affected_DataSNPs*.fam")

print(len(inputfiles))

i = 1
for inputfile in inputfiles:
	jobname = "phase2_"+str(i)+".sh"	
	f = open(jobname,"w")

	f.write("../../MAIN_DIR/minimac/Minimac3Executable/bin/Minimac3 --refHaps ../../RefData/ALL_1000G_phase1integrated_v3_impute/1000GP_Phase3/ALL.chr"+str(i)+".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --vcfReference --rs --haps gwas_Data.chr"+str(i)+".vcf --doseOutput --hapOutput --chr "+str(i)+" --prefix Gwas.Chr"+str(i)+"_Study.Imputed.Output")

	f.close()
	
#	os.system("qsub "+jobname)

	i = i+1





