import glob
import os
import optparse
from __future__ import print_function

parser = optparse.OptionParser()
parser.add_option("-f", "--factor", action="store", dest="fact", help="Factor X")
(options, args) = parser.parse_args()

inputfiles = glob.glob("../Stage2_GenoImpute/Data_SNPfile*.raw")

print(len(inputfiles))


i = 1
for inputfile in inputfiles:
	jobname = "job_"+str(i)+"_factor.sh"	
	f = open(jobname,"w")
	f.write("R --no-save --args ../Stage2_GenoImpute/Data_SNPfile"+str(i)+".raw ../RefData/Complete_DEdata_ForRegression.txt JA Age,Sex Site ResultsMergedJA_"+str(i)+".txt<  Code5_Regression.R")
	f.close()
	
#	os.system("qsub "+jobname)

	i = i+1
