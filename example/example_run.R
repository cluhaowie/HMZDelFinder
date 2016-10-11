############################################################
# Running the HMZDelFinder on 50 samples from 1000 genomes #
############################################################

# define working directory:
workDir <- getwd() 
# set project and data directory
mainDir <- paste0(workDir ,"/HMZDelFinder/"); if (!dir.exists(mainDir)){dir.create(mainDir)}
dataDir <- paste0(mainDir, "data/" , sep=""); if (!dir.exists(dataDir)){dir.create(dataDir)} # data directory
# Load the source code 
source("https://raw.githubusercontent.com/BCM-Lupskilab/HMZDelFinder/master/src/HMZDelFinder.R")

# Install missing packages
list.of.packages <- c("RCurl", "gdata", "data.table", "parallel", "Hmisc", "matrixStats")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
biocLitePackages <- c("DNAcopy", "GenomicRanges")
new.biocLitePackage <- list.of.packages[!(biocLitePackages %in% installed.packages()[,"Package"])]
if(length(new.biocLitePackage)) { source("http://bioconductor.org/biocLite.R"); biocLite(new.biocLitePackage)}

# load packages
library(RCurl)
library(data.table)
library(gdata)
library(parallel)
library(Hmisc)
library(matrixStats)
library(DNAcopy)
library(GenomicRanges)


# download RPKM data for 50 samples from 1000genomes
if (!file.exists(paste0(dataDir, "TGP/"))){
	if (file.exists(paste0(dataDir, "TGP.tar.gz")))file.remove(paste0(dataDir, "TGP.tar.gz"))
	dl_from_dropbox( paste0(dataDir, "TGP.tar.gz"), "6y14wftyhh6r2j0")
	untar(paste0(dataDir, "TGP.tar.gz"), exdir = dataDir)
}
# download BED file
if (!file.exists(paste0(dataDir, "tgp_hg19.bed"))){ 
	if (file.exists(paste0(dataDir, "tgp_hg19.bed.tar.gz"))){file.remove(paste0(dataDir, "tgp_hg19.bed.tar.gz"))}
	dl_from_dropbox( paste0(dataDir, "tgp_hg19.bed.tar.gz"), "1v5jbbm2r809ssy")
	untar(paste0(dataDir, "tgp_hg19.bed.tar.gz"), exdir = dataDir)
}

# set/create other paths and identifiers
bedFile <- paste0(dataDir, "tgp_hg19.bed") # set path to BED file
outputDir <- paste0(mainDir, "out/" , sep=""); if (!dir.exists(outputDir)){dir.create(outputDir)} # create output directory
plotsDir <- paste0(mainDir, "plots/" , sep=""); if (!dir.exists(plotsDir)){dir.create(plotsDir)} # create output plots directory
rpkmFiles <- dir(paste(dataDir, "TGP/",sep=""), "rpkm2.txt$")	# list of RPKM file names
rpkmFids <- gsub(".rpkm2.txt", "", rpkmFiles) 					# list of sample identifiers
rpkmPaths <- paste0(paste0(dataDir, "TGP/"), rpkmFiles) 		# list of paths to RPKM files
aohDir <- paste0(mainDir, "AOH/" , sep=""); if (!dir.exists(aohDir)){dir.create(aohDir)} 
aohRDataOut <- paste(mainDir, "AOH/extAOH_small.RData", sep="")	# temprary file to store AOH data


# In this example, we are not performing AOH filtering and VCF files are not required.
# If one wants to perform AOH filtering than have to provide:
# vcfPaths - the list of paths to VCF files
# vcfFids - the list of sample identifiers that corresponds to VCF files (same order)

#snpFiles <- dir (inputDir,"SNPs_Annotated.vcf.bz2$", recursive=T, include.dirs=FALSE)
#vcfFids <- sapply(strsplit(snpFiles,"[/\\.]"), function(x){x[2]})
#vcfPaths<- paste(inputDir,vcfNames,"/", sapply(strsplit(snpFiles,"/"), function(x){x[2]}), sep="")



########################################
# THRESHOLDS
#
# See description of HMZDelFinder function for details
########################################
is_cmg <- FALSE 		# only for CMG project - otherwhise use FALSE
lowRPKMthreshold <- 0.65# RPKM threshold  
maxFrequency <- 0.05	# max frequncy of HMZ deletion; default =0.005
minAOHsize <- 1000		# min AOH size
minAOHsig <- 0.45		# min AOH signal threshold
mc.cores=4 				# number of cores



# running HMZDelFinder
results <- runHMZDelFinder (NULL,		# vcfPaths - paths to VCF files for AOH analysis (not used for 1000 genomes) 
							NULL,		# vcfFids - sample identifiers corresponding to VCF files  (not used for 1000 genomes) 
							rpkmPaths, 	# paths to RPKM files 
							rpkmFids,	# samples identifiers corresponding to RPKM files
							mc.cores,	# number of CPU cores
							aohRDataOut,# temp file to store AOH data
							bedFile,	# bed file with target 
							lowRPKMthreshold, #  RPKM threshold 
							minAOHsize, # min AOH size
							minAOHsig,	# min AOH signal threshold
							is_cmg)		# flag used for CMG specific annotations; TRUE samples are from BHCMG cohort, FALSE otherwhise
					
					
# saving results in csv files
write.csv(results$filteredCalls, paste0(outputDir,"hmzCalls.csv"), row.names=F )

# plotting deletions
lapply(1:nrow(results$filteredCalls),function(i){
			plotDeletion (results$filteredCalls, i, results$bedOrdered, results$rpkmDtOrdered,  plotsDir, mainText=""  )})
	
## Selected columns from the results$filteredCalls object:					
#					Chr     Start      Stop   Genes Start_idx     FID
#					1:   5 140235634 140236833 PCDHA10    133937 NA11919
#					2:   X  47918257  47919256  ZNF630    167263 NA18856
#					3:  11   7817521   7818489   OR5P2     27561 NA19137
#					4:  11   7817521   7818489   OR5P2     27561 NA19236
#					5:   9 107379729 107380128  OR13C9    161704 NA19473
#					6:   1 196795959 196796135   CFHR1     15101 NA20798
#					7:   5  42629139  42629205     GHR    130161 NA07347
#					8:   5  42629139  42629205     GHR    130161 NA12342
#					9:   5  42629139  42629205     GHR    130161 NA19213
#					10:  16  55866915  55866967    CES1     64208 NA18553
## NOTE: Deletions of CES1, CFHR1 and OR13C9 are located within segmental duplications


					