#!/bin/bash -euo pipefail

function args()
{
	options=$(getopt -o svh --long help --long reads: --long genome: \
	--long cleaned: --long layout: --long saveIntermediateFiles --long saveIntermediateHisat2Bam \
	--long threads: --long maxForks: --long profile: -- "$@")
	[ $? -eq 0 ] || {
		echo "Incorrect option provided"
		exit 1
	}
	
	eval set -- "$options"
	while true; do
		case "$1" in
		-v)
			VERBOASE=1
			;;
		-h)
			HELP=1
			;;
		--help)
			HELP=1
			;;
		--reads)
			shift;
			reads=$1
			;;
		--genome)
			shift;
			genome=$1
			;;
		--cleaned)
			shift;
			cleaned=$1
			;;
		--layout)
			shift;
			layout=$1
			;;
		--threads)
			shift;
			threads=$1
			;;
		--maxForks)
			shift;
			maxForks=$1
			;;
		--profile)
			shift;
			profile=$1
			;;
		--saveIntermediateFiles)
			saveIntermediateFiles=1
			;;
		--saveIntermediateHisat2Bam)
			saveIntermediateHisat2Bam=1
			;;
		--)
			shift
			break
			;;
		esac
		shift
	done
				
}

HELP=0
reads=""
genome=""
maxForks=6
threads=8
saveIntermediateFiles=0
saveIntermediateHisat2Bam=0
layout="paired"
profile="docker"
cleaned="false"
args $0 "$@"


usage(){
echo -n "
 Usage:
 pipeOne.sh [options]* --reads <> --genome <>
 
 Require:
  --reads  <string>	input Fastq files, for example: \"/home/reads/*_{1,2}.fastq.gz\"
  --genome <string>	the genome profile your set in config files: conf/igenomes.config

 Optional:
  --cleaned <booloen> true or false 
  --layout <string>	paired or single. defualt [paired]
  --threads <int>	number of CPU process for each steps
  --maxForks <int>	max forks number of parrallel
  --profile <str>	execution envirenment. defualt [docker]
  --saveIntermediateFiles	save intermediate files defualt [off]
  --saveIntermediateHisat2Bam	save intermediate hisat2 BAM files. defualt [off]
  -h --help	print usage
  
"

}

echo "HELP: $HELP"
echo "genome: $genome"

if [ $HELP -eq 1 ]; then
	usage
	exit 1
fi


baseDir=$(dirname "$0")
workDir=$(pwd)



if [  "$genome" == '' ] || [ "$reads" == '' ]; then
	echo "--genome and --reads are required!"
	exit 1
	
fi


reads_basename=$(basename $reads)
reads="$(dirname $(realpath $reads))/${reads_basename}"
echo $reads



s1_Dir=${workDir}/s1_lncRNA
mkdir -p $s1_Dir; cd $s1_Dir
nextflow run ${baseDir}/s1_lncRNA.nf -resume -profile $profile --genome $genome --reads $reads --cleaned $cleaned

if [ "$cleaned" != "true" ]; then
	reads="../s1_lncRNA/results/fastp/clean/*.R{1,2}.fastp.fq.gz"
	cleaned="true"

fi

s2_Dir=${workDir}/s2_circRNA
mkdir -p $s2_Dir; cd $s2_Dir
#nextflow run ${baseDir}/s2_circRNA.nf -resume -profile $profile --genome $genome --reads $reads --cleaned $cleaned
s3_Dir=${workDir}/s3_APA-3TUR
mkdir -p $s3_Dir; cd $s3_Dir
nextflow run ${baseDir}/s3_APA-3TUR.nf -resume -profile $profile --genome $genome --reads $reads --cleaned $cleaned

s4_Dir=${workDir}/s4_retrotranscriptome
mkdir -p $s4_Dir; cd $s4_Dir
nextflow run ${baseDir}/s4_retrotranscriptome.nf -resume -profile $profile --genome $genome  --reads $reads --cleaned $cleaned

s5_Dir=${workDir}/s5_fusion
mkdir -p $s5_Dir; cd $s5_Dir
nextflow run ${baseDir}/s5_fusion.nf -resume -profile $profile --genome $genome --reads $reads --cleaned $cleaned

s6_Dir=${workDir}/s6_rnaEditing
mkdir -p $s6_Dir; cd $s6_Dir
nextflow run ${baseDir}/s6_rnaEditing.nf -resume -profile $profile --genome $genome --reads $reads --cleaned $cleaned

s7_Dir=${workDir}/s7_alternative_splicing
mkdir -p $s7_Dir; cd $s7_Dir
# nextflow run ${baseDir}/s7_alternative_splicing.nf -resume -profile $profile --genome $genome --reads $reads --cleaned  $cleaned

s8_Dir=${workDir}/s8_SNP
mkdir -p $s8_Dir; cd $s8_Dir
nextflow run ${baseDir}/s8_SNP.nf -resume -profile $profile --genome $genome --bam "../s7_alternative_splicing/results/star2pass/*.bam" 
