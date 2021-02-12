#!/bin/bash -euo pipefail
version=1.0.0
function args()
{
	options=$(getopt -o svh --long help --long version --long update_GTF --long scratch: --long reads: --long genome: \
	--long cleaned: --long layout: --long saveIntermediateFiles: --long saveIntermediateHisat2Bam \
	--long threads: --long maxForks: --long profile: --long library:  -- "$@")
	[ $? -eq 0 ] || {
		echo "Incorrect option provided"
		usage
		exit 1
	}
	
	eval set -- "$options"
	while true; do
		case "$1" in
		-v)
			VERBOASE=1
			;;
		--version)
			VERBOASE=1
			;;
		-h)
			HELP=1
			;;
		--help)
			HELP=1
			;;
		--update_GTF)
			
			update_GTF=1
			;;
		--scratch)
			shift;
		   scratch=$1
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
		--library)
			shift;
			library=$1
			;;
		--profile)
			shift;
			profile=$1
			;;
		--saveIntermediateFiles)
			shift;
			saveIntermediateFiles=$1
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
VERBOASE=0
update_GTF=0
reads=""
genome=""
maxForks=2
threads=8
saveIntermediateFiles=false
saveIntermediateHisat2Bam=0
scratch="false"
layout="paired"
profile="docker"
cleaned="false"
library="polyA"
args $0 "$@"

# color--------------------------------
NOCOLOR='\033[0m'
RED='\033[0;31m'
GREEN='\033[0;32m'
ORANGE='\033[0;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
LIGHTGRAY='\033[0;37m'
DARKGRAY='\033[1;30m'
LIGHTRED='\033[1;31m'
LIGHTGREEN='\033[1;32m'
YELLOW='\033[1;33m'
LIGHTBLUE='\033[1;34m'
LIGHTPURPLE='\033[1;35m'
LIGHTCYAN='\033[1;36m'
WHITE='\033[1;37m'

usage(){
echo -e -n "${RED}Usage:${NOCOLOR}
  pipeOne.sh [options]* --reads <> --genome <>
 
${RED}Require:${NOCOLOR}
  --reads  <string>	input Fastq files, for example: \"/home/reads/*_{1,2}.fastq.gz\"
  --genome <string>	the genome profile you set in config files: conf/igenomes.config

${RED}Optional:${NOCOLOR}
  --cleaned	<boolean>	true or false. defualt [true]
  --layout	<string>	paired or single. defualt [paired]
  --library	<string>	polyA or total. defualt [polyA]
  --threads	<int>	number of CPU process for each steps. default [8]
  --maxForks	<int>	max forks number of parrallel. default [2]
  --profile	<str>	execution envirenment. defualt [docker]
  --saveIntermediateFiles 	<boolean>	true or false, save intermediate files . defualt [false]
  --update_GTF 	<boolean>	<boolean>	true or false, use customize GTF generated in step s1.1_lncRNA.nf instand of GENCODE GTF as input for step: s1.5_fusion.nf and s1.7_alternative_splicing.nf .defualt [false] 
  -h --help	print usage
"
exit 1
}

if [ $HELP -eq 1 ]; then
	usage
fi

if [ $VERBOASE -eq 1 ]; then
	echo "version: $version"
	exit 1
fi

SCRIPT=`realpath $0`
baseDir=$(dirname $SCRIPT)
workDir=$(pwd)

if [  "$genome" = '' ] || [ "$reads" = '' ]; then
	echo -e "${RED}--genome${NOCOLOR} and ${RED}--reads${NOCOLOR} are required!\n"
	usage
fi

#reads_basename=$(basename $reads)
#reads="$(dirname $(realpath $reads))/${reads_basename}"
echo -e "\n${RED}--genome:${NOCOLOR} ${GREEN}$genome${NOCOLOR}"
echo -e "${RED}--reads:${NOCOLOR} ${GREEN}\"$reads\"${NOCOLOR}\n"

bash_input="--profile $profile \
--genome $genome  \
--layout $layout \
--library $library \
--threads $threads  --maxForks $maxForks \
--saveIntermediateFiles $saveIntermediateFiles \
--update_GTF $update_GTF \
--scratch $scratch \
--cleaned $cleaned --reads \"$reads\""

echo "bash $0  $bash_input" >one_command.sh

set_args_base="-resume -profile $profile --genome $genome \
--layout $layout --threads $threads --maxForks $maxForks \
--saveIntermediateFiles $saveIntermediateFiles --scratch $scratch "


s1_Dir=${workDir}/s1.1_lncRNA
mkdir -p $s1_Dir; cd $s1_Dir
set_args="$set_args_base --cleaned $cleaned  --reads $reads"
nextflow run ${baseDir}/s1.1_lncRNA.nf $set_args

if [ "$cleaned" != "true" ]; then
	reads="../s1.1_lncRNA/results/fastp/clean/*.R{1,2}.fastp.fq.gz"
	cleaned="true"
fi

set_args="$set_args_base --cleaned $cleaned  --reads $reads"

if [ "$library" = "total" ]; then
	s2_Dir=${workDir}/s1.2_circRNA
	mkdir -p $s2_Dir; cd $s2_Dir
	nextflow run ${baseDir}/s1.2_circRNA.quant.nf $set_args --hisat2_bam "../s1.1_lncRNA/results/hisat2/bam/*.hisat2.sortbycoordinate.bam" --update_GTF $update_GTF
fi

s3_Dir=${workDir}/s1.3_APA-3TUR
mkdir -p $s3_Dir; cd $s3_Dir
nextflow run ${baseDir}/s1.3_APA-3TUR.nf $set_args

s4_Dir=${workDir}/s1.4_retrotranscriptome
mkdir -p $s4_Dir; cd $s4_Dir
nextflow run ${baseDir}/s1.4_retrotranscriptome.nf $set_args

s5_Dir=${workDir}/s1.5_fusion
mkdir -p $s5_Dir; cd $s5_Dir
nextflow run ${baseDir}/s1.5_fusion.nf $set_args --update_GTF $update_GTF

s6_Dir=${workDir}/s1.6_rnaEditing
mkdir -p $s6_Dir; cd $s6_Dir
nextflow run ${baseDir}/s1.6_rnaEditing.nf $set_args

s7_Dir=${workDir}/s1.7_alternative_splicing 
mkdir -p $s7_Dir; cd $s7_Dir
nextflow run ${baseDir}/s1.7_alternative_splicing.nf $set_args --update_GTF $update_GTF

s8_Dir=${workDir}/s1.8_SNP
mkdir -p $s8_Dir; cd $s8_Dir
nextflow run ${baseDir}/s1.8_SNP.nf  $set_args_base --bam "../s1.7_alternative_splicing/results/star2pass/*.bam"

table_dir="${workDir}/00_tables"
mkdir -p $table_dir; cd $table_dir
conda_base=`conda info --base`
source ${conda_base}/etc/profile.d/conda.sh
set +u; conda activate pipeOne_ml; set -u
python3 ${baseDir}/bin/summary_table.py check_tables $library
python3 ${baseDir}/bin/summary_table.py mark_feature $library
conda deactivate