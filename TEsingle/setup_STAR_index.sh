#!/bin/sh -l
#SBATCH -t 12:0:0
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=3G
#SBATCH -J STAR_genomeGen
#SBATCH -e %x-%j.err
#SBATCH -o %x-%j.out
#SBATCH --export=ALL

if [ -z "$1" ]; then
    echo "Usage: sbatch/sh $0 [output folder]" >&2
    exit 1
fi

FOLDER="$1"

if [ ! -d "${FOLDER}" ]; then
    mkdir ${FOLDER}
fi

cd ${FOLDER}

FAURL="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
GTFURL="https://www.dropbox.com/scl/fi/wqjuuhs3wfiztx491gu83/T2T_TEsingle_gene.gtf.gz?rlkey=0w78ul7p3t789uucg8nsudxyn&st=x3qyntu9&dl=1"
TEURL="https://www.dropbox.com/scl/fi/xzutkue61wtnx2r2yr7y5/T2T_TEsingle_TE.gtf.gz?rlkey=6w14hslzfvd7d7smni9uan683&st=tqhfwea3&dl=1"
WLURL="https://www.dropbox.com/scl/fi/1qjh5i2pqelznz9gyv3dp/barcode_whitelist.txt.gz?rlkey=vi4z117kayv8xb8hnl65rulsv&st=x2wrbe2d&dl=1"

DTOOL="wget"
if ! command -v ${DTOOL} &> /dev/null
then
    DTOOL="curl"
    if ! command -v ${DTOOL} &> /dev/null
    then
       echo "Neither wget or curl are installed. Please ensure a command-line tool for downloading from URL is installed" >&2
       exit 1
    fi
fi

LOG="T2T_STAR_index_generation.log"
if [ -f "${LOG}" ]; then
    rm ${LOG}
fi

FASTA="T2T_CHM13v2.fa"

if [ ! -f "${FASTA}" ]; then
    if [ "${DTOOL}" == "wget" ]; then
	wget -O ${FASTA}.gz "${FAURL}" 2>>${LOG}
    else
	curl -o ${FASTA}.gz "${FAURL}" 2>>${LOG}
    fi
    gunzip ${FASTA}.gz
fi

GTF="T2T_TEsingle_gene.gtf"
if [ ! -f "${GTF}" ]; then
    if [ "${DTOOL}" == "wget" ]; then
	wget -O ${GTF}.gz "${GTFURL}" 2>>${LOG}
    else
	curl -o ${GTF}.gz "${GTFURL}" 2>>${LOG}
    fi
    gunzip ${GTF}.gz
fi

TE="T2T_TEsingle_TE.gtf"
if [ ! -f "${TE}" ]; then
    if [ "${DTOOL}" == "wget" ]; then
	wget -O ${TE}.gz "${TEURL}" 2>>${LOG}
    else
	curl -o ${TE}.gz "${TEURL}" 2>>${LOG}
    fi
    gunzip ${TE}.gz
fi

WL="barcode_whitelist.txt"
if [ ! -f "${WL}" ]; then
    if [ "${DTOOL}" == "wget" ]; then
	wget -O ${WL}.gz "${WLURL}" 2>>${LOG}
    else
	curl -o ${WL}.gz "${WLURL}" 2>>${LOG}
    fi
    gunzip ${WL}.gz
fi


if ! command -v STAR &>/dev/null
then
    echo "STAR could not be found. Please ensure it is installed and in the PATH variable" >&2
    exit 1
fi

CMD="STAR --runMode genomeGenerate --limitGenomeGenerateRAM 39000000000 --runThreadN 10"

CMD="${CMD} --genomeDir ${FOLDER}  --genomeFastaFiles ${FASTA}"

CMD="${CMD} --sjdbGTFfile ${GTF} --sjdbOverhang 100"

${CMD} 2>>${LOG}

if [ $? -ne 0 ]; then
    echo "Error with STAR run. See ${LOG} for details" >&2
    exit 1
else
    echo "Done"
fi
