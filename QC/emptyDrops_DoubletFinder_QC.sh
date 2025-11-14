#!/bin/sh -l
#SBATCH -t 5-0:0:0
#SBATCH --mem-per-cpu=100G
#SBATCH -J emptyDrops_DoubletFinder
#SBATCH -e %x-%j.err
#SBATCH -o %x-%j.out
#SBATCH --export=ALL

if [ -z "$1" ]; then
    echo "sbatch $0 [TEsingle mtx file]" >&2
    echo " Assumes cbcs and annots file in same folder" >&2
    exit 1
fi

SCRIPTDIR=$(dirname $0)
SCRIPTDIR="${SCRIPTDIR}/src"

MTX="$1"
MTX=$(realpath ${MTX})
ANN=$(echo ${MTX} | sed 's/mtx$/annots/')
CBC=$(echo ${MTX} | sed 's/mtx$/cbcs/')
DIR=$(basename ${MTX} \.mtx)

if [ ! -d "${DIR}" ]; then
    mkdir ${DIR}
else
    echo "${DIR} already exists. Files within will be overwritten" >&2
fi

ln -sf ${MTX} ${DIR}/matrix.mtx
ln -sf ${ANN} ${DIR}/genes.tsv
ln -sf ${CBC} ${DIR}/barcodes.tsv

SCRIPT="${SCRIPTDIR}/snRNA/emptyDrops_DoubletFinder_QC.R"

LOG="${DIR}_QC.log"

echo "Rscript ${SCRIPT} ${DIR}"

Rscript ${SCRIPT} ${DIR} > ${LOG}

if [ $? -ne 0 ]; then
    echo "Error with QC step. Check the log file" >&2
    exit 1
else
    echo "Done"
fi

