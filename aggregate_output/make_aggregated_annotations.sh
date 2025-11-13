#!/bin/sh -l
#SBATCH -t 12:0:0
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1G
#SBATCH -J annotation_aggregation
#SBATCH -e %x-%j.err
#SBATCH -o %x-%j.out
#SBATCH --export=ALL

if [ -z "$3" ]; then
    echo"
    Usage: sbatch/sh $0 [aggregated output folder]
                        [Martirosyan library annotations]
                        [QC output file(s)]
"
    exit 1
fi

SCRIPTDIR=$(dirname $0)
SCRIPTDIR="${SCRIPTDIR}/src"

FOLDER="$1"
if [ -f "${SCRIPTDIR}/Martirosyan_library_annotations.txt" ]; then
    ANNOT="${SCRIPTDIR}/Martirosyan_library_annotations.txt"
else
    ANNOT="$2"
fi
shift 2

SCRIPT="${SCRIPTDIR}/annotate_barcode.pl"
OUTFILE="Martirosyan_barcode_annotations.csv"
BC="${FOLDER}/barcodes.tsv.gz"

perl ${SCRIPT} ${ANNOT} ${BC} > ${OUTFILE}

if [ $? -ne 0 ]; then
    echo "Error with annotating barcodes" >&2
    exit 1
fi

SCRIPT="${SCRIPTDIR}/combine_QC_output.pl"
OUTFILE="Martirosyan_ED_DF_output.csv"
LIBORDER="${FOLDER}/libOrder.txt"

perl ${SCRIPT} ${LIBORDER} $@ > ${OUTFILE}

echo "Done"

