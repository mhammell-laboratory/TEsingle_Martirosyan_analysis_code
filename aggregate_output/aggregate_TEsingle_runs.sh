#!/bin/sh -l
#SBATCH -t 12:0:0
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=1G
#SBATCH -J aggregate_TEsingle
#SBATCH -e %x-%j.err
#SBATCH -o %x-%j.out
#SBATCH --export=ALL

if [ -z "$2" ]; then
    echo "
    Usage: sbatch/sh $0 [output folder name] 
                        [TEsingle annots files ...]
           Assumes cbcs and mtx files are in the same folder
"
fi

SCRIPTDIR=$(dirname $0)
SCRIPT="${SCRIPTDIR}/src/TEsingle_aggregate.pl"
OUTDIR="$1"
shift

LOG=${OUTDIR}_aggregation.log

perl ${SCRIPT} ${OUTDIR} $@ 2> ${LOG}

if [ $? -ne 0 ]; then
    echo "Error with aggregating runs" >&2
    exit 1
else
    echo "Done"
fi

   
