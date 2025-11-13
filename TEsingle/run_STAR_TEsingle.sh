#!/bin/sh -l
#SBATCH -J run_STAR_TEsingle
#SBATCH -e %x-%j.err
#SBATCH -o %x-%j.out
#SBATCH --export=ALL
#SBATCH --partition=fn_mediuam
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=10
#SBATCH -t 5-0:0:0

THREADS=10
MAXNUM=100
ANCHOR=200
MISMATCH=999
MISMATCH_LMAX="0.04"

if [ -z "$3" ]; then
    echo "sbatch/sh run_STAR_TEsingle.sh [STAR index] [R1] [R2]" >&2
    exit 1
fi

GENOME="$1"
UMI="$2"
INPUT="$3"

WHITELIST="${GENOME}/barcode_whitelist.txt"
GENEGTF="${GENOME}/T2T_TEsingle_gene.gtf"
TEGTF="${GENOME}/T2T_TEsingle_TE.gtf"

if [ ! -f "${WHITELIST}" ]; then
    echo "Barcode whitelist not found in STAR index folder." >&2
    echo "Please run setup_STAR_index.sh to obtain all necesary files" >&2
    exit 1
fi

CURRDIR=$PWD
ANCHOR=$((MAXNUM + 50))

FILEBASE=`basename $INPUT`
UMIBASE=`basename $UMI`
BASE=`basename $FILEBASE \.gz`
BASE=`basename $BASE \.fastq`
BASE=`basename $BASE \.fq`
BASE=`basename $BASE _R2`
OUTDIR="${CURRDIR}/${BASE}"

if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

ln -s "$INPUT" "${CURRDIR}/${FILEBASE}"
ln -s "$UMI" "${CURRDIR}/${UMIBASE}"

cd $OUTDIR

CMD="STARsolo --genomeLoad NoSharedMemory --outSAMunmapped None --outFilterType BySJout --outSAMstrandField intronMotif --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --sjdbScore 1"

CMD="$CMD --soloType CB_samTagOut --soloCBmatchWLtype 1MM --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outSAMattributes NH HI AS nM CR CY UR UY CB GX GN sS sQ sM --soloFeatures GeneFull"

CMD="$CMD --genomeDir $GENOME --readFilesIn ${CURRDIR}/${FILEBASE} ${CURRDIR}/${UMIBASE} --outFilterMultimapNmax $MAXNUM --winAnchorMultimapNmax $ANCHOR --runThreadN $THREADS --outFilterMismatchNmax $MISMATCH --outFilterMismatchNoverReadLmax $MISMATCH_LMAX"

CMD="$CMD --soloCBwhitelist $WHITELIST --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000"

GZIP=`basename $FILE1 \.gz`

if [ "$GZIP" != "$FILEBASE" ]; then
    CMD="$CMD --readFilesCommand 'zcat'"
fi

CMD="$CMD --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 --limitBAMsortRAM 40000000000"

CMD="$CMD --soloUMIlen 12"

${CMD}

if [ $? -ne 0 ]; then
    echo "Error encountered during mapping" >&2
    exit 1
else
    mv Log.final.out "${BASE}_STAR_mapping.log"
    mv Aligned.*.bam "${BASE}_STAR_10x.bam"
    samtools index "${BASE}_STAR_10x.bam"
    echo "Mapping completed" >&2
fi

cd ..

BAM="${OUTDIR}/${BASE}_STAR_10x.bam"
BASE=$(basename ${BAM} \.bam)
BASE="${BASE}_TEsingle"

STRAND="forward"
MINUMI=100
LOG="${BASE}.log"

CMD="TEsingle -b ${BAM} --GTF $GENEGTF --TE $TEGTF --stranded $STRAND --threads 10 --cutoff ${MINUMI} --project ${BASE}" > ${LOG}

${CMD}

if [ $? -ne 0 ]; then
    echo "Error with TEsingle" >&2
    exit 1
else
    echo "Done"
fi
