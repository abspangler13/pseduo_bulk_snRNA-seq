###
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=32G,h_fsize=100G
#$ -pe local 1
#$ -N bam_split
#$ -o ./logs/split_cellType.txt
#$ -e ./logs/split_cellType.txt
#$ -m e
#$ -t 1-41
#$ -tc 20

module load samtools
module load python/3.6.9

SUB=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/subset-bam-1.0-x86_64-linux/subset-bam

BAMFILE=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/bam_pseudobulk_abby/sample_level_cellType_NAc_map.tsv
BARFILE=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/bam_pseudobulk_abby/barcode_level_cellType_NAc_map.tsv

## get bam and sample info
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $9}' $BAMFILE | awk "NR==${SGE_TASK_ID}")
CELLTYPE=$(awk 'BEGIN {FS="\t"} {print $8}' $BAMFILE | awk "NR==${SGE_TASK_ID}")
BAM=$(awk 'BEGIN {FS="\t"} {print $10}' $BAMFILE | awk "NR==${SGE_TASK_ID}")

mkdir -p /dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/$SAMPLE/CellType/ ##why don't I see the new dir

## get out barcodes
BC_FILE=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/$SAMPLE/CellType/${SAMPLE}_${CELLTYPE}_barcodes.txt
awk -v a=$SAMPLE -v b=$CELLTYPE '$10==a && $9==b {print $1}' $BARFILE > $BC_FILE ##produces empty file

# ## subset
NEWBAM=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/${SAMPLE}/CellType/${SAMPLE}_${CELLTYPE}.bam

$SUB --bam $BAM --cell-barcodes $BC_FILE --cores 6 --out-bam $NEWBAM

# ## index
samtools index $NEWBAM 

# ## dedup
NEWBAM_DEDUP=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/${SAMPLE}/CellType/${SAMPLE}_${CELLTYPE}_dedup.bam

umi_tools dedup --umi-tag=UB --cell-tag=CB --temp-dir=$TMPDIR --method=unique \
	--extract-umi-method=tag --stdin=$NEWBAM --stdout=$NEWBAM_DEDUP
	
samtools index $NEWBAM_DEDUP

## feature counts- genes
module unload conda_R

GTF=/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
OUTGENE=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/${SAMPLE}/CellType/${SAMPLE}_${CELLTYPE}.genes.counts
OUTEXON=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/${SAMPLE}/CellType/${SAMPLE}_${CELLTYPE}.exons.counts

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts \
		-a $GTF -o $OUTGENE $NEWBAM_DEDUP
# exons	
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts  \
	-O -f -a $GTF -o $OUTEXON $NEWBAM_DEDUP


# junctions	
OUTJXN=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/${SAMPLE}/CellType/${SAMPLE}_${CELLTYPE}.junctions.bed
OUTCOUNT=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/${SAMPLE}/CellType/${SAMPLE}_${CELLTYPE}.junctions.count

module load python/2.7.9

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/regtools/build/regtools junctions extract -i 9 -o ${OUTJXN} ${NEWBAM_DEDUP}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/bed_to_juncs_withCount < ${OUTJXN} > ${OUTCOUNT}

module load ucsctools
BW=/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/${SAMPLE}/CellType/${SAMPLE}_${CELLTYPE}
python ~/.local/bin/bam2wig.py -s /dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Annotation/hg38.chrom.sizes.cellRanger.hg38 \
	-i $NEWBAM_DEDUP -t 4000000000 -o $BW
rm $BW.wig
