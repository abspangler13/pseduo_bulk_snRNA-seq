# pseudoBulk

basically, our current workflow involves aligning reads using [cell/space]ranger, which gives you one BAM file per library, and then your gene x barcode matrix of counts
some of our initial exploratory data analysis of chromium data suggested that there were way more reads outside of the 3' end of genes than expected (edited) 
we are also less interested in single cells/nuclei and more interested in specific cell populations, in general
so it would be good to explore how useful the chromium data is outside of the 3' ends of genes
like if it captures transcript-specific exons or junctions
to do this, you can take the cell populations defined by @mattntran in the human data which should be finalized now
and then you can split the BAM files by the barcodes corresponding to each cell population and then filter by the UMI sequence within each population
some code from visium was: 

```/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/get_barcode_lists.R```

and 

```/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/split_bam_layer.sh``` that you can adapt to the chromium data (edited) 
you have to uncomment out some of that code for the script to run each step
then you can merge the counts across samples, adapting this script: 

```/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/prep_counts.R```  
one thing i hadnt done was make coverage / bigwig files of the pseudobulked profiles, which could be useful for visualizing some of the non-3' signal
especially at genes we think are important/interesting
and to assess the consistency of non-uniformity of coverage across multiple samples to see if there are sequence-specific biases driving the non-3' signal
actually the coverage is in there now, sorry...i also changed commenting so that the entire thing would run
@Sang Ho was looking at the output of this process on visium data which is analogous , but this way you could learn more of the command line tools, but then work together a bit on the downstream analyses (edited) 
the first thing he was going to be doing was comparing the pseudobulked gene counts with featureCounts made here to the original cell ranger outputs we are using in the respective papers

Splitting NAc BAM files:

ran this file line by line on the JHPCE cluster in R.
```/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/bam_pseudobulk_abby/get_barcodes_NAc.Rmd```

produced these two output files:
barcode list of cells
```/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/bam_pseudobulk_abby/barcode_level_cellType_NAc_map.tsv```
sample table 
```/dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/bam_pseudobulk_abby/sample_level_cellType_NAc_map.tsv```

Then run this file on JHPCE:
```qsub /dcl01/ajaffe/data/lab/singleCell/10x_pilot/premRNA/bam_pseudobulk_abby/split_bam_cellType.sh```


052920 Moved a bunch of stuff over from the bam_pseudobulk_abby dir to my home dir. Have to make bam files again and edit sh file to contain correct paths. 


