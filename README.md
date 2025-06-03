# OMICS ✨

This repository houses conceptual perspectives, coding practice, assignment/competition solutions based on the materials from a variety of computational biology/bioinformatics courses, workshops, technical manuals, academic articles, and others. 

*Last updated: 3 Jun 2025*

## Features
* [Single cell RNA-seq data analysis](#Analyze-single-cell-RNA-seq-data)
* [Bulk RNA-seq data analysis](#Analyze-bulk-RNA-seq-data)
* [ATAC-seq data analysis](#Analyze-ATAC-seq-data)

## Technical procedures
### Analyze single cell RNA-seq data
- [Introduction to single cell RNA-seq](https://github.com/hbctraining/Intro-to-scRNAseq/blob/master/schedule/links-to-lessons.md)
- [2025 single cell genomics day](https://satijalab.org/scgd25/)
- [*Frontier*](Science_Reading/scRNAseq.md) single-cell RNA-seq analytical progress
- [Single-cell eQTL](https://github.com/XFWuCN/scQTLtools)
- [COVID-19 RNA-seq data resources](https://github.com/ScienceComputing/COVID-19-RNA-Seq-datasets)
- A mini scRNA-seq [pipeline](SingleCellRNASeq/Scanpy/pipeline.py) | [Why apply pipelines?](https://github.com/SciComp8/NGSOmics_Programming/blob/main/SingleCellRNASeq/Scanpy/Why_Pipeline.md)
- If given raw `bcl` files, we [convert them to fastq files](FastQC/bcl_to_fastq.sh)
- As inputs are `fastq` files, we can ...
  - Run FastQC to [evaluate sequence quality and content](FastQC/Run_FastQC.sh)
  - Use Trim Galore to [trim reads](FastQC/Trim_Read.sh) if we spot unexpected low-quality base calls/adaptor contamination
  - Re-run FastQC to [re-evaluate sequence quality and content](FastQC/Run_FastQC.sh)
  - If single-cell RNA-seq data is generated from the plate-based protocol, we can ...
    - Use STAR to perform alignment and FeatureCounts to generate the count matrix
  - Else if single-cell RNA-seq data is generated from the droplet-based protocol, we can ...
    - Use kb-python package to perform [pseudo sequence alignment and generate the count matrix](SingleCellRNASeq/kb-python)
    - Use Cell Ranger pipelines to perform [sequence alignment and generate the count matrix](SingleCellRNASeq/CellRanger/cellranger_count.sh)
- After having the `feature-barcode matrices` at hand, we can ...
  - Use Scanpy workflow to perform [quality assurance, cell clustering, marker gene detection for cell identities](SingleCellRNASeq/Scanpy/PBMC), and [trajectory inference](SingleCellRNASeq/Scanpy/Bone_Marrow)
  - Use Seurat workflow to perform quality assurance, cell clustering, and marker gene detection for cell identities [case 1](SingleCellRNASeq/Seurat/scRNAseq_analysis_full.Rmd) | [case 2](SingleCellRNASeq/Seurat/SkinCell.Rmd)
    - If we observe the factor-specific clustering and want cells of the same cell type cluster together across single/multiple confounding factors, we can use canonical correlation analysis or Harmony (suitable for complicated confounding effects) to [integrate](SingleCellRNASeq/Seurat/scRNAseq_analysis_full.Rmd) cells
    - We can leverage SingleR or ScType to partially or fully automate cell-type identification
      - Other options of automating cell-type identification by mapping to references and then transfering labels: scArches, Symphony
  - Use Bioconductor packages to [perform single cell RNA-Seq data analysis](SingleCellRNASeq/Bioconductor/BioconductorSkinCell.Rmd)
  - Generate [pseudobulk](SingleCellRNASeq/Scanpy/Pseudobulk.py), which aggregates the gene expression levels specific to each cell type within an individual
  - Perform pseudobulk-based differentially gene expression analysis in [edgeR](SingleCellRNASeq/Scanpy/scRNAseq_DE_Part1.ipynb) or [DESeq2](SingleCellRNASeq/Bioconductor/Pseudobulk_DE.Rmd)
  - Use bulk RNAseq-based pathway analysis tools (e.g., clusterProfiler, GSEA, GSVA) or single cell RNAseq-based Pagoda2 to evaluate if a predefined set of genes shows statistically significant and consistent variations between biological conditions
  - Use scGen to model [perturbation responses](SingleCellRNASeq/Perturbation/scGen)  

<hr>

### Analyze bulk RNA-seq data
  - [COVID-19 RNA-seq data resources](https://github.com/ScienceComputing/COVID-19-RNA-Seq-datasets)
  - Run [FastQC](FastQC/Run_FastQC.sh) or [fastp](FastQC/Run_fastp.sh) to evaluate sequence quality and content
  - [Recommend] Use splice-aware genome aligner STAR to [align the reads](BulkRNASeq/STAR_Align.sh)
      - Other splice-aware alignment tool options: Olego, HISAT2, MapSplice, ABMapper, Passion, BLAT, RUM ...
      - Other alignment tools that disregard isoforms: BWA, Bowtie2 ...
  - Use Rsubread to [align the reads](BulkRNASeq/AlignmentCountingTCell.Rmd)
    - **Why align?** To pinpoint the specific location on the human genome from which our reads originated
  - Use Qualimap to perform [quality assurance](BulkRNASeq/Qualimap_QC.sh) on the aligned reads
  - Use [MultiQC](BulkRNASeq/multiqc_QC.sh) to harmonize all QC and alignment metadata from FastQC, STAR, Qualimap, and other [tools](https://multiqc.info/modules/)
  - Use GenomicAlignments for aligned reads to [obtain the gene-level or exon-level quantification](BulkRNASeq/AlignmentCountingTCell.Rmd)
  - Use featureCounts for aligned reads to [count the fragments](BulkRNASeq/featureCounts.sh)
  - [Recommend] Use Salmon for unaligned reads to [obtain the transcript-level quantification](BulkRNASeq/Salmon_quant.sh)
    - **Why unalign?** To speed up the counting process of reads
    - Next step: Use tximport to aggregate transcript-level quantification to the gene level
  - [Perform differential gene expression analysis](BulkRNASeq/DEAnalysisTCell.Rmd)
  - [Perform principal component analysis, heatmap, and clustering](BulkRNASeq/PCAHeatmapClusteringTissue.Rmd)
  - [Perform gene set enrichment analysis](BulkRNASeq/GeneSetTCell.Rmd)
  - Achieve cell-type resolution in bulk RNA-Seq through deconvolution techniques (*Under Active Construction*)

<hr>

### Analyze ATAC-seq data
  - [Practical guide](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3)
  - Run [ENCODE ATAC-seq pipeline](https://github.com/ScienceComputing/atac-seq-pipeline/blob/master/README.md) to perform alignment, quality assurance, peaking calling, and signal track generation
  - If we're interested in inspecting every step in each analytical phase, or even leveraging advanced/unique features of other tools that the current pipeline ignores, ...
    - For alignment and post-alignment phases, we can ...
      - Use Rsubread or Rbowtie2 to [align the fastq files relative to hg19/hg38/hs1](ATACSeq/AlignFASTQ.Rmd)
      - Use GenomicAlignments and GenomicRanges to perform post-alignment processing including [reading properly paired reads, estimating MapQ scores/insert sizes, reconstructing the full-length fragment, and others](ATACSeq/PostAlignment.Rmd)
      - Use ATACseqQC to perform [comprehensive ATAC-seq quality assurance](ATACSeq/ATACseqQC.Rmd)
    - For TSS analysis phase, we can ...
      - Use soGGi to [assess the transcriptional start site signal](ATACSeq/EvaluateTSS.Rmd) in the nucleosome-free open region
    - For peaking calling phase, we can ...
      - Use MACS2 and ChIPQC to [call peaks in the nucleosome-free open region, and perform quality assurance](ATACSeq/CallPeak.Rmd)
      - Or use Genrich to call peaks in the nucleosome-free open region
      - Or use MACS3/MACSr (R wrapper of MACS3) to [call peaks in the nucleosome-free open region](ATACSeq/CallPeak.Rmd)
      - Use ChIPseeker to [annotate peak regions with genomic features](ATACSeq/CallPeak.Rmd)
    - For functional analysis phase, we can ...
      - Use rGREAT to [functionally interpret the peak regions based on the GO database](ATACSeq/FunctionalAnalysis.Rmd) 
      - Use GenomicRanges and GenomicAlignments to [select and count non-redundant peaks](ATACSeq/DifferentialAnalysis.Rmd)
      - Use DESeq2/DESeq2-based DiffBind and ChIPseeker to [analyze differences in peaks with gene annotations across conditions](ATACSeq/DifferentialAnalysis.Rmd)
      - Use clusterProfiler to [perform enrichment analysis of differential peak regions](ATACSeq/DifferentialAnalysis.Rmd)
      - However, functional insights gained by peak annotations can hardly illustrate what key regulators shape the transcription mechanism. 
    - To further infer transcription factors acting in peak regions, we can ...
      - Use MotifDb/JASPAR2022 and seqLogo/ [recommend] ggseqlogo to [search and visualize motifs](ATACSeq/Search_Visualize_Motif.Rmd)
      - Use motifmatchr (R wrapper of MOODS) to [map peaks to motifs](ATACSeq/IdentifyMotif.Rmd), DNA sequences preferred by transcription factors
      - Use chromVAR to [analyze differences in motifs across conditions](ATACSeq/Detect_Difference_Motif.Rmd)
  - [Transfer the cell type labels](ATACSeq/Integration/SingleCell/Integration_Full_v1.qmd) from single-cell RNA-seq data to separately collected single-cell ATAC-seq data

<hr>

### Analyze other omics data
#### Genomics
- [MatrixEQTL](https://github.com/andreyshabalin/MatrixEQTL)

#### Proteomics
- A quick start from [loading an online spectrum, performing peak quality control, annotating peaks, to visualizing the annotated peaks](Proteomics/spectrum_utils/0_Quick_Start.py)

#### Metabolomics
- [MetaboAnalystR](https://www.metaboanalyst.ca/docs/RTutorial.xhtml)


<hr>

## Conceptual lens
- [High level multi-omics idea](HighLevelIdea_MultiOmics.md)
- [Case-control design](CaseControl_Design.md)

