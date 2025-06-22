# OIMCS <img src="https://github.com/SciComp8/NGSOmics_Programming/blob/main/dna.png" width="45" height="45">

This repository is created for helping people find what kind of computational biology/bioinformatics questions they can work on. Unlike well-structured frames in mathematics and statistics, where theorems and definitions are clearly laid out, biology and medicine often involve complex, less organized domains of knowledge. However, this repo sets itself apart by cutting through such complexity, with the careful design to offer you a structured and clear overview of approachable problems and their candidate solutions, helping either newcomers or seasoned researchers easily and effectively navigate the computational omics landscape.

*Last updated: 22 Jun 2025*

## Features
* [Single cell RNA-seq analysis](#Analyze-single-cell-RNA-seq-data)
* [Bulk RNA-seq analysis](#Analyze-bulk-RNA-seq-data)
* [Single cell ATAC-seq analysis](#Analyze-single-cell-ATAC-seq-data)
* [Bulk ATAC-seq analysis](#Analyze-bulk-ATAC-seq-data)
* [Proteomics | metabolomics | spatial transcriptomics analysis](#Analyze-other-omics-data)

## Epigenomics sequencing

*Under construction*

### Analyze ChIP-seq data

### Analyze Cut-And-Run data

### Analyze bulk ATAC-seq data
  - [Introduction to transcription regulation](https://www.youtube.com/watch?v=K2gKdFPipv0)
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

<hr>

### Analyze single cell ATAC-seq data
  - Profile distinct open chromatin regions across the genome at single-cell resolution with [Epi ATAC](https://www.10xgenomics.com/support/epi-atac)
  - [Epigenomic and transcriptomic signatures of aging and cancer at single-cell resolution](https://www.nature.com/articles/s43587-024-00751-8)
  - [Technical Q & A](https://kb.10xgenomics.com/hc/en-us/categories/360001072491)
  - [Correct batch effect](https://www.10xgenomics.com/analysis-guides/batch-effect-correction-in-chromium-single-cell-atac-data)
  - [Transfer cell type labels](ATACSeq/Integration/SingleCell/Integration_Full_v1.qmd) from single-cell RNA-seq data to separately collected single-cell ATAC-seq data
  - Find DNA motifs linked to differences in single-cell or bulk chromatin accessibility with [chromVAR](https://github.com/GreenleafLab/chromVAR)
  - Profile somatic mutations with epigenetic alterations at single-cell resolution with [GoT–ChA](https://www.nature.com/articles/s41586-024-07388-y)

<hr>

## DNA sequencing
### Analyze whole genome sequencing data

- **Solve raw read quality control and preprocessing**
  - Run [fastp](https://github.com/OpenGene/fastp) to remove reads with low average quality score, trim [adapters](https://www.thermofisher.com/us/en/home/life-science/cloning/cloning-learning-center/invitrogen-school-of-molecular-biology/next-generation-sequencing/dna-sequencing-preparation-illumina.html), and eliminate [poly-G tails](https://speciationgenomics.github.io/fastp/) in Illumina NovaSeq/NextSeq data
  - Run [MultiQC](https://seqera.io/multiqc/) to evaluate pre- and post-trimming metrics
  - [Validate](https://www.biorxiv.org/content/10.1101/2024.11.23.624993v1.full) sample identity using genetically inferred markers (e.g., sex chromosomes, SNP fingerprinting) and file hashing to ensure data integrity
  - Check [sequencing coverage](https://www.illumina.com/documents/products/technotes/technote_coverage_calculation.pdf) (e.g., 30–50× for human genomes)/read length uniformity/read quality score distribution/GC content distribution
- **Solve alignment and variant calling**
  - For human study: adopt and adapt [GRCh38](https://www.genome.ucsc.edu/cgi-bin/hgGateway) build by [Genome Reference Consortium](https://www.ncbi.nlm.nih.gov/grc), considering [ambiguous mapping](https://www.illumina.com/science/genomics-research/articles/dragen-demystifying-reference-genomes.html)
  - Run [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) for short-read alignment, while [Minimap2](https://github.com/lh3/minimap2) for long-read alignment
  - Detect **single-nucleotide polymorphism (SNP)/indel** with [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us) for germline DNA
  - Detect **structural variant (SV)** with [VISTA](https://academic.oup.com/bib/article/25/5/bbae462/7761957), which optimizes the F1-score of SV calls by combining different high-performing SV callers; or run multiple SV callers (e.g., Manta, DELLY, GRIDSS), and then infer shared SV calls
  - Detect **copy number gains and losses (CNV)** detection with [CNVKit](https://cnvkit.readthedocs.io)
  - Detect **splice-altering intronic variants** with [spliceAI](https://github.com/Illumina/SpliceAI)
  - Understand [VCF](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format), a common variant report file format, and [MAF](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/), which aggregates mutation information from VCF
  - Translate a VCF file from its current reference genome build to another build version: run [LiftoverVcf (Picard)](https://gatk.broadinstitute.org/hc/en-us/articles/360036363632-LiftoverVcf-Picard)
  - Common variants may less likely cause rare or highly [penetrant](https://link.springer.com/article/10.1007/s00439-013-1331-2) diseases: exclude variants with allele frequency > 1% in [Genome Aggregation Database (gnomAD)](https://gnomad.broadinstitute.org/)
 
- **Solve genome annotation**
  - Annotate and assess how genetic variants affect genes and proteins, including specific changes to amino acids with GATK-compatible [SnpEff](http://pcingola.github.io/SnpEff/snpeff/introduction/)
  - Look out for coding, splice-site, and [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) pathogenic variants
  - In a clinical setting, prioritize summary and report of variants using [Human Phenotype Ontology (HPO) terms](https://www.genomicseducation.hee.nhs.uk/genotes/knowledge-hub/the-human-phenotype-ontology); for rare diseases, to improve certainty about whether a variant is pathogenic, look for the one showing phenotypic features in common with [DECIPHER](https://www.deciphergenomics.org); for cancers, mine the biological consequences and therapeutic/diagnostic/prognostic implications of genetic variants with [OncoKB](https://www.oncokb.org) and [oncokb-annotator](https://github.com/oncokb/oncokb-annotator)

- **Solve downstream analysis for biomedical insights**
  - Assess the relationship between genotype and gene expression with [MatrixEQTL](https://github.com/andreyshabalin/MatrixEQTL), which operates linear regression with [additive genotype effect](https://plato.stanford.edu/entries/heritability/#HeriPopuGene)/[ANOVA genotype effect](https://www.fao.org/4/y4391e/y4391e07.htm)
  - Test whether groups of SNPs, often linked to sets of functionally related genes, show a stronger overall association with a phenotype than would be expected by randomness with [INRICH](https://zzz.bwh.harvard.edu/inrich/)
  - Infer differentially expressed genes and enriched pathyways for the trait-associated SNPs with [GIGSEA](https://github.com/zhushijia/GIGSEA?tab=readme-ov-file)

- **Solve scalable and reproducible processing**
  - Leverage cloud-based computating with [Terra](https://www.broadinstitute.org/videos/broade-introduction-terra-scalable-platform-biomedical-research) and [Galaxy](https://galaxyproject.org)
  - Make WGS analysis workflow reproducible with Python-based [snakemake](https://snakemake.readthedocs.io/en/stable/)
  - Do whole‑genome association analysis at biobank scale (i.e., thousands of phenotypes across hundreds of thousands of samples) for both quantitative and binary traits, while being computationally and statistically efficient with [regenie](https://github.com/rgcgithub/regenie)

<hr>

## RNA sequencing (RNA-seq)
- When was the term 'RNA-seq' [first coined](https://www.cell.com/cell/fulltext/S0092-8674(08)00448-0)? At what point did researchers [start RNA-seq before the term 'RNA-seq' was formally introduced](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-7-246)?
- Besides its common use in understanding gene expression differences, and isoform and splicing patterns across tissues or patient of comparison interest, what [other insights](https://academic.oup.com/bib/article/22/6/bbab259/6330938) can RNA-seq technology uncover?

### Analyze single cell RNA-seq data
- [Introduction to single cell RNA-seq](https://github.com/hbctraining/Intro-to-scRNAseq/blob/master/schedule/links-to-lessons.md)
- Profile single cell transcriptome with [Chromium Single Cell Universal 3' Gene Expression](https://www.10xgenomics.com/support/universal-three-prime-gene-expression)
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
  - Run scGen to model [perturbation responses](SingleCellRNASeq/Perturbation/scGen); for heterogeneous perturbed cell population, run [CellOT](https://github.com/bunnech/cellot)
- Review [2025 single cell genomics day](https://satijalab.org/scgd25/)
- Learn [frontier](Science_Reading/scRNAseq.md) single-cell RNA-seq analytical progress
- Identify [single-cell eQTL](https://github.com/XFWuCN/scQTLtools)
- Define [spatial architecture in single cell data](https://www.10xgenomics.com/analysis-guides/integrating-single-cell-and-visium-spatial-gene-expression-data) | [spacexr](https://www.10xgenomics.com/analysis-guides/integrating-10x-visium-and-chromium-data-with-r)
- Capture [gene expression and chromatin accessibility together in a single cell](https://www.10xgenomics.com/support/epi-multiome)
- Experiment with your data analysis process using [COVID-19 RNA-seq data resources](https://github.com/ScienceComputing/COVID-19-RNA-Seq-datasets)

<hr>

### Analyze bulk RNA-seq data
  - Can we [obtain cell-type-specific gene expression information without using single-cell or single-nucleus RNA-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03016-6), as they are costly in clinical research?
  - How do RNA molecules get [prepared](https://github.com/hbctraining/Intro-to-bulk-RNAseq/blob/main/lessons/01_intro-to-RNAseq.md#illumina-library-preparation) and [sequenced](https://github.com/hbctraining/Intro-to-bulk-RNAseq/blob/main/lessons/01_intro-to-RNAseq.md#illumina-sequencing) using Illumina technology?
    - What types of [library preparation kits](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits.html?tags=&pageSize=20&pageNumber=1&sortBy=title) are offered to prepare a [complementary DNA (cDNA) library](https://plato.stanford.edu/entries/molecular-biology/#GoinGenoPostGeno) for sequencing?
      -  Prior to cDNA library preparation, is RNA extracted sufficiently? How is the [degradation level](https://nanoporetech.com/document/requirements/rna-stability) of extracted RNA? Here is a more thorough guide on [checking RNA integrity](https://www.rnaseqcore.vet.cornell.edu/files/TREx_Guide_to_Fragment_Analyzer_RNA_QC%20_V1_2.pdf)
    - What types of [sequencers](https://www.illumina.com/systems/sequencing-platforms.html) are offered to sequence the stable double-stranded cDNA?
  - How accurate are the transcripts measured by the sequencer?
    - Run [FastQC](FastQC/Run_FastQC.sh) or [fastp](FastQC/Run_fastp.sh) to evaluate sequence quality and content
  - Which genome regions are transcribed? What are the exact genomic coordinates (of the reference genome) our sequencing reads come from?
    - There are many ways to find the reads location by aligning reads with the reference genome, but you can choose a tool that's especially useful for your own scientific design. The splice-aware genome aligner [STAR](BulkRNASeq/STAR_Align.sh) is strongly recommened.
        - Other splice-aware alignment tool options include Olego, HISAT2, MapSplice, ABMapper, Passion, BLAT, RUM ...
        - Other alignment tools disregarding isoforms include BWA, Bowtie2 ...
    - Use Rsubread to [align the reads](BulkRNASeq/AlignmentCountingTCell.Rmd)
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
  - Experiment with your data analysis process using [COVID-19 RNA-seq data resources](https://github.com/ScienceComputing/COVID-19-RNA-Seq-datasets)
    
<hr>

### Analyze other omics data
#### Proteomics
- A quick start from [loading an online spectrum, performing peak quality control, annotating peaks, to visualizing the annotated peaks](Proteomics/spectrum_utils/0_Quick_Start.py)

#### Metabolomics
- [MetaboAnalystR](https://www.metaboanalyst.ca/docs/RTutorial.xhtml)

#### Spatial transcriptomics
- Profile the whole transcriptome from Formalin-Fixed Paraffin-Embedded tissue sections with [Visium](https://www.10xgenomics.com/support/spatial-gene-expression-ffpe)
- Profile subcellular-level RNA targets with [Xenium](https://www.10xgenomics.com/support/in-situ-gene-expression)
- Read spatial transcriptomics data with [SpatialData](Spatial_Transcriptomics/Read_Data.py)
- Analyze spatial DNA/RNA/protein data in subcellular/single cell/multiple cells with [Giotto Suite](https://github.com/drieslab/Giotto)
  
<hr>

## Conceptual lens
- [High level multi-omics idea](HighLevelIdea_MultiOmics.md)
- [Case-control design](StudyDesign/CaseControl_Design.md)

<hr>

## Citation

> Mary Piper, Meeta Mistry, Jihe Liu, William Gammerdinger, & Radhika Khetani. (2022, January 6). hbctraining/scRNA-seq_online: scRNA-seq Lessons from HCBC (first release). Zenodo. https://doi.org/10.5281/zenodo.5826256.

