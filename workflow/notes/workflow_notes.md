# Workflow Notes

## Overview

This document provides an overview of the creation process for the Snakemake workflow. It includes explanations of design decisions, steps taken, and rationale behind choices.

---

## Table of Contents

- [Introduction](#introduction)
- [Workflow Strategy](#workflow-strategy)
- [Rule Dependencies](#rule-dependencies)
- [Configuration](#configuration)
- [Optimizations](#optimizations)
- [Challenges](#challenges)
- [Future Enhancements](#future-enhancements)
- [Resources](#resources)

---

## Introduction

This project is about quantifying the [retrotranscriptome](https://en.wikipedia.org/wiki/Endogenous_retrovirus#Human_endogenous_retroviruses) at a locus specific level of B cells from the [germinal center](https://en.wikipedia.org/wiki/Germinal_center) in [single cells](https://youtu.be/k9VFNLLQP8c?si=WGxfFhNQz_19hPYx).

The purpose of this project is to identify specific HERV loci that are transcribed in specific B cell subtypes from the germinal center.

In this project we analyze publicly available scRNA-seq data from [Holmes et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32603407/) ([GEO Series GSE139891](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139891)) and from [Corinaldesi et al. 2022](https://pubmed.ncbi.nlm.nih.gov/35095922/) ([GEO Series GSE188617](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188617)).
 
We use our newly developed tool Stellarscope to obtain counts matrices for each cell that include both protein-coding transcripts and TE transcripts. Then we analyze the matrices to characterize the retrotranscriptome in B cells from the Light, Intermediate, and Dark zones of the Germinal Center, as well as memory precursors and plasma cells.

### Chemistry 

[10X what is the difference between 3' and 5'](https://kb.10xgenomics.com/hc/en-us/articles/360000939852-What-is-the-difference-between-Single-Cell-3-and-5-Gene-Expression-libraries-)

[10X genomics single cell VDJ](https://mixcr.com/mixcr/reference/overview-built-in-presets/#10xgenomics)

[Comparing 10x Genomics single-cell 3' and 5' assay in short-and long-read sequencing](https://www.biorxiv.org/content/10.1101/2022.10.27.514084v1)

[Cell Ranger V(D)J Algorithms Overview](https://github.com/Dongfang1021/10X_vdj_gex#cell-ranger-vdj-algorithms-overview)

[What is Cell Ranger for Immune Profiling?](https://github.com/Dongfang1021/10X_vdj_gex#what-is-cell-ranger-for-immune-profiling)

[Biostars STARsolo config for 10x Chromium v1, v2, v3](https://www.biostars.org/p/462568/)

[scg_lib_structs](https://teichlab.github.io/scg_lib_structs/)

---

## Workflow Strategy

- Download raw data from GEO
  - Set up and parse SRA run table 
  - Use prefetch to obtain `.sra` files
  - Use fasterqdump to obtain `.fastq` files
- Align raw reads to reference genome
  - create reference genome index
    - download genome sequence
    - download genome annotation
    - build index
  - download whitelist
  - figure out 10x version for parameters
  - run alignment
- Run stellarscope 
  - build stellarscope annotation
    - download gtf
  - run stellarscope pseudobulk
  - run stellarscope celltype
    - download available cell types for cell barcodes
- Merge star and stellarscope matrices
- secondary analysis
  - check that expressed HERVs from LZ and DZ match HERVs from GC clusters
  - integrate all samples
  - HERV load, markers
  - pseudotime
    - pseudotime heatmap
  - networks of selected HERVs


---

## Rule Dependencies

List the dependencies for each rule in the workflow. Include external tools, software, and data sources that each rule relies on.

---

## Configuration

Detail any configuration files used to parameterize the workflow. Explain the settings and how they affect workflow behavior.

---

## Optimizations

Describe any optimizations implemented to speed up the workflow or improve efficiency. This could include parallelization, caching, or other techniques.

---

## Challenges

Discuss any challenges or obstacles encountered during the workflow creation process. Explain how you overcame them.

---

## Future Enhancements

Outline potential improvements or additions that could enhance the workflow in the future. This could include adding new rules, optimizing existing ones, or incorporating new tools.

---

## Resources

List the resources, tutorials, documentation, and references that were helpful in creating the workflow. This can be a valuable section for yourself and others who may want to understand the workflow in more detail.

---

## Conclusion

Sum up the document, reiterating the main decisions made, lessons learned, and the value of the workflow.

---

## Revision History

Record changes and updates to this document over time, along with dates and descriptions of modifications.

---

