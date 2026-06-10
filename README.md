# DNA-Seq Germline Workflows in Nextflow

Nextflow reimplementation of Broad Institute's WARP germline DNA-seq pipelines. These workflows implement GATK Best Practices for variant discovery and analysis in next-generation sequencing data.

There are some changes to the original WDL material. Preemptible and disk settings have been removed for workflow portability, enabling their execution in local or AWS environments. Process output naming prioritizes prefix convention instead of full filenames in order to standardize file extensions. The bioinformatics logic remains faithful to WARP.

The workflows included are:
* **Exome Germline Single Sample**: performs read mapping, small variants calling and quality control in human exome sequencing data.
* **Whole Genome Germline Single Sample**: performs read mapping, small variants calling and quality control in human whole-genome sequencing data. It implements the DRAGEN-GATK mode for functional equivalence to the DRAGEN software.

References:
* [WARP](https://broadinstitute.github.io/warp/docs/get-started)
* [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us)
* [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices)
* [Nextflow](https://docs.seqera.io/nextflow/)
