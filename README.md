# DNA-Seq Germline Workflows in Nextflow

DNA-Seq Germline pipelines from Broad Institute's WARP project rewritten in Nextflow. All workflows implement the GATK Best Practices for variant analysis in Next-Generation Sequencing data.

The WARP workflows included in this project are:
* **Exome Germline Single Sample**: performs read mapping, germline small variants calling and quality check in human exome sequencing data;
* **Whole Genome Germline Single Sample**: performs read mapping, germline small variants calling and quality check in human whole-genome sequencing data. It implements the DRAGEN-GATK mode for functional equivalence to the DRAGEN software.

References:
* [WARP GitHub](https://github.com/broadinstitute/warp)
* [WARP Documentation](https://broadinstitute.github.io/warp/docs/get-started)
* [Genome Analysis Toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us)
* [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894711-About-the-GATK-Best-Practices)
