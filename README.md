# Single-Nucleus RNA Sequencing of Developing and Mature Superior Colliculus Identifies Neuronal Diversity and Candidate Mediators of Circuit Assembly

This repository accompanies the study: **Single-Nucleus RNA Sequencing of Developing and Mature Superior Colliculus Identifies Neuronal Diversity and Candidate Mediators of Circuit Assembly** by Ayupe and Choi, *et al*. Manuscript can be found here: [bioRxiv link]().

Code used for analysis in producing results presented in the paper can be found in the *scripts* folder. If you have any questions regarding this repository, study, data, or code, please contact James Choi (jsc228 at miami dot edu) or Kevin Park, PhD (kpark at med dot miami dot edu, corresponding author, [lab website](https://www.parklabmiami.com/)).

## Data Availability

Raw data are available from the SRA (Sequence Read Archive) database under study accession [SRP\*\*\*\*\*\*](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP295673). Processed scRNAseq data are available at NCBI GEO accession [GSE\*\*\*\*\*\*](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162610).

Direct downloads for the processed data and Seurat objects can be found here: [link]()

-   **sc.rds**:
-   **neuron.rds**:

## Reproducibility

Analysis was run using the files in `scripts` in the following sequence (which more or less follows the sequence of figures in the manuscript):

1.  `preprocessing-quality-control.qmd`
2.  `integrated-cluster-analysis.qmd`
3.  `neuronal-cluster-analysis.qmd`
4.  `development-differential-expression.qmd`
5.  `cadherin-expression.qmd`
6.  `axon-guidance-molecules.qmd`
7.  `cross-reference-comparison.qmd`
8.  `Pax7-validation.qmd`

Additional `.R` files include various analyses performed for the study but whose results were not included.

## superior-colliculus-snRNAseq browser

A ShinyApps browser is available to more easily query the gene expression data. Portal is available at <https://parklabmiami.shinyapps.io/superior-colliculus-snRNAseq/>

Code repository to build the ShinyApp can be found in the **superior-colliculus-snRNAseq** subfolder.

## Acknowledgements

This work was supported by grants from the National Eye Institute (NEI) R01EY022961 (K.K.P.), NEI R01EY032542 (K.K.P), NEI 1U01EY027257 (K.K.P), NEI R21EY031026 (K.K.P), DOD W81XWH-19-1-0736 (K.K.P), The Miami Project to Cure Paralysis, and the Buoniconti Fund (K.K.P) and Glaucoma Research Foundation (K.K.P).
