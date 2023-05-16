# HRSgenes

Our study on assessment of polygenic risk scores (PRS)

Current title: **An Applied Econometric Assessment of Polygenic Indices**

by Weili Ding (Queen's University), Steven F. Lehrer (Queen's University and NBER), Evgeniya Lukinova (NYU Shanghai and University of Nottingham)

To cite this study: TBD

This repository contains both inline code and code in files with the explanation of how PRS were constructed.

## Data

The [HRS](https://hrs.isr.umich.edu/about) has genotyped 2.5 million single nucleotide polymorphisms (SNPs) on
respondents using Illumina’s Human Omni2.5-Quad (Omni2.5) BeadChip. We used two types of genotyped data: 
- phase 1 and 2 (phg000207.v2.CIDR, PLINK format);
- phases 1, 2 and 3 (phg000842.v1.CIDR, PLINK format); and
imputed data (phg000515.v1; phases 1, 2 and 3).

This data is restricted and available from NCBI Database of Genotypes and Phenotypes (dbGaP) as the distribution set phs000428.v2.p2.c1, which is collected from 15620 participants. 

We used GWAS statistics from Okbay et al. (2016) and Lee et al. (2018) for Educational Attainment (EA) to create PRS.

## Figures & Tables

Most of the figures and regression tables are created via adapted code (not provided here) from replication packages available publicly for Papageorge and Thom (2019) and Barth, Papageorge and Thom (2020) papers in STATA.

## PRS score creation

### PRS for EA in HRS data using PRSice
In order to replicate the polygenic score for EA available at HRS we 
- aligned the genetic (genotyped) data to ‘+’ strand (instead of it being aligned to the illumina TOP strand) using `plink` (Chang et al., 2015):

```
plink --bfile [...]/ncbi/public/files/files/untar/phg000207.v2.CIDR_HRS_phase1.genotype-calls-matrixfmt.c1/[...] --flip [...]/data/fliplist.txt --out [...]/PLINK_sets/target --make-bed
``` 

where 'fliplist.txt' is created from an annotation file provided by HRS support team using `omni_flip.py`.

- changed 'kgp#' markers to 'rs' following ('kgprs.txt') created from an annotation file provided to us by HRS support team: 

```
plink --bfile [...]/PLINK_sets/target --update-name [...]/data/kgprs.txt --out [...]/PLINK_sets/targetDR --make-bed
```

- converted the beta measures to positive values and flipped the reference allele to represent phenotype-increasing PRS, if the beta value from the GWAS meta-analysis was negative (using `betaflip_prsice.py` and GWAS statistics excluding 23andMe provided by Aysu Okbay 'EA2_excl_23andMe_HRS.meta'), and
- modified the source code for PRSice (v.1.25, Euesden, Lewis and O’Reilly, 2015), to achieve the summation of the score using all SNPs: PRSice_v1.25_mod.R, not provided).

The score was standardized within the European population (N = 8598 individuals). 

For EA2, we used phase 1 and 2 genetic data and GWAS statistics from Okbay et al. (2016). For EA3, we
used phase 1 and 2 genetic data and GWAS statistics from Lee et al. (2018).

**Example shell script excerpt for EA2 creation**
```
R --file=PRSice_v1.25_mod.R -q --args \
        plink  /usr/local/bin/plink \
        base [...]/data/cleaning/forprsice/EA2_HRS2_bf.txt \
        target [...]/data/PLINK_sets123/targetDR \
        wd [...]/PRSice_out/PRSice_123EA2/ \
        slower 0 \
        sinc 0.5 \
        supper 1.0 \
        report.best.score.only F \
        report.individual.scores T \
        no.regression T \
        covary F \
        pheno.file [...]/data/cleaning/randomphenotypes.txt \
        clump.snps F \
        allow.no.sex T \
        remove.mhc T
```

where 
- 'EA2_HRS2_bf.txt' is the base dataset, i.e. GWAS statistics adjusted to represent phenotype-increasing PRS acquired via `betaflip_prsice.py`.
- 'targetDR' is the target dataset in PLINK binary format (i.e., consisting of three files: .bed, .bim, and a .fam files, where bed contains the compressed genotype data, bim contains the SNP information and fam contains the family information) created via `plink` calls above (aligning and renaming).
- 'randomphenotypes.txt' is the dummy phenotype file, since no regression is run (as suggested [here](https://www.biostars.org/p/206740/))

#### Split scores creation

We split the EA3 score into four components: significant positive ($β ≥ 0$, $p ≤ 5e^−8$), significant negative ($β < 0$, $p ≤ 5e^−8$), nonsignificant positive and nonsignificant negative. 

I.e., for a split score creation the *base* in the shell script above was changed to the one corresponding to a particular component, e.g. "EA2_HRS2_bf_sigpos.txt" corresponding to the significant positive component (created via `betaflip_prsice.py`.

### PRS for EA using LDPred
We also constructed a polygenic score with adjusted weights (just EA3). Using the LDPred software tool (ver. 1.0.8, Vilhjalmsson et al., 2015) the weights were adjusted for linkage disequilibrium. The LD-adjusted univariate GWAS weights were obtained for 1,433,221 SNPs that are common across the genetic data (after aligning to ‘+’ strand and changing ‘kgp#’ markers to ‘rs’) and the GWAS summary statistics for the educational attainment phenotype (EA3, Lee et al., 2018), and that pass the filters imposed by LDpred: (i) the variant has a minor allele frequency (MAF) greater than 1% in the reference data, (ii) the variant does not have ambiguous nucleotides, (iii) there is no mismatch between nucleotides in the summary statistics and reference data, and (iv) there is no high (> 0.1) MAF discrepancy between summary statistics and validation
sample.

The score was created in the PLINK using the LDPred adjusted weights and then standardized within the European population (N = 8530 individuals). 

#### Split scores creation

First, the SNPs were split into four components using raw β, then the adjusted weights corresponding to each SNP were used to create the scores.

### References

Barth, Daniel, Nicholas W. Papageorge, and Kevin Thom. 2020. "Genetic Endowments and Wealth Inequality." Journal of Political Economy, 128(4): 1474–1522. eprint: https://doi.org/10.1086/705415

Chang, Christopher C, Carson C Chow, Laurent CAM Tellier, et al. 2015. "Second-generation PLINK: rising to the challenge of larger and richer datasets." Gigascience, 4(1): s13742–015. Publisher: Oxford University Press.

Euesden, Jack, Cathryn M Lewis, and Paul F O’Reilly. 2015. "PRSice: polygenic risk score software." Bioinformatics, 31(9): 1466–1468. Publisher: Oxford University Press.

Lee, James J., 23andMe Research Team, COGENT (Cognitive Genomics Consortium), et al. 2018. "Gene discovery and polygenic prediction from a genome-wide association study of educational attainment in 1.1 million individuals." Nature Genetics, 50(8): 1112–1121.

Okbay, Aysu, Jonathan P. Beauchamp, Mark Alan Fontana, et al. 2016. "Genome-wide association study identifies 74 loci associated with educational attainment." Nature, 533(7604): 539–542.

Papageorge, Nicholas W, and Kevin Thom. 2019. "Genes, education, and labor market outcomes: evidence from the health and retirement study.” Journal of the European Economic Association."

Vilhjalmsson, Bjarni J, Jian Yang, Hilary K Finucane, et al. 2015. "Modeling linkage disequilibrium increases accuracy of polygenic risk scores." The american journal of human genetics, 97(4): 576–592. Publisher: Elsevier.

