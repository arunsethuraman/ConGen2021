# ConGen2021 - Day 4 - Correcting for Multiple Testing
Arun Sethuraman, PhD
San Diego State University
asethuraman@sdsu.edu
Twitter: @arunsethuraman

# Workshop Objectives
1. Review performing tests of Hardy-Weinberg Equilibrium, Linkage Disequilibrium, and Differentiation
2. Learn how to correct for multiple testing
3. Compare results and inferences between different correction techniques, as against uncorrected tests

# Pre-requisites
1. vcftools - https://vcftools.github.io/index.html
2. Base R or RStudio
3. Input files - accessible on the Box page under Lecture Materials > Arun Sethuraman > Hands-on Files

# Exercise 0 - Reviewing tests of HWE, LD, Differentiation
We will be using the data from Bataillon et al. 2015 (Bataillon, T., Duan, J., Hvilsom, C., Jin, X., Li, Y., Skov, L., Glemin, S., Munch, K., Jiang, T., Qian, Y. and Hobolth, A., 2015. Inference of purifying and positive selection in three subspecies of chimpanzees (Pan troglodytes) from exome sequencing. Genome biology and evolution, 7(4), pp.1122-1132.) for this exercise. You can download the VCF file here: http://datadryad.org/resource/doi:10.5061/dryad.56m2g or the SNP file is available as "SNP_v37_version3.vcf" in the Hands-on Files folder.

Back story - this study analyzes whole exome data from three subspecies of chimpanzees - Pan troglodytes troglodytes (Central), P.t.verus (Western), and P.t.schweinfurthii (Eastern) which have very distinctive geographical ranges separated by large river systems in Africa, and have not been observed to hybridize in the wild with each other, or with Pan paniscus (Bonobos, but genomic evidence states otherwise!). Nonetheless, chimpanzee numbers are dwindling in the wild, with estimates of effective population sizes between 22,000 - 27,000 (Fischer et al., 2004, Sethuraman and Hey 2016, Won and Hey 2005).

![image](https://user-images.githubusercontent.com/5439390/132733181-192de700-0c7f-409c-aedb-621b9e674bed.png)

Further references:
1. Won, Y.J. and Hey, J., 2005. Divergence population genetics of chimpanzees. Molecular biology and evolution, 22(2), pp.297-307.
2. Sethuraman, A. and Hey, J., 2016. IM a2p–parallel MCMC and inference of ancient demography under the Isolation with migration (IM) model. Molecular ecology resources, 16(1), pp.206-215.
3. Kuhlwilm, M., Han, S., Sousa, V.C., Excoffier, L. and Marques-Bonet, T., 2019. Ancient admixture from an extinct ape lineage into bonobos. Nature ecology & evolution, 3(6), pp.957-965.
4. Fischer, A., Wiebe, V., Pääbo, S. and Przeworski, M., 2004. Evidence for a complex demographic history of chimpanzees. Molecular biology and evolution, 21(5), pp.799-808.


1. Create a new directory, and place the SNP VCF file in this folder.
```shell
mkdir chimps
cd chimps
mv ../SNP_v37_version3.vcf .
```

2. Compute heterozygosity across the entire file
```shell
vcftools --vcf SNP_v37_version3.vcf --het --out allchroms_het
R #this will open R
```

3. Visualize this in R, to observe some basic statistics about the individual chimps.What all statistics were obtained? You’ll notice that the O.HOM. column contains the number of observed homozygote sites, and the E.HOM. contains the number of expected homozygote sites (under HWE), N.SITES is the total number of loci that were analyzed. The first column (INDV) contains the names of all individuals, while the last column, F contains the inbreeding coefficient. 
What do you observe? Who is most inbred? Who is least inbred?

 ```R
 allchroms_het<-read.table("allchroms_het.het",header=TRUE)
 allchroms_het
 summary(allchroms_het)
 plot(allchroms_het[,5])
 q()
 ```
 4. Testing for Hardy-Weinberg Equilibrium, computing Diversities, Tajima's D across chromosomes in 10k windows, visualization
 ```Shell
 vcftools --vcf SNP_v37_version3.vcf --hardy --out allchroms_hwe
 vcftools --vcf SNP_v37_version3.vcf --TajimaD 10000 --out allchroms_tajimad
 vcftools --vcf SNP_v37_version3.vcf --window-pi 10000 --window-pi-step 10000 --out allchromspi_10kwindow
```

```R
allchroms_hwe<-read.table("allchroms_hwe.hwe",header=TRUE)
tajimad<-read.table("allchroms_tajimad.Tajima.D",header=TRUE)
tajimad_nomissing<-na.omit(tajimad)
hist(tajimad_nomissing$TajimaD)
#Alternately, if we plot it by chromosome

#Plot of log p-values for sites in HWE
boxplot(log(P_HWE,10)~CHR,data=allchroms_hwe,xlab="Chromosome",ylab="log10 P-value")

#Plot of Tajima's D
plot(tajimad_nomissing$CHROM,tajimaD_nomissing$TajimaD,xlab="Chromosome",ylab="Tajima's D")

```
These should produce rather ugly plots like the ones below:

<img width="674" alt="Screen Shot 2021-09-09 at 11 30 53 AM" src="https://user-images.githubusercontent.com/5439390/132742740-5b788ca2-9e90-4193-8e6c-d246c14502e5.png">

![image](https://user-images.githubusercontent.com/5439390/132740592-346ddaf2-130c-4e5b-891b-928102d05ba3.png)


So let's make these prettier, and visualize as a Manhattan Plot instead - this requires the qqman package in R. Note that qqman and its dependencies (e.g. calibrate) are available for R.3.5.0 and above - so please note that the Manhattan plots will only work with these versions and above).

```R
install.packages("qqman")
library(qqman)
manhattan(allchroms_hwe,chr="CHR",bp="POS",p="P_HWE",snp="POS",logp=TRUE,ylab="log p-values",ylim=c(-3.0,3.0))
manhattan(tajimad_nomissing,chr="CHROM",bp="BIN_START",p="TajimaD",snp="N_SNPS",logp=FALSE,ylab="Tajima’s D",ylim=c(-3.0,3.0))
```
![image](https://user-images.githubusercontent.com/5439390/132736520-2acb3475-bd12-474c-8474-94862eaf6776.png)[tajimad_hist.pdf]


Voila! Much prettier(?) What do you notice about these?


Now we'll do the same thing, but with the diversity estimates.

```R
diversity<-read.table("allchromspi_10kwindow.windowed.pi",header=TRUE)
manhattan(diversity,chr="CHROM",bp="BIN_START",p="PI",snp="N_VARIANTS",logp=FALSE,ylab="Pi",ylim=c(0,0.001))
```

![image](https://user-images.githubusercontent.com/5439390/132736812-51634ad5-f8c6-4f91-8bca-ab811d3aa5a1.png)

Now let's try and compute some differentiation statistics based on these data. THe three files that are provided (schweinfurthii.txt, troglodytes.txt, and verus.txt) contain the names of each of the individuals from each separate population. Make sure that these files are present in the same folder (chimps) prior to running these commands.

```Shell
vcftools --vcf chimps.vcf --weir-fst-pop schweinfurthii.txt --weir-fst-pop troglodytes.txt --weir-fst-pop verus.txt --out allthreepopsfst
R
```

Thereon, let's plot and analyze these in R:

```R
fst<-read.table("allthreepopsfst.weir.fst",header=TRUE)
manhattan(fst,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="POS",logp=FALSE,ylab="Fst")
```

<img width="668" alt="Screen Shot 2021-09-09 at 2 37 14 PM" src="https://user-images.githubusercontent.com/5439390/132765979-97833f8e-d75d-4dbf-a565-74bc00cd4f65.png">

However, the problem with vcftools is that it doesn't compute p-values for estimates of Tajima's D, diversity, or Fst. So there are other packages that might be used for the same purpose - e.g. hierfstat, ARLEQUIN

# Exercise 1 - Correcting for Multiple Testing - Tests of HWE

1. Let's use the same P-values from the tests of HWE - how do we correct for multiple testing in base R?

Turns out, base R has some very useful functions to do this - we will just add new columns to our HWE results, to compare them against each other.

```R
allchroms_hwe$bonferroni=p.adjust(allchroms_hwe$P_HWE,"bonferroni")
allchroms_hwe$holm=p.adjust(allchroms_hwe$P_HWE,"holm")
allchroms_hwe$bh=p.adjust(allchroms_hwe$P_HWE,"BH")
allchroms_hwe$fdr=p.adjust(allchroms_hwe$P_HWE,"fdr")

#Now you can determine which ones are statistically significant at the alpha (FPR) of 0.05
bonferroni<-which(allchroms_hwe$bonferroni < 0.05)
holm<-which(allchroms_hwe$holm < 0.05)
bh<-which(allchroms_hwe$bh < 0.05)
fdr<-which(allchroms_hwe$fdr < 0.05)

#Some examples
head(allchroms_hwe[holm,])
tail(allchroms_hwe[bh,])
```

What do you notice about these loci?

Some thought exercises:

1) Based on your tests of Tajima’s D, do you identify particular loci that have highly positive, and highly negative values across the genome? Make a list of these sites, as well as their coordinates. HINT - you can do this in R by summarizing the Tajima’s D values (using the summary() function), then using the which() function to find out which values are either greater than or less than a “cutoff” value.

2) Repeat the analyses above for identifying so called “outlier” loci/windows that have extremely high levels of diversity or extremely low levels of diversity, and extremely high levels of differentiation, and extremely low levels of differentiation. Are there any overlaps between these lists?

3) Interpret your results - why would you have overlapping regions between these lists at all? Explain.

4) Now do something fun - go to the Ensembl Chimpanzee genome database: https://uswest.ensembl.org/Pan_troglodytes/Info/Index

Now search for some of these loci that you’ve identified as outliers. E.g. The highest Tajima’s D value happens to be 2.97, which is at a window on Chromosome 11, between 56110000 and 56120000. So searching for 11:5611000-5621000 shows up several genes in and around that region. What are some genes that you notice? What would a high degree of positive Tajima’s D at this locus mean? 

Now redo this analysis with some other Tajima’s D outliers. What genes do you identify? What do you hypothesize is happening at this locus across chimpanzees?

# Exercise 2 - Correcting for Multiple Testing - Tests of HWE - Data from Wei and Nielsen 2019

1. Make sure to have downloaded the data (Exercise2.txt) from the Box page: Lecture Materials > Arun Sethuraman > Hands-on Files > Exercise2.txt

2. In Unix shell, you can "explore" this data by using head/less, or other commands - for instance, here you'll notice that the file contains count information for each genotype - the first column contains the name of the SNP, while the rest of the columns contain genotype counts of A/A, A/B, B/B respectively, where A and B are the two alleles at each locus. For our purpose, we will use the "HardyWeinberg" package in R to test, correct, and obtain "outliers".

3. Now in R:

```R
install.packages("HardyWeinberg")
library(HardyWeinberg)
x<-read.table("Exercise2.txt")
#note that this is of dimension 5933 x 4 - you can find this by doing dim(x)
pvals<-c() #Empty array to store p-values from HWE tests
i<-1
while (i < 5933) {
pvals[i]<-HWChisq(c(x[i,2],x[i,3],x[i,4]))$pval
i=i+1
}
summary(pvals)
```

4. Determine how many loci are out of HWE at an alpha of 0.05, then correct the p-values and then determine the same numbers

```R
length(which(pvals<0.05))
# Turns out there are 1477 loci that are out of HWE out of the 5932, if we don't correct for FPR
# Now let's try a bunch of corrections, to see if we can whittle this down to correct for FPR
bonferroni<-p.adjust(pvals,"bonferroni")
bh<-p.adjust(pvals,"BH")
holm<-p.adjust(pvals,"holm")
#Note that from the manuscript, rs62625034 is the deletion that is determined as the delta32 mutation - so you can search for this particular variant in the dataset, see what you obtain as a p-value prior to correction, and after correction
z<-which(x$V1=="rs62625034")
pvals[z]
which(pvals < pvals[z])
```

What do you notice? Do we recapitulate the results from Wei and Nielsen? 

Thought Exercise:

Why or why don't we recapitulate the results from Wei and Nielsen? Can you think of some other ways for us to make this more efficient? 

# Exercise 3 - Tests of differentiation - Fst outliers using OutFLANK, or similar methods

This tutorial has been adapted from the OutFLANK vignette that can be accessed here: https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html

```R
library(OutFLANK)
library(vcfR)
#Read the chimps.vcf file, convert it into OutFLANK format
chimps<-read.vcfR("chimps.vcf")
geno <- extract.gt(chimps) # Character matrix containing the genotypes
position <- getPOS(chimps) # Positions in bp
chromosome <- getCHROM(chimps) # Chromosome information
#Let's specify the population ID's as 1, 2 and 3 for P.t.schweinfurthii, P.t.troglodytes, P.t.verus respectively
pop <- c("1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","3","3","3","3","3","3","3","3")

G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))

G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
G[is.na(G)]<-9 #important! Missing data calls need to be coded as 9 in OutFLANK

table(as.vector(G))

# Now G should be in OutFLANK format
# Calculate Fst - note I'm only doing this for 1000 loci

my_fst<-MakeDiploidFSTMat(t(G[1:1000,]), locusNames = position[1:1000], popNames = pop)
# You can view this
head(fst)

#Now run OutFLANK
out_trim<-OutFLANK(fst,NumberOfSamples=31,qthreshold=0.05,Hmin=0.1)
head(out_trim$results)

#There are several pruning steps, sanity checks that I am skipping here for time - please refer to the original vignette for more details!

#OutFLANK automates p-value correction for multiple testing using the FDR method, and threshold - please see Whitlock and Lotterhos 2015
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
head(P1)
summary(P1$pvalues)
```

Thought Exercise:


What do you determine? Do you find any outlier loci? If you were interpreting these plainly based on p-values, would you identify outliers?
What if you had done the same analyses, but with a different q threshold?



















