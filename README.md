# Analyzing population genomic data with ANGSD  
**updated Feb. 2022**

This repository contains practical training material as part of the ACE bioinformatics program. This lesson focuses on exploring Hardy-Weinberg Equilibrium and performing population structure inference and population assignment from Next-Generation Sequencing data.

*author*: David B. Stern, Ph.D.  
*email*: david.stern@nih.gov  

This tutorial uses the following software  
* [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD)
* [NGSAdmix](http://www.popgen.dk/software/index.php/NgsAdmix)
* [Tidyverse R packages](https://www.tidyverse.org/)

R code for making plots are embedded here and also available in [R_plots.R](R_plots.R) so it can be imported into RStudio.

The data are from the [1000 Genomes Project](https://www.internationalgenome.org/).

Pre-prepared output files are available on in this Github repo (download with `git clone https://github.com/niaid/ACE_2021_HWE/`)

### connect to ACE HPC and start an interactive session
```bash
ssh <username>@biocompace.ace.ac.ug

# I will start an interactive session with: srun --pty bash
# We can all connect to the interactive session using ssh.
# I will tell you the node name once I start the interactive session.
ssh <nodename>

#create variables that point to the ANGSD and NGSAdmix software locations
ANGSD=/etc/ace-data/bcbb_teaching_files/HWE_Structure_Oct2021/software/angsd
NGSADMIX=/etc/ace-data/bcbb_teaching_files/HWE_Structure_Oct2021/software/NGSadmix

```
### Process data
Copy the data from the shared directory to your directory.  
```bash
mkdir angsd_tutorial #make a working directory
cd angsd_tutorial #change to that directory
cp /etc/ace-data/bcbb_teaching_files/HWE_Structure_Oct2021/bamfiles.tar.gz . #copy the data to this directory
tar -xvzf bamfiles.tar.gz #extract the files. This should create a directory called bamfiles
```

The directory should contain 30 .bam files, each with an index generated by Samtools. Each .bam file is the result of mapping low-coverage genome sequencing data from one individual to the human genome reference. The dataset contains 10 individuals from each of 3 putative ancestries, indicated in the file names.

Make a list with the bam file names  
`ls -1 bamfiles/*.bam > bamfiles.txt`  

Create a population information file for later

`paste -d " " <( cut -f 5 -d"." bamfiles.txt ) <(cut -f 1 -d"." bamfiles.txt | xargs -n1 basename) > bamfiles.popinfo.txt`  

### Simulatenously detect SNPs, calculate genotype likelihoods, and test for HWE  

```bash
## analyze all 30 individuals from the 3 putative ancestries
$ANGSD -bam bamfiles.txt -doHWE 1 -domajorminor 1 -GL 2 -doGlf 2 -doMaf 1 -SNP_pval 2e-6 -minMapQ 30 -minQ 20 -minInd 25 -minMaf 0.05 -out bamfiles
```

While that is running, let's look at the command flags
```
-bam bamfiles.txt: tell ANGSD the locations of the bam files to process
-doHWE 1: calculate HWE test p-value based on genotype likelihoods
-domajorminor 1: Infer major and minor alleles from genotype likelihoods
-GL 2: use the GATK model for genotype likelihoods
-doGlf 2: Output genotype likelihoods in beagle format
-doMAF 1: tells ANGSD to calculate minor allele frequencies
-SNP_pval 2e-6: only run the test on sites that are called SNPs (p-value < 2e-6)
-minMapQ 30: Only use reads with a mapping score above 30
-minQ 20: Only count bases with a quality score above 20
-minInd 25: Only consider sites with data from at least 25 individuals
-minMaf 0.05: Only consider SNPs with a minor allele frequency of 0.05
-out bamfiles: The prefix for output files (e.g. bamfiles.beagle.gz)
```

### Inspect the output files
```bash
#view first few rows
gunzip -c bamfiles.hwe.gz | head | column -t
#count the total number of SNPs
gunzip -c bamfiles.hwe.gz | wc -l
#count the number of sites with HWE p-value < 0.05
gunzip -c bamfiles.hwe.gz | awk '{if($9<0.05) total+=1}END{print total}'
```
What proportion of sites are not in Hardy-Weinberg Equilibrium?  

Does that mean all the other sites are in HWE?  

### Let's run the same test with 10 individuals from one of the 3 locations / putative ancestries
```bash
#make a like of bam files for JPT individuals
ls -1 bamfiles/*JPT*.bam > JPT.files
$ANGSD -bam JPT.files -doHWE 1 -domajorminor 1 -GL 2 -doGlf 2 -doMaf 1 -SNP_pval 2e-6 -minMapQ 30 -minQ 20 -minInd 8 -minMaf 0.05 -out JPT
gunzip -c  JPT.hwe.gz | wc -l
gunzip -c JPT.hwe.gz | awk '{if($9<0.05) total+=1}END{print total}'
```
What proportion of sites are not in Hardy-Weinberg Equilibrium?

### View the genotype likelihood file
```bash
gunzip -c bamfiles.beagle.gz | head -n 10 | cut -f 1-10 | column -t
```
Explanation of the beagle genotype likelihood file
```bash
 -first column: the name of the SNP (usually chromosome_position)
 -second and third columns: the major and minor alleles (0=A, 1=C, 2=G, 3=T)
 -fourth, fifth, and sixth columns: the estimated likelihoods for the 3 genotype for individual 1
        - Homozygous major (A/A), Heterozygous (A/B), Homozygous minor (B/B)
```

## Plot genotype vs allele frequencies along with the HWE expectations  
### Adapted from:  
 https://gcbias.org/2011/10/13/population-genetics-course-resources-hardy-weinberg-eq/

**Collect the beagle genotype likelihood file we generated in the above exercise:**  
File name: `bamfiles.beagle.gz`

Open RStudio, set the working directory to the directory with the `bamfiles.beagle.gz` file and run the following R code:  

```R
genotypes.beagle<-read.table('bamfiles.beagle.gz',header=T)

#get the counts of each genotype
AA_count <- round(rowSums(genotypes.beagle[,seq(4,91,3)]))
Aa_count <- round(rowSums(genotypes.beagle[,seq(5,92,3)]))
aa_count <- round(rowSums(genotypes.beagle[,seq(6,93,3)]))
counts<-cbind(AA_count,Aa_count,aa_count)
#calculate the total
tot.counts<-rowSums(counts)
#calculate the frequency for each genotype
geno.freq<-counts/tot.counts
#calculate the frequency for the first allele
#frequency of homozygous genotype plus 1/2 heterozygous
allele.freq<-(geno.freq[,1]+.5*geno.freq[,2])

##alleles are ordered by minor allele, to make this prettier flip 1/2 the alleles around
these.minor<-sample(1:nrow(geno.freq),nrow(geno.freq)/2)
these.major<-sample(1:nrow(geno.freq),nrow(geno.freq)/2)
ss.allele<-c(allele.freq[these.minor],1-allele.freq[these.major])
ss.geno<-rbind(geno.freq[these.minor,],geno.freq[these.major,c(3,2,1)])

#plot the allele frequency vs the 3 genotype frequencies
plot(ss.allele,ss.geno[,1],xlim=c(0,1),ylim=c(0,1),col=adjustcolor("red",0.1),xlab="allele frequency",ylab="genotype frequency")
points(ss.allele,ss.geno[,3],xlim=c(0,1),ylim=c(0,1),col=adjustcolor("blue",0.1))
points(ss.allele,ss.geno[,2],xlim=c(0,1),ylim=c(0,1),col=adjustcolor("green",0.1))

#plot the mean for each plot
smooth=1/5
lines(lowess(ss.geno[,1]~ss.allele,f = smooth),col="black")
lines(lowess(ss.geno[,3]~ss.allele,f = smooth),col="black")
lines(lowess(ss.geno[,2]~ss.allele,f = smooth),col="black")

#calculate and plot the Hardy-Weinberg expectations
x=1:1000/1000 #allele frequency from 0 to 1 (i.e. p)
lines(x,x^2,lty=2) #homozygous AA
lines(x,2*x*(1-x),lty=2) #heterozygous Aa. q = 1-p
lines(x,(1-x)^2,lty=2) #homozygous aa
legend(x=0.3,y=1,col=c("red","blue","green",rep("black",2)),legend=c("Homozygote AA","Homozygote aa","Heterozygote Aa","Mean","Hardy Weinberg Expectation"),pch=c(rep(1,3),rep(NA,2)),lty=c(rep(NA,3),1,2))
```


# Using NGSAdmix to estimate assignment and ancestry proportions

Back on the ACE HPC...  
We will run NGSAdmix on the genotype likelihood file we generated previously.  

```bash
$NGSADMIX -likes bamfiles.beagle.gz -K 3 -minMaf 0.05 -seed 1 -o bamfiles.ngsadmix.k3
```

Explanation of the command flags
```bash
 -likes : the genotype likelihood file
 -K : the number of populations to assign ancestry to
 -minMaf: the minimum minor allele frequency to include a SNP -- to avoid biases from sequencing errors
 -seed : random seed
 -o : output prefix
```

View the output files
```bash
#log file
cat bamfiles.ngsadmix.k3.log
#ancestry proportions
head -10 bamfiles.ngsadmix.k3.qopt | column -t
```

## Plot ancestry proportions  
**Collect pre-prepared output files**  
File names:  
_bamfiles.ngsadmix.k3.qopt_  
_bamfiles.popinfo.txt_  


Open RStudio, set the working directory to the directory with the downloaded files, and run the following R code:   

### Plot inferred admixture proportions
```R
q<-read.table("bamfiles.ngsadmix.k3.qopt")
pop<-read.table("bamfiles.popinfo.txt")$V1
# Plot them (ordered by population)
ord = order(pop)

par(mar=c(7,4,1,1))
barplot(t(q)[,ord], col=c("chocolate1","deepskyblue4","brown3"),
        names=pop[ord],las=2,
        space=0,
        ylab="Admixture proportions",
        xlab="Individuals",
        border=NA,
        cex.names=0.75,main="NGSAdmix K=3")
```

### Estimating the best-fit K  
##### We will use a larger dataset for which genotype likelihoods have already been estimated  
**100 individuals**  
**50,000 sites**  
Individuals are labeled by their geographic location / putative ancestry:   
ASW	- HapMap African Americans from SW US  
CEU	- European individuals  
CHB	- Han Chinese in Beijing  
JPT	- Japanese individuals  
YRI	- Yoruba individuals from Nigeria  
MXL	- Mexican individuals from LA California  

Here is how the analysis was run, but **no need to run it now**.  
It will take quite a long time:  

```bash
for k in {1..15};
do
for rep in {1..10};
do
$NGSADMIX -likes 100ind.GL.gz -K $k -minMaf 0.05 -seed $RANDOM -o 100ind.ngsadmix.k${k}.rep${rep} -P 16
#capture likelihood
grep 'best like' 100ind.ngsadmix.k${k}.rep${rep}.log | cut -d '=' -f 2 | cut -d ' ' -f 1 > 100ind.ngsadmix.k${k}.rep${rep}.loglik
done
done
```

**Collect the files from Github**   
File name: 100ind.tar.gz  

Extract the files by double clicking or using the command line:   
`tar -xvzf 100ind.tar.gz`

In RStudio --  

```R
library(tidyverse)
# create list of all input files
files <- list.files(path="100ind",pattern='loglik',full.names=TRUE)

#collect the log likelihoods for each K and put into a data frame
logliks <- data.frame()
for (file in files){
    #pull the K value from the file name
    k <- substr(strsplit(file,'\\.')[[1]][[3]],2,3)
    #read in the log likelihood
    loglik <- read.table(file)[1,]
    df <- data.frame("K"=k, "logLik"=loglik)
    logliks <- rbind(logliks, df)
}

#plot the log-likelihoods for each K
#these are negative log likelihoods, so higher is better
plot(logliks$K,logliks$logLik, col="blue",
    ylab="Log Likelihood",
    xlab="K",
    main="NGSAdmix K vs Likelihood")

#high K = more parameters, so necessarily higher likelihood. But, we want to avoid over-fitting. Perhaps what we want is the largest change in likelihood

difK_df <- data.frame()
#starting from 4
for (k in 4:15){
    #calculate mean K
    k_mean <- mean(filter(logliks, K==k)$logLik)
    #calculate mean K of the previous K
    km1_mean <- mean(filter(logliks, K==k-1)$logLik)
    #take the difference
    difK <- k_mean - km1_mean
    dfout <- data.frame('K'=k, "difK"=difK)
    difK_df <- rbind(dfout, difK_df)
}

plot(difK_df$K,difK_df$difK, col="dodgerblue4", type="l", lwd=3,
    ylab="L'(K)",
    xlab="K",
    main="L'(K)")

#perhaps K=4 is the best fit due to the largest increase in likelihood at K=4 from K=3

#Let's calculate the deltaK statistic from Evanno et al

dK_df <- data.frame()
for (k in 2:14){
    #1 average the L(K) over the replicates
    k_mean <- mean(filter(logliks, K==k)$logLik)
    #2 estimate from these averages L''(K) as abs( L(K+1) - 2L(K) + L(K-1) )
    kp1_mean <- mean(filter(logliks, K==k+1)$logLik)
    km1_mean <- mean(filter(logliks, K==k-1)$logLik)
    kpp <- abs(kp1_mean - 2*k_mean + km1_mean)
    #3 divide by the standard deviation of L(K) (sd of the different replicates for the same K)
    #some of the standard deviations are unreasonably small, so lets remove those K values
    k_sd <- ifelse(abs((sd(filter(logliks, K==k)$logLik))) < 10,0,abs((sd(filter(logliks, K==k)$logLik))))
    deltaK <- kpp / k_sd
    dfout <- data.frame("K"=k,"deltaK"=deltaK)
    dK_df <- rbind(dK_df,dfout)
}

plot(dK_df$K,dK_df$deltaK, col="dodgerblue4", type="l", lwd=5,
    ylab="DeltaK",
    xlab="K",
    main="NGSAdmix K vs deltaK")

#better option might be to use:
#STUCTURE harvester: http://taylor0.biology.ucla.edu/structureHarvester/
#Clumpak: http://clumpak.tau.ac.il/bestK.html
#POPHELPER: http://www.royfrancis.com/pophelper/


# So K=13 is the highest. Let's plot the admixture proportions for K=13
files <- list.files(path="100ind",pattern='qopt',full.names=TRUE)
pop<-read.table("100ind/poplabel.txt")$V1
k13_files <- files[grepl("k13",files)]

#Let's pick one of the replicates
#For a real analysis, it would make sense to take the average over multiple replicates
data_k13 <- read.table(k13_files[1])
# Plot them (ordered by population)
ord = order(pop)

par(mar=c(7,4,1,1))
cols <- rainbow(13)

barplot(t(data_k13)[,ord],col=cols,names=pop[ord],las=2,
space=0,
ylab="Admixture proportions",
xlab="Individuals",
border=NA,
cex.names=0.75,main="NGSAdmix K=13")

#let's try K=4 to compare and plot below K=13
k4_files <- files[grepl("k4",files)]
data_k4<- read.table(k4_files[1])

#create space for two bar plots
par(mfrow=c(2,1))
#plot K=13
barplot(t(data_k13)[,ord],col=cols,names=pop[ord],las=2,
space=0,
ylab="Admixture proportions",
xlab="Individuals",
border=NA,
cex.names=0.75,main="NGSAdmix K=13")

#plot K=4
barplot(t(data_k4)[,ord],col=cols,names=pop[ord],las=2,
space=0,
ylab="Admixture proportions",
xlab="Individuals",
border=NA,
cex.names=0.75,main="NGSAdmix K=4")
```
