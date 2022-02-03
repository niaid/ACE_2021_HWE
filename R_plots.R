library(tidyverse)

###### HWE exercise #####
# Plot genotype vs allele frequencies along with the HWE expectations
# Adapted from https://gcbias.org/2011/10/13/population-genetics-course-resources-hardy-weinberg-eq/

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


######### NGSAdmix exercise ###########

######### Plot inferred admixture proportions for a single run of K=3
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


######### Larger analysis of 100 individuals

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
#POPGELPER: http://www.royfrancis.com/pophelper/


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
