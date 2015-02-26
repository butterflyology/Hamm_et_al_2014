# Here are the data and code from the Hamm et al. 2014 Molecular Ecology paper on Wolbachia in Drosophila suzukii and D. subpulchrella. Everything in this code that you need to re-create the analyses conducted in the paper and the figures (except the phylogenetics, which were not done in R)

# The code chunks below generate output realted to:
#  - make a map with states coded by date of first *D. suzukii* detection
#  - determine infection frequencies across and among populations
#  - determine rates of vertical transmission
#  - CI assays
#  - MK assays
#  - qPCR Ct values
#  - Desiccation experiment
#  - Starvation experiment

set.seed(1234)
setwd("")
library(exactRankTests)
library(bootstrap)
library(maps)
library(mapdata)
options(warn = -1) #I do this here only because of an encoding issue for Figure 3. 
#save(list = ls(), file = "SWD_Wolbachia.R")
load("SWD_Wolbachia.R") # This R data object contains all data used in the paper. Just load this object rather than import all of the files individually. 

# Here I make the map that acted as Figure 1. I could not figure out how to get Nova Scotia to be filled in, so I just used Photoshop. I know, that is lame. If you fingure out how to color in NS, please send me an email or a tweet and let me know.

data(worldMapEnv)
a <- map("state", plot = FALSE)
# states <- read.csv("", header = TRUE) 
head(states)
states$State %in% a$names # check to ensure states in my file are spelled identically to R. This is vital because you won't get an error if there is a mispelling.

#Here I make a wrapper that colors states by date of detection:
cols <- NA
for(i in 1:length(states[,1])){
  if (states$Year[i] == "2008")cols[i]="#000000"
	if (states$Year[i] == "2009")cols[i]="#2E2E2E"
	if (states$Year[i] == "2010")cols[i]="#5D5D5D"
	if (states$Year[i] == "2011")cols[i]="#8B8B8B"
	if (states$Year[i] == "2012")cols[i]="#D1D1D1"
	if (states$Year[i] == "2013")cols[i]="#E8E8E8"
	if (states$Year[i] == "none")cols[i]="#FFFFFF"
}

states <- cbind(states, cols)
states <- states[order(states$State),]
states$cols <- as.character(states$cols)

#Here are the coordinates 
lat <- c(44.64, 43.12, 38.52, 36.91)
long <- c(-63.57, -77.67, -121.97, -121.77)

map(database="world",ylim=c(25, 50), xlim=c(-130, -60), lwd=2)
map("state", region=states$State, interior=TRUE, boundary=TRUE, fil=TRUE, col=states$cols, lwd=2, add=TRUE)
points(long, lat, pch=21, cex=2, bg="white", lwd=2)
legend("bottomright", legend=c("2013", "2012", "2011", "2010", "2009", "2008"), col=c("#E8E8E8", "#D1D1D1", "#8B8B8B", "#5D5D5D", "#2E2E2E", "#000000"), bty="n", pch=15, pt.cex=2)

# Here I define a few functions for use later:

# To calculate the exact binomial confidence intervals
bn.conf.exact <- function(x,n, conf.level){
  p <- x / n
	alpha <- 1 - conf.level
	alpha <- rep(alpha, length = length(p))
	alpha2 <- 0.5 * alpha
	x1 <- x == 0
	x2 <- x == n
	lb <- ub <- x
	lb[x1] <- 1
	ub[x2] <- n[x2] - 1
	lcl <- 1 - qbeta(1 - alpha2, n + 1 - x, lb)
	ucl <- 1 - qbeta(alpha2, n - ub, x + 1)
	if(any(x1))
	lcl[x1] <- rep(0, sum[x1])
	if(any(x2))
	ucl[x2] <- rep(1, sum[x1])
	out <- matrix(c(p, lcl, ucl), nrow = 1, ncol = 3, dimnames=list("data", c("mean", "lower", "upper")))
	print(out)
}

# To calculate standard error
se <- function(x){
sqrt(var(x, na.rm = TRUE) / length(na.omit(x)))
} 

# To calculate the traditional 95% CI
ci95 <- function(x){
	t.value <- qt(0.975, length(x) - 1)
	standard.error <- se(x)
	ci <- t.value * standard.error
	cat("mean" = mean(x), "95% CI =", mean(x) - ci, "to", mean(x) + ci, "\n")
}

# The first results presented in the paper are the Wolbachia infection frequencies by population, we then compare frequencies among populations.

bn.conf.exact(10, 109, conf.level = 0.95) # Rochester August 2012
bn.conf.exact(12, 178, conf.level = 0.95) # Rochester Sept 2012
bn.conf.exact(7, 38, conf.level = 0.95) # Wolfskill June 2012
bn.conf.exact(41, 71, conf.level = 0.95) # Wolfskil May 2013
bn.conf.exact(32, 192, conf.level = 0.95) # Watsonville September 2012
bn.conf.exact(7, 40,conf.level = 0.95) # Watsonville, October 2012
bn.conf.exact(13, 57, conf.level = 0.95) # Watsonville August 2013
bn.conf.exact(4, 34, conf.level = 0.95) # Nova Scotia September 2013
bn.conf.exact(19, 136, conf.level  = 0.95) # Watsonville October 2013

# the tests for homogeneity of infection frequencies used Chi-squared test with P-values estimated using Monte Carlo sampling. 

swd_x1 <- matrix(c(10, 12, 4, 7, 41, 32, 7, 13, 33, 99, 166, 30, 31, 30, 160, 33, 44, 177), nrow = 9) # Includes all populations
chisq.test(swd_x1, simulate.p.value = TRUE, B = 100000)

swd_East <- matrix(c(10, 12, 4, 99, 166, 30), nrow = 3) # Eastern samples
chisq.test(swd_East, simulate.p.value = TRUE, B = 100000)
bn.conf.exact(26, 321, conf.level = 0.95)

swd_West <- matrix(c(7, 32, 7, 13, 4, 31, 160, 33, 44, 30), nrow = 5) # Western samples, exlcuding the outlier from California 
chisq.test(swd_West, simulate.p.value = TRUE, B = 100000)
bn.conf.exact(78, 463, conf.level = 0.95)

swd_x2 <- matrix(c(10, 12, 4, 7, 32, 7, 13, 33, 99, 166, 30, 31, 160, 33, 44, 177), nrow = 8) #excluding the outlier
chisq.test(swd_x2, simulate.p.value = TRUE, B = 100000)
chisq.test(swd_x2)

swd_Wolf <- matrix(c(7, 41, 31, 30), nrow = 2) # Wolfskill/Winters population only, excluding the outlier
chisq.test(swd_Wolf, simulate.p.value = TRUE, B = 100000)

# 1df Chisq test for East vs. West
east <- c(0.0917, 0.0674, 0.1170)
west <- c(0.184, 0.167, 0.175, 0.228, 0.140, 0.578)
kruskal.test(list(west, east))

# Now we examine the data for vertical transmission and make Figure 3.

tf <- c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.2, 0.4, 0.9, 0.9, 42/44, 24/64, 1, 1, 1, 1) # create a vector of transmission rates

theta <- function(x){
  mean(x)
}

mean(tf)
bcanon(tf, nboot = 10000, theta = theta, alpha = c(0.025, 0.975)) # Bias Corrected and accelerated (BCa) bootstrap.

# The maternal transmission data and Figure 3. 
JJ1 <- seq(0, 1, by = 0.1)
JJ2 <- c(0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 10)
JJ <- as.data.frame(cbind(JJ1, JJ2))
JJ
colnames(JJ)[1] <- "Proportion"
colnames(JJ)[2] <- "F-infect"
UCD_M <- as.data.frame(cbind(1, c(20 / 21), c(13 / 32)))

plot(JJ, pch = 19, cex = 1.5, ylab = "# females", cex.axis = 1.2, las = 1, cex.lab = 1.2, xlab = (expression(paste(italic("Wolbachia"), " transmission rate"))), ylim = c(0, 10), type = "n", bty = "l")
points(0.2, 1, pch = 23, cex = 1.5, bg = "grey")
points(0.4, 1, pch = 23, cex = 1.5, bg = "grey")
points(0.9, 2, pch = 23, cex = 1.5, bg = "grey")
points(1, 10, pch = 23, cex = 1.5, bg = "grey")
points(c(42 / 43), 1, pch = 25, cex = 1.5, bg = "grey")
points(c(11 / 32), 1, pch = 25, cex = 1.5, bg = "grey")
points(1, 4, pch = 25, cex  =1.5, bg = "grey")
points(1, 4, pch = 17, cex = 1.5) 
points(c(20 / 21), 1, pch = 17, cex = 1.5)
points(c(13 / 32), 1, pch = 17, cex = 1.5)
legend("topleft", legend = c("NY \u2640 offspring", "CA \u2640 offspring", "CA \u2642 offspring"), pt.cex = 1.5, pch = c(23, 25, 17), bty = "n", pt.bg = c("grey", "grey", "black"))

#Fisher's exact test asking if Wolbachia is randomly associated with haplotype
htype <- matrix(c(5, 2, 2, 3, 4, 0, 2, 0, 1, 4), ncol = 2, nrow = 5, dimnames = list(c("H1", "H2", "H3", "H4", "H5"), c("U", "I")))
fisher.test(htype, simulate.p.value = TRUE, B = 10000)

# Here are the cytoplasmic incompatibility (CI) assays for *D. suzukii*:

# TurCI <- read.csv("", header = TRUE) # For D. suzukii from California

TurCI$paternal <- as.factor(TurCI$paternal)
TurCI$maternal <- as.factor(TurCI$maternal)

pat1 <- subset(TurCI, paternal == 1, drop = TRUE) # 38 males
mat1 <- subset(TurCI, maternal == 1, drop = TRUE) # 51 females
pat0 <- TurCI[TurCI$paternal == 0, ]
mat0 <- TurCI[TurCI$maternal == 0, ]

p1m1 <- pat1[pat1$maternal == 1, ]
p1m0 <- pat1[pat1$maternal == 0, ]
p0m1 <- mat1[mat1$paternal == 0, ]
p0m0 <- pat0[pat0$maternal == 0, ]

mean((p0m0$cum_total - p0m0$cum_unhatch) / p0m0$cum_total)
se((p0m0$cum_total - p0m0$cum_unhatch) / p0m0$cum_total)
w00 <- ((p0m0$cum_total - p0m0$cum_unhatch) / p0m0$cum_total)

mean((p0m1$cum_total - p0m1$cum_unhatch) / p0m1$cum_total)
se((p0m1$cum_total - p0m1$cum_unhatch) / p0m1$cum_total)
w01 <- ((p0m1$cum_total - p0m1$cum_unhatch) / p0m1$cum_total)

mean((p1m0$cum_total - p1m0$cum_unhatch) / p1m0$cum_total)
se((p1m0$cum_total - p1m0$cum_unhatch) / p1m0$cum_total)
w10 <- ((p1m0$cum_total - p1m0$cum_unhatch) / p1m0$cum_total)

mean((p1m1$cum_total - p1m1$cum_unhatch) / p1m1$cum_total)
se((p1m1$cum_total - p1m1$cum_unhatch) / p1m1$cum_total)
w11 <- ((p1m1$cum_total - p1m1$cum_unhatch) / p1m1$cum_total)

kruskal.test(list(w00, w01, w10, w11))


# Jaen <- read.csv("", header = TRUE) # For D. suzukii from New York

Jaen$Female <- as.factor(Jaen$Female)
Jaen$Male <- as.factor(Jaen$Male)

fem1 <- Jaen[Jaen$Female == 1, ]
mal1 <- Jaen[Jaen$Male ==1, ]
fem0 <- Jaen[Jaen$Female == 0, ]
mal0 <- Jaen[Jaen$Male == 0, ]

jm0f0 <- mal0[mal0$Femal == 0, ]
jm0f1 <- mal0[mal0$Female == 1, ]
jm1f0 <- mal1[mal1$Female == 0, ]
jm1f1 <- mal1[mal1$Female == 1, ]

mean(jm0f0$n_hatched / jm0f0$Total_eggs)
se(jm0f0$n_hatched / jm0f0$Total_eggs)
j00 <- (jm0f0$n_hatched / jm0f0$Total_eggs)

mean(jm0f1$n_hatched / jm0f1$Total_eggs)
se(jm0f1$n_hatched / jm0f1$Total_eggs)
j01 <- (jm0f1$n_hatched / jm0f1$Total_eggs)

mean(jm1f0$n_hatched / jm1f0$Total_eggs)
se(jm1f0$n_hatched / jm1f0$Total_eggs)
j10 <- (jm1f0$n_hatched / jm1f0$Total_eggs)

mean(jm1f1$n_hatched / jm1f1$Total_eggs)
se(jm1f1$n_hatched / jm1f1$Total_eggs)
j11 <- (jm1f1$n_hatched / jm1f1$Total_eggs)

kruskal.test(list(j00, j01, j10, j11)) # statistically significant, but not because of CI. The result is due to the low hatch rate of m0 f1, which is not a CI interaction

# Combining the New York and California data
kruskal.test(list(w00, w01, w10, w11, j00, j01, j10, j11)) # still statistically significant

# Here are the CI assays for *D. subpulchrella*:

# subp <- read.csv("", header = TRUE)
subp$Spoon <- as.factor(subp$Spoon)
subp <- subp[complete.cases(subp), ]

IMIF <- subset(subp, cross == "IMxIF")
dim(IMIF)
IMIF <- IMIF[IMIF$grand_total >= 10, ]
mean(IMIF$grand_total)
se(IMIF$grand_total)

IMUF <- subset(subp, cross == "IMxUF")
IMUF <- IMUF[IMUF$grand_total >= 10, ]
dim(IMUF)
mean(IMUF$grand_total, na.rm = TRUE)
se(IMUF$grand_total)

UMIF <- subset(subp, cross == "UMxIF")
UMIF <- UMIF[UMIF$grand_total >= 10, ]
dim(UMIF)
mean(UMIF$grand_total, na.rm = TRUE)
se(UMIF$grand_total)

UMUF <- subset(subp, cross == "UMxUF")
UMUF <- UMUF[UMUF $grand_total >= 10, ]
dim(UMUF)
mean(UMUF$grand_total, na.rm = TRUE)
se(UMUF$grand_total)

kruskal.test(list(IMIF$grand_total, IMUF$grand_total, UMIF$grand_total, UMUF$grand_total), na.rm = TRUE)

#exact ranks test for CI crosses
wilcox.exact(IMUF$grand_total, UMIF$grand_total, exact = TRUE, alternative = "two.sided")

# Now the data on male / female ratios of *Wolbachia* infected offspring and Figure 4:

# lf1 <- read.csv("", header = TRUE) # First the data from D. suzukii
lf1$Wol <- as.factor(lf1$Wol)
lf1$Total

WPosC <- subset(lf1, lf1$Wol == "1")
sum(WPosC$Total)
mean(WPosC$Total)

WNegC <- subset(lf1, Wol == "0")
sum(WNegC$Total)
mean(WNegC$Total)

WNegM <- WNegC[, c(seq(3, 61, 5))]
WNegM <- rowSums(WNegM, na.rm=TRUE)
WNegF <- WNegC[, c(seq(4, 61, 5))]
WNegF <- rowSums(WNegF, na.rm=TRUE)

WPosM <- WPosC[, c(seq(3, 61, 5))]
WPosM <- rowSums(WPosM, na.rm = TRUE)
WPosF <- WPosC[, c(seq(4, 61, 5))]
WPosF <- rowSums(WPosF, na.rm=TRUE)

binom.test(sum(WPosF), sum(WPosM, WPosF), p = 0.5) # no significant difference between males and females

# subp <- read.csv("", header = TRUE) # The D. subpulchrella data
spINF <- subset(subp, infection == "yes")

infM <- spINF[, c(seq(4, 16, 4))]
infM <- rowSums(infM, na.rm = TRUE)
infF <- spINF[, c(seq(5, 17, 4))]
infF <- rowSums(infF, na.rm = TRUE)

spUN <- subset(subp, infection == "no")
unM <- spUN[, c(seq(4, 16, 4))]
unM <- rowSums(unM, na.rm = TRUE)
unF <- spUN[, c(seq(5, 17, 4))]
unF <- rowSums(unF, na.rm = TRUE)

binom.test(sum(infF), sum(infF, infM), p = 0.5) # Compare number of males and females from infected mothers

plot(jitter(WNegF), jitter(WNegM), pch = 19, las = 1, ylab = '# males', xlab = '# females', ylim = c(0, 80), xlim = c(0, 80), col = 'black', cex = 1.2)
points(jitter(WPosF), jitter(WPosM), col = 'dark grey', pch = 19, cex = 1.2)
points(jitter(unF), jitter(unM), col = "black", pch = 17, cex = 1.2)
points(jitter(infF), jitter(infM), col = "dark grey", pch = 17, cex = 1.2)
abline(0, 1, lwd = 2)
legend("topleft", legend=c(expression(paste(italic("D. suzukii"), " -")), expression(paste(italic("D. suzukii"), " +")), expression(paste(italic("D. subpulchrella"), " -")), expression(paste(italic("D. subpulchrella"), " +"))), bty = "n", pch = c(19, 19, 17, 17), col = c("black", "dark grey"))

# Here we compare the delta Ct values from qPCR with the fecundity

# q1 <- read.csv("", header = TRUE) # read in the data

qWol <- subset(q1, Assay == "Wol")
qRef <- subset(q1, Assay == "Ref")

w28 <- subset(qWol, Sample == "w28f")
r28 <- subset(qRef, Sample == "w28f")
cq28 <- mean(r28$Cq - w28$Cq)

w59 <- subset(qWol, Sample == "w59f")
r59 <- subset(qRef, Sample == "w59f")
cq59 <- mean(r59$Cq - w59$Cq)

w60 <- subset(qWol, Sample == "w60f")
r60 <- subset(qRef, Sample == "w60f")
cq60 <- mean(r60$Cq - w60$Cq)

w94 <- subset(qWol, Sample == "w94f")
r94 <- subset(qRef, Sample == "w94f")
cq94 <- mean(r94$Cq - w94$Cq)

w140 <- subset(qWol, Sample == "w140f")
r140 <- subset(qRef, Sample == "w140f")
cq140 <- mean(r140$Cq - w140$Cq)

w227 <- subset(qWol, Sample == "227")
r227 <- subset(qRef, Sample == "227")
cq227 <- mean(r227$Cq - w227$Cq)

w231 <- subset(qWol, Sample == "231")
r231 <- subset(qRef, Sample == "231")
cq231 <- mean(r231$Cq - w231$Cq)

w208 <- subset(qWol, Sample == "208")
r208 <- subset(qRef, Sample == "208")
cq208 <- mean(r208$Cq - w208 $Cq)

w210 <- subset(qWol, Sample == "210")
r210 <- subset(qRef, Sample == "210")
cq210 <- mean(r210$Cq - w210 $Cq)

w228 <- subset(qWol, Sample == "228")
r228 <- subset(qRef, Sample == "228")
cq228 <- mean(r228$Cq - w228 $Cq)

w229 <- subset(qWol, Sample == "229")
r229 <- subset(qRef, Sample == "229")
cq229 <- mean(r229$Cq - w229$Cq)

w232 <- subset(qWol, Sample == "232")
r232 <- subset(qRef, Sample == "232")
cq232 <- mean(r232$Cq - w232 $Cq)

w237 <- subset(qWol, Sample == "237")
r237 <- subset(qRef, Sample == "237")
cq237 <- mean(r237$Cq - w237 $Cq)

wpos1 <- c(cq28, cq59, cq60, cq94, cq140)
wpos2 <- c(cq208, cq210, cq227, cq231, cq237)

w <- as.data.frame(c(wpos1, wpos2), nrow = 10, byrow = TRUE)
colnames(w) <- "titer"
w$fec <- c(0, 24, 13, 0, 0, 0, 62, 29, 66, 39)

lm1 <- lm(w$titer ~ w$fec)
summary(lm1)

# Here are the data and analyses from the desiccation experiment:

# Desicc <- read.csv("", header = TRUE)

DesWPos <- subset(Desicc, Wol == "A")
DesWNeg <- subset(Desicc, Wol == "B")

ci95(DesWNeg$Total)
ci95(DesWPos$Total)

wilcox.exact(DesWPos$Total, DesWNeg$Total, exact = TRUE)

# Here are the data form the starvation experiment:

# Starve <- read.csv("", header = TRUE)

StarWPos <- subset(Starve, Wol == "A")
StarWPos <- StarWPos[complete.cases(StarWPos), ]

StarWNeg <- subset(Starve, Wol == "B")
StarWNeg <- StarWNeg[complete.cases(StarWNeg), ]

ci95(StarWPos$Hrs_Surv)
ci95(StarWNeg$Hrs_Surv)