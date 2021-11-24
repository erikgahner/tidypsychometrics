# Chapter 4: Item Response Theory

library("tidyverse")

# Page 96-97: dimensionality assessment

library("MPsychoR")
library("mirt")

data("zareki")

zarsub <- zareki |> 
  select(starts_with("subtr"))

library("Gifi")

prinzar <- princals(zarsub)

plot(prinzar, main = "Zareki Loadings")

fitifa1 <- mirt(zarsub, 1, verbose = FALSE)
fitifa2 <- mirt(zarsub, 2, verbose = FALSE, TOL = 0.001) 

anova(fitifa1, fitifa2, verbose = FALSE)


# Page 99-110: unidimensional dichotomous IRT models
## the Rasch model

library("eRm")

fitrasch1 <- RM(zarsub)

fitrasch1

fitrasch1 |> 
  pluck("betapar") |> 
  round(3)

fitrasch1 |> 
  pluck("betapar") |> 
  round(3) |> 
  sort()

timecat <- zareki |> 
  mutate(timecat = ifelse(time <= median(time), "fast", "slow")) |> 
  pull(timecat)

fitLR <- LRtest(fitrasch1, timecat)

fitLR

Waldtest(fitrasch1, timecat)

plotGOF(fitLR, ctrline = list(col = "gray"), conf = list())

fitrasch2 <- zarsub |> 
  select(-subtr5) |> 
  RM()

LRtest(fitrasch2, timecat)

set.seed(123)

T1 <- zarsub |> 
  select(-subtr5) |> 
  as.matrix() |> 
  NPtest(n = 1000, method = "T1")

T1

T11 <- zarsub |> 
  select(-subtr5) |> 
  as.matrix() |> 
  NPtest(n = 1000, method = "T11")

T11


fitrasch2 |> 
  pluck("betapar") |> 
  sort() |> 
  round(2)

plotjointICC(fitrasch2, xlab = "Subtraction Trait", main = "ICCs Subtraction Items")

zarppar <- person.parameter(fitrasch2)


zareki$theta <- zarppar |> 
  pluck("theta.table") |> 
  pull(`Person Parameter`)

aov(theta ~ class, data = zareki) |> 
  summary()

## Two-Parameter Logistic Model

library("ltm")

data("RWDQ")

fit2pl1 <- ltm(RWDQ ~ z1)

fit2pl1 |> 
  coef() |> 
  head()

RWDQ1 <- RWDQ |> 
  dplyr::select(-wdq_22)

fit2pl2 <- ltm(RWDQ1 ~ z1)

fit2pl2 |> 
  coef() |> 
  head()

item.fit(fit2pl2)

plot(fit2pl2, item = 1:5, legend = TRUE)

coef(fit2pl2)[1:5, 2] |> 
  round(3)

ppars <- ltm::factor.scores(fit2pl2, resp.patterns = RWDQ1) |> 
  pluck("score.dat") |> 
  pull(z1)

## Three-Parameter Logistic Model

data("Wilmer")

VPMT <- Wilmer |> 
  dplyr::select(starts_with("vpmt"))

fit3pl <- tpm(VPMT)

fit3pl |> 
  coef() |> 
  head() |> 
  round(3)

plot(fit3pl, item = 1:6, legend = TRUE)


# Page 112-122: unidimensional polytomous IRT models

## Rating Scale Model

data("CEAQ")

itceaq <- CEAQ |> 
  dplyr::select(starts_with("ceaq")) |> 
  mutate_if(is.numeric, ~ .x - 1)

fitrsm <- RSM(itceaq)

ppar <- person.parameter(fitrsm)

ifit0 <- eRm::itemfit(ppar)

ifit0

itceaq1 <- itceaq |> 
  dplyr::select(-ceaq10)

fitrsm1 <- RSM(itceaq1)

ppar1 <- person.parameter(fitrsm1)

ifit1 <- eRm::itemfit(ppar1)

itceaq2 <- itceaq1 |> 
  dplyr::select(-ceaq15)

fitrsm2 <- RSM(itceaq2)

ppar2 <- person.parameter(fitrsm2)

ifit2 <- eRm::itemfit(ppar2)

library("mice")

set.seed(222)

imp <- mice(CEAQ)

gradevec <- complete(imp) |>
  pull(grade)

levels(gradevec) <- c("grade56", "grade56", "grade78", "grade78")

LRtest(fitrsm2, gradevec)

thpar <- thresholds(fitrsm2)

thpar

plotPImap(fitrsm2, latdim = "Empathy",
          main = "Person-Item Map CEAQ")

## Partial Credit Model and Generalizations

data("ASTI")

PGitems <- ASTI |> 
  dplyr::select(ASTI11, ASTI14, ASTI15, ASTI17, ASTI18, ASTI23)

fitpcm <- PCM(PGitems)

thresholds(fitpcm)

plotPImap(fitpcm, latdim = "Presence/Growth", main = "Person-Item Map ASTI")

data("ASTI")

STitems <- ASTI |> 
  dplyr::select(ASTI2, ASTI4, ASTI7, ASTI13, ASTI16, ASTI24, ASTI25)

stpcm <- gpcm(STitems, constraint = "rasch")

stgpcm <- gpcm(STitems)

anova(stpcm, stgpcm)

## Graded Response Model

fitgrm <- grm(STitems)

ppargrm <- ltm::factor.scores(fitgrm)

## Nominal Response Model

library("Gifi")

data("WilPat")

wpit15 <- WilPat |> 
  dplyr::select(Patriotism, Capitalism, Privatization, Nationalism, FreeMarket,
                LowerTaxes, FreeTrade, ChurchAuthority, PrivateHealthcare,
                NuclearEnergy, PrivatePensions, SmallGovernment, Obedience)

wpiprin <- princals(wpit15, ordinal = FALSE)

wpitnew <- wpit15 |> 
  dplyr::select(-c("Nationalism", "Patriotism", "ChurchAuthority", "Obedience"))

wpihom <- homals(wpitnew)

library("mirt")

nrmwp <- mirt(wpitnew, 1, itemtype = "nominal")

# Page 123-125: item and test information

plot(nrmwp, type = "infotrace", main = "Item Information")

plot(nrmwp, type = "info")

# Page 127-129: IRT sample size determination

library("SimDesign")

m <- 20
n <- c(50, 75, 100, 150, 200, 300)

design <- as.data.frame(n)

set.seed(222)

poppars <- rbind(alpha = round(rlnorm(m, 0, 0.25), 2),
                 d = round(rnorm(m), 2))

irtGenerate <- function(condition, fixed_objects = FALSE) {
  n <- condition$n
  a <- fixed_objects['alpha', ]
  d <- fixed_objects['d', ]
  dat <- simdata(a, d, n, itemtype = '2PL')
  return(dat)
}

irtAnalyze <- function(condition, dat, fixed_objects = NULL) {
  mod <- mirt(dat, 1, itemtype = '2PL', verbose = FALSE)
  simpars <- coef(mod, simplify = TRUE, digits = Inf)$items
  irtpars <- c(a = simpars[,1], d = simpars[,2])
  return(irtpars)
}

irtSummarize <- function(condition, results,
                         fixed_objects = NULL) {
  apop <- fixed_objects['alpha', ]
  dpop <- fixed_objects['d', ]
  simrmse <- RMSE(results, c(apop, dpop))
  out <- c(RMSE = simrmse)
  return(out)
}

set.seed(222)

simres <- runSimulation(design, replications = 100,
                        parallel = TRUE,  generate = irtGenerate,
                        analyse = irtAnalyze, summarise = irtSummarize,
                        packages = c('mirt'), fixed_objects = poppars)

simres

sima <- simres |> 
  dplyr::select(contains(".a."))

matplot(simres$n, log(sima), type = "l", col = 1, lty = 1,
        ylab = "log(RMSE)", xlab = "sample size",
        main = "2-PL Monte Carlo", xaxt = "n")
axis(1, at = simres$n)

meanRMSE <- rowMeans(sima)
names(meanRMSE) <- n
round(meanRMSE, 2)

# Page 132-136: differential item functioning

## Logistic Regression DIF Detection

library("lordif")
library("MPsychoR")

data("YouthDep")

cdi <- YouthDep |> 
  dplyr::select(starts_with("CDI"))

cdiDIF <- lordif(cdi, YouthDep$race, criterion = "Chisqr")

cdiDIF |> 
  pluck("stats") |> 
  slice(1:3) |> 
  dplyr::select(item, ncat, chi12, chi13, chi23)

plot(cdiDIF, labels = c("White", "Black", "Asian", "Latino"))

cdiDIF |> 
  pluck("ipar.sparse") |> 
  head(10)

ppar <- cdiDIF |> 
  pluck("calib.sparse") |> 
  pluck("theta")

## Tree-Based DIF Detection

library("psychotree")
library("psychotools")
library("strucchange")

data("MathExam14W")

itmath <- MathExam14W |> 
  pull(solved) |> 
  as.list.data.frame()

covars <- MathExam14W |> 
  dplyr::select(nsolved, tests, gender, study, semester, attempt, group)

mex <- data.frame(solved = itemresp(itmath), covars)

mex <- subset(mex, nsolved > 0 & nsolved < 13)

mex$tests <- ordered(mex$tests)

mex$nsolved <- ordered(mex$nsolved)

mex$attempt <- ordered(mex$attempt)

set.seed(1)

mrt <- raschtree(solved ~ group + tests + nsolved + gender +
                   attempt + study + semester, data = mex, 
                 vcov = "info", minsize = 50, ordinal = "l2", nrep = 1e5)

plot(mrt)

mrt |> 
  itempar() |> 
  as.data.frame() |> 
  dplyr::select(solvedquad, solvedderiv, solvedelasticity, solvedintegral) |> 
  round(2)


# Page 137-144: multidimensional IRT models

## IRT and factor analysis

library("MPsychoR")
library("ltm")

data("RWDQ")

RWDQ1 <- RWDQ |> 
  dplyr::select(-wdq_22)

irtpar <- ltm(RWDQ1 ~ z1)

fapar <- ltm(RWDQ1 ~ z1, IRT.param = FALSE)

cbind(coef(irtpar), coef(fapar)) |> 
  head() |> 
  round(3)

irtppar <- irtpar |> 
  factor.scores() |> 
  pluck("score.dat") |> 
  pluck("z1")

fappar <- fapar |> 
  factor.scores() |> 
  pluck("score.dat") |> 
  pluck("z1")

identical(irtppar, fappar)

## exploratory multidimensional IRT

library("MPsychoR")
library("Gifi")

data("zareki")

itzareki <- zareki |> 
  dplyr::select(starts_with("addit"), starts_with("subtr"))

przar <- princals(itzareki)

plot(przar)

plot(przar, "screeplot")

zar1d <- mirt(itzareki, 1, itemtype = "2PL")

zar2d <- mirt(itzareki, 2, itemtype = "2PL")

anova(zar1d, zar2d)

M2(zar2d)

ifit2D2pl <- mirt::itemfit(zar2d)

ifit2D2pl |> 
  filter(RMSEA.S_X2 < 0.05)

summary(zar2d, rotate = "varimax")

summary(zar2d, rotate = "oblimin")

itemplot(zar2d, 3, main = "ICS addit3", 
         rot = list(xaxis = -70, yaxis = 50, zaxis = 10))

head(MDIFF(zar2d))

head(fscores(zar2d))

class2 <- zareki$class
levels(class2) <- c("second", "thirdfourth", "thirdfourth")

modMG <- multipleGroup(itzareki, model = 2, group = class2,
                       SE = TRUE, verbose = FALSE)

astiDIF <- DIF(modMG, c('a1', 'd'), Wald = TRUE, p.adjust = 'fdr')

astiDIF |> 
  filter(adj_pvals < 0.05) |> 
  pull(adj_pvals) |> 
  round(4)

## confirmatory multidimensional IRT

library("mirt")

data("ASTI")

itasti <- ASTI |> 
  dplyr::select(starts_with("ASTI"))

modASTI <- mirt.model('
    si = 10,19,20,21
    pm = 1,5,9,22
    na = 3,6,8,12
    st = 2,4,7,13,16,24,25
    pg = 11,14,15,17,18,23
    COV = si*pm*na*st*pg
  ')

asti5d <- mirt(itasti, model = modASTI, itemtype = 'graded',
               method = 'MHRM', SE.type = 'MHRM', verbose = FALSE)

astisum <- summary(asti5d, verbose = FALSE)

astisum |> 
  pluck("fcor") |> 
  round(3)

astisum$rotF["ASTI18",] |> 
  round(4)

M2(asti5d, QMC = TRUE)

# Page 145-152: longitudinal IRT Models

## linear logistic models for measuring change

library("MPsychoR")

data("SDOwave")

SDO3 <- SDOwave |> 
  dplyr::select(contains("1996"), contains("1997"), contains("1998")) |> 
  mutate_if(is.numeric, ~ ifelse(.x == 1, 0, 1))

library("eRm")

sdolltm1 <- LLTM(SDO3, mpoints = 3)

sdolltm1 |> 
  pluck("W")

summary(sdolltm1)

W0 <- sdolltm1 |> 
  pluck("W") |> 
  as_tibble() |> 
  select(-`eta 4`, `eta 5`)

sdolltm0 <- LLTM(SDO3, W0)

anova(sdolltm0, sdolltm1)

group <- rep(c(1, 2), each = nrow(SDO3)/2)

sdolltm2 <- LLTM(SDO3, mpoints = 3, group = group)

anova(sdolltm2, sdolltm1)

sdollra <- LLRA(SDO3, mpoints = 3)

summary(sdollra)

## two-tier approach to longitudinal IRT

library("mirt")

SDO2 <- SDOwave |> 
  dplyr::select(contains("1996"), contains("1997"))

SDO2 <- sapply(SDO2, function(co) cut(co, c(0,1,2,3,4,7),
                                      labels = 1:5))
class(SDO2) <- "numeric"

iloads <- rep(1:4, 2)

ttmodel <- mirt.model('
    T1996 = 1-4
    T1997 = 5-8
    COV = T1996*T1997, T1997*T1997
    MEAN = T1997
    CONSTRAIN = (1, 5, d1), (2, 6, d1), (3, 7, d1), (4, 8, d1),
                (1, 5, d2), (2, 6, d2), (3, 7, d2), (4, 8, d2),
                (1, 5, d3), (2, 6, d3), (3, 7, d3), (4, 8, d3),
                (1, 5, d4), (2, 6, d4), (3, 7, d4), (4, 8, d4)')

fitSDO2 <- bfactor(SDO2, iloads, ttmodel, SE = TRUE)

fitSDO2 |> 
  coef() |> 
  pluck("GroupPars") |> 
  data.frame() |> 
  dplyr::select(MEAN_2) |> 
  round(4)

## latent growth IRT models

itloads <- rep(1:4, 2)

modgr <- mirt.model('
      Intercept = 1-8
      Slope = 1-8
      COV = Intercept*Slope, Intercept*Intercept, Slope*Slope
      MEAN = Intercept, Slope
      START = (1-8, a1, 1), (1-4, a2, 0), (5-8, a2, 1)
      FIXED = (1-8, a1), (1-4, a2), (5-8, a2)
      CONSTRAIN =(1, 5, d1), (2, 6, d1), (3, 7, d1), (4, 8, d1),
                 (1, 5, d2), (2, 6, d2), (3, 7, d2), (4, 8, d2),
                 (1, 5, d3), (2, 6, d3), (3, 7, d3), (4, 8, d3),
                 (1, 5, d4), (2, 6, d4), (3, 7, d4), (4, 8, d4)')

fitGIRT <- bfactor(SDO2, itloads, modgr, SE = TRUE)

fitGIRT |> 
  coef() |> 
  pluck("GroupPars") |> 
  data.frame() |> 
  dplyr::select(MEAN_2, COV_21) |> 
  round(4)

# Page 153-156: Bayesian IRT

## Bayesian 2-PL estimation

library("MPsychoR")
library("ltm")

data("RWDQ")

RWDQ1 <- RWDQ |> 
  dplyr::select(wdq_23, wdq_24, wdq_25, wdq_26, wdq_27, wdq_28, wdq_29, wdq_30)

freq2pl <- ltm(RWDQ1 ~ z1)

intstart <- -coef(freq2pl)[,1]

discstart <- coef(freq2pl)[,2]

library("MCMCpack")

chainWDQ1 <- MCMCirt1d(RWDQ1, burnin = 5000, mcmc = 50000,
                       seed = 111, AB0 = 0.15, store.item = TRUE,
                       store.ability = FALSE, verbose = TRUE,
                       alpha.start = intstart, beta.start = discstart)

chainWDQ2 <- MCMCirt1d(RWDQ1, burnin = 5000, mcmc = 50000,
                       seed = 222, AB0 = 0.15, store.item = TRUE,
                       store.ability = FALSE, verbose = TRUE,
                       alpha.start = intstart, beta.start = discstart)

chainWDQ3 <- MCMCirt1d(RWDQ1, burnin = 5000, mcmc = 50000,
                       seed = 333, AB0 = 0.15, store.item = TRUE,
                       store.ability = FALSE, verbose = TRUE,
                       alpha.start = intstart, beta.start = discstart)

WDQlist <- mcmc.list(chainWDQ1, chainWDQ2, chainWDQ3)

plot(WDQlist, auto.layout = FALSE, ask = FALSE)

## dynamic 2-PL model

data("HRB")

HRB1 <- HRB |> 
  filter(if_any(where(is.numeric), ~ .x > 0)) 

rownames(HRB1) <- 1:nrow(HRB1)

time <- rep(1:5, each = 4)

fitdyn <- MCMCdynamicIRT1d(HRB1, item.time.map = time,
                           mcmc = 20000, burnin = 5000,seed = 111, store.ability = TRUE,
                           store.item = FALSE, verbose = TRUE)

dynsum <- summary(fitdyn)

nt <- 5

postmean <- dynsum$statistics[,1]

pertraj <- t(matrix(postmean[1:(nrow(HRB1)*nt)], nrow = nt))

colnames(pertraj) <- paste0("T", 1:5)

pertraj |> 
  head() |> 
  round(3)

matplot(t(pertraj), type = "l", lty = 1, cex = 0.8,
        col = adjustcolor(1, alpha.f = 0.3),
        ylab = "person parameter", xlab = "time points",
        main = "Individual Trajectories")
