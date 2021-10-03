# Chapter 2: Factor Analysis

library("tidyverse")

# Page 19: dichotomize items for 2x2 table

library("MPsychoR")
data("YouthDep")

item1 <- YouthDep |> 
  mutate(CDI1 = fct_recode(CDI1, c("2" = "1"))) |> 
  pull(CDI1) 

item2 <- YouthDep |> 
  mutate(CDI15r = fct_recode(CDI15r, c("2" = "1"))) |> 
  pull(CDI15r) 

table(item1, item2)

# Page 20: tetrachoric correlation

library("psych")

tetcor <- cbind(item1, item2) |> tetrachoric()

tetcor

draw.tetra(0.35, 1.16, 0.36)

# Page 21: polychoric correlation

item1 <- YouthDep |> 
  pull(CDI1) 

item2 <- YouthDep |> 
  pull(CDI15r) 

polcor <- cbind(item1, item2) |> polychoric()
polcor

DepItems <- YouthDep |> 
  select(1:26)
Depnum <- data.matrix(DepItems) - 1
Rdep <- polychoric(Depnum)

# Page 22: tetrachoric correlation

data("Rmotivation")

Rmotivation1 <- Rmotivation |> 
  select(starts_with("ext"), starts_with("int"))

Rmot1 <- tetrachoric(Rmotivation1, smooth = FALSE)

Rmot1 |> 
  pluck("rho") |> 
  eigen() |> 
  pluck("values") |> 
  round(3) |> 
  tail()

Rmot <- tetrachoric(Rmotivation1)
Rmot |> 
  pluck("rho") |> 
  eigen() |> 
  pluck("values") |> 
  round(3) |> 
  tail()

# Page 25-28: Explanatory factor analysis

motFA <- Rmot |> 
  pluck("rho") |> 
  fa(nfactors = 2, rotate = "none", fm = "ml")

motFA |> 
  pluck("loadings") |> 
  print(cutoff = 0.2)

## communalities
motFA |> 
  pluck("communality") |> 
  round(2)

## communalities
motFA |> 
  pluck("communality") |> 
  round(2)

## varimax rotation
motFA2 <- Rmot |> 
  pluck("rho") |> 
  fa(nfactors = 2, rotate = "varimax", fm = "ml")

## non-orthogonal rotation
Rmot2 <- Rmotivation |> 
  select(1:36) |> 
  tetrachoric()

motFA3 <- Rmot2 |> 
  pluck("rho") |> 
  fa(nfactors = 3, rotate = "oblimin", fm = "ml")

motFA3 |> 
  pluck("loadings")

## factor correlation matrix
motFA3 |> 
  pluck("Phi") |> 
  round(3)

## Page 29: factor scores
motFA2 <- Rmotivation1 |> 
  fa(nfactors = 2, rotate = "varimax",
     cor = "tet", fm = "ml", scores = "regression",
     missing = TRUE, impute = "median")

motFA2 |> 
  pluck("scores") |> 
  dim()


# Page 30-31: Scree plot
Rdep <- Depnum |> 
  polychoric() |> 
  pluck("rho")

evals <- Rdep |> 
  eigen() |> 
  pluck("values")

scree(Rdep, factors = FALSE)

## proportion of explained variance (in %)
(evals / sum(evals) * 100) |> 
  head(2)

## parallel analysis,
set.seed(123)

resPA <- fa.parallel(Depnum, fa = "pc", cor = "poly", fm = "ml")

# Page 33: very simple structure
resvss <- vss(Rdep, fm = "ml", n.obs = nrow(Depnum), plot = FALSE)

resvss

## Page 34: Tucker Lewis Index
fadep <- fa(Depnum, 1, cor = "poly", fm = "ml")

summary(fadep)


# Page 34: number of factors with nfactors()
resnf <- nfactors(Depnum, n = 8, fm = "ml", cor = "poly")

resnf

# Page 36-38: Bayesian EFA

library("MPsychoR")
library("corrplot")
library("BayesFM")

data("Privacy")

Privstd <- scale(Privacy)

Privstd |> 
  cor() |> 
  corrplot()

Nid <- 2
pmax <- trunc(ncol(Privstd) / Nid)
pmax

set.seed(123)

Rsim <- simul.R.prior(pmax, nu0 = pmax + c(1, 2, 5, 7, 10))

plot(Rsim)

Ksim <- simul.nfac.prior(nvar = ncol(Privstd), Nid = Nid,
                         Kmax = pmax, kappa = c(.1, .2, .5, 1))

plot(Ksim)

set.seed(222)
fitbefa <- befa(Privstd, Nid = 2, Kmax = pmax, nu0 = 10,
                kappa = 0.2, kappa0 = 0.1, xi0 = 0.1,
                burnin = 5000, iter = 50000)

fitbefa <- post.column.switch(fitbefa)

fitbefa <- post.sign.switch(fitbefa)

sumbefa <- summary(fitbefa)

plot(fitbefa)

# Page 41-44: confirmatory factor analysis
library("MPsychoR")
library("lavaan")

data("Rmotivation")

Rmot <- Rmotivation |> 
  select(starts_with("ext"), starts_with("int")) |> 
  na.omit()

mot_model <- '
    extrinsic  =~ ext1 + ext2 + ext3 + ext4 + ext5 + ext6 +
                  ext7 + ext8 + ext9 + ext10 + ext11 + ext12
    intrinsic =~  int1 + int2 + int3 + int4 + int5 '

fitMot <- lavaan::cfa(mot_model, data = Rmot,
                      ordered = names(Rmot))

library("semPlot")

semPaths(fitMot, what = "est", edge.label.cex = 0.7,
         edge.color = 1, esize = 1, sizeMan = 4.5, asize = 2.5,
         intercepts = FALSE, rotation = 4, thresholdColor = "red",
         mar = c(1, 5, 1.5, 5), fade = FALSE, nCharNodes = 4)


## standardized estimates
inspect(fitMot, what = "est") |> 
  pluck("theta")

## loadings
inspect(fitMot, what = "est") |> 
  pluck("lambda")

inspect(fitMot, what = "std") |> 
  pluck("lambda")

## latent variable covariance matrix
inspect(fitMot, what = "est") |> 
  pluck("psi")

inspect(fitMot, what = "std") |> 
  pluck("psi")

## test whether loadings differ from 0
parameterEstimates(fitMot, standardized = TRUE)

## full model output including goodness-of-fit measures
summary(fitMot, standardized = TRUE, fit.measures = TRUE)

parameterEstimates(fitMot) |> 
  slice(5)

mot_model2 <- '
    extrinsic  =~ ext1 + ext2 + ext3 + ext4 + ext6 + ext7 +
                  ext8 + ext9 + ext10 + ext11 + ext12
    intrinsic =~  int1 + int2 + int3 + int4 + int5 '
fitMot2 <- lavaan::cfa(mot_model2, data = Rmot,
                       ordered = names(Rmot)[-5])

# Page 45: higher-order confirmatory factor analysis models

Rmot2 <- Rmotivation |> 
  select(1:4, 13:16, 32:35) |> 
  na.omit()

mot_model3 <- '
    extrinsic  =~ ext1 + ext2 + ext3 + ext4
    hybrid =~ hyb1 + hyb2 + hyb3 + hyb4
    intrinsic =~  int1 + int2 + int3 + int4
    motivation =~ extrinsic + hybrid + intrinsic '

fitMot3 <- lavaan::cfa(mot_model3, data = Rmot2,
                       ordered = names(Rmot2))

summary(fitMot3, standardized = TRUE, fit.measures = TRUE)


# Page 47: confirmatory factor analysis with multiple indicators multiple independent causes

Rmot3 <- Rmotivation |> 
  select(1:4, 13:16, 32:35, 39:41) |> 
  na.omit()

mot_model4 <- '
    extrinsic  =~ ext1 + ext2 + ext3 + ext4
    hybrid =~ hyb1 + hyb2 + hyb3 + hyb4
    intrinsic =~  int1 + int2 + int3 + int4
    motivation =~ extrinsic + hybrid + intrinsic
    motivation ~ npkgs + phd '

fitMot4 <- lavaan::cfa(mot_model4, data = Rmot3,
                       ordered = names(Rmot3[1:12]))

parameterEstimates(fitMot4) |> 
  slice(16, 17)

# Page 49-52: multigroup confirmatory factor analysis

library("semTools")
library("MPsychoR")
library("lavaan")

data("Bergh")

GP_model <- ' GP =~ EP + HP + DP + SP '

minvfit <- measurementInvariance(model = GP_model, data = Bergh,
                                 group = "gender", estimator = "MLR")

minvfit |> 
  pluck("fit.configural") |> 
  summary(standardized = TRUE, fit.measures = TRUE)

GP_model <-' GP =~ c(v1,v1)*EP + c(v2,v2)*HP + c(v3,v3)*DP + SP '

fitBase <- lavaan::cfa(GP_model, data = Bergh, group = "gender",
                       estimator = "MLR")


GP_model <- ' GP =~ EP + HP + DP + SP '
fitBase <- lavaan::cfa(GP_model,data = Bergh, group = "gender",
                       group.equal = c("loadings"),
                       group.partial = c("GP=~ SP"), estimator = "MLR")

fitMeasures(fitBase)

fitBase1 <- lavaan::cfa(GP_model, data = Bergh,
                        group = "gender", group.equal = c("loadings", "intercepts"),
                        group.partial = c("GP=~SP", "DP~1", "HP~1", "SP~1"),
                        estimator = "MLR")

GP_model2 <- ' GP =~ c(v1,v1)*EP + c(v2,v2)*HP + c(v3,v3)*DP +
                      c(NA, 0)*SP '
fitIO <- lavaan::cfa(GP_model2, data = Bergh, group = "gender",
                     group.equal = c("intercepts"),
                     group.partial = c("DP~1", "HP~1", "SP~1"),
                     estimator = "MLR")

fitMarg <- lavaan::cfa(GP_model, data = Bergh,group = "gender",
                       group.equal = c("loadings", "intercepts"),
                       group.partial = c("DP~1", "HP~1", "SP~1"),
                       estimator = "MLR")

anova(fitMarg, fitBase1)

# Page 54-55: longitudinal confirmatory factor analysis

library("MPsychoR")
library("lavaan")

data("SDOwave")

model_sdo1 <- '
    SDO1996 =~ 1*I1.1996 + a2*I2.1996 + a3*I3.1996 + a4*I4.1996
    SDO1998 =~ 1*I1.1998 + a2*I2.1998 + a3*I3.1998 + a4*I4.1998
    SDO1996 ~~ SDO1998
    ## intercepts
    I1.1996 ~ int1*1; I1.1998 ~ int1*1
    I2.1996 ~ int2*1; I2.1998 ~ int2*1
    I3.1996 ~ int3*1; I3.1998 ~ int3*1
    I4.1996 ~ int4*1; I4.1998 ~ int4*1
    ## residual covariances
    I1.1996 ~~ I1.1998
    I2.1996 ~~ I2.1998
    I3.1996 ~~ I3.1998
    I4.1996 ~~ I4.1998
    ## latent means: 1996 as baseline
    SDO1996 ~ 0*1
    SDO1998 ~ 1 '

fitsdo1 <- cfa(model_sdo1, data = SDOwave, estimator = "MLR")

fitsdo1 |> 
  parameterEstimates() |> 
  slice(22, 23)

model_sdo2 <- '
    ## 1st CFA level, constant loadings across time
    SDOD1996 =~ 1*I1.1996 + d1*I2.1996
    SDOD1998 =~ 1*I1.1998 + d1*I2.1998
    SDOD1999 =~ 1*I1.1999 + d1*I2.1999
    SDOE1996 =~ 1*I3.1996 + a1*I4.1996
    SDOE1998 =~ 1*I3.1998 + a1*I4.1998
    SDOE1999 =~ 1*I3.1999 + a1*I4.1999
    ## 2nd CFA level, constant loadings across time
    SDO1996 =~ 1*SDOD1996 + sd1*SDOE1996
    SDO1998 =~ 1*SDOD1998 + sd1*SDOE1998
    SDO1999 =~ 1*SDOD1999 + sd1*SDOE1999
    ## Constant 1st level intercepts
    I1.1996 ~ iI1*1; I1.1998 ~ iI1*1; I1.1999 ~ iI1*1
    I2.1996 ~ iI2*1; I2.1998 ~ iI2*1; I2.1999 ~ iI2*1
    I3.1996 ~ iI3*1; I3.1998 ~ iI3*1; I3.1999 ~ iI3*1
    I4.1996 ~ iI4*1; I4.1998 ~ iI4*1; I4.1999 ~ iI4*1
    ## residual covariances
    I1.1999 ~~ I1.1998; I1.1996 ~~ I1.1998; I1.1999 ~~ I1.1996
    I2.1999 ~~ I2.1998; I2.1996 ~~ I2.1998; I2.1999 ~~ I2.1996
    I3.1999 ~~ I3.1998; I3.1996 ~~ I3.1998; I3.1999 ~~ I3.1996
    I4.1999 ~~ I4.1998; I4.1996 ~~ I4.1998; I4.1999 ~~ I4.1996
    ## latent means
    SDO1996 ~ 0*1
    SDO1998 ~ 1
    SDO1999 ~ 1
    ## 1996 baseline year
    ## 1998 vs. 1996
    ## 1999 vs. 1996 '

fitsdo2 <- cfa(model_sdo2, data = SDOwave, estimator = "MLR")

fitsdo2 |> 
  parameterEstimates() |> 
  slice(43:45)

# Page 56-57: multilevel confirmatory factor analysis

data("FamilyIQ")

modelIQ <- '
   level: 1
    numeric =~ wordlist + cards + matrices
    perception =~ figures + animals + occupation
   level: 2
    general =~ wordlist + cards +  matrices + figures + animals +
                                                     occupation '

fitIQ <- cfa(modelIQ, data = FamilyIQ, cluster = "family", std.lv = TRUE)
fitIQ

# Page 57-59: Bayesian confirmatory factor analysis

library("blavaan")

dpriors()[c("lambda", "theta", "psi")]

library("MPsychoR")

data("Bergh")

GP_model <- 'GP =~ EP + HP + DP + SP'

set.seed(123)

fitBCFA <- bcfa(GP_model, data = Bergh, burnin = 2000,
                sample = 10000, n.chains = 2)

plot(fitBCFA, pars = 1:2, plot.type = "trace")
plot(fitBCFA, pars = 1:2, plot.type = "autocorr")

summary(fitBCFA)