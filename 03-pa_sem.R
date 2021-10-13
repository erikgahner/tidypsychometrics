# Chapter 3: Path Analysis and Structural Equation Models

library("tidyverse")

# Page 64: MANOVA

library("MPsychoR")

data("Bergh")

## fit two multiple regression models
fitmvreg <- lm(cbind(EP, DP) ~ A1 + A2 + O1 + O2, data = Bergh)

library("car")

Manova(fitmvreg)

# Page 65: path model

library("lavaan")

mvreg_model <- '
    EP ~ b11*A1 + b12*A2 + b13*O1 + b14*O2
    DP ~ b21*A1 + b22*A2 + b23*O1 + b24*O2 '

fitmvreg2 <- sem(mvreg_model, data = Bergh)

library("semPlot")

semPaths(fitmvreg2, what = "est", edge.label.cex = 1,
         layout = "tree", residuals = FALSE, edge.color = 1,
         esize = 1, rotation = 3, sizeMan = 8, asize = 2.5,
         fade = FALSE, optimizeLatRes = TRUE)


# Page 68-69: moderation model

library("MPsychoR")

data("Paskvan")

Paskvan <- Paskvan |> 
  mutate(wintense_c = scale(wintense, scale = FALSE),
         pclimate_c = scale(pclimate, scale = FALSE))

fit_YX <- lm(cogapp ~ wintense_c, data = Paskvan)  ## Y on X

fit_YX |> 
  summary() |> 
  pluck("coefficients") |> 
  round(4)

fit_YZ <- lm(cogapp ~ pclimate_c, data = Paskvan)  ## Y on Z

fit_YZ |> 
  summary() |> 
  pluck("coefficients") |> 
  round(4)

library("QuantPsyc")

fit_mod <- moderate.lm(x = wintense, z = pclimate, y = cogapp, data = Paskvan)

fit_mod |> 
  summary() |> 
  pluck("coefficients") |> 
  round(4)

fit_ss <- Paskvan |> 
  pull(pclimate) |> 
  sim.slopes(mod = fit_mod)

round(fit_ss, 4)

# Page 71-72: mediator model


## model with mediation package

library("mediation")

fit_MX <- lm(cogapp ~ wintense, data = Paskvan)
fit_YXM <- lm(emotion ~ wintense + cogapp, data = Paskvan)

set.seed(123)

fitmed <- mediation::mediate(fit_MX, fit_YXM,
                             treat = "wintense", mediator = "cogapp",
                             sims = 999, boot = TRUE, boot.ci.type = "bca")


## model with lavaan package

library("lavaan")

med_model <- '
    emotion ~ c*wintense + b*cogapp
    cogapp ~ a*wintense
    ind := a*b
    tot := ind+c
    prop := ind/tot '

set.seed(123)

fitmedsem <- lavaan::sem(med_model, Paskvan, se = "bootstrap", bootstrap = 999)

parameterEstimates(fitmedsem, zstat = FALSE, pvalue = FALSE, boot.ci.type = "bca.simple") |> 
  slice(7, 1, 8, 9)


# Page 74-75: combined moderator-mediator model

quantile(Paskvan$pclimate)

medmod_model <- '
  ## set of regressions
  cogapp ~ a1*wintense + a2*pclimate + a3*wintense:pclimate
  emotion ~ c*wintense + b*cogapp
  ## conditional indirect effects 
  cie.q1 := (a1 + a3*2)*b   ## first quartile
  cie.q2 := (a1 + a3*3)*b   ## median
  cie.q3 := (a1 + a3*3.5)*b ## third quartile '

set.seed(123)

fitmedmod <- lavaan::sem(medmod_model, data = Paskvan,
                         se = "bootstrap", bootstrap = 999)

semPaths(fitmedmod, layout = "spring", asize = 2.5,
         sizeMan = 10, residuals = FALSE, nCharNodes = 7,
         edge.label.cex = 1)

parameterEstimates(fitmedmod, zstat = FALSE, pvalue = FALSE,
                   boot.ci.type = "bca.simple") |> 
  slice(3, 4, 14:16) 

# Page 77-78: Structural Equation Models

library("MPsychoR")
library("lavaan")

data("Bergh")

Bergh_model <- ' GP =~ EP + HP + DP + SP
                 Agree =~ A1 + A2 + A3
                 Open =~ O1 + O2 + O3
                 GP ~ Agree + Open '

fitGP <- sem(Bergh_model, data = Bergh, estimator = "MLR")

semPaths(fitGP, what = "std", edge.label.cex = 0.7, esize = 1,
         intercepts = FALSE,rotation = 4, edge.color = 1, asize = 2.5,
         sizeMan = 5, mar = c(1, 1.5, 1.5, 3), fade = FALSE)

summary(fitGP, standardized = TRUE, fit.measures = TRUE)


# Page 79-81: multigroup SEM

fit_free <- sem(Bergh_model, group = "gender",
                group.equal = c("intercepts"),
                group.partial = c("DP~1", "HP~1", "SP~1"),
                data = Bergh, estimator = "MLR")

fit_load <- sem(Bergh_model, group = "gender",
                group.equal = c("loadings", "intercepts"),
                group.partial = c("GP=~SP", "DP~1", "HP~1", "SP~1"),
                data = Bergh, estimator = "MLR")

fit_prestrict <- sem(Bergh_model, group = "gender",
                     group.equal = c("intercepts", "regressions"),
                     group.partial = c("DP~1", "HP~1", "SP~1"),
                     data = Bergh, estimator = "MLR")

anova(fit_free, fit_prestrict)

library("nonnest2")

fit_load1 <- update(fit_load, estimator = "ML")

fit_prestrict1 <- update(fit_prestrict, estimator = "ML")

compIC <- icci(fit_load1, fit_prestrict1)
compIC

## Vuong’s non-nested testing strategy
vuongtest(fit_load1, fit_prestrict1)

# Page 85: latent growth models
library("lavaan")
library("aspect")
library("semPlot")

data("duncan")

model_shape <- '
     inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
     shape =~ 0*CIG_T1 + 1*CIG_T2 + CIG_T3 + CIG_T4 '

fitCig1 <- growth(model_shape, data = duncan, estimator = "WLS")

semPaths(fitCig1, what = "std", edge.label.cex = 0.7, esize = 1,
         edge.color = 1, sizeMan = 6, asize = 2.5, intercepts = FALSE,
         rotation = 4, mar = c(3, 5, 3.5, 5), fade = FALSE)

summary(fitCig1, header = FALSE)

model_lin <- '
    inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
    linear =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4 '
fitCig2 <- growth(model_lin, data = duncan, estimator = "WLS") 

parameterEstimates(fitCig2) |> 
  slice(21)

fitMeasures(fitCig2)[c("rmsea", "cfi", "srmr")] |> 
  round(3)

model_quad <- '
    inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
    linear =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
    quad =~ 0*CIG_T1 + 1*CIG_T2 + 4*CIG_T3 + 9*CIG_T4 '

fitCig3 <- growth(model_quad, data = duncan, estimator = "WLS")

parameterEstimates(fitCig3) |> 
  slice(28, 29)

fitMeasures(fitCig3)[c("rmsea", "cfi", "srmr")] |> 
  round(3)

anova(fitCig1, fitCig2)

# Page 86-90: extended latent growth modeling
model_ord <- '
    inter =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
    linear =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
    CIG_T1 | 0*t1 + t2 + t3 + t4
    CIG_T2 | 0*t1 + t2 + t3 + t4 '

fitCigord <- growth(model_ord, data = duncan,
                    ordered = names(duncan)[5:8])


model_pc <- '
    cint =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
    clin =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
    pint =~ 1*POT_T1 + 1*POT_T2 + 1*POT_T3 + 1*POT_T4
    plin =~ 0*POT_T1 + 1*POT_T2 + 2*POT_T3 + 3*POT_T4
    ## correlated errors
    CIG_T1 ~~ CIG_T2; CIG_T2 ~~ CIG_T3; CIG_T3 ~~ CIG_T4
    POT_T1 ~~ POT_T2; POT_T2 ~~ POT_T3; POT_T3 ~~ POT_T4
    ## fix error variances
    CIG_T1 ~~ rc*CIG_T1
    CIG_T2 ~~ rc*CIG_T2
    CIG_T3 ~~ rc*CIG_T3
    CIG_T4 ~~ rc*CIG_T4
    POT_T1 ~~ rp*POT_T1
    POT_T2 ~~ rp*POT_T2
    POT_T3 ~~ rp*POT_T3
    POT_T4 ~~ rp*POT_T4 '

fitPC1 <- growth(model_pc, data = duncan, estimator = "WLS")

fitMeasures(fitPC1)[c("rmsea", "cfi", "srmr")] |> 
  round(3)

fitPC1 |> 
  inspect("std") |> 
  pluck("psi")

duncan <- duncan |> 
  mutate(ALCavg = rowMeans(dplyr::select(duncan, ALC_T1, ALC_T2, ALC_T3, ALC_T4)))

model_pca1 <- '
    cint =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
    clin =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
    pint =~ 1*POT_T1 + 1*POT_T2 + 1*POT_T3 + 1*POT_T4
    plin =~ 0*POT_T1 + 1*POT_T2 + 2*POT_T3 + 3*POT_T4
    ## effects of alcohol on marijuana
    pint ~ ALCavg
    plin ~ ALCavg '

fitPCA1 <- growth(model_pca1, data = duncan, estimator = "MLR")

fitMeasures(fitPCA1)[c("rmsea", "cfi", "srmr")] |> 
  round(3)

parameterEstimates(fitPCA1) |> 
  slice(17, 18)

model_pca2 <- '   
    cint =~ 1*CIG_T1 + 1*CIG_T2 + 1*CIG_T3 + 1*CIG_T4
    clin =~ 0*CIG_T1 + 1*CIG_T2 + 2*CIG_T3 + 3*CIG_T4
    pint =~ 1*POT_T1 + 1*POT_T2 + 1*POT_T3 + 1*POT_T4
    plin =~ 0*POT_T1 + 1*POT_T2 + 2*POT_T3 + 3*POT_T4
    ## effects of alcohol on marijuana
    POT_T1 ~ ALC_T1
    POT_T2 ~ ALC_T2
    POT_T3 ~ ALC_T3
    POT_T4 ~ ALC_T4 '

fitPCA2 <- growth(model_pca2, data = duncan, estimator = "MLR")

fitMeasures(fitPCA2)[c("rmsea", "cfi", "srmr")] |> 
  round(3)

parameterEstimates(fitPCA2) |> 
  slice(17:20)
