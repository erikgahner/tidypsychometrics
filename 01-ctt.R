# Chapter 1: Classical Test Theory

library("tidyverse")

# Page 4-5: Computation of Cronbach’s α

library("MPsychoR")
library("psych")

data("Rmotivation")

## item selection
HybMotivation <- Rmotivation |> 
  select(starts_with("hyb")) |> 
  na.omit()

## number of items
k <- ncol(HybMotivation) 

## calculate the variance-covariance (VC) matrix 
vcmat <- cov(HybMotivation)

## sum of the main diagonal elements
sigma2_Xi <- tr(vcmat)

## total variance
sigma2_X <- sum(vcmat)

## compute α
cronalpha <- k / (k - 1) * (1 - sigma2_Xi / sigma2_X)
round(cronalpha, 2)

alpha_hyb <- psych::alpha(HybMotivation)

alpha_hyb$total |> 
  pull(raw_alpha) |> 
  round(2)

# Page 6: greatest lower bound reliability measure 

library("GPArotation")

glb(HybMotivation)


# Page 7: ω-coefficients

psych::omega(HybMotivation)


# Page 8: data to long format, fixed-effects ANOVA and approximation of Cronbach’s α

Hyb1 <- HybMotivation |>
  mutate(person = as.factor(row_number()))

Hyblong <- Hyb1 |> 
  pivot_longer(cols = starts_with("hyb"), names_to = "item")

Hyblong |> 
  aov(formula = value ~ person + item) |> 
  summary()

round((0.85 - 0.15) / 0.85, 2)

icchyb <- ICC(HybMotivation)

# Page 9: random-effects ANOVA
sqrt((0.85 - 0.15) / 19)

sqrt((31.88 - 0.15) / 777)

library("lme4")

Hyblong |> 
  lmer(formula = value ~ (1 | person) + (1 | item)) |> 
  VarCorr()


# Page 9: generalizability coefficient
library("gtheory")

gfit <- Hyblong |> 
  gstudy(formula = value ~ (1 | person) + (1 | item))

dfit <- dstudy(gfit, colname.objects = "person",
               colname.scores = "value", data = as.data.frame(Hyblong))

round(dfit$generalizability, 3)

# Page 10: data preparation for multi-facet G-theory application
data("Lakes")

phydat <- Lakes |> 
  filter(subtest == "physical") |> 
  mutate(item = droplevels(item))

head(phydat)

# Page 11: multi-facet G-theory application

formula <- score ~ (1 | personID) + ( 1 | raterID) + (1 | item) + (1 | personID:raterID) + (1 | personID:item) + (1 | raterID:item)

gfit <- gstudy(formula = formula, data = phydat)
gfit

# Page 12-13: D-study variance components

dfit <- dstudy(gfit, colname.objects = "personID", colname.scores = "score", data = phydat)

dfit |> 
  pluck("components")

## absolute error variance
dfit |> 
  pluck("var.error.abs")

## absolute standard error of measurement
dfit |> 
  pluck("sem.abs")

## relative standard error of measurement
dfit |> 
  pluck("var.error.rel")

dfit |> 
  pluck("sem.rel")

## dependability coefficient
dfit |> 
  pluck("dependability")

## generalizability coefficient
dfit |> 
  pluck("generalizability")