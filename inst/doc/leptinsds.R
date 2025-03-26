## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(childsds)
library(dplyr)
data(leptin.ref)

## ----include=TRUE-------------------------------------------------------------
sds_2d(value = 20, age = 10, x2 = 1, sex = "male", item = "leptin_until6", ref = leptin.ref)

## ----include=TRUE-------------------------------------------------------------
df <- data.frame(age = seq(0.5, 6.5, by = 1),
                 sex = sample(c("m","f"),7, replace = T),
                 bmisds = rnorm(7),
                 leptin = runif(7, min = 0.01, max = 5))
df

## ----include=TRUE-------------------------------------------------------------

df$leptin_sds <- sds_2d(value = df$leptin,
                       age = df$age,
                       x2 = df$bmisds,
                       sex = df$sex, male = "m", female = "f",
                       item = "leptin_until6",
                       ref = leptin.ref)
df

## ----include=TRUE-------------------------------------------------------------
sds_pub2d(value = 20,
          pubstat = 2,
          x2 = 1,
          sex = "male",
          item = "lep_pub", ref = leptin.ref)

## ----include=TRUE-------------------------------------------------------------
df <- data.frame(age = seq(0.5, 6.5, by = 1),
                 sex = sample(c("m","f"),7, replace = T),
                 bmisds = rnorm(7),
                 leptin = runif(7, min = 0.01, max = 5))
df

## ----include=TRUE-------------------------------------------------------------
df$leptin_sds <- sds_2d(value = df$leptin,
                       age = df$age,
                       x2 = df$bmisds,
                       sex = df$sex, male = "m", female = "f",
                       item = "leptin_until6",
                       ref = leptin.ref)
df

## ----include=TRUE-------------------------------------------------------------
df$leptin_perc <- sds_2d(value = df$leptin,
                       age = df$age,
                       x2 = df$bmisds,
                       sex = df$sex, male = "m", female = "f",
                       item = "leptin_until6",
                       type = "perc",
                       ref = leptin.ref)
df

## ----include=TRUE-------------------------------------------------------------
sds_2d(value = 20, age = 20, x2 = 25, sex = "male", item = "lep_bmi", ref = leptin.ref)

## ----include=TRUE-------------------------------------------------------------
df <- data.frame(age = seq(20, 80, by = 10),
                 sex = sample(c("M","F"),7, replace = T),
                 bmi = runif(7, 20, 40),
                 leptin = runif(7, min = 0.01, max = 20))
df

## ----include=TRUE-------------------------------------------------------------

df$leptin_sds <- sds_2d(value = df$leptin,
                       age = df$age,
                       x2 = df$bmi,
                       sex = df$sex, male = "M", female = "F",
                       item = "lep_bmi",
                       ref = leptin.ref)
df

