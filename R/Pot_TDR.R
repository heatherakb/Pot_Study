#TDR for Pots####

##task: model the change in salt over season
#1. load and clean data
#2. make models
  #does species play an important role? sep by species or no
  #how do these TDR measurements compare to Belanger tree measures.
    #could plot on the same axes with different colors
    #spots in fosse locations as different shapes, treatments as different shades of blue
    #species as a little letters, treatments different shades of pink
    #
#3. make ggplots
#4. make excel summary

# 1. Set up---------------------------------------
##prep data and libraries####
library(dplyr)
library(emmeans)
library(multcomp)
library(multcompView)
library(officer)
library(flextable)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(glmmTMB)
library(readr)

project_root <- here::here()
data_dir <- file.path(project_root, "data", "raw")
figure_dir <- file.path(project_root, "outputs", "figures")
table_dir <- file.path(project_root, "outputs", "tables")


#path <- file.path(data_dir, "TDR_pots_all_dates.csv")

TDRpot_2025 <- read.csv(file.path(data_dir,"TDR_pots_all_dates.csv"), fileEncoding = "ISO-8859-1")
summary(TDRpot_2025)
TDRpot_2025 <- TDRpot_2025 |>
  mutate(
    Species_latin = ifelse(Species_latin=="na", "No tree", Species_latin),
    Species_latin = factor(trimws(Species_latin)),
    Sp = substr(as.character(CODE2), 1, 3),
    Sp = gsub("-", "", Sp, fixed = TRUE),
    Salt    = factor(Salinity, levels = c("N", "L", "H")),
    Biochar = factor(Soil, levels = c("S", "B")),
    Block.old  = substr(as.character(Bloc), nchar(as.character(Bloc)), nchar(as.character(Bloc))),
    Block = substr(as.character(Bloc), 1, 1),
    CODE2 = as.factor(CODE2),
    ID = paste0(CODE2, " ", Block))




#TDRb_allyrs <- TDRb_allyrs |>   #one sampling effort, had to be split into two days.
  #mutate(Humidity= ifelse(Humidity==330.8, 33.8, Humidity)) #typo in humidity

