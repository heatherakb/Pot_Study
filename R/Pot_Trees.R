# Potted Trees 2025 data analysis
## 1. Load and clean ####
library(readr)
library(lubridate)
library(purrr)
library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)   # p-values for lmer / Type III anova
library(emmeans)
library(multcomp)
library(multcompView)
library(ggplot2)
library(tidyr)
library(here)
library(dplyr)

project_root <- here::here()
data_dir <- file.path(project_root, "data", "raw")
figure_dir <- file.path(project_root, "outputs", "figures")
table_dir <- file.path(project_root, "outputs", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

treepot<-read.csv(file.path(data_dir, "Biomass_and_growth_and_chemistry_feb_26.csv"), fileEncoding = "latin1")
treepot <- treepot[rowSums(is.na(treepot) | treepot == "") < ncol(treepot), ]

treepot <- treepot |>
  mutate(
    Species.latin = factor(trimws(Species.latin)),
    Salt    = factor(Salinity, levels = c("N", "L", "H")),
    Biochar = factor(Soil, levels = c("S", "B")),
    ID      = factor(ID),
    Year    = 2025,
    vigor   = Height.cm / Total.diamm.at.15mm,
    Block = substr(as.character(Block), 1, 1),
    Block.last  = substr(as.character(Block), nchar(as.character(Block)), nchar(as.character(Block))),
    Crown.diam.avg.diff.25.23 = (Crown.diam.1.Diff.25.23- Crown.Diameter.2.Diff.25.23)/2
  )


## Growth traits

##ID vars
V_ID<-c("ID", "Species.latin", "Salt", "Biochar", "Block", "Samples")

##Growth vars
V_growth<-c("Total.above.ground.biomass.g", "Below.ground.biomass.dry.g", "Root.shoot.Ratio", "Total.diamm.at.15mm","Diam.Diff.25.23", "Height.cm", "Height.Difference.25.23","Crown.diameter.average", "Crown.diam.avg.diff.25.23")
treepot |>
  group_by(Species.latin) |>
  summarise(across(V_growth, ~ sum(!is.na(.x)))) |>
  print()

##Chemistry vars
V_chem<-c("NITROGEN_ASIS", "P.", "K.", "MG.", "CA.", "NA.2","CL.2")
treepot |>
  group_by(Species.latin) |>
  summarise(across(all_of(V_chem), ~ sum(!is.na(.x)))) |>
  print()
    #note nitrogen is an incomplete dataset, because we were unable to collect enough leaves for some species. 
#Leaf condition vars
V_lc<-c("Licor.Aug.07","Licor.July.27","LICOR.Avg", "leaf.burns...leaves.on.tree","Leaf.burns...of.leaf.burned", "X..discolored","X..crown.missing")

treepot |>
  group_by(Species.latin) |>
  summarise(across(all_of(V_lc), ~ sum(!is.na(.x)))) |>
  print()
  #one tree dead (Poplus deltoides, thus no LICOR on July 27th)
#tree condition ratings
V_tc <- treepot |>
  dplyr::select(Bris_cl0.5spring:etat_general_cl0.4summer) |>
  names()

treepot<-treepot|>
  mutate(
    across(all_of(V_ID), as.factor),
    across(all_of(V_lc), as.numeric),
    across(all_of(V_chem), as.numeric),
    across(all_of(V_tc), as.integer),
    across(all_of(V_growth), as.numeric)
  )



## Saving cleaned dataset for possible future analysis
treeb_cleaned<-treeb
write.csv(treeb_cleaned, file.path("data", "cleaned", "treeb_cleaned.csv"))
save(treeb_cleaned,
     file = file.path("data", "cleaned", "treeb_cleaned.RData"))
