# Potted Trees 2025 data analysis-----
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
    Block.old  = substr(as.character(Block), nchar(as.character(Block)), nchar(as.character(Block))),
    Block = substr(as.character(Block), 1, 1),
    Crown.diam.avg.diff.25.23 = (Crown.diam.1.Diff.25.23 + Crown.Diameter.2.Diff.25.23)/2,
    Block.first = substr(as.character(Block), 1, 1),
    Sp = substr(as.character(ID), 1, 2),
  )


## Growth traits

##ID vars
V_ID<-c("ID","Sp" , "Salt", "Biochar", "Block", "Samples", "Block.old","Species.latin")

##Growth vars
V_growth<-c("Total.above.ground.biomass.g", "Below.ground.biomass.dry.g", "Root.shoot.Ratio", "Total.diamm.at.15mm","Diam.Diff.25.23", "Height.cm", "Height.Difference.25.23","Crown.diameter.average", "Crown.diam.avg.diff.25.23","vigor")
g<-treepot |>
  group_by(Species.latin) |>
  summarise(across(all_of(V_growth), ~ sum(!is.na(.x))))

  print(g)

##Chemistry vars
V_chem<-c("NITROGEN_ASIS", "P.", "K.", "MG.", "CA.", "NA.2","CL.2")
treepot |>
  group_by(Species.latin) |>
  summarise(across(all_of(V_chem), ~ sum(!is.na(.x)))) |>
  print()
    #note nitrogen is an incomplete dataset, because we were unable to collect enough leaves for some species. 
#Leaf condition vars
V_lc<-c("Licor.Aug.07","Licor.July.27","LICOR.Avg", "leaf.burns...leaves.on.tree","Leaf.burns...of.leaf.burned", "X..discolored","X..crown.missing")

lc<-treepot |>
  group_by(Species.latin) |>
  summarise(across(all_of(V_lc), ~ sum(!is.na(.x))))
  print(lc)
  #one tree dead (Poplus deltoides, during season, and Gelditsia triacanthos has leaves to small to scan with LICOR)
#tree condition ratings
V_tc <- treepot |>
  dplyr::select(Bris_cl0.5spring:etat_general_cl0.4summer) |>
  names()
tc<-treepot |>
  group_by(Species.latin) |>
  summarise(across(all_of(V_tc), ~ sum(!is.na(.x))))
print(tc)

treepot <- treepot |>
  mutate(
    across(all_of(V_ID), as.factor),
    across(all_of(V_lc), as.numeric),
    across(all_of(V_chem), as.numeric),
    across(all_of(V_tc), as.integer),
    across(all_of(V_growth), as.numeric)
  )

treepot_all<-treepot
treepot <- treepot |>
  dplyr::select(all_of(c(V_ID, V_growth, V_chem, V_lc, V_tc)))


## Saving cleaned dataset for possible future analysis
write.csv(treepot, file.path("data", "cleaned", "treepot_cleaned.csv"))
save(treepot,
     file = file.path("data", "cleaned", "treepot_cleaned.RData"))
write.csv(treepot_all, file.path("data", "cleaned", "treepot_all.csv"))
save(treepot_all,
     file = file.path("data", "cleaned", "treepot_all.RData"))


##2. Growth ----
#2 types of growth data, whole, and whole-initial (2023 measurements)
#should we remove one species?

#growth vars are V_growth
#equation
rhs<-"Salt*Biochar"
growth_models<- list()

#lmer model loop for growth models
library(lmerTest)   # lmer with p-values
library(emmeans)
species_list <- split(treepot, droplevels(treepot$Sp))
rhs <- "Salt * Biochar"
growth_models <- list()

# traits where zero/negative values are impossible (not the *Diff* columns)
V_positive <- setdiff(V_growth, grep("Diff", V_growth, value = TRUE))

options(contrasts = c("contr.sum", "contr.poly"))   # for valid Type III tests
alpha <- 0.1
growth_models <- list()
failed <- list()

# helper: print estimates + Tukey pairwise for one emmeans spec
show_emm <- function(m, specs, label) {
  emm <- emmeans(m, specs)
  cat("\n--- emmeans:", label, "---\n")
  print(summary(emm, infer = TRUE))
  print(pairs(emm, adjust = "tukey"))
}

for (sp in names(species_list)) {
  cat("\n\n##################################\n")
  cat("Species:", sp, "\n")
  cat("##################################\n")
  
  sp_dat <- species_list[[sp]]
  growth_models[[sp]] <- list()
  
  for (v in V_growth) {
    cat("\n==================================\n")
    cat(sp, "|", "Growth trait:", v, "\n")
    cat("==================================\n")
    
    tryCatch({
      d <- sp_dat |> filter(!is.na(.data[[v]]))
      if (v %in% V_positive) d <- d |> filter(.data[[v]] > 0)
      
      f <- as.formula(paste0(v, " ~ ", rhs, " + (1|Block)"))
      m <- lmer(f, data = d, REML = TRUE)
      growth_models[[sp]][[v]] <- m
      
      a <- anova(m, type = 3)
      print(a)
      
      p <- a[["Pr(>F)"]]
      names(p) <- rownames(a)
      is_sig <- function(term) isTRUE(unname(p[term]) < alpha)   # NA/missing -> FALSE
      
      if (is_sig("Salt"))         show_emm(m, ~ Salt,           "Salt")
      if (is_sig("Biochar"))      show_emm(m, ~ Biochar,        "Biochar")
      if (is_sig("Salt:Biochar")) {
        show_emm(m, ~ Salt | Biochar, "Salt | Biochar")
        show_emm(m, ~ Biochar | Salt, "Biochar | Salt")
      }
      
      if (!any(sapply(c("Salt", "Biochar", "Salt:Biochar"), is_sig)))
        cat("\n(no term significant at p <", alpha, "- emmeans skipped)\n")
      
    }, error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      failed[[length(failed) + 1]] <<- data.frame(
        Species = sp, Trait = v, Error = conditionMessage(e)
      )
    })
  }
}

do.call(rbind, failed)   # which species/trait combos didn't fit
alpha <- 0.1

for (sp in names(growth_models)) {
  for (v in names(growth_models[[sp]])) {
    a <- anova(growth_models[[sp]][[v]], type = 3)
    p <- a[["Pr(>F)"]]
    
    if (any(p < alpha, na.rm = TRUE)) {
      cat("\n=== ANOVA (Type III, Satterthwaite):", sp, "|", v, "===\n")
      print(a[!is.na(p) & p < alpha, ])
    }
  }
}

anova_tbl <- do.call(rbind, lapply(names(growth_models), function(sp) {
  do.call(rbind, lapply(names(growth_models[[sp]]), function(v) {
    a <- as.data.frame(anova(growth_models[[sp]][[v]], type = 3))
    a$Term    <- rownames(a)
    a$Species <- sp
    a$Trait   <- v
    a
  }))
}))
rownames(anova_tbl) <- NULL

# formatted copy for display/export; numeric Pr(>F) is kept for filtering
anova_tbl$'Pr(>F)' <- ifelse(is.na(anova_tbl$`Pr(>F)`),
                            NA_character_,
                            formatC(anova_tbl$`Pr(>F)`, format = "f", digits = 9))



##3. Chemistry ####
#same as above species list and rhs
#species_list <- split(treepot, droplevels(treepot$Sp))
#rhs <- "Salt * Biochar"
chem_models <- list()

options(contrasts = c("contr.sum", "contr.poly"))   # for valid Type III tests
alpha <- 0.1
chem_models <- list()
failed <- list()

# helper: print estimates + Tukey pairwise for one emmeans spec
show_emm <- function(m, specs, label) {
  emm <- emmeans(m, specs)
  cat("\n--- emmeans:", label, "---\n")
  print(summary(emm, infer = TRUE))
  print(pairs(emm, adjust = "tukey"))
}

for (sp in names(species_list)) {
  cat("\n\n##################################\n")
  cat("Species:", sp, "\n")
  cat("##################################\n")
  
  sp_dat <- species_list[[sp]]
  chem_models[[sp]] <- list()
  
  for (v in V_chem) {
    cat("\n==================================\n")
    cat(sp, "|", "Leaf chemistry %:", v, "\n")
    cat("==================================\n")
    
    tryCatch({
      d <- sp_dat |> filter(!is.na(.data[[v]]))
      if (v %in% V_chem) d <- d |> filter(.data[[v]] > 0)
      
      f <- as.formula(paste0(v, " ~ ", rhs, " + (1|Block)"))
      m <- lmer(f, data = d, REML = TRUE)
      chem_models[[sp]][[v]] <- m
      
      a <- anova(m, type = 3)
      print(a)
      
      p <- a[["Pr(>F)"]]
      names(p) <- rownames(a)
      is_sig <- function(term) isTRUE(unname(p[term]) < alpha)   # NA/missing -> FALSE
      
      if (is_sig("Salt"))         show_emm(m, ~ Salt,           "Salt")
      if (is_sig("Biochar"))      show_emm(m, ~ Biochar,        "Biochar")
      if (is_sig("Salt:Biochar")) {
        show_emm(m, ~ Salt | Biochar, "Salt | Biochar")
        show_emm(m, ~ Biochar | Salt, "Biochar | Salt")
      }
      
      if (!any(sapply(c("Salt", "Biochar", "Salt:Biochar"), is_sig)))
        cat("\n(no term significant at p <", alpha, "- emmeans skipped)\n")
      
    }, error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      failed[[length(failed) + 1]] <<- data.frame(
        Species = sp, Trait = v, Error = conditionMessage(e)
      )
    })
  }
}

do.call(rbind, failed)   # which species/trait combos didn't fit
alpha <- 0.1

for (sp in names(chem_models)) {
  for (v in names(chem_models[[sp]])) {
    a <- anova(chem_models[[sp]][[v]], type = 3)
    p <- a[["Pr(>F)"]]
    
    if (any(p < alpha, na.rm = TRUE)) {
      cat("\n=== ANOVA (Type III, Satterthwaite):", sp, "|", v, "===\n")
      print(a[!is.na(p) & p < alpha, ])
    }
  }
}

anova_tbl_chem <- do.call(rbind, lapply(names(chem_models), function(sp) {
  do.call(rbind, lapply(names(chem_models[[sp]]), function(v) {
    a <- as.data.frame(anova(chem_models[[sp]][[v]], type = 3))
    a$Term    <- rownames(a)
    a$Species <- sp
    a$Trait   <- v
    a
  }))
}))
rownames(anova_tbl_chem) <- NULL

# formatted copy for display/export; numeric Pr(>F) is kept for filtering
anova_tbl_chem$'Pr(>F)' <- ifelse(is.na(anova_tbl_chem$`Pr(>F)`),
                             NA_character_,
                             formatC(anova_tbl_chem$`Pr(>F)`, format = "f", digits = 9))
##4. Leaf health models####
lf_models <- list()

options(contrasts = c("contr.sum", "contr.poly"))   # for valid Type III tests
alpha <- 0.1
lf_models <- list()
failed <- list()

# helper: print estimates + Tukey pairwise for one emmeans spec
show_emm <- function(m, specs, label) {
  emm <- emmeans(m, specs)
  cat("\n--- emmeans:", label, "---\n")
  print(summary(emm, infer = TRUE))
  print(pairs(emm, adjust = "tukey"))
}

for (sp in names(species_list)) {
  cat("\n\n##################################\n")
  cat("Species:", sp, "\n")
  cat("##################################\n")
  
  sp_dat <- species_list[[sp]]
  lf_models[[sp]] <- list()
  
  for (v in V_lc) {
    cat("\n==================================\n")
    cat(sp, "|", "Leaf chemistry %:", v, "\n")
    cat("==================================\n")
    
    tryCatch({
      d <- sp_dat |> filter(!is.na(.data[[v]]))
      if (v %in% V_lc) d <- d |> filter(.data[[v]] > 0)
      
      f <- as.formula(paste0(v, " ~ ", rhs, " + (1|Block)"))
      m <- lmer(f, data = d, REML = TRUE)
      lf_models[[sp]][[v]] <- m
      
      a <- anova(m, type = 3)
      print(a)
      
      p <- a[["Pr(>F)"]]
      names(p) <- rownames(a)
      is_sig <- function(term) isTRUE(unname(p[term]) < alpha)   # NA/missing -> FALSE
      
      if (is_sig("Salt"))         show_emm(m, ~ Salt,           "Salt")
      if (is_sig("Biochar"))      show_emm(m, ~ Biochar,        "Biochar")
      if (is_sig("Salt:Biochar")) {
        show_emm(m, ~ Salt | Biochar, "Salt | Biochar")
        show_emm(m, ~ Biochar | Salt, "Biochar | Salt")
      }
      
      if (!any(sapply(c("Salt", "Biochar", "Salt:Biochar"), is_sig)))
        cat("\n(no term significant at p <", alpha, "- emmeans skipped)\n")
      
    }, error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      failed[[length(failed) + 1]] <<- data.frame(
        Species = sp, Trait = v, Error = conditionMessage(e)
      )
    })
  }
}

do.call(rbind, failed)   # which species/trait combos didn't fit
alpha <- 0.1

for (sp in names(lf_models)) {
  for (v in names(lf_models[[sp]])) {
    a <- anova(lf_models[[sp]][[v]], type = 3)
    p <- a[["Pr(>F)"]]
    
    if (any(p < alpha, na.rm = TRUE)) {
      cat("\n=== ANOVA (Type III, Satterthwaite):", sp, "|", v, "===\n")
      print(a[!is.na(p) & p < alpha, ])
    }
  }
}

anova_tbl_lf <- do.call(rbind, lapply(names(lf_models), function(sp) {
  do.call(rbind, lapply(names(lf_models[[sp]]), function(v) {
    a <- as.data.frame(anova(lf_models[[sp]][[v]], type = 3))
    a$Term    <- rownames(a)
    a$Species <- sp
    a$Trait   <- v
    a
  }))
}))
rownames(anova_tbl_lf) <- NULL

# formatted copy for display/export; numeric Pr(>F) is kept for filtering
anova_tbl_lf$'Pr(>F)' <- ifelse(is.na(anova_tbl_lf$`Pr(>F)`),
                                  NA_character_,
                                  formatC(anova_tbl_lf$`Pr(>F)`, format = "f", digits = 9))
##5. Tree health####

library(lmerTest)
library(emmeans)
library(ARTool)

alpha <- 0.1
tc_models <- list()
failed <- list()

# ART main-effect contrasts (Tukey)
show_art_main <- function(m, term) {
  cat("\n--- art.con:", term, "---\n")
  print(art.con(m, term, adjust = "tukey"))
}

# ART conditional contrasts for the interaction (e.g. Salt | Biochar)
show_art_cond <- function(m, specs, label) {
  cat("\n--- art contrasts:", label, "---\n")
  emm <- emmeans(artlm(m, "Salt:Biochar"), specs)
  print(pairs(emm, adjust = "tukey"))
}

# term names live in a "Term" column or in the rownames, depending on ARTool version
get_terms <- function(a) if ("Term" %in% names(a)) as.character(a$Term) else rownames(a)

for (sp in names(species_list)) {
  cat("\n\n##################################\n")
  cat("Species:", sp, "\n")
  cat("##################################\n")
  
  sp_dat <- species_list[[sp]]
  tc_models[[sp]] <- list()
  
  for (v in V_tc) {
    cat("\n==================================\n")
    cat(sp, "|", "Tree Rating Scale:", v, "\n")
    cat("==================================\n")
    
    tryCatch({
      d <- sp_dat |> filter(!is.na(.data[[v]]))
      d$Block <- droplevels(factor(d$Block))
      
      # informative guard for the "grouping factor levels" error
      if (nlevels(d$Block) >= nrow(d))
        stop("one tree per Block (n = ", nrow(d), "); random intercept not estimable")
      
      m <- art(as.formula(paste0(v, " ~ Salt * Biochar + (1|Block)")), data = d)
      a <- anova(m)
      tc_models[[sp]][[v]] <- m
      tc_anova[[sp]][[v]] <- as.data.frame(a)
     
      p <- a[["Pr(>F)"]]
      names(p) <- get_terms(a)
      is_sig <- function(term) isTRUE(unname(p[term]) < alpha)
      
      if (any(sapply(c("Salt", "Biochar", "Salt:Biochar"), is_sig))) {
        print(a)
        if (is_sig("Salt"))    show_art_main(m, "Salt")
        if (is_sig("Biochar")) show_art_main(m, "Biochar")
        if (is_sig("Salt:Biochar")) {
          show_art_cond(m, ~ Salt | Biochar, "Salt | Biochar")
          show_art_cond(m, ~ Biochar | Salt, "Biochar | Salt")
        }
      } else {
        cat("(no term significant at p <", alpha, ")\n")
      }
      
    }, error = function(e) {
      cat("FAILED:", conditionMessage(e), "\n")
      failed[[length(failed) + 1]] <<- data.frame(
        Species = sp, Trait = v, Error = conditionMessage(e)
      )
    })
  }
}

do.call(rbind, failed)   # which species/trait combos didn't fit

# ---- significant-only summary 
for (sp in names(tc_models)) {
  for (v in names(tc_models[[sp]])) {
    a <- anova(tc_models[[sp]][[v]])
    p <- a[["Pr(>F)"]]
    if (any(p < alpha, na.rm = TRUE)) {
      cat("\n=== ART ANOVA:", sp, "|", v, "===\n")
      print(a[!is.na(p) & p < alpha, ])
    }
  }
}

# ---- full table
anova_tbl_tc <- do.call(rbind, lapply(names(tc_anova), function(sp) {
  do.call(rbind, lapply(names(tc_anova[[sp]]), function(v) {
    a <- tc_anova[[sp]][[v]]
    a$Term    <- get_terms(a)
    a$Species <- sp
    a$Trait   <- v
    a
  }))
}))
rownames(anova_tbl_tc) <- NULL

anova_tbl_tc$p_value <- ifelse(is.na(anova_tbl_tc$`Pr(>F)`), NA_character_,
                               formatC(anova_tbl_tc$`Pr(>F)`, format = "f", digits = 9))
