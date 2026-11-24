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
##4. Leaf health ####
lic_vars <- c("Licor.July.27", "Licor.Aug.07")
V_lc_other <- setdiff(V_lc, c(lic_vars, "LICOR.Avg"))

# (a) your existing loop, but over V_lc_other instead of V_lc


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
  
  for (v in V_lc_other) {
    cat("\n==================================\n")
    cat(sp, "|", "Leaf chemistry %:", v, "\n")
    cat("==================================\n")
    
    tryCatch({
      d <- sp_dat |> filter(!is.na(.data[[v]]))
      if (v %in% V_lc_other) d <- d |> filter(.data[[v]] > 0)
      
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
#     -> change:  for (v in V_lc)  to  for (v in V_lc_other)
#        and:     if (v %in% V_lc)  to  if (v %in% V_lc_other)

# (b) LICOR: one model, both dates, Date as a random effect
for (sp in names(species_list)) {
  tryCatch({
    d <- species_list[[sp]] |>
      dplyr::select(ID, Block, Salt, Biochar, all_of(lic_vars)) |>
      tidyr::pivot_longer(all_of(lic_vars), names_to = "Date", values_to = "LICOR") |>
      filter(!is.na(LICOR), LICOR > 0) |>
      mutate(Date  = factor(Date),
             ID    = droplevels(factor(ID)),
             Block = droplevels(factor(Block)))
    
    m <- lmer(LICOR ~ Salt * Biochar + (1 | Block) + (1 | Date) + (1 | ID),
              data = d, REML = TRUE)
    lf_models[[sp]][["LICOR"]] <- m
    
    cat("\n==================================\n")
    cat(sp, "Stomatal conductance",
        if (isSingular(m)) "| SINGULAR FIT" else "", "\n")
    cat("==================================\n")
    
    a <- anova(m, type = 3)
    p <- a[["Pr(>F)"]]; names(p) <- rownames(a)
    is_sig <- function(term) isTRUE(unname(p[term]) < alpha)
    
    if (any(sapply(c("Salt", "Biochar", "Salt:Biochar"), is_sig))) {
      print(a)
      if (is_sig("Salt"))         show_emm(m, ~ Salt,    "Salt")
      if (is_sig("Biochar"))      show_emm(m, ~ Biochar, "Biochar")
      if (is_sig("Salt:Biochar")) {
        show_emm(m, ~ Salt | Biochar, "Salt | Biochar")
        show_emm(m, ~ Biochar | Salt, "Biochar | Salt")
      }
    } else cat("(no term significant at p <", alpha, ")\n")
    
  }, error = function(e) {
    cat("FAILED:", sp, "| LICOR -", conditionMessage(e), "\n")
    failed[[length(failed) + 1]] <<- data.frame(
      Species = sp, Trait = "LICOR", Error = conditionMessage(e))
  })
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
#Aligned Rank Transformation to treat ordinal data for anova.
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
# 6. plotting results####
short_names <- c(
  # growth
  "Total.above.ground.biomass.g" = "Shoot mass (g)",
  "Below.ground.biomass.dry.g"   = "Root mass (g)",
  "Root.shoot.Ratio"             = "Root:shoot",
  "Total.diamm.at.15mm"          = "Diameter",
  "Diam.Diff.25.23"              = "Diam. change",
  "Height.cm"                    = "Height (cm)",
  "Height.Difference.25.23"      = "Height change",
  "Crown.diameter.average"       = "Crown diam.",
  "Crown.diam.avg.diff.25.23"    = "Crown change",
  "vigor"                        = "Vigor",
  # chemistry
  "NITROGEN_ASIS" = "N", "P." = "P", "K." = "K", "MG." = "Mg",
  "CA." = "Ca", "NA.2" = "Na", "CL.2" = "Cl",
  # leaf health
  "LICOR"                        = "LICOR",
  "leaf.burns...leaves.on.tree"  = "burns % per tree",
  "Leaf.burns...of.leaf.burned"  = "burns % leaf surface",
  "X..discolored"                = "Discolored (%)",
  "X..crown.missing"             = "Crown missing (%)",
  # tree health ratings (French -> English)
  "Bris_cl0.5spring"           = "Breakage (spr)",
  "Defoliation_cl_0.4spring"   = "Defoliation (spr)",
  "Gelure_cl0.4spring"         = "Frost damage (spr)",
  "Gelivures_cl.0.4spring"     = "Frost cracks (spr)",
  "Carence_cl_0.3_spring"      = "Deficiency (spr)",
  "Infection_cl_0.3_spring"    = "Infection (spr)",
  "Etat_general_cl_0.4_spring" = "Gen. cond. (spr)",
  "Bris_cl0.5summer2"          = "Breakage (sum)",
  "Defoliation_cl0.4summer3"   = "Defoliation (sum)",
  "Carence_cl0.4summer"        = "Deficiency (sum)",
  "Infection_cl0.4summer"      = "Infection (sum)",
  "etat_general_cl0.4summer"   = "Gen. cond. (sum)"
)

# unmatched names (e.g. the tree-rating columns) pass through unchanged
shorten <- function(x) {
  x <- as.character(x)
  ifelse(x %in% names(short_names), short_names[x], x)
}

term_short <- c("Salt:Biochar" = "Salt\u00D7Biochar")
shorten_term <- function(x) ifelse(x %in% names(term_short), term_short[x], x)

library(ggplot2)
library(dplyr)

alpha     <- 0.1
trt_terms <- c("Salt", "Biochar", "Salt:Biochar")
fig_dir   <- file.path(project_root, "outputs", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# species code -> latin name
sp_lookup <- unique(as.data.frame(treepot[c("Sp", "Species.latin")]))
sp_label  <- setNames(as.character(sp_lookup$Species.latin), as.character(sp_lookup$Sp))

# ---- 1. significant terms for every model -----------------------------------
get_sig <- function(model_list, type, domain) {
  rows <- list()
  for (sp in names(model_list)) for (v in names(model_list[[sp]])) {
    m <- model_list[[sp]][[v]]
    a <- tryCatch(
      if (type == "lmer") as.data.frame(anova(m, type = 3)) else as.data.frame(anova(m)),
      error = function(e) NULL)
    if (is.null(a)) next
    terms <- if ("Term" %in% names(a)) trimws(as.character(a$Term)) else rownames(a)
    p     <- a[["Pr(>F)"]]
    keep  <- terms %in% trt_terms & !is.na(p) & p < alpha
    if (any(keep))
      rows[[length(rows) + 1]] <- data.frame(Domain = domain, Species = sp, Trait = v,
                                             Term = terms[keep], p = p[keep])
  }
  if (length(rows)) do.call(rbind, rows) else NULL
}

model_sets <- list(
  "Growth"      = list(models = growth_models, type = "lmer"),
  "Chemistry"   = list(models = chem_models,   type = "lmer"),
  "Leaf health" = list(models = lf_models,     type = "lmer"),
  "Tree health" = list(models = tc_models,     type = "art")
)

sig_all <- do.call(rbind, lapply(names(model_sets), function(d)
  get_sig(model_sets[[d]]$models, model_sets[[d]]$type, d)))

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3))
sig_all$lab <- paste0(sub(":", "\u00D7", sig_all$Term), " p=", fmt_p(sig_all$p))

sig_combo <- sig_all |>
  group_by(Domain, Species, Trait) |>
  summarise(terms = paste(lab, collapse = "; "), .groups = "drop")

# ---- 2. data behind each significant combination ----------------------------
get_data <- function(domain, sp, v) {
  m <- model_sets[[domain]]$models[[sp]][[v]]
  if (model_sets[[domain]]$type == "lmer") {
    mf <- model.frame(m)                       # exactly the rows the model used
    data.frame(Salt = mf$Salt, Biochar = mf$Biochar, value = mf[[1]])
  } else {
    d <- species_list[[sp]]; d <- d[!is.na(d[[v]]), ]
    data.frame(Salt = d$Salt, Biochar = d$Biochar, value = d[[v]])
  }
}

plot_df <- do.call(rbind, lapply(seq_len(nrow(sig_combo)), function(i) {
  r <- sig_combo[i, ]
  d <- get_data(r$Domain, r$Species, r$Trait)
  d$Domain   <- r$Domain
  d$Species  <- r$Species
  d$SpLatin  <- ifelse(is.na(sp_label[r$Species]), r$Species, sp_label[r$Species])
  d$Trait    <- r$Trait
  d$terms    <- r$terms
  d
}))
plot_df$Salt    <- factor(plot_df$Salt,    levels = c("N", "L", "H"))
plot_df$Biochar <- factor(plot_df$Biochar, levels = c("S", "B"))
plot_df$Trait <- shorten(plot_df$Trait)

# ---- 3. one plotting function for both layouts -------------------------------
dodge <- position_dodge(width = 0.8)

base_plot <- function(d, facet_var) {
  ggplot(d, aes(Salt, value, fill = Biochar, colour = Biochar)) +
    geom_boxplot(position = dodge, alpha = 0.35, outlier.shape = NA, width = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               alpha = 0.7, size = 1.6) +
    stat_summary(fun = mean, geom = "point", position = dodge,
                 shape = 23, size = 3, fill = "white", colour = "black") +
    facet_wrap(as.formula(paste("~", facet_var)), scales = "free_y") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    labs(x = "Salt treatment", y = NULL) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 8, face = "bold"))
}

page_size <- function(n, ncol = 3) c(w = 4 * min(n, ncol) + 1, h = 3.4 * ceiling(n / ncol) + 1.2)

# ---- A. one figure per VARIABLE, a panel per species -------------------------
plot_df$panel_sp <- paste0(plot_df$SpLatin, "\n", plot_df$terms)

plots_by_trait <- lapply(split(plot_df, plot_df$Trait), function(d) {
  base_plot(d, "panel_sp") +
    labs(title = unique(d$Trait), subtitle = paste("Domain:", unique(d$Domain)),
         y = unique(d$Trait))
})

# ---- B. one figure per SPECIES, a panel per significant variable --------------
plot_df$panel_tr <- paste0(plot_df$Domain, ": ", plot_df$Trait, "\n", plot_df$terms)

plots_by_species <- lapply(split(plot_df, plot_df$Species), function(d) {
  base_plot(d, "panel_tr") +
    labs(title = paste0(unique(d$Species), "  ", unique(d$SpLatin)),
         subtitle = paste0("Significant treatment effects (p < ", alpha, ")"))
})

# ---- 4. save to multi-page PDFs ----------------------------------------------
save_pdf <- function(plot_list, d_list, file) {
  sz <- page_size(max(sapply(d_list, function(d) length(unique(d$panel_sp)))))
  pdf(file, width = sz["w"], height = sz["h"], onefile = TRUE)
  for (p in plot_list) print(p)
  dev.off()
}

# size each page to its own number of panels
pdf_pages <- function(plot_list, d_list, panel_col, file) {
  pdf(file, width = 12, height = 8, onefile = TRUE)   # fixed page, panels scale to fit
  for (nm in names(plot_list)) print(plot_list[[nm]])
  dev.off()
}
pdf_pages(plots_by_trait,   split(plot_df, plot_df$Trait),   "panel_sp",
          file.path(fig_dir, "Significant_by_variable.pdf"))
pdf_pages(plots_by_species, split(plot_df, plot_df$Species), "panel_tr",
          file.path(fig_dir, "Significant_by_species.pdf"))

# or save individual PNGs sized to their panel count
for (nm in names(plots_by_species)) {
  n  <- length(unique(plot_df$panel_tr[plot_df$Species == nm]))
  sz <- page_size(n)
  ggsave(file.path(fig_dir, paste0("species_", nm, ".png")), plots_by_species[[nm]],
         width = sz["w"], height = sz["h"], dpi = 300)
}
for (nm in names(plots_by_trait)) {
  n  <- length(unique(plot_df$panel_sp[plot_df$Trait == nm]))
  sz <- page_size(n)
  ggsave(file.path(fig_dir, paste0("var_", make.names(nm), ".png")), plots_by_trait[[nm]],
         width = sz["w"], height = sz["h"], dpi = 300)
}

# view one on screen
plots_by_species[[3]]
plots_by_trait[["LICOR"]]

# 7. Significant results summary: one sheet per species -----------------------
library(openxlsx)
library(emmeans)
library(ARTool)

alpha     <- 0.1
trt_terms <- c("Salt", "Biochar", "Salt:Biochar")

# ---- helpers ----------------------------------------------------------------
# significant treatment rows of one model's ANOVA table
sig_anova <- function(m, type) {
  a <- if (type == "lmer") as.data.frame(anova(m, type = 3)) else as.data.frame(anova(m))
  terms <- if ("Term" %in% names(a)) trimws(as.character(a$Term)) else rownames(a)
  out <- data.frame(
    Term = terms,
    F    = a[["F value"]],
    df1  = if ("NumDF" %in% names(a)) a$NumDF else a$Df,
    df2  = if ("DenDF" %in% names(a)) a$DenDF else a$Df.res,
    p    = a[["Pr(>F)"]]
  )
  out[out$Term %in% trt_terms & !is.na(out$p) & out$p < alpha, , drop = FALSE]
}

# contrast object -> tidy data frame (same columns whatever the "by" variable)
tidy_con <- function(con, label) {
  df <- as.data.frame(con)
  std <- c("contrast", "estimate", "SE", "df", "t.ratio", "z.ratio",
           "p.value", "lower.CL", "upper.CL", "null")
  by_cols <- setdiff(names(df), std)
  cont <- as.character(df$contrast)
  
  # direction from the sign of the estimate: positive = first level is higher
  dir <- mapply(function(ct, est) {
    lv <- trimws(gsub("[()]", "", strsplit(ct, " - ", fixed = TRUE)[[1]]))
    if (length(lv) != 2 || is.na(est)) return(NA_character_)
    if (est > 0) paste(lv[1], ">", lv[2]) else paste(lv[2], ">", lv[1])
  }, cont, df$estimate, USE.NAMES = FALSE)
  
  if (length(by_cols) > 0) {
    by_txt <- apply(df[by_cols], 1, function(r)
      paste(paste(by_cols, trimws(r), sep = "="), collapse = ", "))
    cont <- paste0(cont, "  [", by_txt, "]")
  }
  data.frame(Comparison = label, Contrast = cont,
             Estimate = df$estimate, SE = df$SE,
             Direction = dir, p = df$p.value)
}
# significant Tukey contrasts for the significant terms of one model
sig_posthoc <- function(m, sig_terms, type) {
  out <- list()
  for (term in sig_terms) {
    if (term %in% c("Salt", "Biochar")) {
      con <- if (type == "lmer") {
        pairs(emmeans(m, as.formula(paste("~", term))), adjust = "tukey")
      } else {
        art.con(m, term, adjust = "tukey")
      }
      out[[term]] <- tidy_con(con, term)
    } else {   # Salt:Biochar -> both conditional versions
      specs <- list("Salt | Biochar" = ~ Salt | Biochar,
                    "Biochar | Salt" = ~ Biochar | Salt)
      for (lab in names(specs)) {
        emm <- if (type == "lmer") emmeans(m, specs[[lab]])
        else emmeans(artlm(m, "Salt:Biochar"), specs[[lab]])
        out[[lab]] <- tidy_con(pairs(emm, adjust = "tukey"), lab)
      }
    }
  }
  res <- do.call(rbind, out)
  res[!is.na(res$p) & res$p < alpha, , drop = FALSE]
}

# one model list (species -> trait -> model) -> species -> list(anova, posthoc)
collect_sig <- function(model_list, type) {
  lapply(model_list, function(sp_models) {
    an <- list(); ph <- list()
    for (v in names(sp_models)) {
      m <- sp_models[[v]]
      a <- tryCatch(sig_anova(m, type), error = function(e) NULL)
      if (is.null(a) || nrow(a) == 0) next
      an[[v]] <- cbind(Trait = shorten(v), a)
      p <- tryCatch(sig_posthoc(m, a$Term, type), error = function(e) NULL)
      if (!is.null(p) && nrow(p) > 0) ph[[v]] <- cbind(Trait = shorten(v), p)
    }
    bind <- function(l) if (length(l)) { x <- do.call(rbind, l); rownames(x) <- NULL; x } else NULL
    list(anova = bind(an), posthoc = bind(ph))
  })
}

# ---- collect ----------------------------------------------------------------
results <- list(
  "Growth"      = collect_sig(growth_models, "lmer"),
  "Chemistry"   = collect_sig(chem_models,   "lmer"),
  "Leaf health" = collect_sig(lf_models,     "lmer"),
  "Tree health" = collect_sig(tc_models,     "art")
)

# species code -> latin name, for the sheet titles
sp_lookup <- unique(as.data.frame(treepot[c("Sp", "Species.latin")]))
sp_label  <- setNames(as.character(sp_lookup$Species.latin), as.character(sp_lookup$Sp))

# ---- write workbook -----------------------------------------------------------
wb <- createWorkbook()
header_style  <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "bottom")
sig05_style   <- createStyle(textDecoration = "bold", fgFill = "#FFF2CC")   # p < 0.05
p_style       <- createStyle(numFmt = "0.000000000")
title_style   <- createStyle(textDecoration = "bold", fontSize = 14)
section_style <- createStyle(textDecoration = "bold", fontSize = 12, fgFill = "#E2EFDA")
label_style   <- createStyle(textDecoration = "bold")

write_block <- function(wb, sheet, df, row) {
  df <- as.data.frame(df)
  for (cn in intersect(c("F", "df1", "df2", "Estimate", "SE"), names(df)))
    df[[cn]] <- round(df[[cn]], 3)
  writeData(wb, sheet, df, startRow = row, headerStyle = header_style)
  rows <- (row + 1):(row + nrow(df))
  addStyle(wb, sheet, p_style, rows = rows, cols = which(names(df) == "p"),
           gridExpand = TRUE, stack = TRUE)
  hi <- which(!is.na(df$p) & df$p < 0.05)
  if (length(hi) > 0)
    addStyle(wb, sheet, sig05_style, rows = row + hi, cols = seq_len(ncol(df)),
             gridExpand = TRUE, stack = TRUE)
  row + nrow(df) + 2
}

all_sp <- sort(unique(unlist(lapply(results, names))))

for (sp in all_sp) {
  sheet <- substr(sp, 1, 31)
  addWorksheet(wb, sheet)
  lat <- unname(sp_label[sp]); if (is.na(lat)) lat <- ""
  writeData(wb, sheet, paste0(sp, "  ", lat, "  (terms with p < ", alpha, ")"), startRow = 1)
  addStyle(wb, sheet, title_style, rows = 1, cols = 1)
  writeData(wb, sheet, "Bold/yellow rows: p < 0.05", startRow = 2)
  row <- 4
  
  for (dom in names(results)) {
    writeData(wb, sheet, dom, startRow = row)
    addStyle(wb, sheet, section_style, rows = row, cols = 1:8, gridExpand = TRUE)
    row <- row + 1
    
    r <- results[[dom]][[sp]]
    if (is.null(r) || is.null(r$anova)) {
      writeData(wb, sheet, paste0("No treatment term significant at p < ", alpha), startRow = row)
      row <- row + 3
      next
    }
    
    writeData(wb, sheet, "ANOVA", startRow = row)
    addStyle(wb, sheet, label_style, rows = row, cols = 1)
    row <- write_block(wb, sheet, r$anova, row + 1)
    
    writeData(wb, sheet, paste0("Tukey-adjusted contrasts (p < ", alpha, ")"), startRow = row)
    addStyle(wb, sheet, label_style, rows = row, cols = 1)
    if (is.null(r$posthoc)) {
      writeData(wb, sheet, "No significant pairwise contrasts", startRow = row + 1)
      row <- row + 4
    } else {
      row <- write_block(wb, sheet, r$posthoc, row + 1)
    }
    row <- row + 1
  }
  setColWidths(wb, sheet, cols = 1:8, widths = c(34, 16, 34, 12, 12, 12, 14, 14))
}

saveWorkbook(wb, file.path(table_dir, "Treepot_r_species.xlsx"), overwrite = TRUE)
