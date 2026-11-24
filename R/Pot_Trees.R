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

###plotting grid####
library(ggplot2)
library(patchwork)
library(dplyr)

V_growth_matrix<-c(
 "Shoot mass (g)",
 "Root mass (g)",
 "Root:shoot",
 "Diam. change",
 "Height change",
 "Crown change",
 "Vigor",
)
V_growth_change <- c("Total.above.ground.biomass.g",   # final biomass, no "change" version exists
                  "Below.ground.biomass.dry.g",
                  "Root.shoot.Ratio",
                  "Diam.Diff.25.23",                # diameter change (replaces Total.diamm.at.15mm)
                  "Height.Difference.25.23",        # height change (replaces Height.cm)
                  "Crown.diam.avg.diff.25.23")      # crown change (replaces Crown.diameter.average)
# ---- 1. matrix definition ----------------------------------------------------
col_def <- bind_rows(
  data.frame(Domain = "Growth",      Trait = V_growth_change),
  data.frame(Domain = "Chemistry",   Trait = V_chem),
  data.frame(Domain = "Leaf health", Trait = c("LICOR", setdiff(V_lc, c(lic_vars, "LICOR.Avg")))),
  data.frame(Domain = "Tree health", Trait = V_tc)
) |> mutate(Trait = shorten(Trait))

sp_codes  <- names(species_list)
sp_latin  <- ifelse(is.na(sp_label[sp_codes]), "", sp_label[sp_codes])
sp_rowlab <- setNames(paste0(sp_codes, "\n", sp_latin), sp_codes)

dom_fill <- c("Growth" = "#C6E0B4", "Chemistry" = "#BDD7EE",
              "Leaf health" = "#FFE699", "Tree health" = "#F4B183")

# ---- 2. building blocks ------------------------------------------------------
dodge <- position_dodge(width = 0.8)

cell_plot <- function(d) {
  ggplot(d, aes(Salt, value, fill = Biochar, colour = Biochar)) +
    geom_boxplot(position = dodge, alpha = 0.35, outlier.shape = NA, width = 0.7,
                 linewidth = 0.3) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               alpha = 0.7, size = 0.7) +
    stat_summary(fun = mean, geom = "point", position = dodge,
                 shape = 23, size = 1.6, fill = "white", colour = "black") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    labs(x = NULL, y = NULL, subtitle = unique(d$terms)) +
    theme_bw(base_size = 7) +
    theme(plot.subtitle = element_text(size = 5.5),
          legend.position = "bottom",
          panel.grid.minor = element_blank())
}

wrap_txt <- function(x, w) paste(strwrap(x, w), collapse = "\n")

text_cell <- function(label, fill, size = 3, face = "bold", wrap = 18) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = wrap_txt(label, wrap),
             size = size, fontface = face, lineheight = 0.9) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = fill, colour = "white", linewidth = 0.5))
}

# ---- 3. one matrix (rows = species, columns = traits) -------------------------
build_grid <- function(sp_vec, cols, title = NULL, txt = 3) {
  nr <- length(sp_vec); nc <- nrow(cols)
  cells <- list()
  
  # header row 1: domain (name shown on the first column of each domain only)
  cells <- c(cells, list(plot_spacer()))
  for (j in seq_len(nc)) {
    first <- j == 1 || cols$Domain[j] != cols$Domain[j - 1]
    cells <- c(cells, list(text_cell(if (first) cols$Domain[j] else "",
                                     dom_fill[[cols$Domain[j]]], size = txt + 0.5)))
  }
  # header row 2: trait names
  cells <- c(cells, list(plot_spacer()))
  for (j in seq_len(nc))
    cells <- c(cells, list(text_cell(cols$Trait[j], dom_fill[[cols$Domain[j]]],
                                     size = txt, wrap = 16)))
  
  # species rows
  for (sp in sp_vec) {
    cells <- c(cells, list(text_cell(sp_rowlab[[sp]], "grey90", size = txt, wrap = 16)))
    for (j in seq_len(nc)) {
      d <- plot_df[plot_df$Species == sp & plot_df$Trait == cols$Trait[j] &
                     plot_df$Domain == cols$Domain[j], ]
      cells <- c(cells, list(if (nrow(d) == 0) plot_spacer() else cell_plot(d)))
    }
  }
  
  g <- wrap_plots(cells, ncol = nc + 1,
                  widths  = c(0.8, rep(1, nc)),
                  heights = c(0.25, 0.45, rep(1, nr)),
                  guides  = "collect")
  if (!is.null(title))
    g <- g + plot_annotation(title = title,
                             theme = theme(plot.title = element_text(size = 12, face = "bold")))
  g & theme(legend.position = "bottom")
}

n_sig <- function(sp_vec, cols)
  sum(plot_df$Species %in% sp_vec & paste(plot_df$Domain, plot_df$Trait) %in%
        paste(cols$Domain, cols$Trait))

# ---- 4a. BIG version: one large page per domain -------------------------------
cell_in <- 2.6                                   # inches per plot
pdf(file.path(fig_dir, "Results_matrix_big.pdf"), onefile = TRUE,
    width = 1.6 + cell_in * 12,  height = 2.2 + cell_in * length(sp_codes))
# note: Width is fixed to the widest domain (12 columns) so all pages are the same size.
for (dom in unique(col_def$Domain)) {
  cols <- col_def[col_def$Domain == dom, ]
  print(build_grid(sp_codes, cols,
                   title = paste0(dom, " - significant treatment effects (p < ", alpha, ")"),
                   txt = 4.5))
}
dev.off()

# ---- 4b. STANDARD version: 20 cells per US Letter page (5 species x 4 variables)
sp_blocks <- split(sp_codes, ceiling(seq_along(sp_codes) / 5))

pdf(file.path(fig_dir, "Results_matrix_letter.pdf"), onefile = TRUE,
    width = 8.5, height = 11)
for (dom in unique(col_def$Domain)) {
  cols_dom   <- col_def[col_def$Domain == dom, ]
  col_blocks <- split(seq_len(nrow(cols_dom)), ceiling(seq_len(nrow(cols_dom)) / 4))
  for (sb in seq_along(sp_blocks)) for (cb in col_blocks) {
    cols <- cols_dom[cb, ]
    if (n_sig(sp_blocks[[sb]], cols) == 0) next          # skip pages that would be empty
    print(build_grid(sp_blocks[[sb]], cols,
                     title = paste0(dom, " (species ", paste(range(match(sp_blocks[[sb]], sp_codes)), collapse = "-"),
                                    " of ", length(sp_codes), ", p < ", alpha, ")"),
                     txt = 2.6))
  }
}
dev.off()


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

###PCA####
library(vegan)      # adonis2, how()
library(ggplot2)
library(dplyr)

lic_vars <- c("Licor.July.27", "Licor.Aug.07")
pca_sets <- list(
  "All"         = c(V_growth, V_chem, "LICOR.Avg",
                    setdiff(V_lc, c(lic_vars, "LICOR.Avg")), V_tc),
  "Growth"      = V_growth,
  "Chemistry"   = V_chem,
  "Leaf health" = c("LICOR.Avg", setdiff(V_lc, c(lic_vars, "LICOR.Avg"))),
  "Tree health" = V_tc
)

run_pca <- function(sp, vars, set_name, min_complete = 0.8) {
  d <- species_list[[sp]]
  X <- as.data.frame(d[, vars])
  X <- X[, colMeans(!is.na(X)) >= min_complete, drop = FALSE]
  X <- X[, sapply(X, function(x) isTRUE(sd(x, na.rm = TRUE) > 0)), drop = FALSE]
  ok   <- complete.cases(X)
  X    <- X[ok, , drop = FALSE]
  meta <- d[ok, c("Salt", "Biochar", "Block")] |> droplevels()
  if (ncol(X) < 3 || nrow(X) < 10) return(NULL)
  
  pca <- prcomp(X, scale. = TRUE)
  ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  
  # PERMANOVA on standardised data, permuting within Block
  pm <- tryCatch(
    adonis2(scale(X) ~ Salt * Biochar, data = meta, method = "euclidean",
            permutations = how(nperm = 999, blocks = meta$Block), by = "terms"),
    error = function(e) NULL)
  pp <- function(term) if (is.null(pm) || !term %in% rownames(pm)) NA_real_ else pm[term, "Pr(>F)"]
  
  sc <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta)
  ld <- data.frame(var = shorten(rownames(pca$rotation)),
                   PC1 = pca$rotation[, 1], PC2 = pca$rotation[, 2])
  ld <- ld[order(-(ld$PC1^2 + ld$PC2^2)), ][seq_len(min(8, nrow(ld))), ]   # 8 strongest
  mult <- 0.8 * min(diff(range(sc$PC1)) / diff(range(ld$PC1)),
                    diff(range(sc$PC2)) / diff(range(ld$PC2)))
  ld[c("PC1", "PC2")] <- ld[c("PC1", "PC2")] * mult
  
  p <- ggplot(sc, aes(PC1, PC2, colour = Salt, shape = Biochar)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    stat_ellipse(aes(group = Salt), level = 0.68, linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 inherit.aes = FALSE, colour = "grey30",
                 arrow = arrow(length = unit(0.12, "cm"))) +
    geom_text(data = ld, aes(PC1, PC2, label = var), inherit.aes = FALSE,
              size = 2.8, colour = "grey20", vjust = -0.4) +
    scale_colour_brewer(palette = "Dark2") +
    labs(title = paste0(sp, " ", ifelse(is.na(sp_label[sp]), "", sp_label[sp]), " - ", set_name),
         subtitle = sprintf("n = %d trees, %d variables | PERMANOVA: Salt p=%.3f, Biochar p=%.3f, Salt\u00D7Biochar p=%.3f",
                            nrow(X), ncol(X), pp("Salt"), pp("Biochar"), pp("Salt:Biochar")),
         x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_bw() + theme(plot.subtitle = element_text(size = 8), legend.position = "bottom")
  
  list(plot = p,
       info = data.frame(Species = sp, Set = set_name, n = nrow(X), p = ncol(X),
                         PC1_pct = round(ve[1], 1), PC2_pct = round(ve[2], 1),
                         p_Salt = pp("Salt"), p_Biochar = pp("Biochar"),
                         p_SaltxBiochar = pp("Salt:Biochar")))
}

# run everything
pca_res <- list()
for (s in names(pca_sets)) for (sp in names(species_list))
  pca_res[[paste(s, sp)]] <- tryCatch(run_pca(sp, pca_sets[[s]], s),
                                      error = function(e) NULL)
pca_res <- Filter(Negate(is.null), pca_res)

pca_summary <- do.call(rbind, lapply(pca_res, `[[`, "info"))
rownames(pca_summary) <- NULL
pca_summary

# one page per species x set
pdf(file.path(fig_dir, "PCA_by_species.pdf"), width = 9, height = 7.5, onefile = TRUE)
for (r in pca_res) print(r$plot)
dev.off()

pca_res[["All MB"]]$plot   # view one
pca_res


library(vegan)
library(ggplot2)
library(dplyr)

lic_vars <- c("Licor.July.27", "Licor.Aug.07")
pca_vars <- c(V_growth_change, V_chem, "LICOR.Avg",
              setdiff(V_lc, c(lic_vars, "LICOR.Avg")), V_tc)

dat <- as.data.frame(treepot)
dat$Species <- droplevels(factor(dat$Sp))

# ---- 1. choose variables so no species is lost ------------------------------
# keep a variable only if every species has >= 50% of it measured
# (otherwise one species with no LICOR or N data would remove all its trees)
comp_by_sp <- sapply(pca_vars, function(v) tapply(!is.na(dat[[v]]), dat$Species, mean))
keep <- colnames(comp_by_sp)[apply(comp_by_sp, 2, min) >= 0.5]
cat("Dropped (too incomplete in at least one species):\n  ",
    paste(shorten(setdiff(pca_vars, keep)), collapse = ", "), "\n")

X0 <- dat[, keep, drop = FALSE]
X0 <- X0[, sapply(X0, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]
ok <- complete.cases(X0)
cat("Trees used:", sum(ok), "of", nrow(dat), "| variables:", ncol(X0), "\n")
print(table(dat$Species[ok]))                      # trees per species that survive

X    <- X0[ok, , drop = FALSE]
meta <- dat[ok, c("Species", "Salt", "Biochar", "Block")] |> droplevels()
meta$SpLatin <- ifelse(is.na(sp_label[as.character(meta$Species)]),
                       as.character(meta$Species), sp_label[as.character(meta$Species)])

# ---- 2. two scalings -----------------------------------------------------------
X_all <- scale(X)                                              # across species
X_within <- as.data.frame(lapply(X, function(x)                # within species
  ave(x, meta$Species, FUN = function(z) {
    s <- sd(z); if (is.na(s) || s == 0) rep(0, length(z)) else (z - mean(z)) / s })))
X_within <- as.matrix(X_within)

# ---- 3. plotting helper ----------------------------------------------------------
make_pca <- function(Xm, colour_var, shape_var = NULL, title, subtitle, ellipse_by) {
  pca <- prcomp(Xm, scale. = FALSE)               # already standardised
  ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  sc  <- cbind(data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2]), meta)
  
  ld <- data.frame(var = shorten(rownames(pca$rotation)),
                   PC1 = pca$rotation[, 1], PC2 = pca$rotation[, 2])
  ld <- ld[order(-(ld$PC1^2 + ld$PC2^2)), ][1:min(10, nrow(ld)), ]
  mult <- 0.8 * min(diff(range(sc$PC1)) / diff(range(ld$PC1)),
                    diff(range(sc$PC2)) / diff(range(ld$PC2)))
  ld[c("PC1", "PC2")] <- ld[c("PC1", "PC2")] * mult
  
  aes_pts <- if (is.null(shape_var)) aes(colour = .data[[colour_var]])
  else aes(colour = .data[[colour_var]], shape = .data[[shape_var]])
  
  p <- ggplot(sc, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    stat_ellipse(aes(group = .data[[ellipse_by]], colour = .data[[colour_var]]),
                 level = 0.68, linewidth = 0.5, show.legend = FALSE) +
    geom_point(aes_pts, size = 2.2, alpha = 0.8) +
    geom_segment(data = ld, aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 inherit.aes = FALSE, colour = "grey30",
                 arrow = arrow(length = unit(0.12, "cm"))) +
    geom_text(data = ld, aes(PC1, PC2, label = var), inherit.aes = FALSE,
              size = 2.8, colour = "grey15", vjust = -0.4) +
    labs(title = title, subtitle = subtitle,
         x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
    theme_bw() + theme(plot.subtitle = element_text(size = 8), legend.position = "right")
  list(plot = p, pca = pca, scores = sc, ve = ve)
}

# ---- 4A. species grouping --------------------------------------------------------
pm_sp <- adonis2(X_all ~ Species, data = meta, method = "euclidean", permutations = 999)
A <- make_pca(X_all, "SpLatin", title = "PCA of all species, all traits (grouped by species)",
              subtitle = sprintf("n = %d trees, %d variables | PERMANOVA Species R\u00B2 = %.2f, p = %.3f",
                                 nrow(X), ncol(X), pm_sp$R2[1], pm_sp$`Pr(>F)`[1]),
              ellipse_by = "SpLatin") 
A$plot <- A$plot + scale_colour_manual(values = scales::hue_pal()(nlevels(meta$Species)),
                                       name = "Species")

# ---- 4B. treatment grouping (species effect removed) -------------------------------
pm_tr <- adonis2(X_within ~ Salt * Biochar, data = meta, method = "euclidean",
                 permutations = how(nperm = 999, blocks = meta$Species), by = "terms")
B <- make_pca(X_within, "Salt", "Biochar",
              title = "PCA of all species, all traits (grouped by treatment; scaled within species)",
              subtitle = sprintf("n = %d trees, %d variables | PERMANOVA (permuted within species): Salt p=%.3f, Biochar p=%.3f, Salt\u00D7Biochar p=%.3f",
                                 nrow(X), ncol(X), pm_tr["Salt", "Pr(>F)"],
                                 pm_tr["Biochar", "Pr(>F)"], pm_tr["Salt:Biochar", "Pr(>F)"]),
              ellipse_by = "Salt")
B$plot <- B$plot + scale_colour_brewer(palette = "Dark2")

# treatment PCA, one panel per species (same axes, so panels are comparable)
B_facet <- ggplot(B$scores, aes(PC1, PC2, colour = Salt, shape = Biochar)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_wrap(~ SpLatin, ncol = 4) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Treatment PCA split by species",
       x = sprintf("PC1 (%.1f%%)", B$ve[1]), y = sprintf("PC2 (%.1f%%)", B$ve[2])) +
  theme_bw() + theme(legend.position = "bottom", strip.text = element_text(size = 7))

# ---- 5. save -----------------------------------------------------------------------
pdf(file.path(fig_dir, "PCA_all_species_final.pdf"), width = 10, height = 7.5, onefile = TRUE)
print(A$plot); print(B$plot); print(B_facet)
dev.off()
ggsave(file.path(fig_dir, "PCA_species.png"),   A$plot, width = 10, height = 7.5, dpi = 300)
ggsave(file.path(fig_dir, "PCA_treatment.png"), B$plot, width = 10, height = 7.5, dpi = 300)

# loadings table, to cite which traits drive each axis
load_tbl <- function(R) {
  l <- as.data.frame(R$pca$rotation[, 1:3]); l$Trait <- shorten(rownames(l)); rownames(l) <- NULL
  l[, c("Trait", "PC1", "PC2", "PC3")]
}
load_species   <- load_tbl(A)
load_treatment <- load_tbl(B)


###excel 2: matrix with arrows####
library(openxlsx)
library(emmeans)
library(ARTool)
library(dplyr)
library(ggplot2)

alpha <- 0.1
require_sig_contrast <- TRUE   # TRUE: arrow only if the N vs H contrast itself is p < alpha
# FALSE: arrow whenever the Salt (or Salt x Biochar) F-test is p < alpha

sets <- list(
  "Growth"      = list(models = growth_models, type = "lmer"),
  "Chemistry"   = list(models = chem_models,   type = "lmer"),
  "Leaf health" = list(models = lf_models,     type = "lmer"),
  "Tree health" = list(models = tc_models,     type = "art")
)

# ---- columns (original trait names, so we can find the models) ----------------
col_x <- bind_rows(
  data.frame(Domain = "Growth",      Trait = V_growth_change),
  data.frame(Domain = "Chemistry",   Trait = V_chem),
  data.frame(Domain = "Leaf health", Trait = c("LICOR", setdiff(V_lc, c(lic_vars, "LICOR.Avg")))),
  data.frame(Domain = "Tree health", Trait = V_tc)
) |> mutate(Short = shorten(Trait), col = row_number())

sp_codes <- names(species_list)
sp_label_full <- setNames(paste0(sp_codes, " ", ifelse(is.na(sp_label[sp_codes]), "", sp_label[sp_codes])), sp_codes)

# ---- helpers ------------------------------------------------------------------
anova_p <- function(m, type) {
  a <- as.data.frame(if (type == "lmer") anova(m, type = 3) else anova(m))
  terms <- if ("Term" %in% names(a)) trimws(as.character(a$Term)) else rownames(a)
  setNames(a[["Pr(>F)"]], terms)
}

get_emm <- function(m, specs, term, type)
  if (type == "lmer") emmeans(m, specs) else emmeans(artlm(m, term), specs)

# N vs H contrast only, with direction (H_higher = TRUE if H > N)
nh_contrast <- function(emm) {
  df <- as.data.frame(pairs(emm, adjust = "tukey"))
  lv <- strsplit(gsub("[()]", "", as.character(df$contrast)), " - ", fixed = TRUE)
  keep <- sapply(lv, function(x) setequal(trimws(x), c("N", "H")))
  df <- df[keep, , drop = FALSE]; lv <- lv[keep]
  if (nrow(df) == 0) return(NULL)
  first <- trimws(sapply(lv, `[`, 1))
  df$H_higher <- ifelse(first == "N", df$estimate < 0, df$estimate > 0)
  df
}

cell_arrows <- function(m, type) {
  p   <- anova_p(m, type)
  sig <- function(t) isTRUE(unname(p[t]) < alpha)
  if (sig("Salt:Biochar")) {                       # interaction: one arrow per Biochar level
    term <- "Salt:Biochar"
    ct <- nh_contrast(get_emm(m, ~ Salt | Biochar, term, type))
    bio <- if (is.null(ct)) NULL else as.character(ct$Biochar)
  } else if (sig("Salt")) {                        # main effect: single arrow
    term <- "Salt"
    ct <- nh_contrast(get_emm(m, ~ Salt, term, type))
    bio <- if (is.null(ct)) NULL else ""
  } else return(NULL)
  if (is.null(ct)) return(NULL)
  out <- data.frame(dir = ifelse(ct$H_higher, "up", "down"), Bio = bio,
                    contrast_p = ct$p.value, term = term, anova_p = unname(p[term]))
  if (require_sig_contrast) out <- out[!is.na(out$contrast_p) & out$contrast_p < alpha, ]
  if (nrow(out) == 0) return(NULL)
  out[order(match(out$Bio, c("", "B", "S"))), ]
}

# ---- collect arrows -------------------------------------------------------------
arrows <- do.call(rbind, lapply(seq_len(nrow(col_x)), function(j) {
  d <- col_x$Domain[j]; tr <- col_x$Trait[j]
  do.call(rbind, lapply(sp_codes, function(sp) {
    m <- sets[[d]]$models[[sp]][[tr]]
    if (is.null(m)) return(NULL)
    r <- tryCatch(cell_arrows(m, sets[[d]]$type), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    cbind(Species = sp, Domain = d, Trait = tr, Short = col_x$Short[j], col = j, r)
  }))
}))
rownames(arrows) <- NULL

arrows <- arrows |>
  group_by(Species, col) |>
  mutate(k = n(), idx = row_number(),
         state = ifelse(all(dir == "up"), "up", ifelse(all(dir == "down"), "down", "mixed"))) |>
  ungroup()

# ---- EXCEL ----------------------------------------------------------------------
dom_fill <- c("Growth" = "#C6E0B4", "Chemistry" = "#BDD7EE",
              "Leaf health" = "#FFE699", "Tree health" = "#F4B183")
state_fill <- c(up = "#F8CBAD", down = "#BDD7EE", mixed = "#E4DFEC")

mat <- matrix("", nrow = length(sp_codes), ncol = nrow(col_x),
              dimnames = list(sp_codes, NULL))
for (i in seq_len(nrow(arrows))) {
  a <- arrows[i, ]
  sym <- if (a$dir == "up") "\u2191" else "\u2193"
  lab <- if (a$Bio == "") sym else paste(sym, a$Bio)
  cur <- mat[a$Species, a$col]
  mat[a$Species, a$col] <- if (cur == "") lab else paste(cur, lab, sep = "\n")
}

wb <- createWorkbook()
addWorksheet(wb, "Matrix")
ctr <- createStyle(halign = "center", valign = "center", wrapText = TRUE)

# headers
writeData(wb, "Matrix", "Species", startCol = 1, startRow = 2)
for (dom in unique(col_x$Domain)) {
  cols <- which(col_x$Domain == dom) + 1
  writeData(wb, "Matrix", dom, startCol = min(cols), startRow = 1)
  if (length(cols) > 1) mergeCells(wb, "Matrix", cols = cols, rows = 1)
  addStyle(wb, "Matrix", createStyle(textDecoration = "bold", fgFill = dom_fill[[dom]],
                                     halign = "center", border = "TopBottomLeftRight"),
           rows = 1, cols = cols, gridExpand = TRUE)
  addStyle(wb, "Matrix", createStyle(textDecoration = "bold", fgFill = dom_fill[[dom]],
                                     halign = "center", valign = "center", wrapText = TRUE,
                                     border = "TopBottomLeftRight"),
           rows = 2, cols = cols, gridExpand = TRUE)
}
writeData(wb, "Matrix", t(col_x$Short), startCol = 2, startRow = 2, colNames = FALSE)

# body
writeData(wb, "Matrix", sp_label_full[sp_codes], startCol = 1, startRow = 3, colNames = FALSE)
writeData(wb, "Matrix", mat, startCol = 2, startRow = 3, colNames = FALSE)
nr <- length(sp_codes); nc <- nrow(col_x)
addStyle(wb, "Matrix", createStyle(fontSize = 16, halign = "center", valign = "center",
                                   wrapText = TRUE, border = "TopBottomLeftRight",
                                   borderColour = "grey70"),
         rows = 3:(nr + 2), cols = 2:(nc + 1), gridExpand = TRUE)
addStyle(wb, "Matrix", createStyle(textDecoration = "bold", fgFill = "grey90", valign = "center",
                                   border = "TopBottomLeftRight", borderColour = "grey70"),
         rows = 3:(nr + 2), cols = 1, gridExpand = TRUE)
for (i in seq_len(nrow(arrows))) {
  a <- arrows[i, ]
  if (a$idx > 1) next                                   # style each cell once
  addStyle(wb, "Matrix",
           createStyle(fontSize = 16, halign = "center", valign = "center", wrapText = TRUE,
                       fgFill = state_fill[[a$state]],
                       fontColour = switch(a$state, up = "#C55A11", down = "#2E75B6", mixed = "#5B4B8A"),
                       textDecoration = "bold", border = "TopBottomLeftRight", borderColour = "grey70"),
           rows = match(a$Species, sp_codes) + 2, cols = a$col + 1, stack = FALSE)
}
setColWidths(wb, "Matrix", cols = 1, widths = 30)
setColWidths(wb, "Matrix", cols = 2:(nc + 1), widths = 13)
setRowHeights(wb, "Matrix", rows = 2, heights = 32)
setRowHeights(wb, "Matrix", rows = 3:(nr + 2), heights = 38)
freezePane(wb, "Matrix", firstActiveRow = 3, firstActiveCol = 2)

leg <- nr + 4
writeData(wb, "Matrix", c(
  paste0("\u2191 = H salt higher than N salt;  \u2193 = H salt lower than N salt (L not shown)."),
  "Two arrows in a cell = Salt x Biochar interaction: B = with biochar, S = soil without biochar.",
  paste0("Blank = no significant Salt effect (p < ", alpha, ")",
         if (require_sig_contrast) " in the Tukey N-H contrast." else " in the ANOVA."),
  "Tree health: direction is on the aligned-rank scale (up = higher rating score)."
), startCol = 1, startRow = leg, colNames = FALSE)

# details sheet
addWorksheet(wb, "Details")
det <- arrows |>
  transmute(Species, Domain, Trait = Short,
            Biochar = ifelse(Bio == "", "(averaged)", Bio),
            Direction = ifelse(dir == "up", "H > N", "H < N"),
            ANOVA_term = term, ANOVA_p = anova_p, NH_contrast_p = contrast_p)
writeData(wb, "Details", det, headerStyle = createStyle(textDecoration = "bold", fgFill = "#D9E1F2"))
addStyle(wb, "Details", createStyle(numFmt = "0.000000000"), rows = 2:(nrow(det) + 1),
         cols = 7:8, gridExpand = TRUE)
setColWidths(wb, "Details", cols = 1:8, widths = c(10, 14, 26, 12, 12, 14, 16, 16))

saveWorkbook(wb, file.path(table_dir, "Arrow_matrix_H_vs_N.xlsx"), overwrite = TRUE)

# ---- PICTURE ---------------------------------------------------------------------
grid <- expand.grid(col = col_x$col, row = seq_along(sp_codes))
grid$Species <- sp_codes[grid$row]
grid$Domain  <- col_x$Domain[grid$col]
grid <- left_join(grid, distinct(arrows, Species, col, state), by = c("Species", "col"))
grid$state[is.na(grid$state)] <- "none"
grid$y <- length(sp_codes) - grid$row + 1

ar <- arrows |>
  mutate(y = length(sp_codes) - match(Species, sp_codes) + 1,
         xx = col + ifelse(k == 1, 0, ifelse(idx == 1, -0.2, 0.2)),
         y0 = ifelse(dir == "up", y - 0.22, y + 0.22),
         y1 = ifelse(dir == "up", y + 0.22, y - 0.22))

col_x$Domain <- factor(col_x$Domain, levels = unique(col_x$Domain))
grid$Domain  <- factor(grid$Domain,  levels = levels(col_x$Domain))
ar$Domain    <- factor(ar$Domain,    levels = levels(col_x$Domain))

p_arrow <- ggplot() +
  geom_tile(data = grid, aes(col, y, fill = state), colour = "grey70", width = 1, height = 1) +
  geom_segment(data = ar, aes(x = xx, xend = xx, y = y0, yend = y1, colour = dir),
               linewidth = 1.1, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_text(data = filter(ar, Bio != ""), aes(x = xx + 0.08, y = y, label = Bio),
            hjust = 0, size = 3.2, fontface = "bold") +
  facet_grid(. ~ Domain, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = col_x$col, labels = col_x$Short, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(sp_codes),
                     labels = rev(sp_label_full[sp_codes]), expand = c(0, 0)) +
  scale_fill_manual(values = c(none = "white", state_fill)) +
  scale_colour_manual(values = c(up = "#C55A11", down = "#2E75B6")) +
  labs(x = NULL, y = NULL,
       title = "Effect of high salt (H) relative to neutral (N)",
       subtitle = paste0("\u2191 H higher than N, \u2193 H lower than N (p < ", alpha,
                         "); two arrows = Salt \u00D7 Biochar interaction (B = biochar, S = soil); blank = not significant")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.3, "cm"))

ggsave(file.path(fig_dir, "Arrow_matrix_H_vs_N.png"), p_arrow,
       width = 0.55 * nrow(col_x) + 4, height = 0.6 * length(sp_codes) + 3, dpi = 250)
p_arrow
  

# ---- Biochar arrow matrix: B vs S -------------------------------------------
# up = biochar (B) higher than no biochar (S); down = B lower than S

bs_contrast <- function(emm) {
  df <- as.data.frame(pairs(emm, adjust = "tukey"))
  lv <- strsplit(gsub("[()]", "", as.character(df$contrast)), " - ", fixed = TRUE)
  first <- trimws(sapply(lv, `[`, 1))
  df$B_higher <- ifelse(first == "B", df$estimate > 0, df$estimate < 0)
  df
}

cell_arrows_bio <- function(m, type) {
  p   <- anova_p(m, type)
  sig <- function(t) isTRUE(unname(p[t]) < alpha)
  if (sig("Salt:Biochar")) {                     # interaction: one arrow per Salt level
    term <- "Salt:Biochar"
    ct  <- bs_contrast(get_emm(m, ~ Biochar | Salt, term, type))
    lab <- as.character(ct$Salt)
  } else if (sig("Biochar")) {                   # main effect: single arrow
    term <- "Biochar"
    ct  <- bs_contrast(get_emm(m, ~ Biochar, term, type))
    lab <- rep("", nrow(ct))
  } else return(NULL)
  out <- data.frame(dir = ifelse(ct$B_higher, "up", "down"), Lab = lab,
                    contrast_p = ct$p.value, term = term, anova_p = unname(p[term]))
  if (require_sig_contrast) out <- out[!is.na(out$contrast_p) & out$contrast_p < alpha, ]
  if (nrow(out) == 0) return(NULL)
  out[order(match(out$Lab, c("", "N", "L", "H"))), ]
}

# ---- collect ------------------------------------------------------------------
arrows_b <- do.call(rbind, lapply(seq_len(nrow(col_x)), function(j) {
  d <- col_x$Domain[j]; tr <- col_x$Trait[j]
  do.call(rbind, lapply(sp_codes, function(sp) {
    m <- sets[[d]]$models[[sp]][[tr]]
    if (is.null(m)) return(NULL)
    r <- tryCatch(cell_arrows_bio(m, sets[[d]]$type), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    cbind(Species = sp, Domain = d, Trait = tr, Short = col_x$Short[j], col = j, r)
  }))
}))
rownames(arrows_b) <- NULL

arrows_b <- arrows_b |>
  group_by(Species, col) |>
  mutate(k = n(), idx = row_number(),
         state = ifelse(all(dir == "up"), "up", ifelse(all(dir == "down"), "down", "mixed"))) |>
  ungroup()

# ---- EXCEL --------------------------------------------------------------------
mat <- matrix("", nrow = length(sp_codes), ncol = nrow(col_x), dimnames = list(sp_codes, NULL))
for (i in seq_len(nrow(arrows_b))) {
  a <- arrows_b[i, ]
  sym <- if (a$dir == "up") "\u2191" else "\u2193"
  lab <- if (a$Lab == "") sym else paste(sym, a$Lab)
  cur <- mat[a$Species, a$col]
  mat[a$Species, a$col] <- if (cur == "") lab else paste(cur, lab, sep = "\n")
}

wb <- createWorkbook()
addWorksheet(wb, "Matrix")
nr <- length(sp_codes); nc <- nrow(col_x)

writeData(wb, "Matrix", "Species", startCol = 1, startRow = 2)
for (dom in unique(col_x$Domain)) {
  cols <- which(col_x$Domain == dom) + 1
  writeData(wb, "Matrix", dom, startCol = min(cols), startRow = 1)
  if (length(cols) > 1) mergeCells(wb, "Matrix", cols = cols, rows = 1)
  addStyle(wb, "Matrix", createStyle(textDecoration = "bold", fgFill = dom_fill[[dom]],
                                     halign = "center", border = "TopBottomLeftRight"),
           rows = 1, cols = cols, gridExpand = TRUE)
  addStyle(wb, "Matrix", createStyle(textDecoration = "bold", fgFill = dom_fill[[dom]],
                                     halign = "center", valign = "center", wrapText = TRUE,
                                     border = "TopBottomLeftRight"),
           rows = 2, cols = cols, gridExpand = TRUE)
}
writeData(wb, "Matrix", t(col_x$Short), startCol = 2, startRow = 2, colNames = FALSE)
writeData(wb, "Matrix", sp_label_full[sp_codes], startCol = 1, startRow = 3, colNames = FALSE)
writeData(wb, "Matrix", mat, startCol = 2, startRow = 3, colNames = FALSE)

addStyle(wb, "Matrix", createStyle(fontSize = 14, halign = "center", valign = "center",
                                   wrapText = TRUE, border = "TopBottomLeftRight",
                                   borderColour = "grey70"),
         rows = 3:(nr + 2), cols = 2:(nc + 1), gridExpand = TRUE)
addStyle(wb, "Matrix", createStyle(textDecoration = "bold", fgFill = "grey90", valign = "center",
                                   border = "TopBottomLeftRight", borderColour = "grey70"),
         rows = 3:(nr + 2), cols = 1, gridExpand = TRUE)
for (i in seq_len(nrow(arrows_b))) {
  a <- arrows_b[i, ]
  if (a$idx > 1) next
  addStyle(wb, "Matrix",
           createStyle(fontSize = 14, halign = "center", valign = "center", wrapText = TRUE,
                       fgFill = state_fill[[a$state]],
                       fontColour = switch(a$state, up = "#C55A11", down = "#2E75B6", mixed = "#5B4B8A"),
                       textDecoration = "bold", border = "TopBottomLeftRight", borderColour = "grey70"),
           rows = match(a$Species, sp_codes) + 2, cols = a$col + 1, stack = FALSE)
}
setColWidths(wb, "Matrix", cols = 1, widths = 30)
setColWidths(wb, "Matrix", cols = 2:(nc + 1), widths = 13)
setRowHeights(wb, "Matrix", rows = 2, heights = 32)
setRowHeights(wb, "Matrix", rows = 3:(nr + 2), heights = 52)   # room for up to 3 arrows
freezePane(wb, "Matrix", firstActiveRow = 3, firstActiveCol = 2)

writeData(wb, "Matrix", c(
  "\u2191 = biochar (B) higher than no biochar (S);  \u2193 = biochar lower than S.",
  "Several arrows in a cell = Salt x Biochar interaction: the label (N, L, H) is the salt level the B vs S comparison was made in.",
  paste0("Blank = no significant Biochar effect (p < ", alpha, ")",
         if (require_sig_contrast) " in the Tukey B-S contrast." else " in the ANOVA."),
  "Tree health: direction is on the aligned-rank scale (up = higher rating score)."
), startCol = 1, startRow = nr + 4, colNames = FALSE)

addWorksheet(wb, "Details")
det <- arrows_b |>
  transmute(Species, Domain, Trait = Short,
            Salt_level = ifelse(Lab == "", "(averaged)", Lab),
            Direction = ifelse(dir == "up", "B > S", "B < S"),
            ANOVA_term = term, ANOVA_p = anova_p, BS_contrast_p = contrast_p)
writeData(wb, "Details", det, headerStyle = createStyle(textDecoration = "bold", fgFill = "#D9E1F2"))
addStyle(wb, "Details", createStyle(numFmt = "0.000000000"), rows = 2:(nrow(det) + 1),
         cols = 7:8, gridExpand = TRUE)
setColWidths(wb, "Details", cols = 1:8, widths = c(10, 14, 26, 12, 12, 14, 16, 16))

saveWorkbook(wb, file.path(table_dir, "Arrow_matrix_B_vs_S.xlsx"), overwrite = TRUE)

# ---- PICTURE ------------------------------------------------------------------
grid <- expand.grid(col = col_x$col, row = seq_along(sp_codes))
grid$Species <- sp_codes[grid$row]
grid$Domain  <- factor(col_x$Domain[grid$col], levels = unique(col_x$Domain))
grid <- left_join(grid, distinct(arrows_b, Species, col, state), by = c("Species", "col"))
grid$state[is.na(grid$state)] <- "none"
grid$y <- length(sp_codes) - grid$row + 1

off <- function(k, idx) ifelse(k == 1, 0,
                               ifelse(k == 2, c(-0.2, 0.2)[idx], c(-0.3, 0, 0.3)[idx]))
ar <- arrows_b |>
  mutate(y  = length(sp_codes) - match(Species, sp_codes) + 1,
         xx = col + mapply(off, k, idx),
         y0 = ifelse(dir == "up", y - 0.22, y + 0.22),
         y1 = ifelse(dir == "up", y + 0.22, y - 0.22),
         Domain = factor(Domain, levels = unique(col_x$Domain)))

p_bio <- ggplot() +
  geom_tile(data = grid, aes(col, y, fill = state), colour = "grey70", width = 1, height = 1) +
  geom_segment(data = ar, aes(x = xx, xend = xx, y = y0, yend = y1, colour = dir),
               linewidth = 1.1, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_text(data = filter(ar, Lab != ""), aes(x = xx + 0.06, y = y, label = Lab),
            hjust = 0, size = 2.8, fontface = "bold") +
  facet_grid(. ~ Domain, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = col_x$col, labels = col_x$Short, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq_along(sp_codes), labels = rev(sp_label_full[sp_codes]),
                     expand = c(0, 0)) +
  scale_fill_manual(values = c(none = "white", state_fill)) +
  scale_colour_manual(values = c(up = "#C55A11", down = "#2E75B6")) +
  labs(x = NULL, y = NULL,
       title = "Effect of biochar (B) relative to no biochar (S)",
       subtitle = paste0("\u2191 B higher than S, \u2193 B lower than S (p < ", alpha,
                         "); several arrows = Salt \u00D7 Biochar interaction (label = salt level); blank = not significant")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.3, "cm"))

ggsave(file.path(fig_dir, "Arrow_matrix_B_vs_S.png"), p_bio,
       width = 0.55 * nrow(col_x) + 4, height = 0.6 * length(sp_codes) + 3, dpi = 250)
p_bio
