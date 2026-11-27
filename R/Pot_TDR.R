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
    Salt = as.factor(recode_values(Salinity,
                                   "N" ~ "None",
                                   "L" ~ "Low",
                                   "H" ~ "High")),
    Salt    = factor(Salt, levels = c("None", "Low", "High")),
    Biochar = factor(Soil, levels = c("S", "B")),
    Block.old  = substr(as.character(Bloc), nchar(as.character(Bloc)), nchar(as.character(Bloc))),
    Block = substr(as.character(Bloc), 1, 1),
    CODE2 = as.factor(CODE2),
    ID = paste0(CODE2, " ", Block),
    Date = as.factor(Date),
    Conductivity_mS.cm<-as.factor(Conductivity_mS.cm),
    Humidity_percent<-as.factor(Humidity_percent))


####Visualizing dataset####
TDRpot_2025$Date_vis <- as.Date(TDRpot_2025$Date)        
fertilization_date <- as.Date("2025-05-09") #can change ec because of added N P K.

#conductivity:
conduct_graph <- ggplot(TDRpot_2025,
                               aes(x = Date_vis, y = Conductivity_mS.cm, colour = Salt)) +
         geom_point(alpha = 0.2, size = 1) +
         geom_vline(xintercept = fert_date,
                    linetype = "dashed", colour = "black", linewidth = 0.6) +
         annotate("text", x = fert_date, y = Inf, label = "Fertilizer added",
                  angle = 90, vjust = -0.4, hjust = 1.1, size = 3) +
         labs(x = "Date", y = "Conductivity (mS/cm)", colour = "Salt") +
         theme_bw()
       
       conduct_graph

#compared to Belanger
load(file.path("data", "cleaned", "TDRb_allyrs_cleaned.RData"))  
TDRb_2025<-TDRb_allyrs_cleaned |>
filter(Year == 2025)|>
  mutate(Date = as.Date(Date_Vis))

trt_levelsb <- c("Open", "Partial", "Closed")
trt_colsb   <- c(Open = "#00BA38", Partial = "#619CFF", Closed = "#F8766D")

trt_levelsp <- c("None", "Low", "High")
trt_colsp   <- c("None" = "turquoise3", "Low" = "gold", "High" = "deeppink3")

# colour-blind friendly
#trt_colsp <- c("None" = "#1B9E77", "Low" = "#E6AB02", "High" = "#D95F04")

# sequential, light to dark (suits an increasing salt gradient)
#trt_colsp <- c("None" = "#9ECAE1", "Low" = "#4292C6", "High" = "#08306B")


     
library(ggnewscale)

conduct_study_comp <- ggplot() +
  geom_point(data = TDRpot_2025,
             aes(x = Date_vis, y = Conductivity_mS.cm, colour = Salt),
             alpha = 0.2, size = 1,
             position = position_jitter(width = 0.8, height = 0)) +
  scale_colour_manual(values = trt_colsp, breaks = trt_levelsp,
                      name = "Salinity (pots)",
                      guide = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  new_scale_colour() +
  geom_point(data = TDRb_2025,
             aes(x = Date, y = Conductivity, colour = Treatment),
             alpha = 0.6, size = 1.5, shape = 17,
             position = position_jitter(width = 0.8, height = 0)) +
  scale_colour_manual(values = trt_colsb, breaks = trt_levelsb,
                      name = "Inflow (field)",
                      guide = guide_legend(order = 2, override.aes = list(alpha = 1))) +
  geom_vline(xintercept = fert_date, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = fert_date, y = Inf, label = "Fertilizer added",
           angle = 90, vjust = -0.4, hjust = 1.1, size = 3) +
  labs(x = "Date", y = "Conductivity (mS/cm)") +
  theme_bw()

conduct_study_comp

#humidity 
hum_graph <- ggplot(TDRpot_2025,
                        aes(x = Date_vis, y = Humidity_percent, colour = Salt)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_vline(xintercept = fert_date,
             linetype = "dashed", colour = "black", linewidth = 0.6) +
  annotate("text", x = fert_date, y = Inf, label = "Fertilizer added",
           angle = 90, vjust = -0.4, hjust = 1.1, size = 3) +
  labs(x = "Date", y = "Conductivity (mS/cm)", colour = "Salt") +
  theme_bw()

hum_graph

plot(Humidity_percent~Date, data= TDRpot_2025)
hum_study_comp <- ggplot() +
  geom_point(data = TDRpot_2025,
             aes(x = Date_vis, y = Humidity_percent, colour = Salt),
             alpha = 0.2, size = 1,
             position = position_jitter(width = 0.8, height = 0)) +
  scale_colour_manual(values = trt_colsp, breaks = trt_levelsp,
                      name = "Salinity (pots)",
                      guide = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2))) +
  new_scale_colour() +
  geom_point(data = TDRb_2025,
             aes(x = Date, y = Humidity, colour = Treatment),
             alpha = 0.6, size = 1.5, shape = 17,
             position = position_jitter(width = 0.8, height = 0)) +
  scale_colour_manual(values = trt_colsb, breaks = trt_levelsb,
                      name = "Inflow (field)",
                      guide = guide_legend(order = 2, override.aes = list(alpha = 1))) +
  geom_vline(xintercept = fert_date, linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = fert_date, y = Inf, label = "Fertilizer added",
           angle = 90, vjust = -0.4, hjust = 1.1, size = 3) +
  labs(x = "Date", y = "Conductivity (mS/cm)") +
  theme_bw()

hum_study_comp

##2. Models: Conductivity with Humidity as a cofactor####
###2.1 Humidity####
#testing if species is a significant part of model:
#
m_HumTDRp<-lmer(Humidity_percent ~ Salt*Biochar*Sp+(1|Date_vis)+(1|ID)+(1|Block), data=TDRpot_2025, REML = TRUE)
plot(m_HumTDRp) #no transformation needed

Anova_Soil_Humidity<-anova(m_HumTDRp, type = 3)
Anova_Soil_Humidity
#species is significant therefore seperate by species


plot(residuals(m_HumTDRp))
qqnorm(ranef(m_HumTDRp)$Date_vis[,1]); qqline(ranef(m_HumTDRp)$Date_vis[,1])
qqnorm(ranef(m_HumTDRp)$Block[,1]); qqline(ranef(m_HumTDRp)$Block[,1])

#pairwise comparisons for the timeline model
salt_mtdrh<-emmeans(m_HumTDRp, pairwise~Salt, adjust="tukey", type="response")
Sp_mtdrh<-emmeans(m_HumTDRp, pairwise~Sp, adjust="tukey", type="response")
Salt_Sp_mtdrh<-emmeans(m_HumTDRp, pairwise~Salt|Sp, adjust="tukey", type="response")
Sp_Salt_mtdrh<-emmeans(m_HumTDRp, pairwise~Sp|Salt, adjust="tukey", type="response")


salt_mtdrh
Sp_mtdrh
Sp_Salt_mtdrh
Salt_Sp_mtdrh
###2.2_Conductivity####
#testing if species is a significant part of model:
#
m_ConTDRp<-lmer(log(Conductivity_mS.cm) ~ log(Humidity_percent)+ Salt*Biochar*Date + (1|Sp), data=TDRpot_2025, REML = TRUE)
plot(m_ConTDRp) #no transformation needed

Anova_Soil_Conductivity<-anova(m_ConTDRp, type = 3)
Anova_Soil_Conductivity
#species is significant therefore seperate by species


plot(residuals(m_ConTDRp))
qqnorm(ranef(m_ConTDRp)$Date_vis[,1]); qqline(ranef(m_ConTDRp)$Date_vis[,1])
qqnorm(ranef(m_ConTDRp)$Block[,1]); qqline(ranef(m_ConTDRp)$Block[,1])

#pairwise comparisons for the timeline model
salt_mtdrc<-emmeans(m_ConTDRp, pairwise~Salt, adjust="tukey", type="response")
Sp_mtdrc<-emmeans(m_ConTDRp, pairwise~Sp, adjust="tukey", type="response")
Biochar_mtdrc<-emmeans(m_ConTDRp, pairwise~Biochar, adjust="tukey", type="response")
Salt_Sp_mtdrc<-emmeans(m_ConTDRp, pairwise~Salt|Sp, adjust="tukey", type="response")
Sp_Salt_mtdrc<-emmeans(m_ConTDRp, pairwise~Sp|Salt, adjust="tukey", type="response")


salt_mtdrc
Biochar_mtdrc
Sp_mtdrc
Sp_Salt_mtdrc
Salt_Sp_mtdrc
###

###improved salt models^\


m_F  <- lmer(log(Conductivity_mS.cm) ~ log(Humidity_percent) + Salt*Biochar*Date + (1|Sp),
             data = TDRpot_2025, REML = FALSE)
m_F2 <- lmer(log(Conductivity_mS.cm) ~ log(Humidity_percent) + Salt*Date + Biochar + (1|Sp),
             data = TDRpot_2025, REML = FALSE)
emmeans(m_F2, pairwise ~ Salt | Date, adjust = "tukey")   # salt differences on each date
emmeans(m_F2, pairwise ~ Date | Salt, adjust = "tukey")   # change over time within each salt level
anova(m_F2, m_F)      # LRT: do the extra Biochar interactions improve fit? no
AIC(m_F, m_F2)
#models seperated by species
##3. Models per species: ----
###3.1 Humidity loop####

#equation

#lmer model loop for growth models
library(lmerTest)   # lmer with p-values
library(emmeans)
species_list <- split(treepot, droplevels(treepot$Sp))
rhs <- "Salt * Biochar"
humidity_model <- list()

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


#4. Plots####
###4.1 conductivity per species####
ggplot(TDRpot_2025, aes(Date_vis, Conductivity_mS.cm, colour = Salt, fill = Salt)) +
  geom_vline(xintercept = fert_date, linetype = "dashed") +
  geom_point(alpha = 0.25, size = 1,
             position = position_jitter(width = 0.8, height = 0)) +
  stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.2, colour = NA) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
  facet_wrap(~ Sp, scales = "fixed") +
  scale_colour_manual(values = trt_colsp, breaks = trt_levelsp) +
  scale_fill_manual(values = trt_colsp, breaks = trt_levelsp) +
  scale_y_log10() +
  theme_bw()
###4.2 conductivity avg per species####
emm <- as.data.frame(emmeans(m_ConTDRp, ~ Salt | Sp|Biochar, type = "response"))
ggplot(emm, aes(Salt, response, colour = Salt, shape = Biochar)) +
  geom_pointrange(aes(ymin = lower.CL, ymax = upper.CL)) +
  facet_wrap(~ Sp, scales = "free_y") +
  scale_colour_manual(values = trt_colsp) +
  theme_bw()

### 4.3 humidity vs conductivity per species####
ggplot(TDRpot_2025, aes(Humidity_percent, Conductivity_mS.cm,
                        colour = Salt, shape = Biochar)) +
  geom_point(alpha = 0.5, size = 1.8) +
  facet_wrap(~ Sp) +
  scale_colour_manual(values = trt_colsp, breaks = trt_levelsp) +
  scale_shape_manual(values = 0:2) +
  scale_y_log10() +
  labs(x = "Humidity (%)", y = "Conductivity (mS/cm)",
       colour = "Salt", shape = "Species") +
  theme_bw()
#5.excel####
library(openxlsx)
library(emmeans)
library(lmerTest)

alpha      <- 0.1
skip_terms <- c("Humidity_percent", "Date_f")   # covariates / nuisance terms: no post-hoc

# name = sheet name (keep under 31 characters, no / \ ? * : [ ])
tdr_models <- list(
  "Humidity"           = m_HumTDRp,
  "Conductivity (log)" = m_ConTDRp     # or your final version, e.g. m1
)

# ---- helpers ------------------------------------------------------------------
get_anova <- function(m) {
  a <- as.data.frame(anova(m, type = 3))
  data.frame(Term = rownames(a), NumDF = a$NumDF, DenDF = a$DenDF,
             F = a[["F value"]], p = a[["Pr(>F)"]], row.names = NULL)
}

# which emmeans specs to run for a term (treatments are compared within species)
specs_for_term <- function(term) {
  v <- strsplit(term, ":", fixed = TRUE)[[1]]
  if (length(v) == 1) return(paste("~", v))
  if ("Sp" %in% v) {
    trt <- setdiff(v, "Sp")
    if (length(trt) == 1) return(paste("~", trt, "| Sp"))
    return(c(paste("~", trt[1], "|", trt[2], "* Sp"),
             paste("~", trt[2], "|", trt[1], "* Sp")))
  }
  if (length(v) == 2) return(c(paste("~", v[1], "|", v[2]), paste("~", v[2], "|", v[1])))
  paste("~", v[1], "|", paste(v[-1], collapse = " * "))
}

# contrasts -> tidy table with a direction column
tidy_con <- function(con, label) {
  df   <- as.data.frame(summary(con, type = "response"))   # ratios if the model has log(y)
  vcol <- intersect(c("ratio", "estimate"), names(df))[1]
  std  <- c("contrast", "estimate", "ratio", "SE", "df", "t.ratio", "z.ratio",
            "p.value", "lower.CL", "upper.CL", "null")
  by_cols <- setdiff(names(df), std)
  val   <- df[[vcol]]
  parts <- strsplit(gsub("[()]", "", as.character(df$contrast)), " [-/] ")
  dir <- mapply(function(p, v) {
    if (length(p) != 2 || is.na(v)) return(NA_character_)
    up <- if (vcol == "ratio") v > 1 else v > 0
    if (up) paste(p[1], ">", p[2]) else paste(p[2], ">", p[1])
  }, parts, val, USE.NAMES = FALSE)
  cont <- as.character(df$contrast)
  if (length(by_cols)) {
    by_txt <- apply(df[by_cols], 1, function(r)
      paste(paste(by_cols, trimws(r), sep = "="), collapse = ", "))
    cont <- paste0(cont, "  [", by_txt, "]")
  }
  data.frame(Comparison = label, Contrast = cont, Estimate = val,
             SE = df$SE, Direction = dir, p = df$p.value)
}

posthoc_term <- function(m, term) {
  do.call(rbind, lapply(specs_for_term(term), function(s)
    tidy_con(pairs(emmeans(m, as.formula(s)), adjust = "tukey"),
             sub("^~\\s*", "", s))))
}

model_info <- function(m) {
  vc <- as.data.frame(VarCorr(m))
  list(
    head = data.frame(
      Item  = c("Formula", "Observations", "Singular fit"),
      Value = c(paste(deparse(formula(m), width.cutoff = 500), collapse = " "),
                nobs(m), ifelse(isSingular(m), "YES", "no"))),
    vc = data.frame(Random_effect = vc$grp, Variance = vc$vcov, SD = vc$sdcor))
}

# ---- collect results ----------------------------------------------------------
collect <- lapply(names(tdr_models), function(nm) {
  m  <- tdr_models[[nm]]
  an <- get_anova(m)
  mf <- model.frame(m)
  sig_terms <- an$Term[!is.na(an$p) & an$p < alpha]
  keep <- vapply(strsplit(sig_terms, ":", fixed = TRUE), function(v)
    !any(v %in% skip_terms) &&
      !any(vapply(mf[intersect(v, names(mf))], is.numeric, logical(1))), logical(1))
  ph <- do.call(rbind, lapply(sig_terms[keep], function(t) {
    r <- tryCatch(posthoc_term(m, t), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    cbind(Term = t, ANOVA_p = an$p[an$Term == t], r)
  }))
  if (!is.null(ph)) ph <- ph[!is.na(ph$p) & ph$p < alpha, , drop = FALSE]
  if (!is.null(ph) && nrow(ph) == 0) ph <- NULL
  list(anova = an, posthoc = ph, info = model_info(m))
})
names(collect) <- names(tdr_models)

# ---- write workbook -----------------------------------------------------------
wb   <- createWorkbook()
hdr  <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2", border = "bottom")
s05  <- createStyle(textDecoration = "bold", fgFill = "#FFE699")   # p < 0.05
s10  <- createStyle(fgFill = "#FFF2CC")                            # 0.05 <= p < 0.1
pfmt <- createStyle(numFmt = "0.000000000")
ttl  <- createStyle(textDecoration = "bold", fontSize = 12)
warn <- createStyle(fontColour = "#C00000", textDecoration = "bold")

write_tbl <- function(sheet, df, row) {
  df <- as.data.frame(df)
  num <- setdiff(names(df)[sapply(df, is.numeric)], c("p", "ANOVA_p"))
  for (cn in num) df[[cn]] <- round(df[[cn]], 4)
  writeData(wb, sheet, df, startRow = row, headerStyle = hdr)
  pc <- which(names(df) %in% c("p", "ANOVA_p"))
  if (length(pc))
    addStyle(wb, sheet, pfmt, rows = row + seq_len(nrow(df)), cols = pc,
             gridExpand = TRUE, stack = TRUE)
  if ("p" %in% names(df)) {
    r10 <- which(!is.na(df$p) & df$p >= 0.05 & df$p < 0.1)
    r05 <- which(!is.na(df$p) & df$p < 0.05)
    if (length(r10)) addStyle(wb, sheet, s10, rows = row + r10, cols = seq_along(df),
                              gridExpand = TRUE, stack = TRUE)
    if (length(r05)) addStyle(wb, sheet, s05, rows = row + r05, cols = seq_along(df),
                              gridExpand = TRUE, stack = TRUE)
  }
  row + nrow(df) + 2
}

# Summary sheet
addWorksheet(wb, "Summary")
writeData(wb, "Summary", paste0("TDR pots: significant model terms (p < ", alpha,
                                "), Type III F-tests, Satterthwaite df"), startRow = 1)
addStyle(wb, "Summary", ttl, rows = 1, cols = 1)
sumtab <- do.call(rbind, lapply(names(collect), function(nm) {
  a <- collect[[nm]]$anova; a <- a[!is.na(a$p) & a$p < alpha, , drop = FALSE]
  if (nrow(a)) cbind(Model = nm, a)
}))
row <- 3
if (is.null(sumtab)) {
  writeData(wb, "Summary", "No term significant in any model.", startRow = row); row <- row + 2
} else row <- write_tbl("Summary", sumtab, row)
singular <- names(collect)[sapply(collect, function(x) x$info$head$Value[3] == "YES")]
if (length(singular)) {
  writeData(wb, "Summary", paste("Singular fit (a random-effect variance is ~0; p-values may be slightly unreliable):",
                                 paste(singular, collapse = ", ")), startRow = row)
  addStyle(wb, "Summary", warn, rows = row, cols = 1); row <- row + 2
}
writeData(wb, "Summary", c(
  "Bold yellow rows: p < 0.05; light yellow: p < 0.1.",
  "Contrasts: Tukey-adjusted within each 'by' group; only contrasts with p < alpha are listed.",
  "Direction: first level > second level (from the sign / ratio of the estimate).",
  "Log-conductivity model: estimates are ratios on the original scale (> 1 means the first level is higher); other models: differences.",
  "Means are averaged over the other factors and evaluated at the mean Humidity_percent."),
  startRow = row, colNames = FALSE)
setColWidths(wb, "Summary", cols = 1:6, widths = c(24, 28, 10, 10, 10, 14))

# Sig_contrasts sheet (all models stacked)
addWorksheet(wb, "Sig_contrasts")
allph <- do.call(rbind, lapply(names(collect), function(nm)
  if (!is.null(collect[[nm]]$posthoc)) cbind(Model = nm, collect[[nm]]$posthoc)))
if (is.null(allph)) writeData(wb, "Sig_contrasts", "No significant Tukey contrasts.") else
  write_tbl("Sig_contrasts", allph, 1)
setColWidths(wb, "Sig_contrasts", cols = 1:9, widths = c(20, 26, 12, 24, 44, 12, 12, 18, 14))

# one sheet per model
for (nm in names(collect)) {
  r  <- collect[[nm]]
  sh <- substr(nm, 1, 31)
  addWorksheet(wb, sh)
  writeData(wb, sh, paste0(nm, ": ANOVA (Type III, Satterthwaite)"), startRow = 1)
  addStyle(wb, sh, ttl, rows = 1, cols = 1)
  row <- write_tbl(sh, r$anova, 2)
  
  writeData(wb, sh, "Model", startRow = row); addStyle(wb, sh, ttl, rows = row, cols = 1)
  row <- write_tbl(sh, r$info$head, row + 1) - 1
  row <- write_tbl(sh, r$info$vc, row + 1)
  if (r$info$head$Value[3] == "YES") {
    writeData(wb, sh, "Singular fit: at least one random-effect variance is ~0. Consider simplifying the random terms.",
              startRow = row - 1)
    addStyle(wb, sh, warn, rows = row - 1, cols = 1); row <- row + 1
  }
  
  writeData(wb, sh, paste0("Tukey-adjusted contrasts for significant terms (p < ", alpha, ")"), startRow = row)
  addStyle(wb, sh, ttl, rows = row, cols = 1)
  if (is.null(r$posthoc)) {
    writeData(wb, sh, "No significant contrasts (or the significant terms are covariates).", startRow = row + 1)
  } else write_tbl(sh, r$posthoc, row + 1)
  setColWidths(wb, sh, cols = 1:8, widths = c(28, 14, 24, 44, 12, 12, 18, 16))
}

saveWorkbook(wb, file.path(table_dir, "TDR_pots_model_summary.xlsx"), overwrite = TRUE)

