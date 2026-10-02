# Shared data preparation for sample.Rmd and replication.Rmd.
# Sourced with source('common.R', local = TRUE) so that the Rmd's `params` are visible.

library(dplyr); library(tidyr); library(readxl); library(readr)
library(forcats); library(stringr); library(ggplot2); library(knitr)

# ---- Czech survey ----
# dataAll.RData is made by preprocess.Rmd: raw answers of the four institutional
# exports with short variable names; one test response (8 Apr 2019) already removed.
load("dataAll.RData")
cz <- data

flow <- tibble(Step = "dataAll.RData (AV, UK, MU, JU merged; test response removed)", N = nrow(cz))
if (params$time_cutoff_min > 0) {
  cz <- filter(cz, as.numeric(cas_sec) >= 60 * params$time_cutoff_min)
  flow <- add_row(flow, Step = paste0("Completion time >= ", params$time_cutoff_min, " min"), N = nrow(cz))
}
if (params$basic_only) {
  cz <- filter(cz, zakl_vyzkum == "Ano")
  flow <- add_row(flow, Step = "Does basic research (Q3 = Ano)", N = nrow(cz))
}
if (params$drop_math) {
  cz <- filter(cz, is.na(spec_hlavni) | spec_hlavni != "Mathematics")
  flow <- add_row(flow, Step = "Mathematicians removed", N = nrow(cz))
}

# ---- Baker (Nature 2016) survey ----
bk_raw <- read_xlsx("../data_and_R_src/dataNature.xlsx", skip = 1, .name_repair = "minimal")
h <- names(bk_raw)
# Columns are selected by position, so check that each position holds the expected question
stopifnot(grepl("mainstream media", h[13]), grepl("scientific journals", h[14]),
          grepl("conferences", h[15]), grepl("colleagues", h[16]), grepl("^No", h[19]),
          grepl("major problem in my field", h[24]), grepl("most often mea", h[50]),
          grepl("rarely detracts", h[51]), grepl("^Fraud", h[52]), grepl("^Bad luck", h[65]),
          grepl("own experiments", h[79]), grepl("someone else", h[80]))

bk_cols <- c(rep_Mainstream = 13, rep_Casopisy = 14, rep_Diskuse = 15, rep_Kolegove = 16,
             rep_Jinde = 17, rep_Ne = 19,
             krize_rep = 20, rep_podil = 22,
             nerep_problem = 24, nerep_problemAll = 25, nerep_chyba = 50, nerep_valid = 51,
             fakt_podvod = 52, fakt_karier = 53, fakt_dohled = 54, fakt_recenz = 55,
             fakt_selekc = 56, fakt_replik = 57, fakt_statis = 58, fakt_odborn = 59,
             fakt_data = 60, fakt_dokume = 61, fakt_metody = 62, fakt_variab = 63,
             fakt_design = 64, fakt_smula = 65,
             nerep_jaMuj = 79, nerep_jaCizi = 80, nerep_oniMuj = 85, spec_hlavni = 91)
bk <- bk_raw[, bk_cols]
names(bk) <- names(bk_cols)

# ---- Harmonised answers (Baker vs Czech) ----
en <- c(
  # Czech
  "Nevím" = "I don't know", "Ano" = "Yes", "Ne" = "No", "Nepamatuji se" = "I can't remember",
  "Nedělám experimenty" = "I don't do experiments",
  "Ano, významná krize" = "Significant crisis", "Ano, nevýznamná krize" = "Slight crisis",
  "Ne, žádná krize není" = "No crisis",
  "Vždy" = "Always", "Velmi často" = "Very often", "Někdy" = "Sometimes",
  "Zřídka" = "Rarely", "Nikdy" = "Never",
  "Rozhodně souhlasím" = "Strongly agree", "Spíše souhlasím" = "Agree",
  "Ani souhlasím/ani nesouhlasím" = "Neither agree nor disagree",
  "Spíše nesouhlasím" = "Disagree", "Rozhodně nesouhlasím" = "Strongly disagree",
  # Baker
  "There is a significant crisis of reproducibility" = "Significant crisis",
  "There is a slight crisis of reproducibility" = "Slight crisis",
  "There is no crisis of reproducibility" = "No crisis",
  "Always contributes" = "Always", "Very often contributes" = "Very often",
  "Sometimes contributes" = "Sometimes", "Rarely contributes" = "Rarely",
  "Never contributes" = "Never")
tr <- function(x) { x <- as.character(x); i <- !is.na(x) & x %in% names(en); x[i] <- en[x[i]]; x }

field <- function(x) case_when(
  x %in% c("Astronomy and Planetary Science", "Astronomy and planetary science") ~ "Astronomy",
  x %in% c("Other", "Other specialization") ~ "Other",
  x == "Earth and Environmental Science" ~ "Earth & environment",
  TRUE ~ x)

compared <- c("krize_rep", "rep_podil", "nerep_problem", "nerep_chyba", "nerep_valid",
              grep("^fakt_", names(bk_cols), value = TRUE),
              "nerep_jaMuj", "nerep_jaCizi", "nerep_oniMuj")

d <- bind_rows(
  cz %>% select(all_of(compared), spec_hlavni) %>% mutate(across(everything(), as.character), Survey = "Czech"),
  bk %>% select(all_of(compared), spec_hlavni) %>% mutate(across(everything(), as.character), Survey = "Baker")
) %>%
  mutate(across(all_of(compared), tr), field = field(spec_hlavni),
         Survey = factor(Survey, levels = c("Baker", "Czech")))

# ---- Field weights for adjusted comparisons ----
# Czech respondents are reweighted so that their field mix (Q50) matches Baker's.
# Czech respondents with no field, or with a field Baker does not have, get no weight
# and are left out of the adjusted figures.
fw <- d %>% filter(!is.na(field)) %>% count(Survey, field) %>%
  group_by(Survey) %>% mutate(share = n / sum(n)) %>% ungroup() %>%
  select(-n) %>% pivot_wider(names_from = Survey, values_from = share) %>%
  mutate(w = Baker / Czech)
d <- d %>% left_join(select(fw, field, w), by = "field") %>%
  mutate(w = if_else(Survey == "Baker", 1, w))
ADJ <- "Czech (field-adjusted)"

# ---- Helper functions ----
DK <- "I don't know"
pal <- c(Baker = "#4C72B0", Czech = "#DD8452", "Czech (field-adjusted)" = "#C44E52")
default_versions <- function(dk) setNames(list(character(0), dk),
                                          c("Including don't know", "Excluding don't know"))

# Several versions of one item, each dropping a different set of answers
make_versions <- function(x, versions) bind_rows(lapply(names(versions), function(v)
  x %>% filter(!value %in% versions[[v]]) %>% mutate(Version = v))) %>%
  mutate(Version = factor(Version, levels = names(versions)))

# Adds the field-adjusted Czech rows (weight w); Baker and unadjusted Czech get weight 1
add_adjusted <- function(x, adjust) {
  x <- mutate(x, Survey = as.character(Survey), wt = 1)
  if (adjust) x <- bind_rows(x, x %>% filter(Survey == "Czech", !is.na(w)) %>% mutate(Survey = ADJ, wt = w))
  mutate(x, Survey = factor(Survey, levels = c("Baker", "Czech", ADJ)))
}

# Chi-square test, Baker vs unadjusted Czech
chi_p <- function(s, v) {
  keep <- s %in% c("Baker", "Czech")
  t <- table(droplevels(factor(s[keep])), droplevels(factor(v[keep])))
  if (min(dim(t)) < 2) NA_real_ else suppressWarnings(chisq.test(t)$p.value)
}
fmt_p <- function(p) format.pval(p, digits = 2, eps = 1e-4)

# One categorical question, Baker vs Czech (and field-adjusted Czech)
compare_cat <- function(var, levels, dk = DK, title = var, versions = default_versions(dk), adjust = TRUE) {
  x <- d %>% select(Survey, w, value = all_of(var)) %>% filter(!is.na(value))
  b <- make_versions(x, versions) %>% add_adjusted(adjust)
  tab <- b %>% group_by(Version, Survey, value) %>% summarise(n = sum(wt), .groups = "drop") %>%
    group_by(Version, Survey) %>% mutate(N = sum(n), pct = 100 * n / N) %>% ungroup() %>%
    mutate(value = factor(value, levels = levels))
  p <- b %>% group_by(Version) %>% summarise(p = chi_p(Survey, value), .groups = "drop")
  Ns <- b %>% filter(Survey != ADJ) %>% count(Version, Survey) %>%
    group_by(Version) %>% summarise(lab = paste0(Survey, " N=", n, collapse = ", ")) %>%
    left_join(p, by = "Version") %>%
    mutate(lab = paste0(Version, ": ", lab, ", chi-square p (unadjusted) = ", fmt_p(p)))
  plot <- ggplot(tab, aes(x = pct, y = fct_rev(value), fill = Survey)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    geom_text(aes(label = sprintf("%.1f", pct)), position = position_dodge(width = 0.85),
              hjust = -0.15, size = 2.6) +
    facet_wrap(~ Version) + scale_fill_manual(values = pal) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(title = title, x = "% of respondents", y = NULL, fill = NULL,
         caption = paste(Ns$lab, collapse = "\n")) +
    theme_minimal() + theme(legend.position = "bottom")
  pw <- tab %>% select(Version, Answer = value, Survey, pct) %>% mutate(pct = round(pct, 1)) %>%
    pivot_wider(names_from = Survey, values_from = pct, values_fill = 0, names_prefix = "% ")
  nw <- tab %>% filter(Survey != ADJ) %>% select(Version, Answer = value, Survey, n) %>%
    pivot_wider(names_from = Survey, values_from = n, values_fill = 0, names_prefix = "n ")
  list(plot = plot, table = left_join(nw, pw, by = c("Version", "Answer")) %>% arrange(Version, Answer))
}

# A battery of items with the same scale, Baker vs Czech (and field-adjusted Czech)
compare_battery <- function(vars, labels, levels, top, dk = DK, title = "",
                            versions = default_versions(dk), adjust = TRUE) {
  x <- d %>% select(Survey, w, all_of(vars)) %>%
    pivot_longer(all_of(vars), names_to = "item", values_to = "value") %>%
    filter(!is.na(value)) %>% mutate(item = unname(labels[item]))
  b <- make_versions(x, versions) %>% add_adjusted(adjust)
  tab <- b %>% group_by(Version, Survey, item, value) %>% summarise(n = sum(wt), .groups = "drop") %>%
    group_by(Version, Survey, item) %>% mutate(N = sum(n), pct = 100 * n / N) %>% ungroup() %>%
    mutate(value = factor(value, levels = levels))
  ord <- tab %>% filter(Version == last(names(versions)), Survey == "Baker", value %in% top) %>%
    group_by(item) %>% summarise(s = sum(pct)) %>% arrange(s) %>% pull(item)
  ord <- c(setdiff(unique(tab$item), ord), ord)
  tab <- mutate(tab, item = factor(item, levels = ord))
  p <- b %>% group_by(Version, item) %>% summarise(p = chi_p(Survey, value), .groups = "drop")
  Ns <- b %>% filter(Survey != ADJ) %>% count(Version, Survey, item) %>%
    pivot_wider(names_from = Survey, values_from = n, names_prefix = "N ")
  summ <- tab %>% group_by(Version, Survey, item) %>%
    summarise(top = round(sum(pct[value %in% top]), 1), .groups = "drop") %>%
    pivot_wider(names_from = Survey, values_from = top, names_prefix = "% ") %>%
    left_join(Ns, by = c("Version", "item")) %>% left_join(p, by = c("Version", "item")) %>%
    mutate(item = factor(item, levels = ord), p = fmt_p(p)) %>% arrange(Version, desc(item))
  plot <- ggplot(tab, aes(x = pct, y = item, fill = fct_rev(value))) +
    geom_col(width = 0.75) + facet_grid(Version ~ Survey) +
    scale_fill_brewer(palette = "RdYlBu", direction = -1, name = NULL) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(title = title, x = "% of respondents", y = NULL,
         caption = paste0("Items ordered by Baker's share of '", paste(top, collapse = "' + '"),
                          "'. p = chi-square test, Baker vs unadjusted Czech.")) +
    theme_minimal() + theme(legend.position = "bottom")
  list(plot = plot, table = summ)
}

# Czech-only frequency table and bar chart
describe <- function(var, title = var, lump = FALSE) {
  x <- cz %>% transmute(value = as.character(.data[[var]]))
  if (lump) x <- mutate(x, value = fct_lump_min(value, 5, other_level = "Other (free text)"))
  tab <- x %>% mutate(value = replace_na(as.character(value), "(no answer)")) %>%
    count(value) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% arrange(desc(n))
  plot <- ggplot(tab, aes(x = pct, y = fct_reorder(value, n))) +
    geom_col(fill = pal["Czech"]) + geom_text(aes(label = n), hjust = -0.2, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = title, x = "% of Czech respondents", y = NULL,
         caption = paste0("N = ", nrow(cz))) + theme_minimal()
  list(plot = plot, table = tab)
}

# Czech-only battery (same scale for all items)
describe_battery <- function(vars, title, levels = NULL) {
  x <- cz %>% select(all_of(vars)) %>% pivot_longer(everything(), names_to = "item") %>%
    mutate(value = replace_na(as.character(value), "(no answer)"))
  tab <- x %>% count(item, value) %>% group_by(item) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup()
  if (!is.null(levels)) tab <- mutate(tab, value = factor(value, levels = c(levels, "(no answer)")))
  plot <- ggplot(tab, aes(x = pct, y = item, fill = value)) + geom_col(width = 0.75) +
    scale_fill_brewer(palette = "Spectral", name = NULL) +
    labs(title = title, x = "% of Czech respondents", y = NULL, caption = paste0("N = ", nrow(cz))) +
    theme_minimal() + theme(legend.position = "bottom")
  list(plot = plot, table = tab %>% select(item, value, pct) %>% pivot_wider(names_from = value, values_from = pct, values_fill = 0))
}
