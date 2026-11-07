# Analysis of de-identified, linkable demographic and prepared data.
# Execute from the folder containing the six country CSVs.
# Exact original DE age imputation value is preserved.

# VERSION: corrected demographic predictors — no responder-age requirement
# =============================================================================
# CROSS-CULTURAL ITCQ / ECCQ ANALYSIS
# Germany (DE), China (CN), and Japan (JP)
# =============================================================================
#
# PURPOSE
# -------
# This script reads the six files produced by prepare_DE.R, prepare_CN.R, and
# prepare_JP.R, combines them, excludes J1-J3 for the current paper, and runs:
#
#   RQ1. Are Final items rated more highly than Excluded and Filler items?
#   RQ2. Do item-set ratings and item profiles differ between countries?
#   RQ3. Which individual items show culturally meaningful differences,
#        especially curiosity consensus combined with desirability divergence?
#
# ITCQ and ECCQ are analysed separately. Curiosity and Desirability are shown
# separately in figures, but modelled jointly so Country x Condition effects and
# curiosity-desirability gaps can be estimated.
#
# IMPORTANT SCORING DECISION
# --------------------------
# Response is coded from 1 to 7. The preparation scripts reverse-score the
# explicitly reverse-worded items with 8 - value: ITCQ SR, Ex2, Ex10; and ECCQ
# N3, N5, Ex10, Ex13. Higher values therefore consistently indicate greater
# attributed curiosity or desirability.
#
# EXPECTED INPUT FILES (place them in DATA_DIR)
# ------------------------------------------------
# ITCQ_ECCQ_DE_prepared.csv    demo_DE.csv
# ITCQ_ECCQ_CN_prepared.csv    demo_CN.csv
# ITCQ_ECCQ_JP_prepared.csv    demo_JP.csv
#
# OUTPUT
# ------
# OUTPUT_DIR/
#   tables/   CSV tables for manuscript and supplementary materials
#   figures/  PNG and PDF versions of figures
#   models/   fitted R model objects (.rds)
#   logs/     data checks, settings, warnings, and session information
#
# =============================================================================

# ---- 0. USER SETTINGS --------------------------------------------------------

DATA_DIR <- "."
OUTPUT_DIR <- "cross_cultural_results"

# A difference smaller than this is treated as practically negligible.
# With EQUIVALENCE_MARGIN = 0.50, curiosity profiles are called equivalent only
# when all relevant pairwise 95% CIs lie fully inside -0.5 to +0.5 points.
EQUIVALENCE_MARGIN <- 0.50

# Item-level p-values are corrected with Benjamini-Hochberg FDR.
ALPHA <- 0.05

# Reference categories for the modelled categorical covariates.
# These values must match the standardized values in the prepared files.
SEX_REFERENCE <- "Girl"
EDUCATION_REFERENCE <- "Undergraduate / bachelor"

# Nominal ages: ITCQ 5-23 months; ECCQ 2-5 years (operationalised as 24-71
# months). The script always reports age-range violations. It retains them by
# default because exclusion should be an explicit research decision.
EXCLUDE_OUT_OF_RANGE_AGES <- FALSE

# Number of labels allowed on the cultural-divergence figures.
MAX_LABELS_PER_DIVERGENCE_PLOT <- 15

# ---- 1. PACKAGES -------------------------------------------------------------

required_packages <- c(
  "tidyverse", "lme4", "lmerTest", "emmeans", "ggrepel",
  "sandwich", "performance", "writexl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running the analysis:\n  ",
    paste(missing_packages, collapse = ", "),
    "\nFor example: install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(ggrepel)
})

# Sum-to-zero contrasts make omnibus tests easier to interpret in models with
# interactions. Save the user's setting and restore it when the script exits.
old_contrasts <- options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
# Contrast options apply for the duration of this stand-alone R session.

theme_set(theme_bw(base_size = 11))

emmeans::emm_options(
  lmer.df = "satterthwaite",
  lmerTest.limit = 25000
)
# ---- 2. OUTPUT FOLDERS AND HELPERS -------------------------------------------

table_dir  <- file.path(OUTPUT_DIR, "tables")
figure_dir <- file.path(OUTPUT_DIR, "figures")
model_dir  <- file.path(OUTPUT_DIR, "models")
log_dir    <- file.path(OUTPUT_DIR, "logs")

walk(c(OUTPUT_DIR, table_dir, figure_dir, model_dir, log_dir),
     dir.create, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, filename) {
  readr::write_csv(x, file.path(table_dir, filename), na = "")
}

save_figure <- function(plot, filename, width, height) {
  ggsave(file.path(figure_dir, paste0(filename, ".png")), plot,
         width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0(filename, ".pdf")), plot,
         width = width, height = height, device = "pdf", bg = "white")
}

# Convert emmeans or contrast output into a consistently named tibble.
tidy_emm <- function(x, infer = c(TRUE, TRUE)) {
  as.data.frame(summary(x, infer = infer)) %>%
    as_tibble() %>%
    rename_with(~str_replace_all(.x, c("lower.CL" = "conf.low",
                                       "upper.CL" = "conf.high",
                                       "asymp.LCL" = "conf.low",
                                       "asymp.UCL" = "conf.high")))
}

# ---- 3. READ AND COMBINE PREPARED DATA --------------------------------------

response_files <- c(
  DE = "ITCQ_ECCQ_DE_prepared.csv",
  CN = "ITCQ_ECCQ_CN_prepared.csv",
  JP = "ITCQ_ECCQ_JP_prepared.csv"
)

demo_files <- c(DE = "demo_DE.csv", CN = "demo_CN.csv", JP = "demo_JP.csv")


all_required_files <- file.path(DATA_DIR, c(response_files, demo_files))
missing_files <- all_required_files[!file.exists(all_required_files)]
if (length(missing_files) > 0) {
  stop("The following prepared files are missing:\n  ",
       paste(missing_files, collapse = "\n  "))
}


responses_raw <- imap_dfr(response_files, function(filename, expected_country) {
  
  dat <- read_csv(
    file.path(DATA_DIR, filename),
    show_col_types = FALSE,
    col_types = cols(Participant_ID = col_character())
  )
  
  required <- c(
    "Country", "Questionnaire", "Condition", "Participant_ID",
    "Age_months", "sex", "Item", "Response"
  )
  
  if (!all(required %in% names(dat))) {
    stop(filename, " is missing required columns.")
  }
  
  if (any(dat$Country != expected_country, na.rm = TRUE)) {
    stop(filename, " contains an incorrect Country value.")
  }
  
  dat
})

# Read demographic files as character data because the country questionnaires
# contain different fields and sometimes encode the same concept differently.
# The standardized fields are converted below; every original field remains in
# the combined supplementary demographic file.
demographics_raw <- imap_dfr(demo_files, function(filename, expected_country) {
  dat <- read_csv(
    file.path(DATA_DIR, filename), show_col_types = FALSE,
    col_types = cols(.default = col_character()), name_repair = "unique"
  )
  if (!"Country" %in% names(dat)) dat$Country <- expected_country

  required <- c(
    "Country", "Questionnaire", "Condition", "Participant_ID",
    "responder_education"
  )

  if (!all(required %in% names(dat))) {
    stop(filename, " is missing one or more required demographic columns: ",
         paste(setdiff(required, names(dat)), collapse = ", "))
  }

  dat
})

demographics_raw <- demographics_raw %>%
  mutate(
    responder_education = as.character(responder_education),
    other_parent_education = as.character(other_parent_education),
    
    responder_education = case_when(
      Country == "CN" ~ coalesce(
        responder_education,
        other_parent_education
      ),
      TRUE ~ responder_education
    ),
    
    responder_education_source = case_when(
      Country == "CN" & !is.na(responder_education) ~
        "Responder",
      TRUE ~ responder_education_source
    )
  )

## replace missing ages with group mean
demographics_raw$Age_months <- as.numeric(demographics_raw$Age_months)
# Reproducibility correction: the ORIGINAL script calculated this mean from
# all German demographic source rows, including incompletely surveyed records
# subsequently excluded from the analysis. Those excluded records must not be
# released. The exact original mean of the 188 records with valid ages is fixed
# here so the published analysis reproduces the original five age imputations.
# In a new analysis, consider recalculating from eligible participants only.
de_mean <- 35.47961170212766

if (!is.finite(de_mean)) {
  stop("A valid German mean child age could not be calculated.")
}

demographics_raw$Age_months[
  demographics_raw$Country == "DE" &
    !is.na(demographics_raw$Age_months) &
    demographics_raw$Age_months < 1
] <- de_mean


demographics_raw$exclusion_reason[
  demographics_raw$Country == "DE" & demographics_raw$Participant_ID == "4"
] <- "Test Response"



# Exclude J1-J3 from the current paper. They remain in the prepared source files
# and can be restored later simply by deleting this filter.
responses <- responses_raw %>%
  filter(!Item %in% c("J1", "J2", "J3")) %>%
  mutate(
    Country = factor(Country, levels = c("DE", "CN", "JP")),
    Questionnaire = factor(Questionnaire, levels = c("ITCQ", "ECCQ")),
    Condition = factor(Condition, levels = c("Curiosity", "Desirability")),
    sex = factor(sex),
    Age_months = as.numeric(Age_months),
    Response = as.numeric(Response),
    # IDs can overlap across countries, questionnaires, or conditions.
    PersonID = interaction(
      Country, Questionnaire, Condition, Participant_ID,
      drop = TRUE
    )
  )


## replace missing ages with group mean

responses$Age_months[
  responses$Country == "DE" &
    !is.na(responses$Age_months) &
    responses$Age_months < 1
] <- de_mean

# exclude test response
responses <- responses[
  !(responses$Country == "DE" & responses$Participant_ID == "4"),
]


# Assign item sets using the canonical labels from the final-item workbook.
responses <- responses %>%
  mutate(
    Item_group = case_when(
      str_detect(Item, "^F([1-9]|1[0-6])$") ~ "Filler",
      str_detect(Item, "^Ex([1-9]|1[0-4])$") ~ "Excluded",
      Questionnaire == "ITCQ" &
        str_detect(Item, "^(S([1-9]|10)|SR|Iv[1-5]|Ia[1-7])$") ~ "Final",
      Questionnaire == "ECCQ" &
        str_detect(Item, "^(N[1-9]|In[1-7]|Ia[1-9])$") ~ "Final",
      TRUE ~ NA_character_
    ),
    Item_group = factor(Item_group, levels = c("Final", "Excluded", "Filler")),
    nominal_age = case_when(
      Questionnaire == "ITCQ" ~ Age_months >= 5 & Age_months <= 23,
      Questionnaire == "ECCQ" ~ Age_months >= 24 & Age_months <= 71,
      TRUE ~ FALSE
    )
  )

if (any(is.na(responses$Item_group))) {
  unknown <- responses %>% filter(is.na(Item_group)) %>% distinct(Questionnaire, Item)
  write_table(unknown, "ERROR_unclassified_items.csv")
  stop("Some item labels could not be classified. See ERROR_unclassified_items.csv")
}

if (any(!is.na(responses$Response) & !responses$Response %in% 1:7)) {
  stop("At least one non-missing response falls outside the expected 1-7 range.")
}

# Record age flags before applying the optional exclusion.
age_range_flags <- responses %>%
  distinct(Country, Questionnaire, Condition, Participant_ID, PersonID,
           Age_months, nominal_age) %>%
  arrange(Questionnaire, Country, Age_months)
write_table(age_range_flags, "00_age_range_flags.csv")


# Models cannot use missing outcomes. Missingness is quantified before removal.
data_quality <- responses %>%
  group_by(Questionnaire, Country, Condition) %>%
  summarise(
    participants = n_distinct(PersonID),
    response_rows = n(),
    missing_responses = sum(is.na(Response)),
    missing_percent = 100 * mean(is.na(Response)),
    out_of_range_age_participants = n_distinct(PersonID[!nominal_age]),
    duplicate_person_item_rows = sum(duplicated(paste(PersonID, Item))),
    .groups = "drop"
  )
write_table(data_quality, "00_data_quality.csv")

# Add the responder's education to each response row using the confirmed raw
# column names. The join includes questionnaire and condition because
# Participant_ID values can overlap across questionnaires and conditions.
expected_education_levels <- c(
  "Lower secondary or below",
  "Upper secondary",
  "Undergraduate / bachelor",
  "Vocational / professional",
  "Postgraduate"
)

responder_demographics <- demographics_raw %>%
  transmute(
    Country = factor(Country, levels = c("DE", "CN", "JP")),
    Questionnaire = factor(Questionnaire, levels = c("ITCQ", "ECCQ")),
    Condition = factor(Condition, levels = c("Curiosity", "Desirability")),
    Participant_ID = as.character(Participant_ID),
    responder_education = stringr::str_squish(responder_education)
  ) %>%
  mutate(
    responder_education = na_if(responder_education, "")
  )

unexpected_education_levels <- responder_demographics %>%
  filter(!is.na(responder_education),
         !responder_education %in% expected_education_levels) %>%
  distinct(responder_education) %>%
  arrange(responder_education)


responder_demographics <- responder_demographics %>%
  mutate(
    responder_education = trimws(responder_education),
    responder_education = case_when(
      responder_education == "College / associate" ~ "Postgraduate",
      responder_education == "Other" ~ "Vocational / professional",
      responder_education %in% c(
        "Lower secondary or below",
        "Upper secondary"
      ) ~ "Upper secondary or below",
      TRUE ~ responder_education
    ),
    responder_education = relevel(
      factor(responder_education),
      ref = "Undergraduate / bachelor"
    )
  )

analysis_data_unfiltered <- responses %>%
  filter(!is.na(Response)) %>%
  mutate(Participant_ID = as.character(Participant_ID)) %>%
  left_join(
    responder_demographics,
    by = c("Country", "Questionnaire", "Condition", "Participant_ID")
  )

covariate_missingness <- analysis_data_unfiltered %>%
  distinct(
    Country, Questionnaire, Condition, PersonID,
    Age_months, sex, responder_education
  ) %>%
  group_by(Questionnaire, Country, Condition) %>%
  summarise(
    participants = n(),
    missing_child_age = sum(is.na(Age_months)),
    missing_child_sex = sum(is.na(sex) | stringr::str_squish(as.character(sex)) == ""),
    missing_responder_education = sum(is.na(responder_education)),
    .groups = "drop"
  )
write_table(covariate_missingness, "00_model_covariate_missingness.csv")


analysis_data <- analysis_data_unfiltered %>%
  mutate(
    sex = stringr::str_squish(as.character(sex)),
    sex = na_if(sex, ""),
    sex_binary = case_when(
      trimws(sex) %in% c("Girl", "Female") ~ "Girl",
      trimws(sex) %in% c("Boy", "Male")    ~ "Boy",
      TRUE                                 ~ NA_character_
    ),
    sex_binary = relevel(
      factor(sex_binary),
      ref = "Girl"
    )
  ) %>%
  filter(
    !is.na(Age_months),
    !is.na(sex),
    !is.na(responder_education)
  ) %>%
  group_by(Questionnaire) %>%
  mutate(child_age_z = as.numeric(scale(Age_months))) %>%
  ungroup() %>%
  mutate(
    sex = factor(sex),
    responder_education = factor(
      responder_education,
      levels = c(
        EDUCATION_REFERENCE,
        "Upper secondary or below",
        "Vocational / professional",
        "Postgraduate"
      )
    )
  )

# Apply treatment coding only to child sex and responder education so their
# coefficients and planned contrasts use the requested reference categories.
# Other factors retain the script's sum-to-zero coding for omnibus tests.
set_treatment_reference <- function(x, reference, variable_name) {
  x <- droplevels(factor(x))

  if (!reference %in% levels(x)) {
    stop(
      "Reference level '", reference, "' is absent from ", variable_name,
      ". Observed levels: ", paste(levels(x), collapse = ", ")
    )
  }

  if (nlevels(x) < 2) {
    stop(variable_name, " has fewer than two observed levels.")
  }

  x <- relevel(x, ref = reference)
  contrasts(x) <- contr.treatment(nlevels(x), base = 1)
  x
}

# Keep a reproducible copy of the exact long dataset used by the models.
write_table(analysis_data, "00_analysis_dataset.csv")
write_table(demographics_raw, "00_combined_demographics_all_fields.csv")

# ---- 4. SAMPLE DESCRIPTION --------------------------------------------------

sample_summary <- analysis_data %>%
  distinct(Country, Questionnaire, Condition, PersonID, Age_months, sex) %>%
  group_by(Questionnaire, Country, Condition) %>%
  summarise(
    N = n(),
    age_mean = mean(Age_months),
    age_sd = sd(Age_months),
    age_min = min(Age_months),
    age_max = max(Age_months),
    n_boy = sum(sex == "Boy"),
    n_girl = sum(sex == "Girl"),
    n_other_unknown = sum(!sex %in% c("Boy", "Girl")),
    .groups = "drop"
  )
write_table(sample_summary, "01_sample_summary.csv")

# Participant-by-set means are used for descriptive means and SDs. Computing an
# SD over all raw item responses would mix between-person and between-item
# variation and would overstate the effective sample size.
participant_set_scores <- analysis_data %>%
  group_by(Questionnaire, Country, Condition, PersonID, sex, Item_group) %>%
  summarise(Set_mean = mean(Response), n_items_answered = n(), .groups = "drop")

set_descriptives <- participant_set_scores %>%
  group_by(Questionnaire, Country, Condition, Item_group) %>%
  summarise(
    N = n(), mean = mean(Set_mean), SD = sd(Set_mean),
    SE = SD / sqrt(N),
    conf.low = mean - qt(.975, pmax(N - 1, 1)) * SE,
    conf.high = mean + qt(.975, pmax(N - 1, 1)) * SE,
    .groups = "drop"
  )
write_table(participant_set_scores, "02_participant_item_set_scores.csv")
write_table(set_descriptives, "02_item_set_descriptives.csv")

# Descriptive summaries by child sex complement the primary model, which tests
# whether child-sex differences vary between Curiosity and Desirability ratings.
sex_descriptives <- participant_set_scores %>%
  group_by(Questionnaire, Country, Condition, Item_group, sex) %>%
  summarise(
    N = n(), mean = mean(Set_mean), SD = sd(Set_mean),
    SE = SD/sqrt(N),
    conf.low = mean-qt(.975,pmax(N-1,1))*SE,
    conf.high = mean+qt(.975,pmax(N-1,1))*SE,
    .groups = "drop"
  )

sex_descriptive_differences <- sex_descriptives %>%
  filter(as.character(sex) %in% c("Boy","Girl")) %>%
  select(Questionnaire,Country,Condition,Item_group,sex,N,mean) %>%
  pivot_wider(names_from=sex,values_from=c(N,mean)) %>%
  mutate(Girl_minus_Boy=mean_Girl-mean_Boy)

write_table(sex_descriptives,"02_item_set_descriptives_by_child_sex.csv")
write_table(sex_descriptive_differences,
            "02_descriptive_girl_minus_boy_item_set_differences.csv")

# ---- 5. RQ1/RQ2: ITEM-SET MIXED MODELS -------------------------------------
#
# Each participant rates many items, and each item is rated by many participants.
# Crossed random intercepts account for both sources of dependence.


fit_group_model <- function(questionnaire_name) {
  
  # Keep only complete observations for the requested questionnaire.
  d <- analysis_data %>%
    filter(
      Questionnaire == questionnaire_name,
      !is.na(Response),
      !is.na(sex_binary),
      !is.na(responder_education),
      !is.na(child_age_z)
    ) %>%
    droplevels()
  
  # Set the reference categories used in the model.
  d$sex_binary <- set_treatment_reference(
    d$sex_binary,
    SEX_REFERENCE,
    "sex_binary"
  )
  
  d$responder_education <- set_treatment_reference(
    d$responder_education,
    EDUCATION_REFERENCE,
    "responder_education"
  )
  
  # Main group-level mixed model.
  model <- lmer(
    Response ~ Country * Condition * Item_group +
      Condition * sex_binary +
      Condition * responder_education +
      child_age_z +
      (1 | PersonID) +
      (1 | Item),
    data = d,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  # Omnibus Type-III tests.
  omnibus <- as.data.frame(
    anova(
      model,
      type = 3,
      ddf = "Satterthwaite"
    )
  ) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # Adjusted means for each Country x Condition x Item-group cell.
  means <- emmeans(
    model,
    ~ Country * Condition * Item_group
  ) %>%
    tidy_emm() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # RQ1: planned item-group comparisons within each country and condition.
  group_emm <- emmeans(
    model,
    ~ Item_group | Country * Condition
  )
  
  planned <- contrast(
    group_emm,
    method = list(
      "Final - Excluded" = c(1, -1, 0),
      "Final - Filler" = c(1, 0, -1),
      "Excluded - Filler" = c(0, 1, -1)
    ),
    adjust = "none"
  ) %>%
    tidy_emm() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # RQ2: pairwise country comparisons within condition and item group.
  country <- pairs(
    emmeans(
      model,
      ~ Country | Condition * Item_group
    ),
    adjust = "none"
  ) %>%
    tidy_emm() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # Standardized child-age slope.
  child_age_effect <- emtrends(
    model,
    ~ 1,
    var = "child_age_z"
  ) %>%
    tidy_emm() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # Adjusted child-sex means within curiosity and desirability.
  child_sex_means <- emmeans(
    model,
    ~ sex_binary | Condition
  ) %>%
    tidy_emm() %>%
    rename(sex = sex_binary) %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # Boy-versus-girl contrast within curiosity and desirability.
  child_sex_effects <- contrast(
    emmeans(
      model,
      ~ sex_binary | Condition
    ),
    method = "trt.vs.ctrl",
    ref = 1,
    adjust = "none"
  ) %>%
    tidy_emm() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # Adjusted education means within curiosity and desirability.
  responder_education_means <- emmeans(
    model,
    ~ responder_education | Condition
  ) %>%
    tidy_emm() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  # Education-category contrasts against Undergraduate / bachelor.
  responder_education_effects <- contrast(
    emmeans(
      model,
      ~ responder_education | Condition
    ),
    method = "trt.vs.ctrl",
    ref = 1,
    adjust = "none"
  ) %>%
    tidy_emm() %>%
    mutate(
      Questionnaire = questionnaire_name,
      .before = 1
    )
  
  list(
    model = model,
    omnibus = omnibus,
    means = means,
    planned = planned,
    country = country,
    child_age_effect = child_age_effect,
    child_sex_means = child_sex_means,
    child_sex_effects = child_sex_effects,
    responder_education_means = responder_education_means,
    responder_education_effects = responder_education_effects
  )
}
group_results <- map(c("ITCQ", "ECCQ"), fit_group_model)
names(group_results) <- c("ITCQ", "ECCQ")

walk(names(group_results), function(q) {
  saveRDS(group_results[[q]]$model,
          file.path(model_dir, paste0("item_set_model_", q, ".rds")))
})

group_omnibus <- map_dfr(group_results, "omnibus")
group_means <- map_dfr(group_results, "means")

child_age_effects <- map_dfr(group_results, "child_age_effect")
child_sex_means <- map_dfr(group_results, "child_sex_means")
child_sex_effects <- map_dfr(group_results, "child_sex_effects") %>%
  group_by(Questionnaire, Condition) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()
responder_education_means <- map_dfr(
  group_results,
  "responder_education_means"
)
responder_education_effects <- map_dfr(
  group_results,
  "responder_education_effects"
) %>%
  group_by(Questionnaire, Condition) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write_table(child_age_effects, "03_child_age_effect.csv")
write_table(
  child_sex_means,
  "03_child_sex_adjusted_means_by_condition.csv"
)
write_table(
  child_sex_effects,
  "03_child_sex_effects_by_condition.csv"
)
write_table(
  responder_education_means,
  "03_responder_education_adjusted_means_by_condition.csv"
)
write_table(
  responder_education_effects,
  "03_responder_education_effects_by_condition.csv"
)

group_planned <- map_dfr(group_results, "planned") %>%
  group_by(Questionnaire) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()
group_country <- map_dfr(group_results, "country") %>%
  group_by(Questionnaire) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write_table(group_omnibus, "03_RQ1_RQ2_item_set_omnibus_tests.csv")
write_table(group_means, "03_RQ1_RQ2_item_set_adjusted_means.csv")
write_table(group_planned, "03_RQ1_planned_item_set_contrasts.csv")
write_table(group_country, "03_RQ2_country_contrasts_by_item_set.csv")

# Figure 1: adjusted item-set means from the mixed models.
p_group <- ggplot(group_means,
                  aes(Item_group, emmean, colour = Country, group = Country)) +
  geom_hline(yintercept = 4, linetype = 3, colour = "grey60") +
  geom_line(position = position_dodge(width = .25), linewidth = .55) +
  geom_point(position = position_dodge(width = .25), size = 2) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(width = .25), width = .12) +
  facet_grid(Questionnaire ~ Condition) +
  scale_y_continuous(breaks = 1:7, limits = c(1, 7)) +
  scale_colour_manual(values = c(DE = "#0072B2", CN = "#D55E00", JP = "#009E73")) +
  labs(x = NULL, y = "Adjusted mean rating (95% CI)",
       title = "Ratings of final, excluded, and filler items",
       subtitle = "Horizontal line marks the neutral midpoint (4)") +
  theme(legend.position = "top")
save_figure(p_group, "Figure_1_item_set_adjusted_means", 9, 6.5)

# ---- 6. RQ2: CROSS-COUNTRY ITEM-PROFILE SIMILARITY --------------------------

item_raw_means <- analysis_data %>%
  group_by(Questionnaire, Country, Condition, Item_group, Item) %>%
  summarise(mean = mean(Response), SD = sd(Response), N = n(), .groups = "drop")
write_table(item_raw_means, "04_RQ2_raw_item_means_and_SDs.csv")

# Per-item within-group consensus complements the mean comparisons. A standard
# inter-rater ICC cannot be calculated meaningfully for one isolated item because
# there is no between-target variance to separate from rater/error variance.
# The indices below instead describe how concentrated the 1-7 responses are:
#
#   SD / IQR:              ordinary spread (lower = more consensus)
#   modal proportion:     share selecting the most common response
#   normalized consensus: 1 - normalized Shannon entropy (0 diffuse, 1 identical)
#   r_wg:                 observed variance relative to a uniform 1-7 null;
#                         expected uniform variance = (7^2 - 1) / 12 = 4.
#
# r_wg can be negative when disagreement exceeds the uniform-null variance. Do
# not silently truncate it to zero; negative values are diagnostically useful.
item_consensus <- analysis_data %>%
  group_by(Questionnaire, Country, Condition, Item_group, Item) %>%
  group_modify(~{
    counts <- tabulate(as.integer(.x$Response), nbins = 7)
    probabilities <- counts / sum(counts)
    nonzero <- probabilities[probabilities > 0]
    normalized_entropy <- -sum(nonzero * log(nonzero)) / log(7)
    observed_variance <- var(.x$Response)
    tibble(
      N = nrow(.x),
      mean = mean(.x$Response),
      SD = observed_variance^.5,
      IQR = IQR(.x$Response),
      modal_response = which.max(counts),
      modal_proportion = max(probabilities),
      normalized_entropy = normalized_entropy,
      normalized_consensus = 1 - normalized_entropy,
      r_wg_uniform = 1 - observed_variance / 4
    )
  }) %>%
  ungroup()
write_table(item_consensus, "04_RQ2_item_level_within_group_consensus.csv")

# Spearman correlations ask whether countries rank the behaviours similarly.
# They do not demonstrate equality of absolute levels, which is why the mixed-
# model country contrasts are also reported.
profile_correlations <- item_raw_means %>%
  select(Questionnaire, Condition, Country, Item, mean) %>%
  pivot_wider(names_from = Country, values_from = mean) %>%
  group_by(Questionnaire, Condition) %>%
  group_modify(~{
    pairs <- list(c("DE", "CN"), c("DE", "JP"), c("CN", "JP"))
    map_dfr(pairs, function(pair) {
      complete <- .x %>% select(all_of(pair)) %>% drop_na()
      test <- cor.test(complete[[pair[1]]], complete[[pair[2]]],
                       method = "spearman", exact = FALSE)
      tibble(country_1 = pair[1], country_2 = pair[2],
             items = nrow(complete),
             spearman_rho = unname(test$estimate), p.value = test$p.value)
    })
  }) %>%
  ungroup() %>%
  group_by(Questionnaire, Condition) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()
write_table(profile_correlations, "04_RQ2_item_profile_correlations.csv")

# Heatmap of item means. Items remain grouped for interpretation.
p_heat <- item_raw_means %>%
  mutate(Item_display = paste(Item_group, Item, sep = ": ")) %>%
  ggplot(aes(interaction(Country, Condition, sep = "\n"),
             fct_rev(factor(Item_display)), fill = mean)) +
  geom_tile(colour = "white", linewidth = .15) +
  facet_wrap(~Questionnaire, scales = "free_y", ncol = 2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 4, limits = c(1, 7)) +
  labs(x = NULL, y = NULL, fill = "Mean",
       title = "Item-level curiosity and desirability profiles") +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8),
        legend.position = "top")
save_figure(p_heat, "Figure_2_item_profile_heatmap", 10, 13)

# ---- FIGURE 2: ORDERED ITEM-PROFILE HEATMAP -------------------------------
#
# Layout:
#   ITCQ Final       | ECCQ Final
#   ITCQ Excluded    | ECCQ Excluded
#   ITCQ Filler      | ECCQ Filler
#
# Within every panel:
#   DE Curiosity | DE Desirability |
#   CN Curiosity | CN Desirability |
#   JP Curiosity | JP Desirability
#
# Thick vertical lines separate the countries. The two conditions belonging
# to the same country therefore appear as a visually connected pair.


# Define the intended item order separately for each questionnaire.
item_order <- bind_rows(
  
  tibble(
    Questionnaire = "ITCQ",
    Item = c(
      # Final items
      paste0("S", 1:10),
      "SR",
      paste0("Iv", 1:5),
      paste0("Ia", 1:7),
      
      # Excluded items
      paste0("Ex", 1:11),
      
      # Filler items
      paste0("F", 1:16)
    )
  ),
  
  tibble(
    Questionnaire = "ECCQ",
    Item = c(
      # Final items
      paste0("N", 1:9),
      paste0("In", 1:7),
      paste0("Ia", 1:9),
      
      # Excluded items
      paste0("Ex", 1:14),
      
      # Filler items
      paste0("F", 1:16)
    )
  )
) %>%
  group_by(Questionnaire) %>%
  mutate(Item_order = row_number()) %>%
  ungroup()


# Prepare the heatmap data.
heatmap_data <- item_raw_means %>%
  
  left_join(
    item_order,
    by = c("Questionnaire", "Item")
  ) %>%
  
  mutate(
    
    # Arrange the three item groups from top to bottom.
    Item_group = factor(
      Item_group,
      levels = c("Final", "Excluded", "Filler")
    ),
    
    # Place Curiosity and Desirability beside each other within country.
    Country_condition = factor(
      paste(Country, Condition, sep = "_"),
      
      levels = c(
        "DE_Curiosity",
        "DE_Desirability",
        "CN_Curiosity",
        "CN_Desirability",
        "JP_Curiosity",
        "JP_Desirability"
      ),
      
      labels = c(
        "DE\nCuriosity",
        "DE\nDesirability",
        "CN\nCuriosity",
        "CN\nDesirability",
        "JP\nCuriosity",
        "JP\nDesirability"
      )
    ),
    
    # Questionnaire is included in the internal item identifier because
    # item labels such as Ia1, Ex1, and F1 occur in both questionnaires.
    Item_key = paste(
      Questionnaire,
      Item,
      sep = "__"
    ),
    
    # Explicitly define the order of the six panels.
    Panel = factor(
      paste(
        Questionnaire,
        Item_group,
        sep = " — "
      ),
      
      levels = c(
        "ITCQ — Final",
        "ECCQ — Final",
        "ITCQ — Excluded",
        "ECCQ — Excluded",
        "ITCQ — Filler",
        "ECCQ — Filler"
      )
    )
  )


# Create the questionnaire-specific ordering of the y-axis items.
#
# The levels are reversed because ggplot places the first factor level at the
# bottom of the y-axis. Reversing them places the first item at the top.
item_levels <- item_order %>%
  mutate(
    Item_key = paste(
      Questionnaire,
      Item,
      sep = "__"
    )
  ) %>%
  pull(Item_key)

heatmap_data <- heatmap_data %>%
  mutate(
    Item_key = factor(
      Item_key,
      levels = rev(item_levels)
    )
  )


# Check that every item received a defined position.
unpositioned_items <- heatmap_data %>%
  filter(is.na(Item_order)) %>%
  distinct(Questionnaire, Item)

if (nrow(unpositioned_items) > 0) {
  
  write_table(
    unpositioned_items,
    "ERROR_Figure_2_unpositioned_items.csv"
  )
  
  stop(
    paste(
      "Some items could not be positioned in Figure 2.",
      "See ERROR_Figure_2_unpositioned_items.csv."
    )
  )
}


# Construct the heatmap.
p_heat <- ggplot(
  heatmap_data,
  aes(
    x = Country_condition,
    y = Item_key,
    fill = mean
  )
) +
  
  # Draw the individual heatmap cells.
  geom_tile(
    colour = "white",
    linewidth = 0.35
  ) +
  
  # Thick lines separate Germany, China, and Japan.
  # The two columns within each country remain visually grouped.
  geom_vline(
    xintercept = c(2.5, 4.5),
    colour = "black",
    linewidth = 1.2
  ) +
  
  # Print the raw item mean inside each cell.
  #
  # Numbers far from the midpoint are displayed in white because the
  # corresponding cells have a darker fill colour.
  geom_text(
    aes(
      label = sprintf("%.1f", mean),
      colour = abs(mean - 4) > 1.5
    ),
    size = 2.4
  ) +
  
  # Arrange the six panels in three rows and two columns.
  #
  # The filler panels are directly beside each other. Because both
  # questionnaires use F1-F16, corresponding filler rows will align.
  facet_wrap(
    ~Panel,
    ncol = 2,
    scales = "free_y"
  ) +
  
  # Blue indicates lower ratings, white the neutral midpoint, and red higher
  # ratings.
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 4,
    limits = c(1, 7),
    breaks = 1:7,
    name = "Mean rating"
  ) +
  
  # Remove the internal questionnaire prefix from the displayed item labels.
  scale_y_discrete(
    labels = function(x) {
      str_remove(x, "^[^_]+__")
    }
      ) +
        
        # Choose text colour according to the darkness of the heatmap cell.
        scale_colour_manual(
          values = c(
            `FALSE` = "grey15",
            `TRUE` = "white"
          ),
          guide = "none"
        ) +
        
        labs(
          x = NULL,
          y = NULL,
          
          title = "Item-level curiosity and desirability profiles",
          
          subtitle = paste(
            "Curiosity and desirability are adjacent within each country;",
            "cell values are raw item means"
          )
        ) +
        
        theme_bw(base_size = 10) +
        
        theme(
          
          # Remove background grid lines because the heatmap cells provide the grid.
          panel.grid = element_blank(),
          
          # Space between the six questionnaire/item-group panels.
          panel.spacing.x = grid::unit(1, "lines"),
          panel.spacing.y = grid::unit(1, "lines"),
          
          # Country and condition labels.
          axis.text.x = element_text(
            size = 8,
            lineheight = 0.9,
            colour = "black"
          ),
          
          # Item labels.
          axis.text.y = element_text(
            size = 7,
            colour = "black"
          ),
          
          axis.ticks = element_blank(),
          
          # Panel headings.
          strip.background = element_rect(
            fill = "#DCE6F1",
            colour = "grey50"
          ),
          
          strip.text = element_text(
            face = "bold",
            size = 10
          ),
          
          # Legend and title formatting.
          legend.position = "top",
          
          plot.title = element_text(
            face = "bold",
            size = 14
          ),
          
          plot.subtitle = element_text(
            size = 10
          )
        )
      
      
      # Save both high-resolution PNG and vector PDF versions.
      save_figure(
        p_heat,
        "Figure_2_ordered_item_profile_heatmap",
        width = 12,
        height = 15
      )
      
# ---- 7. RQ3: ITEM-LEVEL MODELS ---------------------------------------------
#
# Each item is modelled separately. The model adjusts for child age and tests
# Condition interactions with child sex and responder education while
# estimating all Country x Condition cells. Because each participant belongs to
# only one condition and contributes one response per item, an ordinary linear
# model is appropriate at this stage; participant random effects are unnecessary
# within a single item.

# Build and test arbitrary linear contrasts from an emmeans grid. This avoids
# relying on implicit row ordering when calculating difference-in-differences.
make_contrast_table <- function(emm, contrast_vectors) {
  contrast(emm, method = contrast_vectors, adjust = "none") %>% tidy_emm()
}

analyse_one_item <- function(d) {
  q <- as.character(unique(d$Questionnaire))
  item_name <- unique(d$Item)
  d <- droplevels(d)

  d$sex_binary <- set_treatment_reference(
    d$sex_binary,
    SEX_REFERENCE,
    "sex_binary"
  )
  d$responder_education <- set_treatment_reference(
    d$responder_education,
    EDUCATION_REFERENCE,
    "responder_education"
  )

  item_formula <- Response ~ Country * Condition +
    Condition * sex_binary +
    Condition * responder_education +
    child_age_z

  fit <- lm(item_formula, data = d)

  # HC3 covariance makes item-level uncertainty robust to unequal residual
  # variances and unequal country sample sizes.
  robust_vcov <- sandwich::vcovHC(fit,type="HC3")
  emm <- emmeans(fit, ~ Country * Condition, vcov.=robust_vcov)
  grid <- as.data.frame(emm) %>% select(Country, Condition)

  # Model-adjusted country-condition means and CIs.
  means <- tidy_emm(emm) %>%
    mutate(Questionnaire = q, Item = item_name, .before = 1)

  # Helper that returns the row of the emmeans grid for one cell.
  cell <- function(country, condition) {
    which(as.character(grid$Country) == country &
            as.character(grid$Condition) == condition)
  }

  countries <- c("DE", "CN", "JP")
  country_pairs <- list(c("DE", "CN"), c("DE", "JP"), c("CN", "JP"))
  conditions <- c("Curiosity", "Desirability")
  k <- nrow(grid)

  # Country differences separately for Curiosity and Desirability.
  country_vectors <- list()
  country_meta <- list()
  counter <- 0
  for (cond in conditions) {
    for (pair in country_pairs) {
      counter <- counter + 1
      v <- rep(0, k)
      v[cell(pair[1], cond)] <- 1
      v[cell(pair[2], cond)] <- -1
      nm <- paste(cond, pair[1], "-", pair[2])
      country_vectors[[nm]] <- v
      country_meta[[counter]] <- tibble(
        contrast = nm, Condition = cond,
        country_1 = pair[1], country_2 = pair[2]
      )
    }
  }
  country_tests <- make_contrast_table(emm, country_vectors) %>%
    left_join(bind_rows(country_meta), by = "contrast") %>%
    mutate(Questionnaire = q, Item = item_name, .before = 1)

  # Within-country desirability-minus-curiosity gaps.
  gap_vectors <- list()
  for (country in countries) {
    v <- rep(0, k)
    v[cell(country, "Desirability")] <- 1
    v[cell(country, "Curiosity")] <- -1
    gap_vectors[[country]] <- v
  }
  gaps <- make_contrast_table(emm, gap_vectors) %>%
    rename(Country = contrast) %>%
    mutate(Questionnaire = q, Item = item_name, .before = 1)

  # Difference-in-differences: compare the D-C gap between each country pair.
  gap_difference_vectors <- list()
  gap_meta <- list()
  counter <- 0
  for (pair in country_pairs) {
    counter <- counter + 1
    nm <- paste(pair[1], "gap -", pair[2], "gap")
    gap_difference_vectors[[nm]] <- gap_vectors[[pair[1]]] - gap_vectors[[pair[2]]]
    gap_meta[[counter]] <- tibble(contrast = nm,
                                  country_1 = pair[1], country_2 = pair[2])
  }
  gap_differences <- make_contrast_table(emm, gap_difference_vectors) %>%
    left_join(bind_rows(gap_meta), by = "contrast") %>%
    mutate(Questionnaire = q, Item = item_name, .before = 1)

  list(means = means, country_tests = country_tests,
       gaps = gaps, gap_differences = gap_differences)
}

item_results <- analysis_data %>%
  group_by(Questionnaire, Item) %>%
  group_split() %>%
  map(analyse_one_item)

item_means <- map_dfr(item_results, "means")

item_country_tests <- map_dfr(item_results, "country_tests") %>%
  group_by(Questionnaire, Condition) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

item_gaps <- map_dfr(item_results, "gaps") %>%
  group_by(Questionnaire) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

item_gap_differences <- map_dfr(item_results, "gap_differences") %>%
  group_by(Questionnaire) %>%
  mutate(p.value.FDR = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write_table(item_means, "05_RQ3_item_adjusted_means.csv")
write_table(item_country_tests, "05_RQ3_item_country_contrasts.csv")
write_table(item_gaps, "05_RQ3_desirability_minus_curiosity_gaps.csv")
write_table(item_gap_differences, "05_RQ3_country_differences_in_gaps.csv")

# ---- 8. RQ3 CLASSIFICATION OF THEORETICALLY INTERESTING ITEMS ---------------

# Curiosity is classified as:
#
# 1. Equivalent:
#    All three pairwise 95% CIs lie inside the equivalence margin.
#
# 2. Meaningfully different:
#    At least one pair has an absolute difference >= the margin and an
#    FDR-adjusted p-value below alpha.
#
# 3. Inconclusive:
#    Neither equivalence nor a meaningful difference was demonstrated.

curiosity_classification <- item_country_tests %>%
  filter(Condition == "Curiosity") %>%
  group_by(Questionnaire, Item) %>%
  summarise(
    max_abs_curiosity_difference = max(abs(estimate)),
    
    curiosity_equivalent = all(
      conf.low > -EQUIVALENCE_MARGIN &
        conf.high < EQUIVALENCE_MARGIN
    ),
    
    curiosity_divergent = any(
      abs(estimate) >= EQUIVALENCE_MARGIN &
        p.value.FDR < ALPHA
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    curiosity_status = case_when(
      curiosity_equivalent ~
        "Equivalent across countries",
      
      curiosity_divergent ~
        "Meaningful country difference",
      
      TRUE ~
        "Inconclusive"
    )
  )


# Apply the same distinction to desirability.

desirability_classification <- item_country_tests %>%
  filter(Condition == "Desirability") %>%
  group_by(Questionnaire, Item) %>%
  summarise(
    max_abs_desirability_difference = max(abs(estimate)),
    
    desirability_equivalent = all(
      conf.low > -EQUIVALENCE_MARGIN &
        conf.high < EQUIVALENCE_MARGIN
    ),
    
    desirability_divergent = any(
      abs(estimate) >= EQUIVALENCE_MARGIN &
        p.value.FDR < ALPHA
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    desirability_status = case_when(
      desirability_equivalent ~
        "Equivalent across countries",
      
      desirability_divergent ~
        "Meaningful country difference",
      
      TRUE ~
        "Inconclusive"
    )
  )


# Retain the country differences in desirability-curiosity gaps.

gap_classification <- item_gap_differences %>%
  group_by(Questionnaire, Item) %>%
  summarise(
    gap_divergent = any(
      abs(estimate) >= EQUIVALENCE_MARGIN &
        p.value.FDR < ALPHA
    ),
    
    max_abs_gap_difference = max(abs(estimate)),
    
    .groups = "drop"
  )


# Combine all three pieces of information.

item_classification <- curiosity_classification %>%
  
  full_join(
    desirability_classification,
    by = c("Questionnaire", "Item")
  ) %>%
  
  full_join(
    gap_classification,
    by = c("Questionnaire", "Item")
  ) %>%
  
  left_join(
    analysis_data %>%
      distinct(Questionnaire, Item, Item_group),
    
    by = c("Questionnaire", "Item")
  ) %>%
  
  mutate(
    theoretical_pattern = case_when(
      
      curiosity_equivalent & desirability_divergent ~
        "Curiosity equivalent + desirability differs",
      
      curiosity_divergent & desirability_divergent ~
        "Curiosity and desirability both differ",
      
      curiosity_divergent & !desirability_divergent ~
        "Curiosity differs; desirability does not clearly differ",
      
      !curiosity_divergent &
        curiosity_status == "Inconclusive" &
        desirability_divergent ~
        "Desirability differs; curiosity inconclusive",
      
      curiosity_equivalent & gap_divergent ~
        "Curiosity equivalent + country-specific D-C gap",
      
      gap_divergent ~
        "Country-specific D-C gap; other evidence inconclusive",
      
      curiosity_equivalent & desirability_equivalent ~
        "Curiosity and desirability equivalent",
      
      TRUE ~
        "No clear meaningful cultural difference"
    )
  ) %>%
  
  arrange(
    Questionnaire,
    Item_group,
    desc(desirability_divergent),
    desc(max_abs_desirability_difference)
  )



# Determine whether the adjusted curiosity mean is above, around, or below
# the neutral midpoint of 4 in each country.

curiosity_midpoint <- item_means %>%
  filter(Condition == "Curiosity") %>%
  mutate(
    curiosity_midpoint_status = case_when(
      conf.low > 4 ~ "Clearly above 4",
      conf.high < 4 ~ "Clearly below 4",
      TRUE ~ "CI includes 4"
    )
  ) %>%
  select(
    Questionnaire,
    Item,
    Country,
    emmean,
    conf.low,
    conf.high,
    curiosity_midpoint_status
  )

write_table(
  curiosity_midpoint,
  "06_RQ3_curiosity_vs_midpoint_by_country.csv"
)


# Save the corrected item classification.

write_table(
  item_classification,
  "06_RQ3_item_classification.csv"
)
# ---- 9. RQ3 FIGURES ---------------------------------------------------------

# Cultural-divergence plot: upper-left items have small curiosity differences
# but large desirability differences and are especially relevant to the paper.
divergence_plot_data <- item_classification %>%
  group_by(Questionnaire) %>%
  mutate(
    label_priority = curiosity_equivalent & desirability_divergent,
    label_rank = rank(-max_abs_desirability_difference, ties.method = "first"),
    plot_label = if_else(label_priority |
                           label_rank <= MAX_LABELS_PER_DIVERGENCE_PLOT,
                         Item, NA_character_)
  ) %>%
  ungroup()

p_divergence <- ggplot(
  divergence_plot_data,
  aes(max_abs_curiosity_difference, max_abs_desirability_difference,
      colour = Item_group)
) +
  geom_vline(xintercept = EQUIVALENCE_MARGIN, linetype = 2, colour = "grey50") +
  geom_hline(yintercept = EQUIVALENCE_MARGIN, linetype = 2, colour = "grey50") +
  geom_point(aes(shape = desirability_divergent), size = 2.4, alpha = .85) +
  geom_text_repel(aes(label = plot_label), size = 3, max.overlaps = Inf,
                  box.padding = .25, seed = 2026) +
  facet_wrap(~Questionnaire, scales = "free") +
  scale_colour_manual(values = c(Final = "#0072B2", Excluded = "#CC79A7",
                                 Filler = "#999999")) +
  labs(
    x = "Largest absolute country difference in Curiosity",
    y = "Largest absolute country difference in Desirability",
    colour = "Item group", shape = "FDR-significant\nand meaningful",
    title = "Where is cultural divergence specific to desirability?",
    subtitle = paste0("Dashed lines show the ", EQUIVALENCE_MARGIN,
                      "-point practical-difference threshold")
  ) +
  theme(legend.position = "top")
save_figure(p_divergence, "Figure_3_cultural_divergence_by_item", 10, 6)

# Curiosity-desirability maps. Points below the diagonal are rated as less
# desirable than curiosity-expressing. Labels highlight the largest absolute
# gaps within each questionnaire-country panel.
map_data <- item_means %>%
  select(Questionnaire, Item, Country, Condition, emmean) %>%
  pivot_wider(names_from = Condition, values_from = emmean) %>%
  left_join(analysis_data %>% distinct(Questionnaire, Item, Item_group),
            by = c("Questionnaire", "Item")) %>%
  mutate(gap = Desirability - Curiosity) %>%
  group_by(Questionnaire, Country) %>%
  mutate(label_rank = rank(-abs(gap), ties.method = "first"),
         plot_label = if_else(label_rank <= 8, Item, NA_character_)) %>%
  ungroup()

p_map <- ggplot(map_data, aes(Curiosity, Desirability, colour = Item_group)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_hline(yintercept = 4, linetype = 3, colour = "grey75") +
  geom_vline(xintercept = 4, linetype = 3, colour = "grey75") +
  geom_point(alpha = .75, size = 1.8) +
  geom_text_repel(aes(label = plot_label), size = 2.5, max.overlaps = Inf,
                  seed = 2026, box.padding = .2) +
  facet_grid(Questionnaire ~ Country) +
  coord_equal(xlim = c(1, 7), ylim = c(1, 7)) +
  scale_x_continuous(breaks = 1:7) +
  scale_y_continuous(breaks = 1:7) +
  scale_colour_manual(values = c(Final = "#0072B2", Excluded = "#CC79A7",
                                 Filler = "#999999")) +
  labs(title = "Item-level relationship between curiosity and desirability",
       subtitle = "Items below the diagonal are rated as less desirable than curiosity-expressing",
       colour = "Item group") +
  theme(legend.position = "top")
save_figure(p_map, "Figure_4_curiosity_desirability_item_maps", 11, 7.5)

# ---- 10. ANALYSIS SETTINGS AND SESSION INFORMATION --------------------------

settings <- tribble(
  ~setting, ~value,
  "J1-J3 included", "No",
  "Questionnaires pooled", "No; ITCQ and ECCQ modelled separately",
  "Curiosity and Desirability jointly modelled", "Yes",
  "Response scale", "1-7; reverse-worded items scored as 8 - response",
  "Reverse-scored ITCQ items", "SR, Ex2, Ex10",
  "Reverse-scored ECCQ items", "N3, N5, Ex10, Ex13",
  "Equivalence/practical-difference margin", as.character(EQUIVALENCE_MARGIN),
  "Alpha", as.character(ALPHA),
  "Multiple-testing correction", "Benjamini-Hochberg FDR",
  "Exclude ages outside nominal range", as.character(EXCLUDE_OUT_OF_RANGE_AGES),
  "Group model", "Linear mixed model: participant and item random intercepts; child age covariate; Condition x child sex and Condition x responder education",
  "Child-sex reference", SEX_REFERENCE,
  "Responder-education reference", EDUCATION_REFERENCE,
  "Collapsed education category", "Lower secondary or below + Upper secondary = Upper secondary or below",
  "Item model", "Linear model with child age covariate, Condition x child sex, and Condition x responder education; HC3 robust covariance"
)
write_table(settings, "99_analysis_settings.csv")

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
capture.output(
  list(
    group_model_summaries = map(group_results, ~summary(.x$model)),
    singular_models = map_lgl(group_results, ~isSingular(.x$model, tol = 1e-4))
  ),
  file = file.path(log_dir, "group_model_summaries.txt")
)

message("Analysis complete. Results are in: ", normalizePath(OUTPUT_DIR))


# Interactive-only View()/view() calls removed so this file can be run
# non-interactively with Rscript. Tables/figures/models already exported above.

# ---- EXPORT MAIN REPORTING TABLES TO EXCEL ---------------------------------

#install.packages("writexl", type = "binary")

dir.create(
  "reporting_results",
  showWarnings = FALSE
)

reporting_tables <- list(
  
  # Sample information
  "Sample description" =
    sample_summary,
  
  # Raw means and SDs for the three item sets
  "Set descriptives" =
    set_descriptives,
  
  # RQ1: final items compared with excluded and filler items
  "RQ1 set comparisons" =
    group_planned %>%
    filter(
      contrast %in% c(
        "Final - Excluded",
        "Final - Filler"
      )
    ) %>%
    arrange(
      Questionnaire,
      Condition,
      Country,
      contrast
    ),
  
  # RQ2: model-adjusted means
  "RQ2 adjusted means" =
    group_means %>%
    arrange(
      Questionnaire,
      Condition,
      Item_group,
      Country
    ),
  
  # Model covariates
  "Child age effect" =
    child_age_effects,

  "Child sex means" =
    child_sex_means %>%
    arrange(Questionnaire, Condition, sex),

  "Child sex effects" =
    child_sex_effects %>%
    arrange(Questionnaire, Condition, p.value.FDR),

  "Education means" =
    responder_education_means %>%
    arrange(Questionnaire, Condition, responder_education),

  "Education effects" =
    responder_education_effects %>%
    arrange(Questionnaire, Condition, p.value.FDR),

  # RQ2: pairwise country comparisons
  "RQ2 country contrasts" =
    group_country %>%
    arrange(
      Questionnaire,
      Condition,
      Item_group,
      p.value.FDR
    ),
  
  # RQ2: similarity of item rankings between countries
  "RQ2 profile correlations" =
    profile_correlations,
  
  # RQ3: classified theoretically interesting items
  "RQ3 item classification" =
    item_classification %>%
    arrange(
      Questionnaire,
      theoretical_pattern,
      desc(max_abs_desirability_difference)
    ),
  
  "RQ3 item means" =
    item_means %>%
    arrange(
      Questionnaire,
      Item,
      Condition,
      Country
    ),
  
  "RQ3 curiosity contrasts" =
    item_country_tests %>%
    filter(Condition == "Curiosity") %>%
    arrange(
      Questionnaire,
      Item,
      p.value.FDR
    ),
  
  "RQ3 desirability contrasts" =
    item_country_tests %>%
    filter(Condition == "Desirability") %>%
    arrange(
      Questionnaire,
      Item,
      p.value.FDR
    ),
  
  # RQ3: significant and meaningful country differences
  "RQ3 meaningful contrasts" =
    item_country_tests %>%
    filter(
      p.value.FDR < ALPHA,
      abs(estimate) >= EQUIVALENCE_MARGIN
    ) %>%
    arrange(
      Questionnaire,
      Condition,
      desc(abs(estimate))
    ),
  
  # RQ3: differences between countries in desirability-curiosity gaps
  "RQ3 gap differences" =
    item_gap_differences %>%
    filter(
      p.value.FDR < ALPHA,
      abs(estimate) >= EQUIVALENCE_MARGIN
    ) %>%
    arrange(
      Questionnaire,
      desc(abs(estimate))
    ),
  
  # Whether curiosity ratings are above or below the midpoint
  "Curiosity midpoint" =
    curiosity_midpoint %>%
    arrange(
      Questionnaire,
      Item,
      Country
    )
)

writexl::write_xlsx(
  reporting_tables,
  path = file.path(
    "reporting_results",
    "Cross_cultural_main_results.xlsx"
  )
)

message(
  "Main reporting workbook written to: ",
  file.path(
    "reporting_results",
    "Cross_cultural_main_results.xlsx"
  )
)











