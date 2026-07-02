
# Required columns:
#   HRS:
#     id, age, sex, education, bmi, depression, self_rated_health,
#     arthritis, pain_any, walk_blocks_difficulty, walk_room_difficulty,
#     walking_time_mean_sec, falls_any, adl_any, iadl_any, adl_count,
#     iadl_count, cognition_z, memory_rating
#
#   Auxiliary cohorts:
#     id, wave, age, sex, education, depression, self_rated_health,
#     arthritis, functional_any_limitation, adl_any, iadl_any,
#     functional_score, cognition_z
#
# Notes:
#   pain_any is optional in the auxiliary cohorts. If absent or entirely missing,
#   pain models are skipped for that cohort.
################################################################################

rm(list = ls())

options(stringsAsFactors = FALSE)

required_packages <- c(
  "data.table", "dplyr", "tidyr", "purrr", "stringr", "broom",
  "ggplot2", "metafor", "patchwork", "scales", "readr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install required packages before running Stage 7: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(broom)
  library(ggplot2)
  library(metafor)
  library(patchwork)
  library(scales)
  library(readr)
})

input_dir <- Sys.getenv(
  "KOA_WM_STAGE7_INPUT_DIR",
  unset = file.path("data", "stage7_aging_cohorts")
)
output_dir <- Sys.getenv(
  "KOA_WM_STAGE7_RESULTS_DIR",
  unset = file.path("results", "stage7_older_adult_phenotype_context")
)
figure_dir <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

theme_stage7 <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(color = "#66737f"),
      axis.text = element_text(color = "#26313b"),
      plot.title = element_text(face = "bold", color = "#18232e", size = base_size + 2),
      plot.subtitle = element_text(color = "#66737f", size = base_size - 1),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.major.x = element_line(color = "#e7edf2", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      plot.margin = margin(8, 10, 8, 10)
    )
}

stop_if_missing <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

read_stage7 <- function(file_name, required_cols, label) {
  path <- file.path(input_dir, file_name)
  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }
  dat <- fread(path)
  stop_if_missing(dat, required_cols, label)
  dat
}

format_ci <- function(est, low, high, digits = 2) {
  paste0(
    formatC(est, format = "f", digits = digits),
    " (",
    formatC(low, format = "f", digits = digits),
    " to ",
    formatC(high, format = "f", digits = digits),
    ")"
  )
}

empty_model_result <- function() {
  tibble(
    outcome = character(),
    exposure = character(),
    outcome_type = character(),
    n = integer(),
    beta = numeric(),
    se = numeric(),
    OR = numeric(),
    CI_low = numeric(),
    CI_high = numeric(),
    p_value = numeric(),
    effect_metric = character(),
    estimate_95_CI = character()
  )
}

empty_model_grid_result <- function() {
  empty_model_result() %>%
    mutate(model = character())
}

fit_model <- function(data, outcome, exposure, covariates, outcome_type) {
  model_vars <- c(outcome, exposure, covariates)
  model_data <- data[, model_vars, with = FALSE] %>% tidyr::drop_na()
  if (nrow(model_data) < 50 || length(unique(model_data[[exposure]])) < 2) {
    return(empty_model_result())
  }

  form <- as.formula(paste(outcome, "~", paste(c(exposure, covariates), collapse = " + ")))

  if (outcome_type == "binary") {
    if (length(unique(model_data[[outcome]])) < 2) return(empty_model_result())
    fit <- suppressWarnings(
      tryCatch(
        glm(form, data = model_data, family = binomial()),
        error = function(e) NULL
      )
    )
    if (is.null(fit)) return(empty_model_result())
    out <- broom::tidy(fit) %>% filter(term == exposure)
    if (nrow(out) == 0) return(empty_model_result())
    beta <- out$estimate[1]
    se <- out$std.error[1]
    tibble(
      outcome = outcome,
      exposure = exposure,
      outcome_type = outcome_type,
      n = nrow(model_data),
      beta = beta,
      se = se,
      OR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      p_value = out$p.value[1],
      effect_metric = "OR",
      estimate_95_CI = format_ci(exp(beta), exp(beta - 1.96 * se), exp(beta + 1.96 * se))
    )
  } else {
    fit <- tryCatch(lm(form, data = model_data), error = function(e) NULL)
    if (is.null(fit)) return(empty_model_result())
    out <- broom::tidy(fit) %>% filter(term == exposure)
    if (nrow(out) == 0) return(empty_model_result())
    beta <- out$estimate[1]
    se <- out$std.error[1]
    tibble(
      outcome = outcome,
      exposure = exposure,
      outcome_type = outcome_type,
      n = nrow(model_data),
      beta = beta,
      se = se,
      OR = NA_real_,
      CI_low = beta - 1.96 * se,
      CI_high = beta + 1.96 * se,
      p_value = out$p.value[1],
      effect_metric = "beta",
      estimate_95_CI = format_ci(beta, beta - 1.96 * se, beta + 1.96 * se, digits = 3)
    )
  }
}

fit_binary_or_continuous_grid <- function(data, outcomes, exposures, covariates, model_label) {
  out <- purrr::map_dfr(names(outcomes), function(outcome) {
    if (!outcome %in% names(data) || all(is.na(data[[outcome]]))) return(empty_model_grid_result())
    purrr::map_dfr(exposures, function(exposure) {
      if (!exposure %in% names(data) || all(is.na(data[[exposure]]))) return(empty_model_grid_result())
      fit_model(
        data = data,
        outcome = outcome,
        exposure = exposure,
        covariates = covariates,
        outcome_type = outcomes[[outcome]]
      ) %>%
        mutate(model = model_label)
    })
  })
  if (!"model" %in% names(out)) {
    out$model <- character()
  }
  out
}

random_effects_meta <- function(model_results, analysis_label) {
  if (nrow(model_results) == 0) {
    return(tibble())
  }
  model_results %>%
    filter(effect_metric == "OR", is.finite(beta), is.finite(se), se > 0) %>%
    group_by(exposure, outcome) %>%
    group_modify(function(dat, key) {
      if (nrow(dat) < 2) return(tibble())
      fit <- metafor::rma(yi = beta, sei = se, data = dat, method = "DL")
      tibble(
        analysis = analysis_label,
        k = nrow(dat),
        logOR = as.numeric(fit$b[1]),
        se = fit$se,
        OR = exp(as.numeric(fit$b[1])),
        CI_low = exp(fit$ci.lb),
        CI_high = exp(fit$ci.ub),
        tau2 = fit$tau2,
        I2_percent = fit$I2,
        Q = fit$QE,
        cohorts = paste(dat$cohort, collapse = ";")
      )
    }) %>%
    ungroup()
}

leave_one_cohort_out <- function(model_results, analysis_label, target_outcome = "adl_any") {
  if (nrow(model_results) == 0) {
    return(tibble())
  }
  model_results %>%
    filter(effect_metric == "OR", outcome == target_outcome, is.finite(beta), is.finite(se), se > 0) %>%
    group_by(exposure) %>%
    group_modify(function(dat, key) {
      if (n_distinct(dat$cohort) < 3) return(tibble())
      purrr::map_dfr(sort(unique(dat$cohort)), function(excluded) {
        tmp <- dat %>% filter(cohort != excluded)
        fit <- metafor::rma(yi = beta, sei = se, data = tmp, method = "DL")
        tibble(
          analysis = analysis_label,
          outcome = target_outcome,
          excluded_cohort = excluded,
          k = nrow(tmp),
          OR = exp(as.numeric(fit$b[1])),
          CI_low = exp(fit$ci.lb),
          CI_high = exp(fit$ci.ub),
          I2_percent = fit$I2
        )
      })
    }) %>%
    ungroup()
}

################################################################################
# 1. Load harmonized cohort data
################################################################################

hrs_required <- c(
  "id", "age", "sex", "education", "bmi", "depression", "self_rated_health",
  "arthritis", "pain_any", "walk_blocks_difficulty", "walk_room_difficulty",
  "walking_time_mean_sec", "falls_any", "adl_any", "iadl_any",
  "adl_count", "iadl_count", "cognition_z", "memory_rating"
)

aux_required <- c(
  "id", "wave", "age", "sex", "education", "depression", "self_rated_health",
  "arthritis", "functional_any_limitation", "adl_any", "iadl_any",
  "functional_score", "cognition_z"
)

hrs <- read_stage7("HRS_2018_stage7.csv", hrs_required, "HRS")

aux_files <- c(
  ELSA = "ELSA_stage7_long.csv",
  MHAS = "MHAS_stage7_long.csv",
  CHARLS = "CHARLS_stage7_long.csv",
  SHARE = "SHARE_stage7_long.csv"
)

aux_data <- imap(aux_files, function(file_name, cohort_name) {
  dat <- read_stage7(file_name, aux_required, cohort_name)
  if (!"pain_any" %in% names(dat)) dat$pain_any <- NA_real_
  dat %>% mutate(cohort = cohort_name)
}) %>%
  bind_rows()

################################################################################
# 2. Variable coverage
################################################################################

coverage_hrs <- tibble(variable = names(hrs)) %>%
  mutate(
    cohort = "HRS",
    non_missing = map_int(variable, ~ sum(!is.na(hrs[[.x]]))),
    mean = map_dbl(variable, ~ suppressWarnings(mean(as.numeric(hrs[[.x]]), na.rm = TRUE)))
  )

coverage_aux <- aux_data %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    arthritis_non_missing = sum(!is.na(arthritis)),
    pain_non_missing = sum(!is.na(pain_any)),
    adl_non_missing = sum(!is.na(adl_any)),
    iadl_non_missing = sum(!is.na(iadl_any)),
    functional_non_missing = sum(!is.na(functional_any_limitation)),
    cognition_non_missing = sum(!is.na(cognition_z)),
    .groups = "drop"
  )

write_csv(coverage_hrs, file.path(table_dir, "Stage7_HRS_variable_coverage.csv"))
write_csv(coverage_aux, file.path(table_dir, "Stage7_auxiliary_variable_coverage.csv"))

################################################################################
# 3. HRS primary phenotype models and pain attenuation
################################################################################

hrs_outcomes <- c(
  walk_blocks_difficulty = "binary",
  walk_room_difficulty = "binary",
  walking_time_mean_sec = "continuous",
  falls_any = "binary",
  adl_any = "binary",
  iadl_any = "binary",
  adl_count = "continuous",
  iadl_count = "continuous",
  cognition_z = "continuous",
  memory_rating = "continuous"
)

hrs_covariates <- c("age", "sex", "education", "bmi", "depression", "self_rated_health")

hrs_base <- fit_binary_or_continuous_grid(
  data = hrs,
  outcomes = hrs_outcomes,
  exposures = c("arthritis", "pain_any"),
  covariates = hrs_covariates,
  model_label = "base"
)

hrs_plus_pain <- map_dfr(names(hrs_outcomes), function(outcome) {
  fit_model(
    data = hrs,
    outcome = outcome,
    exposure = "arthritis",
    covariates = c(hrs_covariates, "pain_any"),
    outcome_type = hrs_outcomes[[outcome]]
  ) %>%
    mutate(model = "plus_pain")
})

hrs_models <- bind_rows(hrs_base, hrs_plus_pain) %>%
  mutate(
    cohort = "HRS",
    year = 2018,
    estimate_95_CI = if_else(
      effect_metric == "OR",
      format_ci(OR, CI_low, CI_high),
      estimate_95_CI
    )
  ) %>%
  select(cohort, year, exposure, outcome, outcome_type, model, effect_metric, n,
         beta, se, OR, CI_low, CI_high, estimate_95_CI, p_value)

write_csv(hrs_models, file.path(table_dir, "Stage7_HRS_primary_models.csv"))

attenuation_outcomes <- c(
  "walk_blocks_difficulty", "walk_room_difficulty", "falls_any", "adl_any", "iadl_any"
)

hrs_attenuation <- hrs_models %>%
  filter(exposure == "arthritis", outcome %in% attenuation_outcomes, effect_metric == "OR") %>%
  select(outcome, model, n, OR, CI_low, CI_high, beta) %>%
  pivot_wider(
    names_from = model,
    values_from = c(n, OR, CI_low, CI_high, beta),
    names_sep = "_"
  ) %>%
  mutate(
    attenuation_percent_logOR =
      100 * (beta_base - beta_plus_pain) / abs(beta_base),
    crossed_null_after_pain = CI_low_plus_pain <= 1 & CI_high_plus_pain >= 1,
    CI_full = format_ci(OR_base, CI_low_base, CI_high_base),
    CI_plus_pain = format_ci(OR_plus_pain, CI_low_plus_pain, CI_high_plus_pain)
  ) %>%
  transmute(
    outcome,
    n_full = n_base,
    OR_full = OR_base,
    CI_full,
    n_plus_pain,
    OR_plus_pain,
    CI_plus_pain,
    attenuation_percent_logOR,
    crossed_null_after_pain
  )

write_csv(hrs_attenuation, file.path(table_dir, "Stage7_HRS_pain_attenuation_gradient.csv"))

################################################################################
# 4. Auxiliary cross sectional and lagged cohort models
################################################################################

aux_outcomes <- c(
  functional_any_limitation = "binary",
  adl_any = "binary",
  iadl_any = "binary",
  functional_score = "continuous",
  cognition_z = "continuous"
)

aux_covariates <- c("age", "sex", "education", "depression", "self_rated_health", "wave")

aux_cross <- aux_data %>%
  group_split(cohort) %>%
  map_dfr(function(dat) {
    cohort_name <- unique(dat$cohort)
    exposures <- c("arthritis", "pain_any")
    if (all(is.na(dat$pain_any))) exposures <- "arthritis"
    fit_binary_or_continuous_grid(
      data = as.data.table(dat),
      outcomes = aux_outcomes,
      exposures = exposures,
      covariates = aux_covariates,
      model_label = "cross_sectional"
    ) %>%
      mutate(cohort = cohort_name)
  }) %>%
  select(cohort, exposure, outcome, outcome_type, model, effect_metric, n,
         beta, se, OR, CI_low, CI_high, estimate_95_CI, p_value)

write_csv(aux_cross, file.path(table_dir, "Stage7_auxiliary_cross_sectional_models.csv"))

make_lagged_data <- function(dat) {
  dat %>%
    arrange(id, wave) %>%
    group_by(id) %>%
    mutate(
      next_functional_any_limitation = lead(functional_any_limitation),
      next_adl_any = lead(adl_any),
      next_iadl_any = lead(iadl_any),
      next_functional_score = lead(functional_score),
      next_cognition_z = lead(cognition_z)
    ) %>%
    ungroup()
}

lagged_outcomes <- c(
  next_functional_any_limitation = "binary",
  next_adl_any = "binary",
  next_iadl_any = "binary",
  next_functional_score = "continuous",
  next_cognition_z = "continuous"
)

lagged_baseline_map <- c(
  next_functional_any_limitation = "functional_any_limitation",
  next_adl_any = "adl_any",
  next_iadl_any = "iadl_any",
  next_functional_score = "functional_score",
  next_cognition_z = "cognition_z"
)

aux_lagged <- aux_data %>%
  group_split(cohort) %>%
  map_dfr(function(dat) {
    cohort_name <- unique(dat$cohort)
    lagged <- make_lagged_data(dat)
    exposures <- c("arthritis", "pain_any")
    if (all(is.na(lagged$pain_any))) exposures <- "arthritis"

    map_dfr(names(lagged_outcomes), function(outcome) {
      map_dfr(exposures, function(exposure) {
        baseline_outcome <- lagged_baseline_map[[outcome]]
        fit_model(
          data = as.data.table(lagged),
          outcome = outcome,
          exposure = exposure,
          covariates = c(aux_covariates, baseline_outcome),
          outcome_type = lagged_outcomes[[outcome]]
        )
      })
    }) %>%
      mutate(
        cohort = cohort_name,
        model = "lagged_followup_adjusted_for_baseline_outcome",
        outcome = str_replace(outcome, "^next_", "")
      )
  }) %>%
  select(cohort, exposure, outcome, outcome_type, model, effect_metric, n,
         beta, se, OR, CI_low, CI_high, estimate_95_CI, p_value)

write_csv(aux_lagged, file.path(table_dir, "Stage7_auxiliary_lagged_models.csv"))

################################################################################
# 5. Random effects meta analyses, outcome profile, leave one cohort out
################################################################################

meta_cross <- random_effects_meta(aux_cross, "cross_sectional")
meta_lagged <- random_effects_meta(aux_lagged, "lagged")
meta_all <- bind_rows(meta_cross, meta_lagged)

write_csv(meta_all, file.path(table_dir, "Stage7_cross_cohort_random_effects_meta.csv"))

loo_cross <- leave_one_cohort_out(aux_cross, "cross_sectional")
loo_lagged <- leave_one_cohort_out(aux_lagged, "lagged")
loo_all <- bind_rows(loo_cross, loo_lagged)

write_csv(loo_all, file.path(table_dir, "Stage7_leave_one_cohort_out_meta.csv"))

outcome_profile <- bind_rows(aux_cross, aux_lagged) %>%
  mutate(
    endpoint_family = case_when(
      outcome %in% c("functional_any_limitation", "adl_any", "iadl_any", "functional_score") ~ "Functional tests",
      outcome == "cognition_z" ~ "Cognition tests",
      TRUE ~ "Other"
    ),
    adverse_direction = case_when(
      endpoint_family == "Functional tests" & beta > 0 ~ TRUE,
      endpoint_family == "Cognition tests" & beta < 0 ~ TRUE,
      TRUE ~ FALSE
    ),
    nominal = p_value < 0.05
  ) %>%
  group_by(model, exposure, endpoint_family) %>%
  summarise(
    n_tests = n(),
    pct_adverse_direction = 100 * mean(adverse_direction, na.rm = TRUE),
    pct_nominal = 100 * mean(nominal, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(outcome_profile, file.path(table_dir, "Stage7_outcome_specificity_summary.csv"))

################################################################################
# 6. Figure 5 and supplementary Stage 7 visuals
################################################################################

outcome_labels <- c(
  walk_blocks_difficulty = "Walking several blocks",
  walk_room_difficulty = "Walking across room",
  falls_any = "Falls",
  adl_any = "ADL limitation",
  iadl_any = "IADL limitation",
  functional_any_limitation = "Any functional limitation"
)

p_a <- hrs_attenuation %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_labels),
    outcome_label = factor(outcome_label, levels = rev(recode(attenuation_outcomes, !!!outcome_labels)))
  ) %>%
  ggplot(aes(x = attenuation_percent_logOR, y = outcome_label)) +
  geom_col(fill = "#2f74a8", width = 0.55) +
  geom_text(aes(label = paste0(round(attenuation_percent_logOR), "%")),
            hjust = -0.15, size = 3.1, color = "#26313b") +
  coord_cartesian(xlim = c(0, max(hrs_attenuation$attenuation_percent_logOR, na.rm = TRUE) * 1.18)) +
  labs(
    title = "A  Pain attenuation in HRS",
    subtitle = "Reduction in arthritis log odds ratio after adding pain.",
    x = "Attenuation of log odds ratio, %",
    y = NULL
  ) +
  theme_stage7()

hrs_pair <- hrs_models %>%
  filter(exposure == "arthritis", outcome %in% attenuation_outcomes, model %in% c("base", "plus_pain")) %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_labels),
    model_label = recode(model, base = "Base model", plus_pain = "Pain adjusted")
  )

p_b <- hrs_pair %>%
  ggplot(aes(x = OR, y = outcome_label, color = model_label)) +
  geom_vline(xintercept = 1, color = "#b7c3cf", linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.75,
                 position = position_dodge(width = 0.45)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.45)) +
  scale_color_manual(values = c("Base model" = "#2f74a8", "Pain adjusted" = "#cf6d36")) +
  labs(
    title = "B  HRS estimates before and after pain adjustment",
    subtitle = "Paired odds ratios for arthritis models.",
    x = "Odds ratio",
    y = NULL
  ) +
  theme_stage7()

meta_plot_data <- meta_all %>%
  filter(outcome %in% c("adl_any", "functional_any_limitation", "iadl_any")) %>%
  mutate(
    outcome_label = recode(outcome, !!!outcome_labels),
    exposure_label = recode(exposure, arthritis = "Arthritis", pain_any = "Pain"),
    analysis_label = recode(analysis, cross_sectional = "Cross cohort phenotype context", lagged = "Lagged robustness"),
    label = paste0(formatC(OR, format = "f", digits = 2), " (",
                   formatC(CI_low, format = "f", digits = 2), " to ",
                   formatC(CI_high, format = "f", digits = 2), "); I2 ",
                   round(I2_percent), "%")
  )

p_c <- meta_plot_data %>%
  filter(analysis == "cross_sectional") %>%
  mutate(row_label = paste(exposure_label, outcome_label, sep = ": ")) %>%
  ggplot(aes(x = OR, y = row_label, color = exposure_label)) +
  geom_vline(xintercept = 1, color = "#b7c3cf", linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.75) +
  geom_point(size = 2.4) +
  geom_text(aes(label = label), x = Inf, hjust = 1, size = 2.8, color = "#66737f") +
  scale_color_manual(values = c("Arthritis" = "#2a9d8f", "Pain" = "#7a5fb3")) +
  labs(
    title = "C  Cross cohort phenotype context",
    subtitle = "Random effects estimates across harmonized aging cohorts.",
    x = "Random effects odds ratio",
    y = NULL
  ) +
  theme_stage7()

p_d <- meta_plot_data %>%
  filter(analysis == "lagged") %>%
  mutate(row_label = paste(exposure_label, outcome_label, sep = ": ")) %>%
  ggplot(aes(x = OR, y = row_label, color = exposure_label)) +
  geom_vline(xintercept = 1, color = "#b7c3cf", linetype = "dashed") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.75) +
  geom_point(size = 2.4) +
  geom_text(aes(label = label), x = Inf, hjust = 1, size = 2.8, color = "#66737f") +
  scale_color_manual(values = c("Arthritis" = "#2a9d8f", "Pain" = "#7a5fb3")) +
  labs(
    title = "D  Lagged robustness",
    subtitle = "Baseline exposure with next wave functional outcomes.",
    x = "Random effects odds ratio",
    y = NULL
  ) +
  theme_stage7()

p_e <- outcome_profile %>%
  mutate(
    model_label = recode(model,
                         cross_sectional = "Cross sectional",
                         lagged_followup_adjusted_for_baseline_outcome = "Lagged"),
    exposure_label = recode(exposure, arthritis = "Arthritis", pain_any = "Pain"),
    row_label = paste(model_label, exposure_label, sep = ", ")
  ) %>%
  ggplot(aes(x = pct_adverse_direction, y = row_label, color = endpoint_family)) +
  geom_segment(aes(x = 0, xend = pct_adverse_direction, yend = row_label), linewidth = 1.2) +
  geom_point(size = 2.6) +
  geom_text(aes(label = paste0(round(pct_adverse_direction), "%")),
            hjust = -0.15, size = 2.8) +
  scale_color_manual(values = c("Functional tests" = "#4c956c", "Cognition tests" = "#bc4749")) +
  coord_cartesian(xlim = c(0, 105)) +
  labs(
    title = "E  Outcome profile",
    subtitle = "Functional endpoints and cognition sensitivity outcomes.",
    x = "Tests with adverse direction, %",
    y = NULL
  ) +
  theme_stage7()

p_f <- loo_all %>%
  filter(outcome == "adl_any") %>%
  group_by(analysis, exposure) %>%
  summarise(
    OR_min = min(OR, na.rm = TRUE),
    OR_max = max(OR, na.rm = TRUE),
    OR_mid = median(OR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    exposure_label = recode(exposure, arthritis = "arthritis", pain_any = "pain"),
    analysis_label = recode(analysis, cross_sectional = "Cross sectional", lagged = "Lagged"),
    row_label = paste(analysis_label, exposure_label)
  ) %>%
  ggplot(aes(y = row_label, color = exposure_label)) +
  geom_vline(xintercept = 1, color = "#b7c3cf", linetype = "dashed") +
  geom_segment(aes(x = OR_min, xend = OR_max, yend = row_label), linewidth = 4, lineend = "round") +
  geom_point(aes(x = OR_mid), color = "white", size = 2.1) +
  geom_text(aes(x = OR_max, label = paste0(formatC(OR_min, digits = 2, format = "f"),
                                           " to ",
                                           formatC(OR_max, digits = 2, format = "f"))),
            hjust = -0.12, color = "#66737f", size = 2.8) +
  scale_color_manual(values = c("arthritis" = "#2a9d8f", "pain" = "#7a5fb3")) +
  coord_cartesian(xlim = c(0.8, max(loo_all$CI_high, na.rm = TRUE) * 1.08)) +
  labs(
    title = "F  Leave one cohort out ADL stability",
    subtitle = "Random effects estimates after excluding one cohort at a time.",
    x = "Odds ratio range after exclusion",
    y = NULL
  ) +
  theme_stage7()

figure5 <- (p_a | p_b | p_c) / (p_d | p_e | p_f) +
  plot_annotation(
    caption = "HRS, Health and Retirement Study; ADL, activities of daily living; IADL, instrumental activities of daily living."
  )

ggsave(file.path(figure_dir, "Figure_5_older_adult_phenotype_context.pdf"),
       figure5, width = 18, height = 12, units = "in")
ggsave(file.path(figure_dir, "Figure_5_older_adult_phenotype_context.png"),
       figure5, width = 18, height = 12, units = "in", dpi = 400)
if (requireNamespace("svglite", quietly = TRUE) &&
    requireNamespace("systemfonts", quietly = TRUE)) {
  ggsave(file.path(figure_dir, "Figure_5_older_adult_phenotype_context.svg"),
         figure5, width = 18, height = 12, units = "in")
} else {
  message("Skipping SVG export because svglite/systemfonts is not installed.")
}

################################################################################
# 7. Session information
################################################################################

writeLines(capture.output(sessionInfo()), file.path(output_dir, "Stage7_sessionInfo.txt"))

message("Stage 7 older adult phenotype context analysis completed.")
message("Tables: ", normalizePath(table_dir, winslash = "/"))
message("Figures: ", normalizePath(figure_dir, winslash = "/"))
