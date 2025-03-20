# To run this script, you will need to download the data from IPUMs.

# Download the data as a single csv file, ensuring that it contains all these variables:
# YEAR, AGE, BPL, EDUC, EMPSTAT, GQ, HISPAN, INCTOT, OCC, PERWT,
# RACE, RELATE, SERIAL, SEX, STATEICP, WKSWORK2

# Copy the unzipped csv file named data.csv to the same location as this script.

# You can include as many years as you want, as long as they are listed here:

Years <- c(1950, 1960, 1970, 1980, 1990, 2000:2023)

# OCC is the default variable used for fixed effects, but you can change it here:

Occupation_var <- "OCC"

# When you are ready, you can run the script by pressing "Source" in RStudio.

# Text files with the results will be saved in the same folder as this script.

# Set the working directory to the script's location
setwd(getSrcDirectory(function(dummy) {dummy}))
cat("\014")

# Load required packages
library(tidyverse)
library(fixest)
library(modelsummary)
library(Hmisc)
library(knitr)


# Function to clean income variable
clean_income <- function(df, year) {
  # Handle top and bottom codes
  handled_income <- handle_income_codes(df, year)
  
  df <- handled_income$data %>%
    mutate(
      income_clean = case_when(
        INCTOT == 9999999 ~ NA_real_,
        INCTOT == 9999998 ~ NA_real_,
        is_topcode == TRUE ~ INCTOT_adjusted,  # Use adjusted value for topcoded observations
        TRUE ~ as.numeric(INCTOT)  # Use original INCTOT for non-topcoded observations
      ),
      weight = PERWT/100
    )
  
  return(df)
}

# Function to handle top-coded and bottom-coded incomes for a single year
handle_income_codes <- function(data, year) {
  # First filter out missing values
  data <- data %>%
    mutate(
      INCTOT = case_when(
        INCTOT == 9999999 ~ NA_real_,
        INCTOT == 9999998 ~ NA_real_,
        INCTOT < 0 ~ 0,
        TRUE ~ INCTOT
      )
    )
  
  # Define top-code value for each year
  topcode_values <- list(
    "1950" = 10000,
    "1960" = 25000,
    "1970" = 50000,
    "1980" = 75000,
    "1990" = 400000,
    "2000" = 999998
  )
  
  year_str <- as.character(year)
  
  # For years after 2000, return data without topcode adjustment
  if (!year_str %in% names(topcode_values)) {
    return(list(
      data = data %>% mutate(
        is_topcode = FALSE,
        INCTOT_adjusted = INCTOT
      ),
      alpha = NA,
      mean_above_topcode = NA,
      n_topcoded = 0,
      pct_topcoded = 0
    ))
  }
  
  # For years with topcoding, proceed with existing logic
  topcode_value <- topcode_values[[year_str]]
  
  # Identify topcoded observations
  data <- data %>%
    mutate(
      is_topcode = !is.na(INCTOT) & INCTOT >= topcode_value
    )
  
  valid_obs <- sum(!is.na(data$INCTOT))
  n_topcoded <- sum(data$is_topcode, na.rm = TRUE)
  pct_topcoded <- (n_topcoded / valid_obs) * 100
  
  # Add weight for Pareto calculations
  data <- data %>%
    mutate(
      weight = PERWT/100
    )
  
  # Get non-topcoded incomes for Pareto estimation
  income_data <- data %>%
    filter(!is_topcode, 
           !is.na(INCTOT),
           INCTOT > 0)
  
  # If we have topcoded observations, proceed with Pareto adjustment
  if (n_topcoded > 0) {
    # Calculate percentiles for the top 20% of non-topcoded incomes
    cutoff_80th <- wtd.quantile(income_data$INCTOT, income_data$weight, probs = 0.8)
    
    pareto_data <- income_data %>%
      filter(INCTOT >= cutoff_80th) %>%
      arrange(INCTOT) %>%
      mutate(
        log_income = log(INCTOT),
        # Calculate survival function (1 - CDF)
        survival = 1 - (cumsum(weight) / sum(weight))
      ) %>%
      filter(survival > 0) %>%  # Remove any zero survival probabilities
      mutate(
        log_survival = log(survival)
      )
    
    # Estimate Pareto parameter using regression
    pareto_model <- lm(log_survival ~ log_income, data = pareto_data)
    alpha <- -coef(pareto_model)["log_income"]
    
    # Generate random draws from Pareto distribution for topcoded values
    set.seed(505)  # For reproducibility
    u <- runif(n_topcoded)
    adjusted_values <- topcode_value * (1 - u)^(-1/alpha)
    mean_above_topcode <- mean(adjusted_values)
    
    # Assign adjusted values to topcoded observations
    data <- data %>%
      mutate(
        INCTOT_adjusted = case_when(
          is_topcode ~ adjusted_values[row_number() %% n_topcoded + 1],
          TRUE ~ INCTOT
        )
      )
  } else {
    # If no topcoded observations, use original values
    data <- data %>%
      mutate(INCTOT_adjusted = INCTOT)
    alpha <- NA
    mean_above_topcode <- NA
  }
  
  return(list(
    data = data,
    alpha = alpha,
    mean_above_topcode = mean_above_topcode,
    n_topcoded = n_topcoded,
    pct_topcoded = pct_topcoded
  ))
}

# Function to create broad demographic group summary
create_broad_demo_summary <- function(data) {
  # Remove any NA incomes
  data <- data %>% 
    filter(!is.na(income_clean))
  
  # Retrieve statistics from attributes
  national_mean <- attr(data, "national_mean")
  national_median <- attr(data, "national_median")
  
  # Calculate total income across all groups
  total_income <- sum(data$income_clean * data$weight, na.rm = TRUE)
  
  # Add broad group indicators
  data <- data %>%
    mutate(
      # Broad demographic groups
      white_men = ifelse(RACE == 1 & SEX == 1 & HISPAN == 0, 1, 0),
      white_women = ifelse(RACE == 1 & SEX == 2 & HISPAN == 0, 1, 0),
      nonwhite_men = ifelse((RACE != 1 | HISPAN > 0) & SEX == 1, 1, 0),
      nonwhite_women = ifelse((RACE != 1 | HISPAN > 0) & SEX == 2, 1, 0)
    )
  
  # List of broad groups to analyze
  broad_groups <- c("white_men", "white_women", "nonwhite_men", "nonwhite_women")
  
  # Calculate statistics for each broad group
  broad_income_table <- lapply(broad_groups, function(group_var) {
    group_data <- data[data[[group_var]] == 1, ]
    tibble(
      group = case_when(
        group_var == "white_men" ~ "White men",
        group_var == "white_women" ~ "White women",
        group_var == "nonwhite_men" ~ "Non-white men",
        group_var == "nonwhite_women" ~ "Non-white women"
      ),
      n_obs = nrow(group_data),
      pct_pop = nrow(group_data) / nrow(data) * 100,
      mean_income = weighted.mean(group_data$income_clean, group_data$weight, na.rm = TRUE),
      median_income = wtd.quantile(group_data$income_clean, group_data$weight, probs = 0.5, na.rm = TRUE),
      income_share = sum(group_data$income_clean * group_data$weight, na.rm = TRUE) / total_income * 100,
      pct_of_national_mean = (mean_income / national_mean) * 100,
      pct_of_national_median = (median_income / national_median) * 100
    )
  }) %>%
    bind_rows() %>%
    arrange(desc(mean_income))
  
  # Calculate characteristics for each broad group
  broad_characteristics_table <- lapply(broad_groups, function(group_var) {
    group_data <- data[data[[group_var]] == 1, ]
    tibble(
      group = case_when(
        group_var == "white_men" ~ "White men",
        group_var == "white_women" ~ "White women",
        group_var == "nonwhite_men" ~ "Non-white men",
        group_var == "nonwhite_women" ~ "Non-white women"
      ),
      pct_hs_complete = weighted.mean(group_data$high_school_completion, group_data$weight, na.rm = TRUE) * 100,
      labor_force_part = weighted.mean(group_data$in_labor_force, group_data$weight, na.rm = TRUE) * 100,
      unemployment_rate = weighted.mean(group_data$unemployment_rate, group_data$weight, na.rm = TRUE) * 100,
      pct_foreign_born = weighted.mean(group_data$foreign_born, group_data$weight, na.rm = TRUE) * 100,
      avg_weeks_worked = weighted.mean(group_data$weeks_worked, group_data$weight, na.rm = TRUE),
      pct_college = weighted.mean(group_data$college_graduate, group_data$weight, na.rm = TRUE) * 100,
      pct_inst = weighted.mean(group_data$institutionalized, group_data$weight, na.rm = TRUE) * 100
    )
  }) %>%
    bind_rows() %>%
    right_join(broad_income_table %>% select(group), by = "group")
  
  return(list(
    income_stats = broad_income_table,
    characteristics = broad_characteristics_table
  ))
}

# Prepare_data function
prepare_data <- function(df, year) {
  # First filter for valid income and age
  prepared_df <- df %>%
    filter(
      AGE >= 18,
      !is.na(INCTOT),
      INCTOT != 9999999,
      INCTOT != 9999998
    )
  
  # Then create demographic indicators
  prepared_df <- if(year >= 2000) {
    df %>%
      mutate(
        # White non-Hispanic
        white_male = case_when(
          RACE == 1 & SEX == 1 & HISPAN == 0 ~ 1,
          TRUE ~ 0
        ),
        white_female = case_when(
          RACE == 1 & SEX == 2 & HISPAN == 0 ~ 1,
          TRUE ~ 0
        ),
        
        # Black (regardless of Hispanic status)
        black_male = case_when(
          RACE == 2 & SEX == 1 ~ 1,
          TRUE ~ 0
        ),
        black_female = case_when(
          RACE == 2 & SEX == 2 ~ 1,
          TRUE ~ 0
        ),
        
        # Hispanic (any race except Black)
        hispanic_male = case_when(
          HISPAN > 0 & SEX == 1 & RACE != 2 ~ 1,
          TRUE ~ 0
        ),
        hispanic_female = case_when(
          HISPAN > 0 & SEX == 2 & RACE != 2 ~ 1,
          TRUE ~ 0
        ),
        
        # Other: not Black, not Hispanic, and not white
        other_male = case_when(
          SEX == 1 & RACE != 2 & HISPAN == 0 & RACE != 1 ~ 1,
          TRUE ~ 0
        ),
        other_female = case_when(
          SEX == 2 & RACE != 2 & HISPAN == 0 & RACE != 1 ~ 1,
          TRUE ~ 0
        ),
        
        # Group assignment
        group = case_when(
          white_male == 1 ~ "White male",
          white_female == 1 ~ "White female",
          black_male == 1 ~ "Black male",
          black_female == 1 ~ "Black female",
          hispanic_male == 1 ~ "Hispanic male",
          hispanic_female == 1 ~ "Hispanic female",
          other_male == 1 ~ "Other male",
          other_female == 1 ~ "Other female"
        )
      )
  } else {
    df %>%
      mutate(
        # White non-Hispanic
        white_male = case_when(
          RACE == 1 & SEX == 1 & (is.na(HISPAN) | HISPAN == 0) ~ 1,
          TRUE ~ 0
        ),
        white_female = case_when(
          RACE == 1 & SEX == 2 & (is.na(HISPAN) | HISPAN == 0) ~ 1,
          TRUE ~ 0
        ),
        
        # Black (regardless of Hispanic status)
        black_male = case_when(
          RACE == 2 & SEX == 1 ~ 1,
          TRUE ~ 0
        ),
        black_female = case_when(
          RACE == 2 & SEX == 2 ~ 1,
          TRUE ~ 0
        ),
        
        # Hispanic (any race except Black)
        hispanic_male = case_when(
          !is.na(HISPAN) & HISPAN > 0 & SEX == 1 & RACE != 2 ~ 1,
          TRUE ~ 0
        ),
        hispanic_female = case_when(
          !is.na(HISPAN) & HISPAN > 0 & SEX == 2 & RACE != 2 ~ 1,
          TRUE ~ 0
        ),
        
        # Other: not Black, not Hispanic, and not white
        other_male = case_when(
          SEX == 1 & RACE != 2 & (is.na(HISPAN) | HISPAN == 0) & RACE != 1 ~ 1,
          TRUE ~ 0
        ),
        other_female = case_when(
          SEX == 2 & RACE != 2 & (is.na(HISPAN) | HISPAN == 0) & RACE != 1 ~ 1,
          TRUE ~ 0
        ),
        
        # Group assignment
        group = case_when(
          white_male == 1 ~ "White male",
          white_female == 1 ~ "White female",
          black_male == 1 ~ "Black male",
          black_female == 1 ~ "Black female",
          hispanic_male == 1 ~ "Hispanic male",
          hispanic_female == 1 ~ "Hispanic female",
          other_male == 1 ~ "Other male",
          other_female == 1 ~ "Other female"
        )
      )
  }
  
  # Clean income and add core variables
  prepared_df <- prepared_df %>%
    clean_income(year) %>%
    filter(!is.na(income_clean)) %>%
    mutate(
      # Education variables
      high_school_completion = as.numeric(EDUC >= 6),
      college_graduate = as.numeric(EDUC >= 10),
      
      # Labor force variables
      foreign_born = as.numeric(BPL >= 100),
      weeks_worked = case_when(
        WKSWORK2 == 0 ~ 0,
        WKSWORK2 == 1 ~ 12,
        WKSWORK2 == 2 ~ 20,
        WKSWORK2 == 3 ~ 33,
        WKSWORK2 == 4 ~ 43.5,
        WKSWORK2 == 5 ~ 48.5,
        WKSWORK2 == 6 ~ 51
      ),
      in_labor_force = as.numeric(EMPSTAT == 1 | EMPSTAT == 2),
      unemployed = as.numeric(EMPSTAT == 2),
      unemployment_rate = ifelse(EMPSTAT == 1 | EMPSTAT == 2, 
                                 as.numeric(EMPSTAT == 2), NA),
      institutionalized = as.numeric(GQ == 3)
    )
  
  # Calculate white male dropout statistics
  wm_dropout_stats <- prepared_df %>%
    filter(
      AGE >= 18,
      white_male == 1,
      EDUC < 6
    ) %>%
    summarise(
      mean_income = weighted.mean(income_clean, weight, na.rm = TRUE),
      median_income = wtd.quantile(income_clean, weight, probs = 0.5, na.rm = TRUE)
    )
  
  # Calculate national mean and median
  national_stats <- prepared_df %>%
    filter(AGE >= 18) %>%
    summarise(
      mean_income = weighted.mean(income_clean, weight, na.rm = TRUE),
      median_income = wtd.quantile(income_clean, weight, probs = 0.5, na.rm = TRUE)
    )
  
  # Add normalized income and fixed effects
  prepared_df <- prepared_df %>%
    filter(AGE >= 18) %>%
    mutate(
      income_normalized = (income_clean / wm_dropout_stats$mean_income) * 100,
      
      # Fixed effects
      age = factor(AGE),
      state = factor(STATEICP),
      occupation = factor(get(Occupation_var))
    )
  
  # Add attributes for later use
  attr(prepared_df, "wm_dropout_mean") <- wm_dropout_stats$mean_income
  attr(prepared_df, "wm_dropout_median") <- wm_dropout_stats$median_income
  attr(prepared_df, "national_mean") <- national_stats$mean_income
  attr(prepared_df, "national_median") <- national_stats$median_income
  
  return(prepared_df)
}

# Modified function to create demographic summary table with national comparisons
create_demo_summary <- function(data) {
  # Remove any NA incomes
  data <- data %>% 
    filter(!is.na(income_clean))
  
  # Retrieve statistics from attributes
  national_mean <- attr(data, "national_mean")
  national_median <- attr(data, "national_median")
  
  # Calculate total income across all groups
  total_income <- sum(data$income_clean * data$weight, na.rm = TRUE)
  
  # Add aggregate group indicators
  data <- data %>%
    mutate(
      # Race groups
      black = ifelse(RACE == 2, 1, 0),
      hispanic = ifelse(HISPAN > 0, 1, 0),
      white = ifelse(RACE == 1 & HISPAN == 0, 1, 0),
      other = ifelse(RACE != 1 & RACE != 2 & HISPAN == 0, 1, 0),
      
      # Gender groups
      male = ifelse(SEX == 1, 1, 0),
      female = ifelse(SEX == 2, 1, 0)
    )
  
  # List of all groups to analyze
  group_vars <- c("white_male", "black_male", "black_female", "white_female", 
                  "hispanic_male", "hispanic_female", "other_male", "other_female",
                  "black", "hispanic", "white", "other", "male", "female")
  
  # Calculate statistics for each group
  income_table <- lapply(group_vars, function(group_var) {
    if(group_var %in% names(data)) {
      group_data <- data[data[[group_var]] == 1, ]
      tibble(
        group = str_replace_all(group_var, "_", " ") %>% str_to_title(),
        n_obs = nrow(group_data),
        pct_pop = nrow(group_data) / nrow(data) * 100,
        mean_income = weighted.mean(group_data$income_clean, group_data$weight, na.rm = TRUE),
        median_income = wtd.quantile(group_data$income_clean, group_data$weight, probs = 0.5, na.rm = TRUE),
        income_share = sum(group_data$income_clean * group_data$weight, na.rm = TRUE) / total_income * 100,
        pct_of_national_mean = (mean_income / national_mean) * 100,
        pct_of_national_median = (median_income / national_median) * 100
      )
    }
  }) %>%
    bind_rows() %>%
    arrange(desc(mean_income))
  
  # Calculate characteristics for each group
  characteristics_table <- lapply(group_vars, function(group_var) {
    if(group_var %in% names(data)) {
      group_data <- data[data[[group_var]] == 1, ]
      tibble(
        group = str_replace_all(group_var, "_", " ") %>% str_to_title(),
        pct_hs_complete = weighted.mean(group_data$high_school_completion, group_data$weight, na.rm = TRUE) * 100,
        labor_force_part = weighted.mean(group_data$in_labor_force, group_data$weight, na.rm = TRUE) * 100,
        unemployment_rate = weighted.mean(group_data$unemployment_rate, group_data$weight, na.rm = TRUE) * 100,
        pct_foreign_born = weighted.mean(group_data$foreign_born, group_data$weight, na.rm = TRUE) * 100,
        avg_weeks_worked = weighted.mean(group_data$weeks_worked, group_data$weight, na.rm = TRUE),
        pct_college = weighted.mean(group_data$college_graduate, group_data$weight, na.rm = TRUE) * 100,
        pct_inst = weighted.mean(group_data$institutionalized, group_data$weight, na.rm = TRUE) * 100
      )
    }
  }) %>%
    bind_rows() %>%
    right_join(income_table %>% select(group), by = "group")
  
  return(list(
    income_stats = income_table,
    characteristics = characteristics_table
  ))
}

# Modified function to create education summary with column order switched
create_education_summary <- function(data) {
  # Add aggregate group indicators
  data <- data %>%
    mutate(
      # Race groups
      black = ifelse(RACE == 2, 1, 0),
      hispanic = ifelse(HISPAN > 0, 1, 0),
      white = ifelse(RACE == 1 & HISPAN == 0, 1, 0),
      other = ifelse(RACE != 1 & RACE != 2 & HISPAN == 0, 1, 0),
      
      # Gender groups
      male = ifelse(SEX == 1, 1, 0),
      female = ifelse(SEX == 2, 1, 0)
    )
  
  # Access statistics from attributes
  wm_dropout_mean <- attr(data, "wm_dropout_mean")
  wm_dropout_median <- attr(data, "wm_dropout_median")
  national_mean <- attr(data, "national_mean")
  national_median <- attr(data, "national_median")
  
  # List of all groups to analyze
  group_vars <- c("white_male", "black_male", "black_female", "white_female", 
                  "hispanic_male", "hispanic_female", "other_male", "other_female",
                  "black", "hispanic", "white", "other", "male", "female")
  
  # Calculate education statistics for all groups
  results <- lapply(group_vars, function(group_var) {
    if(group_var %in% names(data)) {
      data %>%
        filter(!!sym(group_var) == 1) %>%
        mutate(
          education = case_when(
            college_graduate == 1 ~ "College Graduate",
            high_school_completion == 1 ~ "High School Graduate",
            TRUE ~ "Dropout"
          )
        ) %>%
        group_by(education) %>%
        summarise(
          group = str_replace_all(group_var, "_", " ") %>% str_to_title(),
          n_obs = n(),
          mean_income = weighted.mean(income_clean, weight, na.rm = TRUE),
          median_income = wtd.quantile(income_clean, weight, probs = 0.5, na.rm = TRUE),
          pct_inst = weighted.mean(institutionalized, weight, na.rm = TRUE) * 100,
          pct_of_national_mean = (mean_income / national_mean) * 100,
          pct_of_national_median = (median_income / national_median) * 100,
          .groups = 'drop'
        )
    }
  }) %>%
    bind_rows() %>%
    # Switch education and group columns by reordering them
    select(education, group, everything()) %>%
    arrange(education, group)
  
  return(results)
}

# Function to run models
run_models <- function(data) {
  list(
    # Model 1: Base model - captures all sources of inequality
    m1 = feols(income_normalized ~ black_female + black_male + white_female + 
                 hispanic_female + hispanic_male + other_female + other_male +
                 foreign_born + institutionalized +
                 high_school_completion + college_graduate +
                 black_female:high_school_completion + black_male:high_school_completion +
                 white_female:high_school_completion + hispanic_female:high_school_completion +
                 hispanic_male:high_school_completion + other_female:high_school_completion +
                 other_male:high_school_completion +
                 black_female:college_graduate + black_male:college_graduate +
                 white_female:college_graduate + hispanic_female:college_graduate +
                 hispanic_male:college_graduate + other_female:college_graduate +
                 other_male:college_graduate | 
                 age + state,
               weights = ~weight, cluster = ~state, data = data),
    
    # Model 2: Add labor force controls - inequality net of labor force participation
    m2 = feols(income_normalized ~ black_female + black_male + white_female + 
                 hispanic_female + hispanic_male + other_female + other_male +
                 foreign_born + weeks_worked + institutionalized +
                 high_school_completion + college_graduate +
                 black_female:high_school_completion + black_male:high_school_completion +
                 white_female:high_school_completion + hispanic_female:high_school_completion +
                 hispanic_male:high_school_completion + other_female:high_school_completion +
                 other_male:high_school_completion +
                 black_female:college_graduate + black_male:college_graduate +
                 white_female:college_graduate + hispanic_female:college_graduate +
                 hispanic_male:college_graduate + other_female:college_graduate +
                 other_male:college_graduate | 
                 age + state,
               weights = ~weight, cluster = ~state, data = data),
    
    # Model 3: Add occupation controls - inequality net of both labor force and occupation
    m3 = feols(income_normalized ~ black_female + black_male + white_female + 
                 hispanic_female + hispanic_male + other_female + other_male +
                 foreign_born + weeks_worked + institutionalized +
                 high_school_completion + college_graduate +
                 black_female:high_school_completion + black_male:high_school_completion +
                 white_female:high_school_completion + hispanic_female:high_school_completion +
                 hispanic_male:high_school_completion + other_female:high_school_completion +
                 other_male:high_school_completion +
                 black_female:college_graduate + black_male:college_graduate +
                 white_female:college_graduate + hispanic_female:college_graduate +
                 hispanic_male:college_graduate + other_female:college_graduate +
                 other_male:college_graduate | 
                 age + state + occupation,
               weights = ~weight, cluster = ~state, data = data)
  )
}

# Function to print model results
print_model_results <- function(model, model_num) {
  cat(sprintf("\nModel %d:\n", model_num))
  cat("Number of observations:", model$nobs, "\n")
  cat("R-squared:", round(r2(model)["r2"], 4), "\n")
  cat("Within R-squared:", round(r2(model)["wr2"], 4), "\n\n")
  
  coef_table <- data.frame(
    Estimate = coef(model),
    `Std. Error` = se(model),
    `t-value` = coef(model) / se(model),
    `P-value` = pvalue(model)
  )
  print(round(coef_table, 4))
  cat("\n")
}

# Function to calculate group effects
calculate_group_effects <- function(models, group_var, hs = TRUE, college = TRUE) {
  # Helper to extract coefficients
  get_coef <- function(model, var) {
    coef_val <- coef(model)[var]
    return(if(is.na(coef_val)) 0 else coef_val)
  }
  
  # Special handling for white males
  if(group_var == "white_male") {
    # For white males, base effect is 0 (reference group)
    base <- 0
    # Only include education main effects, no interactions
    hs_main <- if(hs) get_coef(models$m1, "high_school_completion") else 0
    hs_int <- 0  # No interaction term for reference group
    college_main <- if(college) get_coef(models$m1, "college_graduate") else 0
    college_int <- 0  # No interaction term for reference group
  } else {
    # Regular calculation for other groups
    base <- get_coef(models$m1, group_var)
    hs_main <- if(hs) get_coef(models$m1, "high_school_completion") else 0
    hs_int <- if(hs) get_coef(models$m1, paste0(group_var, ":high_school_completion")) else 0
    college_main <- if(college) get_coef(models$m1, "college_graduate") else 0
    college_int <- if(college) get_coef(models$m1, paste0(group_var, ":college_graduate")) else 0
  }
  
  # Calculate total for each model
  calculate_model_total <- function(model) {
    if(group_var == "white_male") {
      # For white males, only include education main effects
      total <- if(hs) get_coef(model, "high_school_completion") else 0
      total <- total + if(college) get_coef(model, "college_graduate") else 0
    } else {
      # Regular calculation for other groups
      base <- get_coef(model, group_var)
      hs_main <- if(hs) get_coef(model, "high_school_completion") else 0
      hs_int <- if(hs) get_coef(model, paste0(group_var, ":high_school_completion")) else 0
      college_main <- if(college) get_coef(model, "college_graduate") else 0
      college_int <- if(college) get_coef(model, paste0(group_var, ":college_graduate")) else 0
      total <- base + hs_main + hs_int + college_main + college_int
    }
    return(total)
  }
  
  # Get totals for each model
  m1_total <- calculate_model_total(models$m1)
  m2_total <- calculate_model_total(models$m2)
  m3_total <- calculate_model_total(models$m3)
  
  # Calculate decomposition
  labor_force <- m1_total - m2_total
  occupation <- m2_total - m3_total
  within_occ <- m3_total
  
  # Create results data frame
  results <- data.frame(
    component = c("Base Effect", "High School Main", "HS Interaction",
                  "College Main", "College Interaction", "Total Effect",
                  "Labor Force Effect", "Occupation Effect", "Within-Occ Effect"),
    value = c(base, hs_main, hs_int, college_main, college_int,
              m1_total, labor_force, occupation, within_occ)
  )
  
  return(results)
}

# Function to calculate weighted Gini coefficient
calculate_weighted_gini <- function(x, w) {
  # Remove NA values
  valid <- !is.na(x) & !is.na(w)
  x <- x[valid]
  w <- w[valid]
  
  # Only calculate if we have valid data
  if(length(x) == 0) return(NA)
  
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  p <- cumsum(w) / sum(w)
  p <- c(0, p)
  L <- cumsum(w * x) / sum(w * x)
  L <- c(0, L)
  gini <- 1 - sum((p[-1] - p[-length(p)]) * (L[-1] + L[-length(L)]))
  return(gini)
}

# Simplified function to calculate Gini coefficients for demographic groups
calculate_gini_by_group <- function(data) {
  # List of groups to analyze
  groups <- c("all", "white_male", "black_male", "black_female", "white_female", 
              "hispanic_male", "hispanic_female", "other_male", "other_female")
  
  # Calculate Gini coefficient for each group
  gini_summary <- lapply(groups, function(group) {
    # Filter data for the group
    group_data <- if(group == "all") {
      data
    } else {
      data %>% filter(!!sym(tolower(group)) == 1)
    }
    
    # Calculate Gini coefficient
    gini <- calculate_weighted_gini(group_data$income_clean, group_data$weight)
    
    # Return the summary
    data.frame(
      Group = gsub("_", " ", toupper(group)),
      Gini = round(gini, 3),
      N = nrow(group_data)
    )
  })
  
  # Combine into a single data frame
  do.call(rbind, gini_summary)
}

# Function for analyzing demographic composition of all deciles
analyze_decile_composition <- function(data) {
  # Add small random noise to break ties
  set.seed(505)
  data <- data %>%
    mutate(
      income_rand = income_clean + runif(n(), -0.0001, 0.0001)
    )
  
  # Calculate decile thresholds using randomized income
  decile_thresholds <- c(
    wtd.quantile(data$income_rand, 
                 weights = data$weight, 
                 probs = seq(0.1, 0.9, by = 0.1), 
                 na.rm = TRUE)
  )
  
  # Create decile indicators
  data <- data %>%
    mutate(
      income_decile = case_when(
        income_rand <= decile_thresholds[1] ~ "0-10",
        income_rand <= decile_thresholds[2] ~ "10-20",
        income_rand <= decile_thresholds[3] ~ "20-30",
        income_rand <= decile_thresholds[4] ~ "30-40",
        income_rand <= decile_thresholds[5] ~ "40-50",
        income_rand <= decile_thresholds[6] ~ "50-60",
        income_rand <= decile_thresholds[7] ~ "60-70",
        income_rand <= decile_thresholds[8] ~ "70-80",
        income_rand <= decile_thresholds[9] ~ "80-90",
        TRUE ~ "90-100"
      )
    )
  
  # Calculate demographic shares within each decile
  decile_composition <- data %>%
    group_by(income_decile) %>%
    summarise(
      `White male` = weighted.mean(white_male, weight) * 100,
      `Black male` = weighted.mean(black_male, weight) * 100,
      `Black female` = weighted.mean(black_female, weight) * 100,
      `White female` = weighted.mean(white_female, weight) * 100,
      `Hispanic male` = weighted.mean(hispanic_male, weight) * 100,
      `Hispanic female` = weighted.mean(hispanic_female, weight) * 100,
      `Other male` = weighted.mean(other_male, weight) * 100,
      `Other female` = weighted.mean(other_female, weight) * 100,
      n = n(),
      .groups = 'drop'
    )
  
  # Calculate overall shares
  overall_shares <- data %>%
    summarise(
      income_decile = "Overall",
      `White male` = weighted.mean(white_male, weight) * 100,
      `Black male` = weighted.mean(black_male, weight) * 100,
      `Black female` = weighted.mean(black_female, weight) * 100,
      `White female` = weighted.mean(white_female, weight) * 100,
      `Hispanic male` = weighted.mean(hispanic_male, weight) * 100,
      `Hispanic female` = weighted.mean(hispanic_female, weight) * 100,
      `Other male` = weighted.mean(other_male, weight) * 100,
      `Other female` = weighted.mean(other_female, weight) * 100,
      n = n()
    )
  
  # Combine and return results
  all_composition <- bind_rows(
    decile_composition,
    overall_shares
  )
  
  all_composition$income_decile <- factor(
    all_composition$income_decile,
    levels = c("0-10", "10-20", "20-30", "30-40", "40-50",
               "50-60", "60-70", "70-80", "80-90", "90-100", 
               "Overall") 
  )
  
  return(arrange(all_composition, income_decile) %>% select(-n))
}

# Simplified function to calculate household Gini coefficient
calculate_household_gini <- function(data) {
  # Use only records with valid income
  data <- data %>%
    filter(!is.na(income_clean))
  
  # Aggregate incomes to household level
  household_income <- data %>%
    group_by(SERIAL) %>%
    # Check if there's a household head (RELATE == 1)
    mutate(has_head = any(RELATE == 1)) %>%
    filter(has_head) %>%  # Exclude households without a head
    summarise(
      household_income = sum(income_clean, na.rm = TRUE),
      household_weight = first(weight[RELATE == 1]),
      .groups = 'drop'
    ) %>%
    filter(!is.na(household_weight))  # Ensure valid weights
  
  # Calculate overall household Gini
  overall_gini <- calculate_weighted_gini(household_income$household_income, 
                                          household_income$household_weight)
  
  return(overall_gini)
}

# Read in the single data file
cat("Reading data.csv file...\n")
full_data <- read.csv("data.csv")

# Main analysis loop for each year
for(Year in Years) {
  cat(sprintf("\n\nAnalyzing Year %d...\n", Year))
  
  # Create output filename for this year
  output_filename <- sprintf("inequality_%d_%s.txt", Year, Occupation_var)
  cat(sprintf("Results will be saved to: %s\n", output_filename))
  
  # Filter the full dataset for this year
  data <- full_data %>% 
    filter(YEAR == Year)
  
  if(nrow(data) == 0) {
    cat(sprintf("WARNING: No data found for year %d\n", Year))
    next
  }
  
  # Start output for this year
  sink(output_filename)
  
  # Prepare data
  data_prep <- prepare_data(data, Year)
  handled_income <- handle_income_codes(data, Year)
  
  # Print header information
  cat(sprintf("Year: %d\n", Year))
  cat(sprintf("Occupational variable: %s\n", Occupation_var))
  cat(sprintf("Total observations: %d\n", nrow(data)))
  cat(sprintf("Observations after filtering: %d\n", nrow(data_prep)))
  cat(sprintf("Observations with negative income: %d (%.2f%%)\n", 
              sum(data$INCTOT < 0, na.rm = TRUE),
              mean(data$INCTOT < 0, na.rm = TRUE) * 100))
  cat(sprintf("Top-coded observations: %d (%.2f%%)\n", 
              handled_income$n_topcoded, 
              handled_income$pct_topcoded))
  cat(sprintf("Pareto alpha: %.2f\n", handled_income$alpha))
  cat(sprintf("Mean income above top-code: %.2f\n", handled_income$mean_above_topcode))
  cat(sprintf("Maximum income after adjustment: %.2f\n", max(data_prep$income_clean, na.rm=TRUE)))
  
  # Print broad demographic summary 
  cat("\nBroad Demographic Group Income Statistics:\n\n")
  broad_demo_summary <- create_broad_demo_summary(data_prep)
  print(kable(broad_demo_summary$income_stats %>% 
                select(group, n_obs, pct_pop, mean_income, median_income, 
                       income_share, pct_of_national_mean, pct_of_national_median), 
              digits = 2,
              col.names = c("Group", "N", "% of Pop", "Mean Income", "Median Income", 
                            "% of Total Income", "% of National Mean", "% of National Median")))
  
  cat("\nBroad Demographic Group Characteristics:\n\n")
  print(kable(broad_demo_summary$characteristics, digits = 2,
              col.names = c("Group", "% HS Complete", "% Labor Force", "% Unemployed", 
                            "% Foreign Born", "Weeks Worked", "% College", 
                            "% Institutionalized")))
  
  # Print demographic summary
  cat("\nDemographic Group Income Statistics:\n\n")
  demo_summary <- create_demo_summary(data_prep)
  print(kable(demo_summary$income_stats %>% 
                select(group, n_obs, pct_pop, mean_income, median_income, 
                       income_share, pct_of_national_mean, pct_of_national_median), 
              digits = 2,
              col.names = c("Group", "N", "% of Pop", "Mean Income", "Median Income", 
                            "% of Total Income", "% of National Mean", "% of National Median")))
  
  cat("\nDemographic Group Characteristics:\n\n")
  print(kable(demo_summary$characteristics, digits = 2,
              col.names = c("Group", "% HS Complete", "% Labor Force", "% Unemployed", 
                            "% Foreign Born", "Weeks Worked", "% College", 
                            "% Institutionalized")))
  
  # Print education-specific summary with switched columns 
  cat("\nDemographic Group by Education Summary Statistics:\n\n")
  edu_summary <- create_education_summary(data_prep)
  print(kable(edu_summary %>% 
                select(education, group, n_obs, mean_income, median_income,
                       pct_inst, pct_of_national_mean, pct_of_national_median), 
              digits = 2,
              col.names = c("Education", "Group", "N", "Mean Income", "Median Income",
                            "% Institutionalized", "% of National Mean", "% of National Median")))
  
  # Run models
  models <- run_models(data_prep)
  
  # Print model results
  for(i in seq_along(models)) {
    print_model_results(models[[i]], i)
  }

  # Calculate and print Gini coefficients by demographic group
  cat("\nGini Coefficients by Demographic Group:\n")
  gini_summary <- calculate_gini_by_group(data_prep)
  print(kable(gini_summary, 
              col.names = c("Demographic Group", "Gini Coefficient", "Number of Observations"),
              digits = 3))
  
  # Analyze and print demographic composition across all deciles
  cat("\nDemographic Composition by Income Decile (%):\n")
  decile_results <- analyze_decile_composition(data_prep)
  print(kable(decile_results,
              col.names = c("Income Decile", "White male", "Black male", "Black female",
                            "White female", "Hispanic male", "Hispanic female",
                            "Other male", "Other female"),
              digits = 2))
  
  # Calculate and print household-level Gini coefficient
  household_gini <- calculate_household_gini(data_prep)
  cat("\nHousehold-Level Income Inequality:\n")
  cat("Overall Household Gini Coefficient:", round(household_gini, 3), "\n\n")
  
  # Close the output file
  sink()
  
  # Clear data to save memory
  rm(data)
  rm(data_prep)
  rm(broad_demo_summary)
  rm(demo_summary)
  rm(edu_summary)
  rm(models)
  rm(gini_summary)
  rm(decile_results)
  rm(household_gini)
  rm(handled_income)
  
  # Force garbage collection
  gc()
  
  cat(sprintf("Analysis for Year %d completed.\n", Year))
}

cat("\nAll analyses completed!\n")
    