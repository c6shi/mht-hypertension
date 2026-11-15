# Data Cleaning

library(here)
library(dplyr)
library(tidyverse)
library(data.table)
library(caret)
library(zoo)


##### Load Data #####
all_icpsr <- list.dirs(here("data"), full.names=F, recursive=F)
all_icpsr <- all_icpsr[sapply(all_icpsr, function(x) startsWith(x, "ICPSR"))]

for (icpsr in all_icpsr) {
  rda <- here("data", icpsr, "DS0001", sprintf("%s-0001-Data.rda", substr(icpsr, 7, 12)))
  load(rda)
}


##### Define Variables #####

# Exposure
baseline_A <- c("HORMPIL0")
A_prefix <- c("COMBIN1", "ESTROG1", "PROGES1", "ESTRNJ1")

# Outcome
Y_prefix <- c("HIGHBP", "SYSBP1", "DIABP1", "BP1")

# Baseline covariates
depress_prefix <- c("BOTHER", "APPETIT", "BLUES", "GOOD",
                    "KEEPMIN", "DEPRESS", "EFFORT", "HOPEFUL",
                    "FAILURE", "FEARFUL", "RESTLES", "HAPPY",
                    "TALKLES", "LONELY", "UNFRNDL", "ENJOY",
                    "CRYING", "SAD", "DISLIKE", "GETGOIN")
anxiety_prefix <- c("IRRITAB", "NRVOUS", "HARTRAC", "FEARFULA")

baseline_W0_inscreener <- c("DIABETE",
                            "PHY_ACT",
                            "DEGREE")

baseline_W0_inbaseline <- c("RACE", 
                            "AGE0", 
                            "HEIGHT0", 
                            "WEIGHT0", 
                            "STATUS0", 
                            "SMOKERE0",
                            "HOTFLAS0",
                            "ALCHSRV0",
                            sapply(depress_prefix, function(x) paste0(x, 0)),
                            sapply(anxiety_prefix, function(x) paste0(x, 0)))
baseline_W0_inbaseline <- stack(baseline_W0_inbaseline)$values

# Time-varying covariates
visit_W <- c("HEIGHT",
             "WEIGHT",
             "STATUS", 
             "DIABETE",
             "SMOKERE",
             "DRNKBEE",
             "PHYSACT",
             "HOTFLAS",
             depress_prefix,
             anxiety_prefix)

##### Construct Cohort from Baseline #####
# We start with N=3302 participants at baseline.

# Add screener variables to baseline dataframe
screener_df <- da04368.0001[da04368.0001$SWANID %in% da28762.0001$SWANID, 
                            c("SWANID", baseline_W0_inscreener)]

# Obtain participants that are 
# 1) pre-menopausal or early peri-menopausal, and 
# 2) have not taken hormone medications in the last month. (N=3243)
baseline_df <- da28762.0001 %>%
  filter(as.integer(STATUS0) %in% c(4, 5)) %>%
  filter(as.integer(HORMPIL0) == 1)

# Remove women with 
# 1) unknown age at menarche, or 
# 2) age at menarche > 18 years. (N=3217)
baseline_df <- baseline_df %>%
  filter(!is.na(as.integer(STARTAG0) <= 18))

# Remove women with 
# 1) NA hormone usage in visit 1 (cannot assess use of hormones between baseline and visit 1), and 
# 2) missing baseline covariates. (N=2574)
visit1_A <- sapply(A_prefix, function(x) paste0(x, 1))
visit1_hormone_df <- da29221.0001 %>% 
  filter_at(vars(visit1_A), all_vars(!is.na(.)))

baseline_df <- baseline_df[baseline_df$SWANID %in% visit1_hormone_df$SWANID, ]

baseline_df <- baseline_df %>%
  filter_at(vars(baseline_W0_inbaseline), all_vars(!is.na(.)))

# Remove women with 
# 1) NA BP0, SYSBP10, or DIABP10 (N=2566)
baseline_df <- baseline_df %>%
  filter_at(vars(BP0, SYSBP10, DIABP10), all_vars(!is.na(.)))

##### Feature Engineering: Baseline #####

# Change ALCHSRV0 to a 0/1 indicator of alcohol use
baseline_df$ALCHSRV0 <- as.numeric(baseline_df$ALCHSRV0 > 0)

# Create depression and anxiety variables
baseline_depress0 <- stack(sapply(depress_prefix, function(x) paste0(x, 0)))$values
baseline_anxiety0 <- stack(sapply(anxiety_prefix, function(x) paste0(x, 0)))$values

# Convert factor variables to numerics - 1
baseline_Y0 <- c("BP0", "SYSBP10", "DIABP10")
baseline_factor_minus1 <- c(baseline_W0_inbaseline[!baseline_W0_inbaseline %in% c("RACE", "AGE0", "HEIGHT0", "WEIGHT0", "STATUS0", "ALCHSRV0")], "HORMUSER0", "BP0")
baseline_df <- baseline_df %>%
  rename("HORMUSER0" = "HORMPIL0") %>%
  select(all_of(c("SWANID", baseline_W0_inbaseline, "HORMUSER0", baseline_Y0))) %>%
  mutate(across(all_of(baseline_factor_minus1), as.numeric)) %>%
  mutate(across(all_of(baseline_factor_minus1), ~ . -1)) %>%
  rowwise() %>%
  mutate("DEPRESSION0" = sum(c_across(all_of(baseline_depress0))),
         "ANXIETY0" = sum(c_across(all_of(baseline_anxiety0)))) %>%
  ungroup()

# Remove depression and anxiety questions
baseline_df <- baseline_df %>%
  select(-all_of(c(baseline_depress0, baseline_anxiety0)))

# One hot encode RACE
race_one_hot <- model.matrix(~ RACE-1, baseline_df)[, c(1:3, 5)]
baseline_df <- baseline_df %>%
  select(-"RACE") %>%
  cbind(race_one_hot)

# Reorder columns for LTMLE (L's, A's, Y's)
baseline_W0 <- names(select(baseline_df, -all_of(c("SWANID", "HORMUSER0", baseline_Y0))))
baseline_W0 <- c(tail(baseline_W0, 3), head(baseline_W0, -3))
baseline_df <- baseline_df %>%
  select(all_of(c("SWANID", baseline_W0, "HORMUSER0", baseline_Y0)))

# Get screener variables for individuals in baseline_df 
# Also, remove women with NA baseline covariates in screener 
screener_df <- screener_df %>%
  filter_at(vars(baseline_W0_inscreener), all_vars(!is.na(.))) %>%
  select(c(SWANID, baseline_W0_inscreener)) 
# mutate(across(where(is.factor), as.numeric))

# Merge screener dataframe and filtered baseline dataframe (N=2473)
clean_df <- merge(screener_df, baseline_df, by="SWANID")

# Ordinal encoding and fix DIABETE
ordinal_cols <- c("PHY_ACT", "DEGREE", "STATUS0")
clean_df <- clean_df %>%
  mutate(across(all_of(c("DIABETE", ordinal_cols)), as.numeric)) %>%
  mutate(across(all_of(c("DIABETE", ordinal_cols)), ~ . -1)) %>%
  rowwise()

# Rename some baseline L's to match visit L's
clean_df <- clean_df %>%
  rename("AGE" = "AGE0",
         "DRNKBEE0" ="ALCHSRV0",
         "PHYSACT0" = "PHY_ACT", 
         "DIABETE0" = "DIABETE",
         "MHT0" = "HORMUSER0")

# Make outcome indicator
clean_df <- clean_df %>%
  mutate(HTN0 = as.integer(BP0 == 1 | SYSBP10 >= 130 | DIABP10 >= 80)) %>%
  select(-all_of(baseline_Y0))

# Final reorder
W <- c("AGE", colnames(race_one_hot), "DEGREE")
L0 <- colnames(clean_df)[!colnames(clean_df) %in% W][-c(1, (ncol(clean_df) - length(W) - 1):(ncol(clean_df) - length(W)))]
clean_df <- clean_df %>%
  select(all_of(c("SWANID", W, L0, "MHT0", "HTN0")))

##### Feature Engineering: Visit T #####
for (i in 3:length(all_icpsr)) {
  icpsr <- all_icpsr[i]
  code <- substr(icpsr, 7, 12)
  visit <- i-2
  
  print(c(code, visit))
  
  visit_i_W <- paste(visit_W, visit, sep="")
  visit_i_A <- paste(A_prefix, visit, sep="")
  visit_i_Y <- paste(Y_prefix, visit, sep="")
  visit_i_df <- eval(parse(text=paste(sprintf("da%s.0001", code))))
  visit_i_df <- tryCatch(
    visit_i_df <- visit_i_df %>%
      select(all_of(c("SWANID", visit_i_W, visit_i_A, visit_i_Y))),
    error = function(cond) {
      message(conditionMessage(cond))
      visit_i_df <- visit_i_df %>%
        select(any_of(c("SWANID", visit_i_W, visit_i_A, visit_i_Y)))
      visit_i_missing_W <- setdiff(c(visit_i_W, visit_i_A, visit_i_Y), 
                                   colnames(visit_i_df))
      for (missing_col in visit_i_missing_W) {
        visit_i_df[missing_col] <- NA
      }
      message("Imputed missing columns.")
      return(visit_i_df)
    },
    warning = function(cond) {
      message(conditionMessage(cond))
      NULL
    },
    finally = {
      message("Working...")
    }
  )
  
  
  # Create depression and anxiety variables
  visit_i_depress <- paste(depress_prefix, visit, sep="")
  visit_i_anxiety <- paste(anxiety_prefix, visit, sep="")
  
  # For the HTN outcome indicator
  bp_col <- paste0("BP1", visit)
  sysbp_col  <- paste0("SYSBP1", visit) 
  diabp_col  <- paste0("DIABP1", visit)
  highbp_col <- paste0("HIGHBP", visit)
  
  # Convert factor variables to numerics - 1
  visit_W_factor <- visit_W[!visit_W %in% c("HEIGHT", "WEIGHT", "STATUS", "PHYSACT")]
  visit_Y_factor <- Y_prefix[!Y_prefix %in% c("SYSBP1", "DIABP1")]
  visit_i_W_factor <- paste(visit_W_factor, visit, sep="")
  visit_i_Y_factor <- paste(visit_Y_factor, visit, sep="")
  visit_i_W_factor_minus1 <- c(visit_i_W_factor, visit_i_A, visit_i_Y_factor)
  visit_i_df <- visit_i_df %>%
    mutate(across(all_of(visit_i_W_factor_minus1), as.numeric)) %>%
    mutate(across(all_of(visit_i_W_factor_minus1), ~ . -1)) %>%
    rowwise() %>%
    mutate("DEPRESSION#" = sum(c_across(all_of(visit_i_depress))),
           "ANXIETY#" = sum(c_across(all_of(visit_i_anxiety))),
           "MHT#" = case_when(
             sum(c_across(all_of(visit_i_A))) == 0 ~ 0,
             sum(c_across(all_of(visit_i_A))) > 0 ~ 1,
             TRUE ~ NA
           ),
           "HTN#" = as.integer(
               .data[[bp_col]] == 1 |
               .data[[sysbp_col]] >= 130 |
               .data[[diabp_col]] >= 80 |
               .data[[highbp_col]] == 1
           )) %>%
    ungroup() %>%
    rename_with(~gsub("#", visit, .x, fixed=T))
  
  # Remove depression, anxiety, and exposure questions
  visit_i_df <- visit_i_df %>%
    select(-all_of(c(visit_i_depress, visit_i_anxiety, visit_i_A, visit_i_Y)))
  
  # Ordinal encoding
  ordinal_cols <- setdiff(c(visit_i_W, visit_i_A), visit_i_W_factor_minus1)[-c(1, 2)]
  visit_i_df <- visit_i_df %>%
    mutate(across(all_of(ordinal_cols), as.numeric)) %>%
    mutate(across(all_of(ordinal_cols), ~ . -1)) %>%
    rowwise()
  
  # Create missingness indicators for covariates and exposure separately (these are part of the covariates)
  # Create censoring (truly right-censoring, LTFU) for outcomes
  visit_i_cols <- colnames(visit_i_df)
  visit_i_W <- visit_i_cols[!visit_i_cols %in% visit_i_Y & visit_i_cols != paste0("MHT", visit) & visit_i_cols != "SWANID"]
  visit_i_df[[paste0("D_TVCOV", visit)]] <- as.numeric(apply(visit_i_df[visit_i_W], 1, function(x) any(is.na(x))))
  visit_i_df[[paste0("D_MHT", visit)]] <- as.numeric(is.na(visit_i_df[[paste0("MHT", visit)]]))
  visit_i_df[[paste0("C_HTN", visit)]] <- as.numeric(is.na(visit_i_df[[paste0("HTN", visit)]]))
  
  # Reorder
  Li <- paste(substr(L0, 1, nchar(L0) - 1), visit, sep="")
  visit_i_df <- visit_i_df %>%
    select(all_of(c("SWANID", Li, 
                    paste0("D_TVCOV", visit),
                    paste0("D_MHT", visit), paste0("MHT", visit), 
                    paste0("C_HTN", visit), paste0("HTN", visit))))
  
  # Merge
  clean_df <- merge(clean_df, visit_i_df, by="SWANID", all.x=T)
  
  # Refill D_TVCOV, D_MHT, C_HTN
  censor_cols <- c(paste0("D_TVCOV", visit), paste0("D_MHT", visit), paste0("C_HTN", visit))
  clean_df <- clean_df %>%
    mutate(across(all_of(censor_cols), ~replace_na(.x, 1)))
  
  # LOCF for intermediate covariates (Comment this  whole section out for eda-data.csv)
  int_cov <- colnames(visit_i_df)[!colnames(visit_i_df) %in% c("SWANID", censor_cols)]
  
  for (col in int_cov) {
    subtract <- ifelse(visit == 10, 2, 1)
    prev_col <- paste0(substr(col, 1, nchar(col) - subtract), visit-1)
    ind <- which(is.na(clean_df[col]))
    clean_df[col][ind, ] <- clean_df[prev_col][ind, ]
  }
}

# Filter out women with HTN0 == 1!
clean_df_htnfree <- clean_df %>%
  filter(HTN0 == 0)

# Write to csv
write.csv(clean_df_htnfree, here("data", "mht_final.csv"), row.names=FALSE)

# Consider if we need to shift everything back:
# to ensure the time-ordering of the A's and L's is correct
shift_clean_df <- copy(clean_df_htnfree)
shift_cols <- paste0(c("D_MHT", "MHT"), rep(c(0:9), each=2))
for (col in shift_cols) {
  current_visit <- as.numeric(substr(col, nchar(col), nchar(col)))
  next_col <- paste0(substr(col, 1, nchar(col) - 1), current_visit + 1)
  shift_clean_df[col] <- shift_clean_df[next_col]
}

shift_clean_df <- shift_clean_df %>%
  select(-c(D_MHT0, D_MHT10, MHT10))

write.csv(shift_clean_df, here("data", "mht_final_shiftexp.csv"), row.names=FALSE)
