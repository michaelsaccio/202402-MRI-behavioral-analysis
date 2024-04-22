args = commandArgs(TRUE)
sub_num = args[1]


library(dplyr)
library(rio)
library(ade4)
library(broom)
library(reshape2)


# Create a list of the scanner files for every subject
scanner_files_path <- ""
scan_files <- list.files(scanner_files_path)

# Behavioral questions
behavioral_questions_path <- ""
behavioral <- import(behavioral_questions_path)

# Meants files
meants_files_path <- ""

# Exported CSV files
export_path <- ""

# Define number of parcels in parcellation scheme
parcs <- 400

# Create tibble for all Mantel tests results 
mantle_output_full <- tibble()



# GLOBAL FUNCTIONS
# Given a dataframe column, fills in NA values using randomly generated imputed that follow mean and SD
impute_na <- function(data_frame, column_name) {
  # Calculate the mean and standard deviation of non-NA values in the specified column
  mean_value <- mean(data_frame[[column_name]], na.rm = TRUE)
  sd_value <- sd(data_frame[[column_name]], na.rm = TRUE)
  
  # Get the indices of NA values
  na_indices <- which(is.na(data_frame[[column_name]]))
  
  # Impute NA values with random values from a normal distribution
  set.seed(123)  # Set seed for reproducibility
  imputed_values <- rnorm(length(na_indices), mean = mean_value, sd = sd_value)
  
  # Fill in NA values in the specified column with imputed values
  data_frame[[column_name]][na_indices] <- imputed_values
  
  # Return the modified data frame
  return(data_frame)
}

# Function to calculate dissimilarity matrix
calculate_dissimilarity_matrix <- function(data, columns) {
  vec <- vector()
  for (i in 1:nrow(data)) {
    anchor <- data[i, columns]
    for (j in 1:nrow(data)) {
      compare <- data[j, columns]
      dist_vec <- rbind(anchor, compare)
      dist <- stats::dist(dist_vec, method = "euclidean")
      vec <- append(vec, dist)
    }
  }
  matrix(vec, nrow = nrow(data), ncol = nrow(data))
}

# Function to update beh_scan_p with social relationship ratings
update_beh_scan_p <- function(scan_p, sub_behavioral, personality_names) {
  for (i in 1:6) {
    beh <- sub_behavioral %>%
      filter(target == i)
    for (j in 1:length(personality_names)) {
      temp <- beh %>%
        select(personality_names[j])
      scan_p[[personality_names[j]]] <- ifelse(scan_p$target == i, temp[[1]], scan_p[[personality_names[j]]])
    }
  }
  scan_p
}

# Check if the dissimilarity matrix has 0-only diagonals
check_diagonal_zero <- function(matrix) {
  all(diag(matrix) == 0)
}

# Masks DSM by setting the lower or upper triangle to NA and optionally melts the matrix into a long format
dsm_mask <- function (dsm, lower = FALSE, melt = TRUE) {
  require(reshape2)
  
  if (melt == FALSE) {
    if (lower == FALSE) {
      mat <- as.matrix(dsm)
      mat[lower.tri(mat, diag = TRUE)] <- NA
      mat
    }
    else {
      mat <- as.matrix(dsm)
      mat[upper.tri(mat, diag = TRUE)] <- NA
      mat
    }
  }
  else {
    if (lower == FALSE) {
      mat <- as.matrix(dsm)
      mat[lower.tri(mat, diag = TRUE)] <- NA
      mat_m <- melt(mat)
      mat_m$Var1 <- as.factor(mat_m$Var1)
      mat_m$Var2 <- as.factor(mat_m$Var2)
      mat_m$Var2 <- factor(mat_m$Var2, levels = rev(levels(mat_m$Var2)))
      mat_m
    }
    else {
      mat <- as.matrix(dsm)
      mat[upper.tri(mat, diag = TRUE)] <- NA
      mat_m <- melt(mat)
      mat_m$Var1 <- as.factor(mat_m$Var1)
      mat_m$Var2 <- as.factor(mat_m$Var2)
      mat_m$Var2 <- factor(mat_m$Var2, levels = rev(levels(mat_m$Var2)))
      mat_m
    }
  }
}

# Get rid of half of the matrix that contains duplicate information and put the values into vector form
vectorize_dsm_matrix <- function(matrix, col_name, expected_length = 64620) {
  dsm_vec <- dsm_mask(matrix) %>%
    filter(!is.na(value)) %>%
    select(value)
  
  # Ensure the vector has the expected length
  dsm_vec <- dsm_vec[1:expected_length, , drop = FALSE]
  
  # Label the vector as neural data
  colnames(dsm_vec) <- col_name
  
  return(dsm_vec)
}

# Function to run linear models and save output
run_linear_model <- function(response_variable, predictor_variables, data_frame) {
  formula_str <- paste(response_variable, "~", paste(predictor_variables, collapse = " + "))
  lm_model <- lm(formula_str, data = data_frame)
  return(lm_model)
}

# Function to save linear model output and create a tibble of results
save_model_output <- function(model, output_name, sub_num, parcel_num) {
  output <- tidy(model)
  output <- output[-1,] %>%
    mutate(subject = sub_num) %>%
    mutate(parcel = parcel_num)
  
  # Check if a tibble with the desired output name already exists in the global environment
  output_tibble_name <- paste(output_name, "_linear_output", sep = "")
  if (exists(output_tibble_name, envir = .GlobalEnv)) {
    # If it exists, retrieve it and bind the new rows
    existing_output_tibble <- get(output_tibble_name, envir = .GlobalEnv)
    output_tibble <- rbind(existing_output_tibble, output)
  } else {
    # If it doesn't exist, use the newly created output as the tibble
    output_tibble <- output
  }
  
  # Assign the updated or new tibble to the dynamically named variable in the global environment
  assign(output_tibble_name, output_tibble, envir = .GlobalEnv)
}

# Function to run Mantel test and save output, now correctly returning the updated tibble
run_mantel_test <- function(predictor_matrix, predictor_name, output_tibble, neural_dist = NULL, nrepet = 10000) {
  predictor_dist <- dist(predictor_matrix)
  mantel_result <- mantel.rtest(neural_dist, predictor_dist, nrepet = nrepet)
  
  # Assuming mantel_result is a list and contains a summary where the p-value and effect can be directly indexed
  p_value <- mantel_result[[5]]  # Check if this directly accesses the p-value
  effect <- mantel_result[[1]]  # Check if this directly accesses the effect size
  
  # Ensure both are numeric for consistent data type in the tibble
  p_value <- as.numeric(p_value)
  effect <- as.numeric(effect)
  
  result_df <- tibble(
    predictor_name = predictor_name,
    p_value = p_value,
    effect = effect
  )
  
  # Update individual output for the predictor in the global environment
  assign(paste0(predictor_name, "_mantel_output"), result_df, envir = .GlobalEnv)
  
  # Append the result to output_tibble and return it
  return(bind_rows(output_tibble, result_df))
}



# RSA DSM CREATION SECTION
# Save the subject number
sub <- sub_num

# Import and order the scanner file
rsa_scan <- import(file.path(scanner_files_path, current_file)) %>%
  arrange(run, target, word)
# Fill in NA values
rsa_scan <- impute_na(rsa_scan, "response.rt")

# Extract unique targets, valence, words, and responses
tars <- unique(rsa_scan$target)
tars <- sort(tars)
valence <- unique(rsa_scan$condition)
words <- unique(rsa_scan$word)
responses <- unique(rsa_scan$response)

# Dummy code target, valence, and word information
for (i in 1:length(tars))
  rsa_scan <- rsa_scan %>% mutate(!! tars[i] := ifelse(target == tars[i], 1, 0))
for (i in 1:length(valence))
  rsa_scan <- rsa_scan %>% mutate(!! valence[i] := ifelse(condition == valence[i], 1, 0))
for (i in 1:length(words))
  rsa_scan <- rsa_scan %>% mutate(!! words[i] := ifelse(word == words[i], 1, 0))
rsa_scan <- rsa_scan %>% mutate(response = ifelse(response.keys == 6, 1, 0))

# Calculate dissimilarity matrices
target_matrix <- calculate_dissimilarity_matrix(rsa_scan, 7:12)
valence_matrix <- calculate_dissimilarity_matrix(rsa_scan, 13:14)
word_matrix <- calculate_dissimilarity_matrix(rsa_scan, 15:74)
rt_matrix <- calculate_dissimilarity_matrix(rsa_scan, 6)

# BEH DSM CREATION SECTION
# Import the scanner file for the subject that the for loop is on
beh_scan <- import(file.path(scanner_files_path, current_file)) %>%
  separate(target, c("dis", "target"), sep = 3) %>% 
  select(!dis)
beh_scan$target <- as.numeric(beh_scan$target)

# Order the scanner output alphabetically by run then target then word so that it matches the order as the merged zstat files (important)
beh_scan <- beh_scan[with(beh_scan, order(run,target, word)), ]
# Isolate the group (g) and sub (s) digits
g <- scan_files[1] %>%
  substr(3, 3)
g <- as.numeric(g)
s <- scan_files[1] %>%
  substr(6, 6)
s <- as.numeric(s)

# Use those digits to filter the behavioral data set to be only the subject that we are on in the loop
sub_behavioral <- behavioral %>%
  filter(group == g) %>%
  filter(subject == s)
# Isolate the social relationship variables in the behavioral data
social <- sub_behavioral[, c(1:3, 55:63)]
# Create a list of the names of the social questions (social_names) as strings
social_names <- colnames(social)
social_names <- social_names[4:12]
# Create a new column in the scanner output dataframe for each social relationship question
for (i in 1:length(social_names)) {
  beh_scan <- beh_scan %>%
    mutate(!! social_names[i] := 0)
}
# Loop through the filtered behavioral data (sub_behavioral) one target at a time
# Loop through the social relationship questions and bring the ratings for those questions into the scan dataset
# Each row should now have the behavioral ratings that the current sub indicated for the target listed on that row in their correct social relationship question columns
for (i in 1:6) {
  beh <- sub_behavioral %>%
    filter(target == i)
  for (j in 1:length(social_names)) {
    temp <- beh %>%
      select(social_names[j])
    beh_scan[[social_names[j]]] <- ifelse(beh_scan$target == i, temp[[1]], beh_scan[[social_names[j]]])
  }
}

# Load scan data
know_matrix <- calculate_dissimilarity_matrix(beh_scan[, c(1:3, 7, 10:12)], c(4:7))
similarity_matrix <- calculate_dissimilarity_matrix(beh_scan[, c(1:3, 8:9, 13)], c(4:6))
friend_matrix <- calculate_dissimilarity_matrix(beh_scan[, c(1:3, 14)], c(4))
like_matrix <- calculate_dissimilarity_matrix(beh_scan[, c(1:3, 15)], c(4))

# Load personality data
beh_scan_p <- import(sprintf("%s/%s", scanner_files_path, current_file)) %>% 
  separate(target, c("dis", "target"), sep = 3) %>% 
  select(!dis)

beh_scan_p$target <- as.numeric(beh_scan_p$target)
beh_scan_p <- beh_scan_p[with(beh_scan_p, order(run, target, word)), ]

# Extract social relationship variables
social_relationship_variables <- sub_behavioral[, c(1:37)]

# Create personality columns in beh_scan_p
personality_names <- colnames(social_relationship_variables)[4:37]
for (i in 1:length(personality_names)) {
  beh_scan_p <- beh_scan_p %>% mutate(!!personality_names[i] := 0)
}

# Update beh_scan_p with social relationship ratings
beh_scan_p <- update_beh_scan_p(beh_scan_p, sub_behavioral, personality_names)

# Calculate dissimilarity matrices for personality traits
personality_col <- beh_scan_p[, 7:40]
personality_matrix <- calculate_dissimilarity_matrix(personality_col)
extroversion_matrix <- calculate_dissimilarity_matrix(personality_col, c(1:4, 9, 14, 19, 24, 29))
agreeable_matrix <- calculate_dissimilarity_matrix(personality_col, c(1:3, 5, 10, 15, 20, 25, 30))
conscientious_matrix <- calculate_dissimilarity_matrix(personality_col, c(1:3, 6, 11, 16, 21, 26, 31))
neuroticism_matrix <- calculate_dissimilarity_matrix(personality_col, c(1:3, 7, 12, 17, 22, 27, 32))
openness_matrix <- calculate_dissimilarity_matrix(personality_col, c(1:3, 8, 13, 18, 23, 28, 33))

# Check the diagonal values for each matrix
check_target <- check_diagonal_zero(target_matrix)
check_valence <- check_diagonal_zero(valence_matrix)
check_word <- check_diagonal_zero(word_matrix)
check_rt <- check_diagonal_zero(rt_matrix)

check_know <- check_diagonal_zero(know_matrix)
check_similarity <- check_diagonal_zero(similarity_matrix)
check_friend <- check_diagonal_zero(friend_matrix)
check_like <- check_diagonal_zero(like_matrix)
check_personality <- check_diagonal_zero(personality_matrix)
check_extroversion <- check_diagonal_zero(extroversion_matrix)
check_agreeable <- check_diagonal_zero(agreeable_matrix)
check_conscientious <- check_diagonal_zero(conscientious_matrix)
check_neuroticism <- check_diagonal_zero(neuroticism_matrix)
check_openness <- check_diagonal_zero(openness_matrix)



# Modeling initialization
groups <- c("G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10",
            "G11", "G12", "G13", "G14", "G15", "G16", "G17", "G18", "G19", "G20")
subjects <- 1:6
files_per_subject <- 100

# Loop through groups, subjects, and meant files
for (group in groups) {
  for (sub in subjects) {
    for (i in 1:files_per_subject) {
      
      # File Processing
      file_path <- file.path(meants_files_path, paste0(group, "S0", sub), paste0("meants_parc", i))
      zstats <- read.table(file_path, quote = "\"", comment.char = "", skip = 3)
      # Turn the zstats into a dissimilarity matrix
      dis_mat <- 1 - cor(t(zstats), method = 'spearman')
      
      # Vectorize dissimilarity matrices
      neural_dsm_vec <- vectorize_dsm_matrix(dis_mat, "neural")
      target_dsm_vec <- vectorize_dsm_matrix(target_matrix, "target")
      valence_dsm_vec <- vectorize_dsm_matrix(valence_matrix, "valence")
      word_dsm_vec <- vectorize_dsm_matrix(word_matrix, "word")
      rt_dsm_vec <- vectorize_dsm_matrix(rt_matrix, "response")
      
      know_dsm_vec <- vectorize_dsm_matrix(know_matrix, "know")
      similarity_dsm_vec <- vectorize_dsm_matrix(similarity_matrix, "similarity")
      friend_dsm_vec <- vectorize_dsm_matrix(friend_matrix, "friend")
      like_dsm_vec <- vectorize_dsm_matrix(like_matrix, "like")
      personality_dsm_vec <- vectorize_dsm_matrix(personality_matrix, "personality")
      extroversion_dsm_vec <- vectorize_dsm_matrix(extroversion_matrix, "extroversion")
      agreeable_dsm_vec <- vectorize_dsm_matrix(agreeable_matrix, "agreeable")
      conscientious_dsm_vec <- vectorize_dsm_matrix(conscientious_matrix, "conscientious")
      neuroticism_dsm_vec <- vectorize_dsm_matrix(neuroticism_matrix, "neuroticism")
      openness_dsm_vec <- vectorize_dsm_matrix(openness_matrix, "openness")
      
      # Combine all dissimilarity matrices into dataframes
      rsa_full_dsm_df <- cbind(neural_dsm_vec, target_dsm_vec, valence_dsm_vec, word_dsm_vec, rt_dsm_vec)
      beh_full_dsm_df <- cbind(neural_dsm_vec, know_dsm_vec, similarity_dsm_vec, friend_dsm_vec, like_dsm_vec, personality_dsm_vec, extroversion_dsm_vec, agreeable_dsm_vec, conscientious_dsm_vec, neuroticism_dsm_vec, openness_dsm_vec)
      
      # Calling the function to run linear models predicting neural dissim from variable dissims
      rsa_fit_full <- run_linear_model("neural", c("target", "valence", "word", "response"), rsa_full_dsm_df)
      fit_tar_val <- run_linear_model("neural", c("target", "valence"), rsa_full_dsm_df)
      fit_tar <- run_linear_model("neural", c("target"), rsa_full_dsm_df)
      fit_val <- run_linear_model("neural", c("valence"), rsa_full_dsm_df)
      fit_word <- run_linear_model("neural", c("word"), rsa_full_dsm_df)
      fit_response <- run_linear_model("neural", c("response"), rsa_full_dsm_df)
      
      beh_fit_full <- run_linear_model("neural", c("know", "similarity", "friend", "like", "personality", "extroversion", "agreeable", "conscientious", "neuroticism", "openness"), beh_full_dsm_df)
      fit_know <- run_linear_model("neural", c("know"), beh_full_dsm_df)
      fit_sim <- run_linear_model("neural", c("similarity"), beh_full_dsm_df)
      fit_friend <- run_linear_model("neural", c("friend"), beh_full_dsm_df)
      fit_like <- run_linear_model("neural", c("like"), beh_full_dsm_df)
      fit_personality <- run_linear_model("neural", c("personality"), beh_full_dsm_df)
      fit_extroversion <- run_linear_model("neural", c("extroversion"), beh_full_dsm_df)
      fit_agreeable <- run_linear_model("neural", c("agreeable"), beh_full_dsm_df)
      fit_conscientious <- run_linear_model("neural", c("conscientious"), beh_full_dsm_df)
      fit_neuroticism <- run_linear_model("neural", c("neuroticism"), beh_full_dsm_df)
      fit_openness <- run_linear_model("neural", c("openness"), beh_full_dsm_df)
      
      # Function calls for each model output
      save_model_output(rsa_fit_full, "rsa_full", sub, i)
      save_model_output(fit_tar_val, "tar_val", sub, i)
      save_model_output(fit_tar, "tar", sub, i)
      save_model_output(fit_val, "val", sub, i)
      save_model_output(fit_word, "word", sub, i)
      save_model_output(fit_response, "rt", sub, i)
      
      save_model_output(beh_fit_full, "beh_full", sub, i)
      save_model_output(fit_know, "know", sub, i)
      save_model_output(fit_sim, "sim", sub, i)
      save_model_output(fit_friend, "friend", sub, i)
      save_model_output(fit_like, "like", sub, i)
      save_model_output(fit_personality, "personality", sub, i)
      save_model_output(fit_extroversion, "extroversion", sub, i)
      save_model_output(fit_agreeable, "agreeable", sub, i)
      save_model_output(fit_conscientious, "conscientious", sub, i)
      save_model_output(fit_neuroticism, "neuroticism", sub, i)
      save_model_output(fit_openness, "openness", sub, i)
      
      # Turn all matrices into distance objects for Mantel tests
      neural_dist <- dist(dis_mat)
      target_dist <- dist(target_matrix)
      val_dist <- dist(valence_matrix)
      word_dist <- dist(word_matrix)
      rt_dist <- dist(rt_matrix)
      
      know_dist <- dist(know_matrix)
      sim_dist <- dist(similarity_matrix)
      friend_dist <- dist(friend_matrix)
      like_dist <- dist(like_matrix)
      personality_dist <- dist(personality_matrix)
      extroversion_dist <- dist(extroversion_matrix)
      agreeable_dist <- dist(agreeable_matrix)
      conscientious_dist <- dist(conscientious_matrix)
      neuroticism_dist <- dist(neuroticism_matrix)
      openness_dist <- dist(openness_matrix)
      
      # Run Mantel tests for each predictor separately
      mantel_output_full <- run_mantel_test(target_matrix, "target", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(valence_matrix, "valence", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(word_matrix, "word", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(rt_matrix, "response", mantel_output_full, neural_dist)
      
      mantel_output_full <- run_mantel_test(know_dist, "know", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(sim_dist, "sim", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(friend_dist, "friend", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(like_dist, "like", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(personality_dist, "personality", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(extroversion_dist, "extroversion", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(agreeable_dist, "agreeable", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(conscientious_dist, "conscientious", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(neuroticism_dist, "neuroticism", mantel_output_full, neural_dist)
      mantel_output_full <- run_mantel_test(openness_dist, "openness", mantel_output_full, neural_dist)
    }
  }
}

# Export all LMs and the full Mantel output as CSVs
export(rsa_full_linear_output, sprintf("%s/rsa_lm_full/%s.csv", export_path, sub))
export(tar_val_linear_output, sprintf("%s/lm_tar_val/%s.csv", export_path, sub))
export(tar_linear_output, sprintf("%s/lm_tar/%s.csv", export_path, sub))
export(val_linear_output, sprintf("%s/lm_val/%s.csv", export_path, sub))
export(word_linear_output, sprintf("%s/lm_word/%s.csv", export_path, sub))
export(rt_linear_output, sprintf("%s/lm_response_time/%s.csv", export_path, sub))

export(beh_full_linear_output, sprintf("%s/beh_lm_full/%s.csv", export_path, sub))
export(know_linear_output, sprintf("%s/lm_know/%s.csv", export_path, sub))
export(sim_linear_output, sprintf("%s/lm_sim/%s.csv", export_path, sub))
export(friend_linear_output, sprintf("%s/lm_friend/%s.csv", export_path, sub))
export(like_linear_output, sprintf("%s/lm_like/%s.csv", export_path, sub))
export(personality_linear_output, sprintf("%s/lm_personality/%s.csv", export_path, sub))
export(extroversion_linear_output, sprintf("%s/lm_extroversion/%s.csv", export_path, sub))
export(agreeable_linear_output, sprintf("%s/lm_agreeable/%s.csv", export_path, sub))
export(conscientious_linear_output, sprintf("%s/lm_conscientious/%s.csv", export_path, sub))
export(neuroticism_linear_output, sprintf("%s/lm_neuroticism/%s.csv", export_path, sub))
export(openness_linear_output, sprintf("%s/lm_openness/%s.csv", export_path, sub))

export(mantel_output_full, sprintf("%s/mantel/%s.csv", export_path, sub))




