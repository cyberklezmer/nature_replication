# Function to get comprehensive variable importance
get_variable_importance <- function(multinom_model, glmnet_model, X, y, original_vars) {
  all_vars <- colnames(X)
  
  # 1. Coefficient-based importance (multinom)
  coeffs <- coef(multinom_model)
  if(is.matrix(coeffs)) {
    coeff_importance <- apply(abs(coeffs), 2, mean)
  } else {
    coeff_importance <- abs(coeffs)
  }
  
  # 2. Penalized coefficients importance (glmnet)
  glmnet_coeffs <- coef(glmnet_model)
  if(is.list(glmnet_coeffs)) {
    # Multinomial case
    glmnet_importance <- apply(sapply(glmnet_coeffs, function(x) abs(x[-1])), 1, mean)
    names(glmnet_importance) <- all_vars
  } else {
    glmnet_importance <- abs(glmnet_coeffs[-1])  # Remove intercept
  }
  
  # 3. Random Forest importance (as reference)
  rf_model <- randomForest(x = X, y = as.factor(y), importance = TRUE)
  rf_importance <- importance(rf_model, type = 1)[, 1]
  
  # Combine all importance measures
  importance_combined <- data.frame(
    Variable = all_vars,
    Multinom_Importance = coeff_importance[all_vars],
    Glmnet_Importance = glmnet_importance[all_vars],
    RF_Importance = rf_importance[all_vars]
  )
  
  # Calculate composite importance (average of standardized measures)
  importance_combined <- importance_combined %>%
    mutate(
      Multinom_Std = scale(Multinom_Importance)[, 1],
      Glmnet_Std = scale(Glmnet_Importance)[, 1],
      RF_Std = scale(RF_Importance)[, 1],
      Composite_Importance = (Multinom_Std + Glmnet_Std + RF_Std) / 3
    ) %>%
    arrange(desc(Composite_Importance))
  
  return(importance_combined)
}

# Function to group dummy variables back to original variables
group_variable_importance <- function(importance_df, original_vars) {
  
  # Initialize results
  grouped_importance <- data.frame(
    Original_Variable = character(),
    Variable_Type = character(),
    Total_Importance = numeric(),
    Max_Importance = numeric(),
    Mean_Importance = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(orig_var in original_vars) {
    # Find all dummy variables that belong to this original variable
    pattern <- paste0("^", orig_var, "($|\\.|[^a-zA-Z0-9_])")
    related_vars <- grep(pattern, importance_df$Variable, value = TRUE)
    
    if(length(related_vars) > 0) {
      related_importance <- importance_df[importance_df$Variable %in% related_vars, ]
      
      # Determine variable type
      var_type <- if(length(related_vars) == 1) {
        "Continuous/Binary"
      } else if(any(grepl("\\.L$|\\.Q$|\\.C$|\\^[0-9]$", related_vars))) {
        "Ordinal"
      } else {
        "Categorical"
      }
      
      # Calculate summary statistics
      total_imp <- sum(related_importance$Importance)
      max_imp <- max(related_importance$Importance)
      mean_imp <- mean(related_importance$Importance)
      
      grouped_importance <- rbind(grouped_importance, data.frame(
        Original_Variable = orig_var,
        Variable_Type = var_type,
        Total_Importance = total_imp,
        Max_Importance = max_imp,
        Mean_Importance = mean_imp
      ))
    }
  }
  
  return(grouped_importance)
}

# print confusion matrix
print.cfM <- function(cfM,file="confusion-matrix.pdf", limits=c(0,max(cfM)), title=""){
  dat2 <- as_tibble(reshape2::melt(cfM)) # converting data to tibble
    ggplot(dat2, aes(Prediction, Reference)) + 
      ggtitle(title) +
      geom_tile(aes(fill = value)) +
      geom_text(aes(label = round(value))) +
      # scale_fill_gradient(low = "white", high = "red", transform = "pseudo_log",limits=limits) +
      scale_fill_gradient(low = "white", high = "red", limits=limits) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), legend.text = element_text(size=5))  
    ggsave(file)
}



# Variable Importance Analysis for Multinom Model
# This script extracts coefficients, calculates p-values, and creates importance plots

library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)

# Function to extract variable importance from multinom model
extract_variable_importance <- function(multinom_model) {
  
  # Extract coefficients and standard errors
  coeffs <- summary(multinom_model)$coefficients
  std_errors <- summary(multinom_model)$standard.errors
  
  # Calculate z-values and p-values
  z_values <- coeffs / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_values)))
  
  # Convert to long format for easier manipulation
  coeff_long <- as.data.frame(coeffs) %>%
    mutate(outcome = rownames(coeffs)) %>%
    pivot_longer(-outcome, names_to = "variable", values_to = "coefficient")
  
  std_err_long <- as.data.frame(std_errors) %>%
    mutate(outcome = rownames(std_errors)) %>%
    pivot_longer(-outcome, names_to = "variable", values_to = "std_error")
  
  p_val_long <- as.data.frame(p_values) %>%
    mutate(outcome = rownames(p_values)) %>%
    pivot_longer(-outcome, names_to = "variable", values_to = "p_value")
  
  # Combine all information
  results <- coeff_long %>%
    left_join(std_err_long, by = c("outcome", "variable")) %>%
    left_join(p_val_long, by = c("outcome", "variable"))
  
  return(results)
}

# Function to aggregate polynomial contrasts by base variable
aggregate_polynomial_importance <- function(results_df) {
  
  # Extract base variable names (remove .L, .Q, .C, ^4, etc.)
  results_df <- results_df %>%
    mutate(
      base_variable = gsub("\\.(L|Q|C)$|\\^[0-9]+$|`([^`]+)\\^[0-9]+`", "\\2", variable),
      base_variable = gsub("`", "", base_variable),
      contrast_type = case_when(
        grepl("\\.L$", variable) ~ "Linear",
        grepl("\\.Q$", variable) ~ "Quadratic", 
        grepl("\\.C$", variable) ~ "Cubic",
        grepl("\\^4", variable) ~ "Quartic",
        grepl("\\^[0-9]+", variable) ~ "Higher-order",
        TRUE ~ "Main"
      )
    )
  
  # Calculate importance metrics by base variable
  variable_importance <- results_df %>%
    filter(variable != "(Intercept)") %>%
    group_by(base_variable) %>%
    summarise(
      # Overall importance metrics
      max_abs_coeff = max(abs(coefficient), na.rm = TRUE),
      mean_abs_coeff = mean(abs(coefficient), na.rm = TRUE),
      sum_abs_coeff = sum(abs(coefficient), na.rm = TRUE),
      
      # P-value metrics
      min_p_value = min(p_value, na.rm = TRUE),
      mean_p_value = mean(p_value, na.rm = TRUE),
      prop_significant = mean(p_value < 0.05, na.rm = TRUE),
      
      # Count of terms
      n_terms = n(),
      n_significant = sum(p_value < 0.05, na.rm = TRUE),
      
      # Contrast information
      contrast_types = paste(unique(contrast_type), collapse = ", "),
      
      .groups = 'drop'
    ) %>%
    arrange(desc(max_abs_coeff))
  
  return(variable_importance)
}

# Function to create detailed coefficient plot
plot_detailed_coefficients <- function(results_df, top_n = 15) {
  
  # Add base variable names
  plot_data <- results_df %>%
    filter(variable != "(Intercept)") %>%
    mutate(
      base_variable = gsub("\\.(L|Q|C)$|\\^[0-9]+$|`([^`]+)\\^[0-9]+`", "\\2", variable),
      base_variable = gsub("`", "", base_variable),
      significance = ifelse(p_value < 0.001, "p < 0.001",
                            ifelse(p_value < 0.01, "p < 0.01", 
                                   ifelse(p_value < 0.05, "p < 0.05", "ns")))
    )
  
  # Get top variables by maximum absolute coefficient
  top_vars <- plot_data %>%
    group_by(base_variable) %>%
    summarise(max_coeff = max(abs(coefficient)), .groups = 'drop') %>%
    top_n(top_n, max_coeff) %>%
    pull(base_variable)
  
  # Filter to top variables
  plot_data <- plot_data %>%
    filter(base_variable %in% top_vars)
  
  # Create the plot
  p1 <- ggplot(plot_data, aes(x = reorder(variable, abs(coefficient)), 
                              y = coefficient, 
                              color = outcome,
                              shape = significance)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    coord_flip() +
    facet_wrap(~outcome, scales = "free_x") +
    scale_shape_manual(values = c("ns" = 1, "p < 0.05" = 16, "p < 0.01" = 17, "p < 0.001" = 18)) +
    labs(
      title = "Detailed Coefficient Plot by Outcome",
      subtitle = paste("Top", top_n, "variables by maximum absolute coefficient"),
      x = "Variables (with polynomial contrasts)",
      y = "Coefficient Value",
      color = "Outcome",
      shape = "Significance"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9)
    )
  
  return(p1)
}

# Function to create variable importance summary plot
plot_variable_importance <- function(importance_df, top_n = 20) {
  
  # Prepare data for plotting
  plot_data <- importance_df %>%
    top_n(top_n, max_abs_coeff) %>%
    mutate(
      base_variable = reorder(base_variable, max_abs_coeff),
      significance_category = case_when(
        min_p_value < 0.001 ~ "Highly significant (p < 0.001)",
        min_p_value < 0.01 ~ "Very significant (p < 0.01)",
        min_p_value < 0.05 ~ "Significant (p < 0.05)",
        TRUE ~ "Not significant"
      )
    )
  
  # Create the importance plot
  p2 <- ggplot(plot_data, aes(x = base_variable, y = max_abs_coeff, 
                              fill = significance_category)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = paste0("n=", n_terms)), 
              hjust = -0.1, size = 3, color = "black") +
    coord_flip() +
    scale_fill_manual(values = c("Highly significant (p < 0.001)" = "#d73027",
                                 "Very significant (p < 0.01)" = "#f46d43", 
                                 "Significant (p < 0.05)" = "#fdae61",
                                 "Not significant" = "#abd9e9")) +
    labs(
      title = "Variable Importance Summary",
      subtitle = paste("Top", top_n, "variables by maximum absolute coefficient"),
      x = "Base Variables",
      y = "Maximum Absolute Coefficient",
      fill = "Significance Level",
      caption = "n = number of polynomial terms per variable"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 10)
    )
  
  return(p2)
}

# Main analysis function
analyze_multinom_importance <- function(multinom_model, top_n_detailed = 15, top_n_summary = 20) {
  
  cat("Extracting coefficients and calculating p-values...\n")
  results <- extract_variable_importance(multinom_model)
  
  cat("Aggregating polynomial contrasts...\n")
  importance <- aggregate_polynomial_importance(results)
  
  cat("Creating plots...\n")
  detailed_plot <- plot_detailed_coefficients(results, top_n_detailed)
  summary_plot <- plot_variable_importance(importance, top_n_summary)
  
  # Print summary table
  cat("\nTop 10 Most Important Variables:\n")
  print(importance %>% 
          dplyr::select(base_variable, max_abs_coeff, min_p_value, n_terms, 
                        n_significant, contrast_types) %>%
          head(10))
  
  return(list(
    results = results,
    importance = importance,
    detailed_plot = detailed_plot,
    summary_plot = summary_plot
  ))
}

# Usage example:
# analysis <- analyze_multinom_importance(model_10f$finalModel)
# 
# # Display plots
# print(analysis$summary_plot)
# print(analysis$detailed_plot)
# 
# # Access results
# View(analysis$importance)
# View(analysis$results)
