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

