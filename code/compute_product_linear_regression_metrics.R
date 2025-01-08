library(Metrics) # For RMSE and MAPE calculations

current_path <- getActiveDocumentContext()$path
setwd(dirname(current_path))

compute_regression_metrics <- function(data, formulas, outcome, k = 10) {
  
  # Helper function for cross-validation
  cross_validate <- function(data, formula, outcome, k) {
    set.seed(123)  # For reproducibility
    n <- nrow(data)
    fold_ids <- sample(rep(1:k, length.out = n))
    
    all_preds <- numeric(n)  # Store predictions
    all_actuals <- numeric(n)  # Store actual values
    all_pred_intervals <- data.frame(lwr = numeric(n), upr = numeric(n))  # Store intervals
    
    for (fold in 1:k) {
      # Split into training and testing folds
      test_indices <- which(fold_ids == fold)
      train_data <- data[-test_indices, ]
      test_data <- data[test_indices, ]
      
      # Fit model on training data
      model <- lm(formula, data = train_data)
      
      # Predict on test data
      preds <- predict(model, test_data, interval = "prediction", level = 0.95)
      
      # Collect predictions and intervals
      
      all_preds[test_indices] <- preds[, "fit"]
      all_actuals[test_indices] <- test_data[[outcome]]
      all_pred_intervals[test_indices,] <- preds[, c("lwr", "upr")]
    }
    
    # Return combined predictions, actuals, and prediction intervals
    return(list(
      predictions = all_preds,
      actuals = all_actuals,
      pred_intervals = all_pred_intervals
    ))
  }
  
  results <- list()
  
  for (i in seq_along(formulas)) {
    formula <- formulas[[i]]
    model_name <- names(formulas)[[i]]
    
    # Fit the model on the entire dataset
    model <- lm(formula, data = data)
    
    # Training predictions
    train_preds <- predict(model, data)
    train_actuals <- data[[outcome]]
    
    # Training metrics
    train_rmse <- rmse(train_actuals, train_preds)
    train_mape <- mape(train_actuals, train_preds) * 100  # Convert to percentage
    
    # Compute 95% prediction intervals
    train_pred_int <- predict(model, data, interval = "prediction", level = 0.95)
    train_interval_width <- mean(train_pred_int[, "upr"] - train_pred_int[, "lwr"])
    train_coverage <- mean(train_actuals >= train_pred_int[, "lwr"] & 
                             train_actuals <= train_pred_int[, "upr"])
    
    # Perform 10-fold cross-validation
    cv_results <- cross_validate(data, formula, outcome, k)
    
    # Compute CV metrics from combined results
    cv_rmse <- rmse(cv_results$actuals, cv_results$predictions)
    cv_mape <- mape(cv_results$actuals, cv_results$predictions) * 100  # Convert to percentage
    cv_interval_width <- mean(cv_results$pred_intervals$upr - cv_results$pred_intervals$lwr)
    cv_coverage <- mean(cv_results$actuals >= cv_results$pred_intervals$lwr & 
                          cv_results$actuals <= cv_results$pred_intervals$upr)
    
    # Save results
    results[[model_name]] <- list(
      Training = list(
        RMSE = train_rmse,
        MAPE = train_mape,
        Coverage = train_coverage,
        Interval_Width = train_interval_width
      ),
      CrossValidation = list(
        RMSE = cv_rmse,
        MAPE = cv_mape,
        Coverage = cv_coverage,
        Interval_Width = cv_interval_width
      )
    )
  }
  
  return(results)
}

print_model_metrics <- function(metrics) {
  for (model_name in names(metrics)) {
    cat("\n", model_name, "\n")
    cat(rep("-", 40), "\n")
    
    # Extract training and CV metrics
    train_metrics <- metrics[[model_name]]$Training
    cv_metrics <- metrics[[model_name]]$CrossValidation
    
    # Combine into a data frame
    results_table <- data.frame(
      RMSE = c(train_metrics$RMSE, cv_metrics$RMSE),
      MAPE = c(train_metrics$MAPE, cv_metrics$MAPE),
      Coverage = c(train_metrics$Coverage, cv_metrics$Coverage),
      Interval_Width = c(train_metrics$Interval_Width, cv_metrics$Interval_Width),
      row.names = c("Training", "CV")
    )
    
    # Print the table
    print(results_table)
  }
}

load('../data/census_ACS_PEP_WP_wDensity_and2018_01022024.RData')
model_formulas <- list(
  FULL = census ~ acs + pep + wp, 
  ACS = census ~ acs, 
  PEP = census ~ pep,
  WP = census ~ wp,
  ACS_WP = census ~ acs + wp
)

metrics <- compute_regression_metrics(df, model_formulas, outcome = 'census')

#print(metrics)
print_model_metrics(metrics)

