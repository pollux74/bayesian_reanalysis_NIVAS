# =============================================
# Project : Bayesian reanalysis of NIVAS
# Script : functions.R
# Author : A. Naudet--Lasserre
# Date last modification : 2025-09-01
# ============================================

# Priors ------------------------------------------------------------------

# Function

#' Title Calcul l'écart type d'une loi normale en fonction de sa moyenne et de la probabilité X>seuil
#'
#' @param mean Moyenne de la loi normale
#' @param target_p Probabilité de X > seuil
#' @param threshold Seuil choisi (défault = 0)
#'
#' @returns sd de la loi normale
calibrate_prior_sd <- function(mean, target_p, threshold = 0) {
  z <- qnorm(target_p)
  sd <- (threshold - mean) / z
  return(sd)
}
fig_path <- function(filename) {
  here::here("outputs", "figures", filename)
}
#' Plot prior distributions with specific color mapping
#'
#' @param priors_df Data.frame with prior_type, prior_mean, prior_sd columns
#' @param x_limits Numeric vector of length 2 for x-axis limits
#' @param title Character string for plot title (optional)
#' @param show_legend Logical, whether to show legend
#' @return ggplot object
#' @export
plot_priors_distributions <- function(priors_df,
                                      x_limits = c(-5, 5),
                                      title = NULL,
                                      show_legend = TRUE) {
  # Input validation
  required_cols <- c("prior_type", "prior_mean", "prior_sd")
  if (!all(required_cols %in% names(priors_df))) {
    stop("Dataframe must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Specific color mapping for prior types
  color_mapping <- c(
    "Enthusiastic" = "#04a5e5", # Peach
    "Pessimistic" = "#8839ef", # Teal
    "Minimally informative" = "#6c6f85", # Catppuccin text
    "Skeptical" = "#1e66f5" # Green
  )

  # Create x values grid for curves
  x_values <- seq(x_limits[1], x_limits[2], length.out = 1000)

  # For second x axis
  # Breaks pour l'axe log-OR (5 valeurs équidistantes)
  log_or_breaks <- c(-4, -2, 0, 2, 4)

  # Breaks pour l'axe OR (valeurs claires et parlantes)
  or_breaks_clear <- c(0.1, 0.5, 1, 2, 10)


  # Calculate densities for each prior
  density_data <- priors_df %>%
    mutate(prior_id = row_number()) %>%
    rowwise() %>%
    do({
      tibble(
        x = x_values,
        density = dnorm(x_values, mean = .$prior_mean, sd = .$prior_sd),
        prior_type = .$prior_type,
        prior_id = .$prior_id
      )
    }) %>%
    ungroup()

  # Build the plot
  p <- ggplot(density_data, aes(x = x, y = density, color = prior_type)) +

    # Density lines
    geom_line(linewidth = 0.8) +

    # Manual color scale with specific mapping
    scale_color_manual(
      values = color_mapping,
      name = "Prior type"
    ) +
    scale_x_continuous(
      limits = x_limits,
      expand = c(0, 0),
      name = "Log-odds ratio",
      breaks = log_or_breaks,
      sec.axis = sec_axis(
        transform = ~ exp(.),
        name = "Odds ratio",
        breaks = or_breaks_clear,
        labels = c("0.1", "0.5", "1", "2", "10")
      )
    ) +

    # Y-axis configuration
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      name = "Density"
    ) +

    # Publication theme
    theme_minimal(base_size = 12) +
    theme(
      # All text in black
      text = element_text(color = "black"),

      # Axes
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 11, color = "black", face = "bold"),

      # Grid (minimal)
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#e6e9ef", linewidth = 0.3),
      panel.grid.major.x = element_blank(),

      # Margins
      plot.margin = margin(15, 15, 15, 15),

      # Legend in top-right corner
      legend.position.inside  = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.background = element_rect(
        fill = "white",
        color = "gray80",
        linewidth = 0.3
      ),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.4, "cm"),

      # Plot border
      panel.border = element_rect(
        color = "gray80",
        fill = NA,
        linewidth = 0.5
      )
    )

  # Hide legend if requested
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }

  # Add title if specified
  if (!is.null(title)) {
    p <- p +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
  }

  return(p)
}
calc_equivalent_sample_size <- function(log_or_mean,
                                        log_or_sd = NULL,
                                        p_control,
                                        p_treatment,
                                        prob_benefit = NULL) {
  # If prob_benefit is provided, calculate log_or_sd
  if (!is.null(prob_benefit)) {
    if (prob_benefit <= 0 | prob_benefit >= 1) {
      stop("prob_benefit must be between 0 and 1")
    }
    z_score <- qnorm(1 - prob_benefit)
    log_or_sd <- abs(log_or_mean) / z_score
  }

  # Check inputs
  if (is.null(log_or_sd)) {
    stop("Either log_or_sd or prob_benefit must be provided")
  }

  # Calculate variance per patient (from theoretical formula)
  variance_per_patient <- (1 / (p_control * (1 - p_control))) +
    (1 / (p_treatment * (1 - p_treatment)))

  # Prior variance
  prior_variance <- log_or_sd^2

  # Equivalent sample size (per group)
  n_per_group <- variance_per_patient / prior_variance
  n_total <- 2 * n_per_group

  # Calculate actual P(log-OR > 0) for verification
  actual_prob_benefit <- pnorm(0,
    mean = log_or_mean, sd = log_or_sd,
    lower.tail = FALSE
  )

  # Return results
  return(list(
    n_per_group = round(n_per_group, 1),
    n_total = round(n_total, 1),
    prior_mean = round(log_or_mean, 3),
    prior_sd = round(log_or_sd, 3),
    prior_variance = round(prior_variance, 3),
    prob_benefit = round(actual_prob_benefit, 3),
    variance_per_patient = round(variance_per_patient, 3)
  ))
}



# Analyse -----------------------------------------------------------------

#' Fit bayesian logistic regression model on a list of priors
#'
#' @param data = dataframe
#' @param outcome = Variable to be explained
#' @param explic = predictors
#' @param intercept = dataframe that contain colon mean and sd specifying the intercept's prior, can be left empty
#' @param verbose  Diagnostic messages
#' @returns model_list -> list of models by priors
#' @export
#'
#' @examples model_list_reint_inter2 <- fit_brm_opti(data = df, outcome = "reintubation", explic = "group", intercept = P_reintubation_inter2)
fit_brm_opti <- function(data, outcome, explic, intercept = NULL) {
  # Récupération dynamique des priors
  prior_var_name <- paste0("priors_", outcome)
  P <- get(prior_var_name)

  if (missing(P)) {
    stop("Prior variable not found in the global environment.")
  }

  cat("Using prior :", prior_var_name, "\n")

  model_list <- vector("list", nrow(P))

  for (i in seq_len(nrow(P))) {
    # Define priors
    if (!is.null(intercept)) {
      print("Intercept detected")
      # VERIFY intercept dataframe structure
      if (!all(c("prior_mean", "prior_sd") %in% names(intercept))) {
        stop("Intercept dataframe must contain 'prior_mean' and 'prior_sd' columns")
      }

      if (any(is.na(c(intercept$prior_mean, intercept$prior_sd)))) {
        stop("Intercept prior_mean and prior_sd cannot be NA")
      }
      prior_l <- c(
        set_prior(
          paste0("normal(", P[i, "prior_mean"], ", ", P[i, "prior_sd"], ")"),
          class = "b"
        ),
        set_prior(
          paste0("normal(", intercept$prior_mean, ", ", intercept$prior_sd, ")"),
          class = "Intercept"
        )
      )
      cat("Prior defined with:")
      print(prior_l[2, ])
    } else {
      print("No prior on intercept using default")
      prior_l <- set_prior(
        paste0("normal(", P[i, "prior_mean"], ", ", P[i, "prior_sd"], ")"),
        class = "b"
      )
      print(prior_l[1, ])
    }

    formula <- as.formula(paste0(outcome, " ~ ", explic))
    
    if (is.factor(data[[outcome]])){
      print("outcome is a factor")
      fam <- bernoulli(link = 'logit')
    }
    if (is.numeric(data[[outcome]])){
      print("outcome is a numeric")
      fam <- gaussian()
    }
    # Fit model
    model_list[[i]] <- tryCatch(
      {
        model <- brm(
          formula,
          data = data,
          family = fam,
          prior = prior_l,
          iter = 4000,
          warmup = 1000,
          seed = 42,
          refresh = 0,
          chains = 4,
          cores = 4,
          backend = "cmdstanr" # required for threading
        )
      },
      error = function(e) {
        cat("Error fitting model with prior index ", i, ": ", e$message, "\n")
        return(NULL)
      }
    )
  }
  print(paste0(">> Model fit finished", outcome, explic))
  return(model_list)
}

estimate_ARR <- function(draws, logOR, model = NULL, df = NULL) {
  if (!is.null(model)) {
    # ARR marginal avec les données d'entraînement
    cat("Using marginal ARR calculation...\n")
    
    model_data <- model$data
    group_levels <- levels(model_data$group)
    
    # Tous les patients en contrôle
    df_control <- model_data
    df_control$group <- factor(group_levels[1], levels = group_levels)
    
    # Tous les patients en traitement
    df_treatment <- model_data
    df_treatment$group <- factor(group_levels[2], levels = group_levels)
    
    # Prédictions
    pred_control <- posterior_epred(model, newdata = df_control)
    pred_treatment <- posterior_epred(model, newdata = df_treatment)
    
    # ARR marginal
    risk_control <- rowMeans(pred_control)
    risk_treatment <- rowMeans(pred_treatment)
    arr <- risk_control - risk_treatment
    arr_hdi <- hdi(arr, .width = 0.95)
    
    # Diagnostics
    cat("Marginal risk control:", round(mean(risk_control) * 100, 1), "%\n")
    cat("Marginal risk treatment:", round(mean(risk_treatment) * 100, 1), "%\n")
    cat("Marginal ARR:", round(mean(arr) * 100, 1), "%\n")
    
    return(list(arr = arr, arr_hdi = arr_hdi))
  } else {
    # ARR original pour modèles bivariés
    cat("Using original ARR calculation...\n")
    intercept <- draws[["b_Intercept"]]
    p_control <- plogis(intercept)
    p_treatment <- plogis(intercept + logOR)
    arr <- p_control - p_treatment
    arr_hdi <- hdi(arr, .width = 0.95)
    return(list(arr = arr, arr_hdi = arr_hdi))
  }
}


#' Extract draws from a single brmsfit model
#'
#' @param model brmsfit object
#' @param group_name Group variable name (without "b_" prefix)
#' @param prior_name Name to identify this set of draws
#' @return data.frame with draws (log_or, arr) and identification columns

get_draws <- function(model, group_name, prior_name = NULL) {
  # Input validation
  stopifnot(
    "model must be brmsfit" = inherits(model, "brmsfit"),
    "group_name must be character" = is.character(group_name)
  )

  # Extract posterior draws
  draws <- as_draws_df(model)
  var_name <- paste0("b_", group_name)

  # Check variable exists
  if (!var_name %in% names(draws)) {
    stop(glue::glue("Variable '{var_name}' not found in model"))
  }

  # Get logOR and compute ARR
  logOR <- draws[[var_name]]
  draws_data <- estimate_ARR(draws, logOR)

  # Return structured draws
  result <- data.frame(
    draw_id = seq_along(draws_data$arr),
    log_or = logOR,
    arr = draws_data$arr,
    stringsAsFactors = FALSE
  )

  # Add prior identification if provided
  if (!is.null(prior_name)) {
    result$prior <- prior_name
  }

  return(result)
}

#' Compute summary metrics from draws
#'
#' @param draws_df data.frame from get_draws()
#' @param prior_name Name to identify this summary
#' @return data.frame with summary metrics

compute_summary_metrics <- function(draws_df, prior_name = NULL) {
  # Extract vectors for calculations
  logOR <- draws_df$log_or
  arr_values <- draws_df$arr

  # Compute HDI for ARR
  arr_hdi <- tidybayes::hdi(arr_values, .width = 0.95)
  cat(mean(arr_values > 0.02))
  # Summary metrics
  result <- data.frame(
    median_log_or = median(logOR),
    median_or = exp(median(logOR)),
    median_arr = median(arr_values),
    ci_lower = exp(tidybayes::hdi(logOR, .width = 0.95)[[1]]),
    ci_upper = exp(tidybayes::hdi(logOR, .width = 0.95)[[2]]),
    p_beneficial = mean(logOR < 0),
    arr_ci_lower = arr_hdi[[1]],
    arr_ci_upper = arr_hdi[[2]],
    arr_gt_2pct = mean(arr_values > 0.02),
    arr_gt_5pct = mean(arr_values > 0.05),
    arr_gt_10pct = mean(arr_values > 0.10),
    stringsAsFactors = FALSE
  )

  # Add prior identification if provided
  if (!is.null(prior_name)) {
    result$prior <- prior_name
  }

  return(result)
}

#' Main function: Compute posterior metrics using modular approach
#'
#' @param model_list List of brmsfit objects
#' @param prior_names Vector of corresponding prior names
#' @param group_name Group variable name (without "b_" prefix)
#' #' @return List with 'summary' and 'draws' dataframes

compute_posterior_metrics <- function(model_list, prior_names, group_name) {
  # Input validation
  stopifnot(
    "model_list must be a list" = is.list(model_list),
    "prior_names and model_list must have same length" =
      length(prior_names) == length(model_list),
    "group_name must be character" = is.character(group_name)
  )

  # Extract all draws using lapply
  all_draws <- lapply(seq_along(model_list), function(i) {
    tryCatch(
      {
        get_draws(model_list[[i]], group_name, prior_names[i])
      },
      error = function(e) {
        warning(glue::glue("Error extracting draws for prior '{prior_names[i]}': {e$message}"))
        return(NULL)
      }
    )
  })

  # Filter out failed extractions
  all_draws <- all_draws[!sapply(all_draws, is.null)]

  if (length(all_draws) == 0) {
    stop("No draws could be extracted from any model")
  }

  # Compute summary metrics using lapply
  all_summaries <- lapply(seq_along(all_draws), function(i) {
    draws_df <- all_draws[[i]]
    prior_name <- unique(draws_df$prior)[1] # Get prior name from draws
    compute_summary_metrics(draws_df, prior_name)
  })

  # Combine results into single dataframes
  draws_combined <- do.call(rbind, all_draws)
  summaries_combined <- do.call(rbind, all_summaries)

  return(list(
    summary = summaries_combined,
    draws = draws_combined,
    metadata = list(
      n_models = length(all_draws),
      group_variable = group_name,
      timestamp = Sys.time()
    )
  ))
}


###### Plot Functions ######


plot_posterior_OR <- function(posterior_draws, outcome) {
  posterior_draws$prior_name <- factor(posterior_draws$prior_name,
    levels = levels(posterior_draws$prior_name)
  )

  p <- ggplot(posterior_draws, aes(x = exp(group_name), y = prior_name)) +
    ggdist::stat_halfeye(
      aes(fill = after_stat(x < 1)),
      point_interval = "mode_hdi",
      .width = 0.95,
      slab_alpha = 0.5,
      interval_size = 1.2,
      slab_linewidth = 0.6,
      linewidth = 3,
      point_size = 2,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("#81c8be", "#99d1db"),
      labels = c("OR ≥ 1", "OR < 1"),
      name = "Direction de l'effet"
    ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "gray20",
      linewidth = 0.6
    ) +
    labs(
      x = "Odds Ratio",
      y = NULL
    ) +
    scale_x_log10() +
    theme_scientific_publication() +
    coord_cartesian(clip = "off")

  return(p)
}

plot_risk_reduction <- function(rr_data, type = "ARR", outcome = "reintubation") {
  # Threshold determination based on measure type
  threshold <- ifelse(type == "RR", 1, 0)

  # Appropriate X-axis label
  x_label <- ifelse(type == "RR",
    "Risk Ratio (RR)",
    "Absolute risk reduction (ARR)"
  )

  # Get unique Prior levels for color generation
  prior_levels <- unique(rr_data$Prior)

  # Create plot with standardized theme
  p <- ggplot(rr_data, aes(x = .data[[type]], y = Prior, fill = Prior)) +
    ggdist::stat_halfeye(
      aes(fill = after_stat(x > 0)),
      point_interval = "mode_hdi",
      .width = 0.95,
      slab_alpha = 0.5,
      interval_size = 1.2,
      slab_linewidth = 0.6,
      linewidth = 3,
      point_size = 2,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("#81c8be", "#99d1db"),
      labels = c("OR ≥ 1", "OR < 1"),
      name = "Direction de l'effet"
    ) +
    # Reference line
    geom_vline(
      xintercept = threshold,
      linetype = "dashed",
      color = "gray20",
      linewidth = 0.6
    ) +
    # Labels
    labs(
      x = x_label,
      y = NULL
    ) +
    # Apply standardized theme
    theme_scientific_publication() +
    coord_cartesian(clip = "off")

  return(p)
}



# Main analysis function with corrected marginal ARR
analyse_outcome <- function(outcome, df, predicteur, m_list, ARR = TRUE, verbose = FALSE) {
  # Verbose control system
  cat_verbose <- function(...) {
    if (verbose) cat(...)
  }

  # Count words in predicteur argument to determine model complexity
  predicteur_words <- trimws(unlist(strsplit(predicteur, "\\s+")))
  is_bivariate <- length(predicteur_words) == 1

  cat_verbose("=== Analysis of", outcome, "===\n")
  cat_verbose("Data source:", deparse(substitute(df)), "\n")
  cat_verbose("Model type:", ifelse(is_bivariate, "BIVARIATE", "MULTIVARIATE"), "\n")
  cat_verbose("ARR calculation:", ifelse(ARR, "ENABLED", "DISABLED"), "\n")

  # Extract metadata
  p_names <- get(paste0("priors_", outcome))$prior_type
  group_name <- rownames(fixef(m_list[[1]]))[2]

  cat_verbose("Priors:", paste(p_names, collapse = ", "), "\n")
  cat_verbose("Groups:", paste(rownames(fixef(m_list[[1]])), collapse = ", "), "\n")
  # Extract posterior draws from each model
  cat_verbose("Extracting posterior draws...\n")
  draws_list <- lapply(seq_along(m_list), function(i) {
    get_draws(m_list[[i]],
      group_name = group_name,
      prior_name = p_names[i]
    )
  })

  # Calculate summary metrics for each model (still using original method for OR calculations)
  cat_verbose("Computing summary metrics...\n")
  posterior_metrics <- lapply(seq_along(draws_list), function(i) {
    compute_summary_metrics(draws_list[[i]], p_names[i])
  })

  # Table formatting based on ARR parameter only
  cat_verbose("Formatting publication table...\n")

  if (ARR) {
    # Format table WITH ARR columns (will be updated later with correct values)
    formatted_metrics <- format_publication_table(
      posterior_metrics,
      p_names,
      paste0(outcome)
    )

    # HTML table with ARR columns
    table_obj <- formatted_metrics %>%
      kable(
        format = "html",
        escape = FALSE,
        col.names = c(
          "Prior Belief",
          "Posterior Median OR (95% CrI)",
          "P(OR < 1)%",
          "Posterior Median ARR, % (95% CrI)",
          "P(ARR > 2%)%",
          "P(ARR ≥5%)%",
          "P(ARR ≥10%)%"
        )
      ) %>%
      kable_styling(
        bootstrap_options = c("striped", "hover"),
        full_width = FALSE,
        position = "center"
      )
  } else {
    # Format table WITHOUT ARR columns
    formatted_metrics <- format_publication_table_OR_only(
      posterior_metrics,
      p_names,
      paste0(outcome)
    )

    # HTML table without ARR columns
    table_obj <- formatted_metrics %>%
      kable(
        format = "html",
        escape = FALSE,
        col.names = c(
          "Prior Belief",
          "Posterior Median OR (95% CrI)",
          "P(OR < 1)%"
        )
      ) %>%
      kable_styling(
        bootstrap_options = c("striped", "hover"),
        full_width = FALSE,
        position = "center"
      )
  }

  # ARR calculation with marginal approach when requested
  if (ARR) {
    cat_verbose("Computing risk reduction estimates...\n")
    arr_data <- lapply(seq_along(m_list), function(i) {
      original_draws <- as_draws_df(m_list[[i]])
      logOR <- original_draws[[paste0("b_", group_name)]]

      # Pass model and dataframe for marginal calculation
      cat_verbose("ici")
      estimate_ARR(
        draws = original_draws,
        logOR = logOR,
        model = m_list[[i]], # Model object for marginal calculation
        df = df # Dataframe (though we'll use model$data internally)
      )
    })

    # Update the formatted_metrics table with correct marginal ARR values
    cat_verbose("Updating table with marginal ARR values...\n")
    for (i in seq_along(arr_data)) {
      # Extract ARR statistics
      arr_values <- arr_data[[i]]$arr * 100 # Convert to percentage
      arr_median <- median(arr_values)
      arr_hdi <- quantile(arr_values, c(0.025, 0.975))

      # Update ARR column in formatted table
      formatted_metrics[i, "ARR_with_CI"] <- sprintf(
        "%.1f%% (%.1f%% to %.1f%%)",
        arr_median, arr_hdi[1], arr_hdi[2]
      )

      # Update probability columns
      formatted_metrics[i, "P_ARR_gte_2pct"] <- paste0(round(mean(arr_values >= 2) * 100), "%")
      formatted_metrics[i, "P_ARR_gte_5pct"] <- paste0(round(mean(arr_values >= 5) * 100), "%")
      formatted_metrics[i, "P_ARR_gte_10pct"] <- paste0(round(mean(arr_values >= 10) * 100), "%")
    }

    # Update HTML table with corrected values
    table_obj <- formatted_metrics %>%
      kable(
        format = "html",
        escape = FALSE,
        col.names = c(
          "Prior Belief",
          "Posterior Median OR (95% CrI)",
          "P(OR < 1)%",
          "Posterior Median ARR, % (95% CrI)",
          "P(ARR > 2%)%",
          "P(ARR ≥5%)%",
          "P(ARR ≥10%)%"
        )
      ) %>%
      kable_styling(
        bootstrap_options = c("striped", "hover"),
        full_width = FALSE,
        position = "center"
      )

    # Prepare data for ARR visualization
    risk_data_combined <- bind_rows(
      lapply(seq_along(arr_data), function(i) {
        data.frame(
          ARR = arr_data[[i]]$arr,
          Prior = p_names[i]
        )
      })
    )
  } else {
    cat_verbose("Skipping ARR calculation (ARR = FALSE)\n")
    arr_data <- NULL
    risk_data_combined <- NULL
  }

  # Prepare data for OR visualization (always needed)
  df_all <- bind_rows(
    lapply(seq_along(m_list), function(i) {
      draws <- as_draws_df(m_list[[i]])
      data.frame(
        group_name = draws[[paste0("b_", group_name)]],
        prior_name = p_names[i]
      )
    })
  )

  df_all$prior_name <- factor(df_all$prior_name, levels = p_names)

  # Plot generation based on ARR parameter
  cat_verbose("Generating visualizations...\n")

  # Always create OR plot
  fig <- plot_posterior_OR(df_all, paste0(outcome))

  if (ARR && !is.null(risk_data_combined)) {
    # Create combined plot (OR + ARR) when ARR is requested
    fig <- fig + theme(plot.margin = margin(r = 0))

    f2.2 <- plot_risk_reduction(risk_data_combined, "ARR") +
      theme(
        plot.margin = margin(l = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )

    combined_plot <- fig + f2.2 +
      plot_layout(ncol = 2, widths = c(1, 1))

    cat_verbose("Created combined plot (OR + ARR)\n")
  } else {
    # Use only OR plot when ARR is not requested
    combined_plot <- fig
    cat_verbose("Created OR-only plot (ARR = FALSE)\n")
  }

  # Clean console output of main results
  cat("\n=== RESULTS FOR", toupper(outcome), "===\n")
  cat("Model type:", ifelse(is_bivariate, "BIVARIATE", "MULTIVARIATE"), "\n")
  cat("ARR calculated:", ifelse(ARR, "YES (marginal approach)", "NO"), "\n")
  print(formatted_metrics)
  cat("\n")

  # Return results with updated structure
  result <- list(
    plot = combined_plot,
    table = table_obj,
    raw = list(
      formatted_metrics = formatted_metrics,
      model_type = ifelse(is_bivariate, "bivariate", "multivariate"),
      arr_calculated = ARR,
      arr_approach = ifelse(ARR, ifelse(is_bivariate, "simple", "marginal"), "none"),
      arr_data = if (ARR) arr_data else NULL
    )
  )

  cat_verbose("Analysis completed.\n")

  # Save outputs
  ggsave(
    plot = combined_plot,
    filename = fig_path(paste0(outcome, deparse(substitute(m_list)), "_combined.png.tiff")),
    width = 7, # Standard journal width in inches
    height = 5, # Height in inches (ratio ~1.4)
    dpi = 600, # 600 dpi for print quality
    compression = "lzw", # Lossless compression
    device = "tiff" # Explicit format specification
  )

  save_kable(table_obj, here("outputs", "tables", paste0(outcome, deparse(substitute(m_list)), ".html")))

  return(result)
}



#' Format publication table for OR-only (multivariate case)
#'
#' @param posterior_metrics list of metrics from compute_summary_metrics
#' @param prior_names vector of prior names
#' @param outcome_name name of outcome variable
#' @return data.frame with formatted OR results only

format_publication_table_OR_only <- function(posterior_metrics, prior_names, outcome_name) {
  # Input validation
  if (length(posterior_metrics) != length(prior_names)) {
    stop("Number of elements in posterior_metrics must match number of priors")
  }

  # Combine metrics into single dataframe if it's a list
  if (is.list(posterior_metrics) && !is.data.frame(posterior_metrics)) {
    metrics_df <- do.call(rbind, posterior_metrics)
  } else {
    metrics_df <- posterior_metrics
  }

  # Create formatted table WITHOUT ARR columns
  formatted_table <- data.frame(
    "Prior_Belief" = prior_names,
    "OR_with_CI" = sprintf(
      "%.2f (%.2f to %.2f)",
      metrics_df$median_or,
      metrics_df$ci_lower,
      metrics_df$ci_upper
    ),
    "P_OR_less_than_1" = sprintf(
      "%.0f%%",
      metrics_df$p_beneficial * 100
    )
  )

  # Add outcome name as attribute
  attr(formatted_table, "outcome") <- outcome_name
  attr(formatted_table, "type") <- "OR_only"
  return(formatted_table)
}



# Bivariate analysis (will compute ARR)
# result_biv <- analyse_outcome("reintubation", df, "treatment_group", models_list)

# Multivariate analysis (will NOT compute ARR)
# result_multi <- analyse_outcome("reintubation", df, "treatment_group age comorbidities", models_list)

# Check what was calculated
# result_biv$raw$arr_calculated    # Should be TRUE
# result_multi$raw$arr_calculated  # Should be FALSE

#' Plot trace plot and autocorrelogram for brm model list
#'
#' @param model_list : the brm model list
#' @param title
#'
#' @returns plot with autocorelograms and traceplot for a model list
#'
#' @examples check_convergence(model_list_reint)
check_convergence <- function(model_list) {
  print("Checking convergence of models\n")
  group_name <- rownames(fixef(model_list[[1]]))[2] # name of outcome in brm model
  pars_name <- paste0("b_", group_name) # Proper format to call

  print(paste(pars_name, "pars"))
  plots_traces <- list()
  plots_acf <- list()

  print("Creating traces plots & autocorrelogram")
  # for loop, generating plots for model in the list
  for (i in seq_along(model_list)) {
    trace_plot <- mcmc_trace(as_draws_df(model_list[[i]]), pars = pars_name) +
      ggtitle(paste("Trace plot "))
    acf_plot <- mcmc_acf(as_draws_df(model_list[[i]]), pars = pars_name) +
      ggtitle(paste("ACF"))
    plots_traces[[i]] <- trace_plot
    plots_acf[[i]] <- acf_plot
  }
  print("Ploting trace plotes")
  print("Loop successful, plots created")
  # Combine all plots (traces first, then ACF)
  all_plots <- c(plots_traces, plots_acf)
  total_plots <- length(model_list)

  ncol_grid <- total_plots # 1 trace, 1 ACF side by side
  nrow_grid <- 2

  converg <- grid.arrange(grobs = all_plots, nrow = nrow_grid, ncol = ncol_grid)
  print("Arranged plot created")
  path <- fig_path(paste0(deparse(substitute(model_list)), "_convergence.tiff"))
  ggsave(
    plot = converg,
    filename = path,
    width = 9, # Width in inches (standard journal)
    height = 5, # Height in inches (ratio ~1.4)
    dpi = 600 # 600 dpi for print quality
    # compression = "lzw", # Lossless compression
    # device = "tiff" # Explicitly specify the format
  )
  print(paste0("Plot saved in", path))
  return(converg)
}

# Theming -----------------------------------------------------------------

# Création d'un thème personnalisé pour publications scientifiques
theme_scientific_publication <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      # Configuration du texte pour impression optimale
      text = element_text(color = "black"),

      # Configuration des axes
      axis.text = element_text(size = base_size - 2, color = "black"),
      axis.title.x = element_text(
        size = base_size - 1,
        face = "bold",
        margin = margin(t = 10)
      ),
      axis.text.y = element_text(face = "bold", hjust = 0),

      # Élimination des éléments de grille non essentiels
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.ticks.x = element_line(linewidth = 0.4, color = "grey80"),
      axis.line = element_line(linewidth = 0.4, color = "grey80"),

      # Ajustement des marges pour optimiser l'espace
      plot.margin = margin(15, 15, 15, 15),

      # Suppression de la légende
      legend.position = "none",


      # Titre et sous-titre
      plot.title = element_text(
        size = base_size + 2,
        face = "bold",
        hjust = 0,
        margin = margin(b = 10)
      ),
      plot.subtitle = element_text(
        size = base_size,
        color = "gray30",
        margin = margin(b = 15)
      ),

      # Paramètres additionnels pour une meilleure lisibilité
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = NA)
    )
}

# Création d'une fonction pour les paramètres standard de stat_halfeye
halfeye_scientific_style <- function() {
  list(
    point_interval = "mode_hdi",
    .width = 0.95,
    slab_alpha = 0.5,
    interval_size = 1.2,
    slab_linewidth = 0.6,
    show.legend = FALSE,
    linewidth = 3,
    point_size = 2
  )
}

# Création d'une fonction pour les échelles logarithmiques standard
scale_log_scientific <- function(limits = c(0.2, 1.5), custom_breaks = NULL) {
  if (is.null(custom_breaks)) {
    custom_breaks <- c(0.2, 0.3, 0.5, 0.8, 1)
  }

  custom_labels <- as.character(custom_breaks)

  scale_x_log10(
    limits = limits,
    breaks = custom_breaks,
    labels = custom_labels
  )
}
#' Format publication table with ARR columns
#'
#' @param posterior_metrics List of data.frames or single data.frame with metrics
#' @param prior_names Character vector of prior names
#' @param outcome_name Character string for outcome name
#' @return Data.frame with formatted results
#' @export
format_publication_table <- function(posterior_metrics, prior_names, outcome_name) {
  # Combine data.frames if input is a list
  if (is.list(posterior_metrics)) {
    df <- do.call(rbind, posterior_metrics)
  } else {
    df <- posterior_metrics
  }

  # Create formatted table with ARR columns
  data.frame(
    "Prior_Belief" = prior_names,
    "OR_with_CI" = sprintf("%.2f (%.2f to %.2f)", df$median_or, df$ci_lower, df$ci_upper),
    "P_OR_less_than_1" = sprintf("%.0f%%", df$p_beneficial * 100),
    "ARR_with_CI" = sprintf(
      "%.1f%% (%.1f%% to %.1f%%)",
      df$median_arr * 100, df$arr_ci_lower * 100, df$arr_ci_upper * 100
    ),
    "P_ARR_gte_2pct" = sprintf("%.0f%%", df$arr_gt_2pct * 100),
    "P_ARR_gte_5pct" = sprintf("%.0f%%", df$arr_gt_5pct * 100),
    "P_ARR_gte_10pct" = sprintf("%.0f%%", df$arr_gt_10pct * 100)
  )
}

#' Format publication table without ARR columns
#'
#' @param posterior_metrics List of data.frames or single data.frame with metrics
#' @param prior_names Character vector of prior names
#' @param outcome_name Character string for outcome name
#' @return Data.frame with formatted results (OR only)
#' @export
format_publication_table_OR_only <- function(posterior_metrics, prior_names, outcome_name) {
  # Combine data.frames if input is a list
  if (is.list(posterior_metrics)) {
    df <- do.call(rbind, posterior_metrics)
  } else {
    df <- posterior_metrics
  }

  # Create formatted table without ARR columns
  data.frame(
    "Prior_Belief" = prior_names,
    "OR_with_CI" = sprintf("%.2f (%.2f to %.2f)", df$median_or, df$ci_lower, df$ci_upper),
    "P_OR_less_than_1" = sprintf("%.0f%%", df$p_beneficial * 100)
  )
}

