library(shiny)
library(DT)
library(nphRshiny)
library(tidyverse)
library(mvtnorm)
library(gsDesign)
library(plotly)
library(ggplot2)
library(SAMprior)
library(shinyjs)
library(BayesianHybridDesign)
library(future)
library(promises)

plan(multisession)

shinyServer(function(input, output, session) {

# ===========================================================================
# Term Definitions
# ===========================================================================

  # Terms & Definitions for the Dynamic Power Prior (DPP) Library tab
  dpp_term_definitions <- list(
    "pt" = "Response rate for experimental arm in current study.",
    "nt" = "Number of patients in experimental arm in current study.",
    "pc" = "Response rate for control arm in current study.",
    "nc" = "Number of patients in control arm in current study.",
    "pc.calib" = "Required for calibration if tau is not provided. Response rate for control arm in current study for calibration. Usually, <code>pc.calib = pch</code>.",
    "pch" = "Response rate for control treatment in historical study.",
    "nche" = "Equivalent number of patients borrowed from historical study.",
    "nch" = "Total number of patients in historical control.",
    "alpha" = "A scalar. One sided type I error rate. Required for calibration if tau is not provided.",
    "tau" = "Calibrated threshold for statistical significance. If tau is not provided, it will be calculated by calibration to type I error alpha.",
    "a0c" = "Hyperprior for control response rate beta(a0c, b0c).",
    "b0c" = "Hyperprior for control response rate beta(a0c, b0c).",
    "a0t" = "Hyperprior for experimental response rate beta(a0t, b0t).",
    "b0t" = "Hyperprior for experimental response rate beta(a0t, b0t).",
    "delta_threshold" = "Borrow when <code>abs(pc_hat (current study) - pch) <= delta_threshold</code>.",
    "method" = "Method for dynamic borrowing, 'Empirical Bayes', 'Bayesian p', 'Generalized BC', 'JSD'.",
    "theta" = "A parameter with a range of (0, 1), and applicable to method: 'Generalized BC'.",
    "eta" = "A parameter with a range of (0, infty), and applicable to method: 'Bayesian p', 'Generalized BC', 'JSD'. 'Generalized BC' method requires two parameters theta and eta.",
    "datamat" = "A matrix with dimension <code>nsim * 2</code> containing the pre-simulated data for the study treatment (1st column) and control (1st column) groups, respectively. If not supplied, binomial random Monte Carlo samples will be generated in the function.",
    "w0" = "Prior power parameters <code>w</code>. If not specified (default), <code>w_d</code> is calculated by the specified method for dynamic borrowing.",
    "nsim" = "Number of replications to calculate power.",
    "seed" = "Seed for simulations."
  )

  # DPP Analysis terms Library
  dpp_analysis_term_definitions <- list(
    "w" = "Borrowing weight. A number between 0 and 1 that shows how much information from historical data is being used.",
    "phat_pt_larger_pc" = "The posterior probability that the experimental arm's response rate is truly higher than the control arm's response rate.",
    "apost_c_trial" = "The posterior alpha parameter for the control arm's response rate, using only data from the current study.",
    "bpost_c_trial" = "The posterior beta parameter for the control arm's response rate, using only data from the current study.",
    "apost_c_hca" = "The posterior alpha parameter for the hybrid control arm's response rate after combining current and historical data.",
    "bpost_c_hca" = "The posterior beta parameter for the hybrid control arm's response rate after combining current and historical data.",
    "apost_t" = "The posterior alpha parameter for the experimental arm's response rate.",
    "bpost_t" = "The posterior beta parameter for the experimental arm's response rate."
  )

  # Terms & Definitions for the Fisher Testing Library tab
  fisher_term_definitions <- list(
    "Yc (fisher)" = "Number of subjects with response in experimental arm for <code>fisher</code> function.",
    "nc (fisher)" = "Number of subjects in control arm for <code>fisher</code> function.",
    "Yt (fisher)" = "Number of subjects with response in control arm for <code>fisher</code> function.",
    "nt (fisher)" = "Number of subjects in experimental arm for <code>fisher</code> function.",
    "alternative (fisher)" = "Type of alternative hypothesis (e.g., 'greater') for <code>fisher</code> function.",
    "pc (fisher.bound)" = "Response rate for control arm for <code>fisher.bound</code> function.",
    "nc (fisher.bound)" = "Number of patients in control arm for <code>fisher.bound</code> function.",
    "nt (fisher.bound)" = "Number of patients in experimental arm for <code>fisher.bound</code> function.",
    "alpha (fisher.bound)" = "P-value threshold for significance. Alpha must be one-sided. Default 0.1 for <code>fisher.bound</code> function.",
    "pt (fisher.power)" = "Probability of success in experimental arm for <code>fisher.power</code> function.",
    "nt (fisher.power)" = "Number of subject in experimental arm for <code>fisher.power</code> function.",
    "pc (fisher.power)" = "Probability of success in control arm for <code>fisher.power</code> function.",
    "nc (fisher.power)" = "Number of subject in control arm for <code>fisher.power</code> function.",
    "alpha (fisher.power)" = "One sided type I error rate for <code>fisher.power</code> function.",
    "nsim (fisher.power)" = "Number of replications to calculate power. Default 100,000 for <code>fisher.power</code> function.",
    "seed (fisher.power)" = "Seed for simulations. Default 2000 for <code>fisher.power</code> function."
  )

  # Terms & Definitions for the SAM Prior Plot Library tab
  sam_plot_term_definitions <- list(
    "alpha_hist" = "Alpha (a) for Beta(a, b) in the informative prior (historical) for SAM Prior Plot.",
    "beta_hist" = "Beta (b) for Beta(a, b) in the informative prior (historical) for SAM Prior Plot.",
    "n_control" = "Number of patients in the control arm data for SAM Prior Plot.",
    "p_control" = "Simulated control rate for SAM Prior Plot.",
    "nf_alpha" = "Alpha for Beta(a, b) in the non-informative prior for SAM Prior Plot.",
    "nf_beta" = "Beta (b) for Beta(a, b) in the non-informative prior for SAM Prior Plot.",
    "delta" = "Clinically meaningful difference (delta) for SAM Prior Plot."
  )

  # Terms & Definitions for the SAM Weight Library tab
  sam_weight_term_definitions <- list(
    "alpha_hist" = "Alpha (a) for Beta(a, b) in the informative prior (historical) for SAM Weight calculation.",
    "beta_hist" = "Beta (b) for Beta(a, b) in the informative prior (historical) for SAM Weight calculation.",
    "n_control" = "Number of patients in the control arm summary data for SAM Weight calculation.",
    "r_control" = "Number of responses in the control arm summary data for SAM Weight calculation.",
    "delta" = "Clinically meaningful difference (delta) for SAM Weight calculation.",
    "method" = "Weight method ('LRT' or 'PPR') for SAM Weight calculation.",
    "prior_odds" = "Prior Odds H0 vs H1 (PPR only) for SAM Weight calculation."
  )

  # ===========================================================================
  # Observe Events for Term Definitions
  # ===========================================================================

  # When pressed on terms -> Outputs DPP Library tab
  observeEvent(input$term_pt, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pt"]])) })
  observeEvent(input$term_nt, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nt"]])) })
  observeEvent(input$term_pc, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pc"]])) })
  observeEvent(input$term_nc, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nc"]])) })
  observeEvent(input$term_p_calib, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pc.calib"]])) })
  observeEvent(input$term_pch, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["pch"]])) })
  observeEvent(input$term_nche, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nche"]])) })
  observeEvent(input$term_nch, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nch"]])) })
  observeEvent(input$term_alpha, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["alpha"]])) })
  observeEvent(input$term_tau, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["tau"]])) })
  observeEvent(input$term_a0c, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["a0c"]])) })
  observeEvent(input$term_b0c, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["b0c"]])) })
  observeEvent(input$term_a0t, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["a0t"]])) })
  observeEvent(input$term_b0t, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["b0t"]])) })
  observeEvent(input$term_delta_threshold, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["delta_threshold"]])) })
  observeEvent(input$term_method, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["method"]])) })
  observeEvent(input$term_theta, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["theta"]])) })
  observeEvent(input$term_eta, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["eta"]])) })
  observeEvent(input$term_datamat, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["datamat"]])) })
  observeEvent(input$term_w0, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["w0"]])) })
  observeEvent(input$term_nsim, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["nsim"]])) })
  observeEvent(input$term_seed, { output$definition_output <- renderUI(HTML(dpp_term_definitions[["seed"]])) })

  # When pressed on terms -> Outputs for DPP Analysis terms
  observeEvent(input$dpp_analysis_term_w, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["w"]])) })
  observeEvent(input$dpp_analysis_term_phat, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["phat_pt_larger_pc"]])) })
  observeEvent(input$dpp_analysis_term_apost_c_trial, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["apost_c_trial"]])) })
  observeEvent(input$dpp_analysis_term_bpost_c_trial, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["bpost_c_trial"]])) })
  observeEvent(input$dpp_analysis_term_apost_c_hca, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["apost_c_hca"]])) })
  observeEvent(input$dpp_analysis_term_bpost_c_hca, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["bpost_c_hca"]])) })
  observeEvent(input$dpp_analysis_term_apost_t, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["apost_t"]])) })
  observeEvent(input$dpp_analysis_term_bpost_t, { output$definition_output <- renderUI(HTML(dpp_analysis_term_definitions[["bpost_t"]])) })

  # When pressed on terms -> Outputs for Fisher Testing Library tab
  observeEvent(input$fisher_term_Yc, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["Yc (fisher)"]])) })
  observeEvent(input$fisher_term_nc_fisher, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nc (fisher)"]])) })
  observeEvent(input$fisher_term_Yt, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["Yt (fisher)"]])) })
  observeEvent(input$fisher_term_nt_fisher, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nt (fisher)"]])) })
  observeEvent(input$fisher_term_alternative, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["alternative (fisher)"]])) })
  observeEvent(input$fisher_term_pc_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["pc (fisher.bound)"]])) })
  observeEvent(input$fisher_term_nc_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nc (fisher.bound)"]])) })
  observeEvent(input$fisher_term_nt_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nt (fisher.bound)"]])) })
  observeEvent(input$fisher_term_alpha_bound, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["alpha (fisher.bound)"]])) })
  observeEvent(input$fisher_term_pt_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["pt (fisher.power)"]])) })
  observeEvent(input$fisher_term_nt_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nt (fisher.power)"]])) })
  observeEvent(input$fisher_term_pc_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["pc (fisher.power)"]])) })
  observeEvent(input$fisher_term_nc_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nc (fisher.power)"]])) })
  observeEvent(input$fisher_term_alpha_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["alpha (fisher.power)"]])) })
  observeEvent(input$fisher_term_nsim_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["nsim (fisher.power)"]])) })
  observeEvent(input$fisher_term_seed_power, { output$fisher_definition_output <- renderUI(HTML(fisher_term_definitions[["seed (fisher.power)"]])) })

  # When pressed on terms -> Outputs for SAM Prior Plot Library tab
  observeEvent(input$sam_plot_library_term_alpha_hist, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["alpha_hist"]])) })
  observeEvent(input$sam_plot_library_term_beta_hist, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["beta_hist"]])) })
  observeEvent(input$sam_plot_library_term_n_control, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["n_control"]])) })
  observeEvent(input$sam_plot_library_term_p_control, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["p_control"]])) })
  observeEvent(input$sam_plot_library_term_nf_alpha, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["nf_alpha"]])) })
  observeEvent(input$sam_plot_library_term_nf_beta, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["nf_beta"]])) })
  observeEvent(input$sam_plot_library_term_delta, { output$sam_plot_definition_output <- renderUI(HTML(sam_plot_term_definitions[["delta"]])) })

  # When pressed on terms -> Outputs for SAM Weight Library tab
  observeEvent(input$sam_weight_library_term_alpha_hist_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["alpha_hist"]])) })
  observeEvent(input$sam_weight_library_term_beta_hist_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["beta_hist"]])) })
  observeEvent(input$sam_weight_library_term_n_control_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["n_control"]])) })
  observeEvent(input$sam_weight_library_term_r_control_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["r_control"]])) })
  observeEvent(input$sam_weight_library_term_delta_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["delta"]])) })
  observeEvent(input$sam_weight_library_term_method_w, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["method"]])) })
  observeEvent(input$sam_weight_library_term_prior_odds, { output$sam_plot_definition_output <- renderUI(HTML(sam_weight_term_definitions[["prior_odds"]])) })


  # ===========================================================================
  # Conditional UI
  # ===========================================================================

  # Conditional UI for Theta for DPP method
  output$theta_ui <- renderUI({
    if (input$method == "Generalized BC") {
      numericInput("theta", "Theta (0-0.5)", value = 0.5, min = 0.001, max = 0.5, step = 0.001)
    }
  })

  # Conditional UI for ETA for DPP methods
  output$eta_ui <- renderUI({
    if (input$method %in% c("Bayesian p", "Generalized BC", "JSD")) {
      numericInput("eta", "Eta (0-Inf)", value = 1, min = 0, step = 0.1)
    }
  })

  # New Conditional UI for Theta in the DPP Analysis tab
  output$dpp_analysis_theta_ui <- renderUI({
    if (input$dpp_analysis_method == "Generalized BC") {
      numericInput("dpp_analysis_theta", "Theta (0-1)", value = 0.5, min = 0.001, max = 0.999, step = 0.001)
    }
  })

  # New Conditional UI for ETA in the DPP Analysis tab
  output$dpp_analysis_eta_ui <- renderUI({
    if (input$dpp_analysis_method %in% c("Bayesian p", "Generalized BC", "JSD")) {
      numericInput("dpp_analysis_eta", "Eta (0-Inf)", value = 1, min = 0, step = 0.1)
    }
  })

  # ===========================================================================
  # Dynamic Power Prior (DPP) Study Design
  # ===========================================================================

  # Create a reactive value to store the DPP output
  dpp_results <- reactiveVal(NULL)

  # Observe the "run_dpp" button click
  observeEvent(input$run_dpp, {
    shinyjs::html("analysis_status_message", "")

    # Input validation
    errors <- c()
    if (is.null(input$pt) || input$pt < 0 || input$pt > 1) errors <- c(errors, "Experimental Arm Response Rate (pt) must be between 0 and 1.")
    if (is.null(input$pc) || input$pc < 0 || input$pc > 1) errors <- c(errors, "Control Arm Response Rate (pc) must be between 0 and 1.")
    if (is.null(input$p_calib) || input$p_calib < 0 || input$p_calib > 1) errors <- c(errors, "Response Rate for Calibration (pc.calib) must be between 0 and 1.")
    if (is.null(input$pch) || input$pch < 0 || input$pch > 1) errors <- c(errors, "Historical Control Response Rate (pch) must be between 0 and 1.")
    if (is.null(input$nt) || input$nt <= 0 || is.null(input$nc) || input$nc <= 0 || is.null(input$nch) || input$nch <= 0 || is.null(input$nche) || input$nche <= 0) {
      errors <- c(errors, "Sample sizes (nt, nc, nch, nche) must be positive integers.")
    }
    if (!is.null(input$nche) && input$nche > 50) {
      errors <- c(errors, "Maximum Number of Patients from Historical Study (nche) cannot be greater than 50.")
    }
    if (is.null(input$alpha) || input$alpha <= 0 || input$alpha > 1) errors <- c(errors, "Type I Error (alpha) must be between 0 and 1.")
    if (!is.null(input$nche) && !is.null(input$nch) && input$nche > input$nch) errors <- c(errors, "Equivalent number of patients borrowed (nche) cannot exceed total historical control patients (nch).")
    if (is.null(input$delta_threshold) || input$delta_threshold < 0 || input$delta_threshold > 1) errors <- c(errors, "Delta Threshold must be between 0 and 1.")

    if (input$method == "Generalized BC" && (is.null(input$theta) || input$theta <= 0 || input$theta >= 1)) {
      errors <- c(errors, "Theta must be between 0 and 1 for 'Generalized BC' method.")
    }
    if (input$method %in% c("Bayesian p", "Generalized BC", "JSD") && (is.null(input$eta) || input$eta <= 0)) {
      errors <- c(errors, "Eta must be positive for selected methods.")
    }

    seed_val <- suppressWarnings(as.integer(input$seed))
    if (is.na(seed_val) || seed_val < 1) {
      errors <- c(errors, "Seed must be a positive integer.")
    }

    if (length(errors) > 0) {
      shinyjs::html("analysis_status_message", paste0("<p style='color: red;'>", paste(errors, collapse = "<br>"), "</p>"))
      return() # Stop execution if there are errors
    }

    if (length(errors) > 0) {
      shinyjs::html("analysis_status_message", paste0("<p style='color: red;'>", paste(errors, collapse = "<br>"), "</p>"))
      return()
    }

    shinyjs::html("analysis_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_dpp")

    power.DPP_args <- list(
      pt = input$pt,
      nt = input$nt,
      pc = input$pc,
      nc = input$nc,
      pc.calib = input$p_calib,
      pch = input$pch,
      nche = input$nche,
      nch = input$nch,
      alpha = input$alpha,
      a0c = input$a0c,
      b0c = input$b0c,
      a0t = input$a0t,
      b0t = input$b0t,
      delta_threshold = input$delta_threshold,
      method = input$method,
      nsim = input$nsim,
      seed = seed_val
    )
    if (input$method == "Generalized BC") {
      power.DPP_args$theta = input$theta
      power.DPP_args$eta = input$eta
    } else if (input$method %in% c("Bayesian p", "JSD")) {
      power.DPP_args$eta = input$eta
    }

    p <- future({
      result <- do.call(BayesianHybridDesign::power.DPP, power.DPP_args)
      result
    }, seed = TRUE)

    then(p, onFulfilled = function(result) {
      dpp_results(result)
      shinyjs::html("analysis_status_message", "<p style='color: green;'>Study Design calculation ran successfully!</p>")
      updateTabsetPanel(session, "bha_tabset_panel", selected = "Study Design Results")
      shinyjs::enable("run_dpp")
    }, onRejected = function(e) {
      shinyjs::html("analysis_status_message", paste0("<p style='color: red;'>Study Design calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
      dpp_results(NULL)
      shinyjs::enable("run_dpp")
    })
  })

  # Outputs in the "Study Design Results" tab
  output$dppPower <- renderPrint({
    req(dpp_results()$power)
    cat(round(dpp_results()$power, 4))
  })

  output$dppTau <- renderPrint({
    req(dpp_results()$tau)
    cat(round(dpp_results()$tau, 4))
  })

  output$dppDeltaBound <- renderPrint({
    req(dpp_results()$delta.bound)
    cat(round(dpp_results()$delta.bound, 4))
  })

  # Updated to reflect the user's request
  output$dppPcPMD <- renderPrint({
    req(dpp_results()$pc.PMD)
    cat(round(dpp_results()$pc.PMD, 4))
  })

  output$dppPcSdPMD <- renderPrint({
    req(dpp_results()$pc.sd)
    cat(round(dpp_results()$pc.sd, 4))
  })

  # Plotting $mean.hca histogram
  output$dpp_mean_hca_hist <- renderPlotly({
    req(dpp_results()$mean_hca)
    plot_data <- data.frame(mean_hca = dpp_results()$mean_hca)
    if (length(plot_data$mean_hca) == 0) {
      return(ggplotly(ggplot() + labs(title = "No data available.")))
    }
    p <- ggplot(plot_data, aes(x = mean_hca)) +
      geom_histogram(binwidth = 0.01, fill = "skyblue", color = "black") +
      labs(title = "Posterior mean ORR for hybrid control by simulated trials",
           x = "Posterior Mean ORR (Hybrid Control)",
           y = "Frequency") +
      theme_minimal()
    ggplotly(p)
  })

  # Plotting $mean.c histogram
  output$dpp_mean_c_hist <- renderPlotly({
    req(dpp_results()$mean_c)
    plot_data <- data.frame(mean_c = dpp_results()$mean_c)
    if (length(plot_data$mean_c) == 0) {
      return(ggplotly(ggplot() + labs(title = "No data available.")))
    }
    p <- ggplot(plot_data, aes(x = mean_c)) +
      geom_histogram(binwidth = 0.01, fill = "lightgreen", color = "black") +
      labs(title = "Posterior mean ORR for concurrent control by simulated trials",
           x = "Posterior Mean ORR (Concurrent Control)",
           y = "Frequency") +
      theme_minimal()
    ggplotly(p)
  })

  output$dpp_plot_pmd <- renderPlotly({
    req(dpp_results()$mean_hca)
    req(dpp_results()$mean_c)
    pmd_distribution <- dpp_results()$mean_hca - dpp_results()$mean_c

    if (length(pmd_distribution) > 1) {
      plot_data <- data.frame(pmd = pmd_distribution)
      p <- ggplot(plot_data, aes(x = pmd)) +
        geom_histogram(binwidth = 0.01, fill = "coral", color = "black") +
        labs(title = "Posterior Mean Difference (PMD) Distribution",
             x = "Posterior Mean Difference",
             y = "Frequency") +
        theme_minimal()
      ggplotly(p)
    } else {
      return(ggplotly(ggplot() + labs(title = "PMD plot data not available or insufficient points.")))
    }
  })

  # ===========================================================================
  # Dynamic Power Prior (DPP) Analysis
  # ===========================================================================

  # Create a reactive value to store the DPP analysis output
  dpp_analysis_results <- reactiveVal(NULL)

  # Observe the "run_dpp_analysis" button click
  observeEvent(input$run_dpp_analysis, {
    shinyjs::html("dpp_analysis_status_message", "")

    # Input validation for DPP.analysis specific inputs
    errors <- c()
    if (is.null(input$dpp_analysis_rt) || input$dpp_analysis_rt < 0) errors <- c(errors, "Observed Responders in Experimental Arm (rt) must be non-negative.")
    if (is.null(input$dpp_analysis_rc) || input$dpp_analysis_rc < 0) errors <- c(errors, "Observed Responders in Control Arm (rc) must be non-negative.")
    if (!is.null(input$dpp_analysis_rt) && !is.null(input$dpp_analysis_nt) && input$dpp_analysis_rt > input$dpp_analysis_nt) errors <- c(errors, "Observed responders in experimental arm (rt) cannot exceed total experimental sample size (nt).")
    if (!is.null(input$dpp_analysis_rc) && !is.null(input$dpp_analysis_nc) && input$dpp_analysis_rc > input$dpp_analysis_nc) errors <- c(errors, "Observed responders in control arm (rc) cannot exceed total control sample size (nc).")

    if (length(errors) > 0) {
      shinyjs::html("dpp_analysis_status_message", paste0("<p style='color: red;'>", paste(errors, collapse = "<br>"), "</p>"))
      return()
    }

    shinyjs::html("dpp_analysis_status_message", "<p style='color: blue;'>DPP Analysis started. Please wait...</p>")
    shinyjs::disable("run_dpp_analysis")

    Ych_calculated <- round(input$dpp_analysis_pch * input$dpp_analysis_nch)

    analysis_args <- list(
      Yt = input$dpp_analysis_rt,
      nt = input$dpp_analysis_nt,
      Yc = input$dpp_analysis_rc,
      nc = input$dpp_analysis_nc,
      Ych = Ych_calculated,
      nch = input$dpp_analysis_nch,
      nche = input$dpp_analysis_nche,
      a0c = input$dpp_analysis_a0c,
      b0c = input$dpp_analysis_b0c,
      a0t = input$dpp_analysis_a0t,
      b0t = input$dpp_analysis_b0t,
      delta_threshold = input$dpp_analysis_delta_threshold,
      method = input$dpp_analysis_method
    )

    if (input$dpp_analysis_method == "Generalized BC") {
      analysis_args$theta = input$dpp_analysis_theta
    }
    if (input$dpp_analysis_method %in% c("Bayesian p", "Generalized BC", "JSD")) {
      analysis_args$eta = input$dpp_analysis_eta
    }

    p <- future({
      result <- do.call(BayesianHybridDesign::DPP.analysis, analysis_args)
      result
    })

    then(p, onFulfilled = function(result) {
      dpp_analysis_results(result)
      shinyjs::html("dpp_analysis_status_message", "<p style='color: green;'>DPP Analysis ran successfully!</p>")
      updateTabsetPanel(session, "bha_tabset_panel", selected = "Statistical Analysis Results (DPP)")
      shinyjs::enable("run_dpp_analysis")
    }, onRejected = function(e) {
      shinyjs::html("dpp_analysis_status_message", paste0("<p style='color: red;'>DPP Analysis calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
      dpp_analysis_results(NULL)
      shinyjs::enable("run_dpp_analysis")
    })
  })

  # Output for DPP.analysis results
  output$dppAnalysisResult <- renderPrint({
    req(dpp_analysis_results())
    print(dpp_analysis_results())
  })

  # Reactive expression to calculate pc
  calculated_pc <- reactive({
    req(input$dpp_analysis_rc, input$dpp_analysis_nc)
    if (input$dpp_analysis_nc > 0) {
      return(input$dpp_analysis_rc / input$dpp_analysis_nc)
    } else {
      return(0)
    }
  })

  # Reactive expression to calculate pt
  calculated_pt <- reactive({
    req(input$dpp_analysis_rt, input$dpp_analysis_nt)
    if (input$dpp_analysis_nt > 0) {
      return(input$dpp_analysis_rt / input$dpp_analysis_nt)
    } else {
      return(0)
    }
  })

  # Display the calculated pc
  output$dpp_analysis_pc_display <- renderPrint({
    cat(sprintf("%.4f", calculated_pc()))
  })

  # Display the calculated pt
  output$dpp_analysis_pt_display <- renderPrint({
    cat(sprintf("%.4f", calculated_pt()))
  })

  # Plot for plotDPP() function
  output$plotDPP <- renderPlot({
    req(dpp_analysis_results())
    BayesianHybridDesign::plotDPP(DPP = dpp_analysis_results())
  })

  # New outputs for the Summary section
  output$concurrent_summary <- renderPrint({
    req(dpp_analysis_results()$apost_c_trial)
    req(dpp_analysis_results()$bpost_c_trial)

    a_post <- dpp_analysis_results()$apost_c_trial
    b_post <- dpp_analysis_results()$bpost_c_trial

    median <- qbeta(p = 0.5, shape1 = a_post, shape2 = b_post)
    ci_lower <- qbeta(p = 0.025, shape1 = a_post, shape2 = b_post)
    ci_upper <- qbeta(p = 0.975, shape1 = a_post, shape2 = b_post)

    cat("Parameters: \n")
    cat(sprintf("  apost_c_trial: %.4f\n", a_post))
    cat(sprintf("  bpost_c_trial: %.4f\n", b_post))
    cat("-----------------------------------\n")
    cat(sprintf("Median: %.4f\n", median))
    cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n", ci_lower, ci_upper))
  })

  output$hybrid_control_summary <- renderPrint({
    req(dpp_analysis_results()$apost_c_hca)
    req(dpp_analysis_results()$bpost_c_hca)

    a_post <- dpp_analysis_results()$apost_c_hca
    b_post <- dpp_analysis_results()$bpost_c_hca

    median <- qbeta(p = 0.5, shape1 = a_post, shape2 = b_post)
    ci_lower <- qbeta(p = 0.025, shape1 = a_post, shape2 = b_post)
    ci_upper <- qbeta(p = 0.975, shape1 = a_post, shape2 = b_post)

    cat("Parameters: \n")
    cat(sprintf("  apost_c_hca: %.4f\n", a_post))
    cat(sprintf("  bpost_c_hca: %.4f\n", b_post))
    cat("-----------------------------------\n")
    cat(sprintf("Median: %.4f\n", median))
    cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n", ci_lower, ci_upper))
  })

  output$experimental_summary <- renderPrint({
    req(dpp_analysis_results()$apost_t)
    req(dpp_analysis_results()$bpost_t)

    a_post <- dpp_analysis_results()$apost_t
    b_post <- dpp_analysis_results()$bpost_t

    median <- qbeta(p = 0.5, shape1 = a_post, shape2 = b_post)
    ci_lower <- qbeta(p = 0.025, shape1 = a_post, shape2 = b_post)
    ci_upper <- qbeta(p = 0.975, shape1 = a_post, shape2 = b_post)

    cat("Parameters: \n")
    cat(sprintf("  apost_t: %.4f\n", a_post))
    cat(sprintf("  bpost_t: %.4f\n", b_post))
    cat("-----------------------------------\n")
    cat(sprintf("Median: %.4f\n", median))
    cat(sprintf("95%% Credible Interval: (%.4f, %.4f)\n", ci_lower, ci_upper))
  })

  # Updated output for Statistical Conclusion
  output$statistical_conclusion <- renderPrint({
    result <- dpp_analysis_results()

    if (is.null(result)) {
      cat("No DPP analysis results available yet. Please run the analysis.")
      return()
    }

    # Extract phat
    phat <- result$phat_pt_larger_pc

    # Try to get tau from result, from design, or from user input
    tau <- NULL

    # 1) From DPP.analysis output
    if (!is.null(result$tau)) {
      tau <- result$tau
    }

    # 2) From study design results (if design was run first)
    if (is.null(tau) && !is.null(dpp_results()$tau)) {
      tau <- dpp_results()$tau
    }

    # 3) From a tau input in the DPP Analysis tab (if you add one in UI)
    if (is.null(tau) && !is.null(input$dpp_analysis_tau)) {
      tau <- input$dpp_analysis_tau
    }

    # 4) Last resort — set a default threshold
    if (is.null(tau)) {
      tau <- 0.95
    }

    # If phat is still missing, we can't continue
    if (is.null(phat)) {
      cat("Analysis did not return 'phat_pt_larger_pc'. Cannot compute conclusion.")
      return()
    }

    # Conclusion text
    cat(sprintf(
      "Statistical significance is defined as P(ORR Treatment > ORR Control | Hybrid Data) >= tau.\n"
    ))

    if (phat >= tau) {
      cat(sprintf(
        "Since P(ORR_t > ORR_c) = %.4f is GREATER THAN OR EQUAL to the significance threshold (tau = %.4f),\n",
        phat, tau
      ))
      cat("The experimental treatment is statistically significant compared to the control treatment in ORR.")
    } else {
      cat(sprintf(
        "Since P(ORR_t > ORR_c) = %.4f is LESS THAN the significance threshold (tau = %.4f),\n",
        phat, tau
      ))
      cat("The experimental treatment is NOT statistically significant compared to the control treatment in ORR.")
    }
  })


  # ===========================================================================
  # SAM Prior Plot Generation
  # ===========================================================================

  # SAM Prior Plot Generation
  observeEvent(input$run_sam, {
    shinyjs::html("sam_plot_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_sam")
    tryCatch({
      set.seed(123)
      prior_hist <- RBesT::mixbeta(c(1, input$alpha_hist, input$beta_hist))
      control_data <- rbinom(input$n_control, size = 1, prob = input$p_control)
      r_control_for_weight <- sum(control_data)
      weight <- SAMprior::SAM_weight(
        if.prior = prior_hist,
        delta = input$delta,
        n = input$n_control,
        r = r_control_for_weight
      )
      nf_prior <- RBesT::mixbeta(c(1, input$nf_alpha, input$nf_beta))

      output$samPlot <- renderPlot({
        sam_prior <- SAMprior::SAM_prior(
          if.prior = prior_hist,
          nf.prior = nf_prior,
          weight = weight
        )
        plot(sam_prior)
      })


      shinyjs::html("sam_plot_status_message", "<p style='color: green;'>SAM Prior Plot generation ran successfully!</p>")
      updateTabsetPanel(session, "sam_design_tabset_panel", selected = "Results and Conclusion")
    }, error = function(e){
      output$samPlot <- renderPlot(NULL)
      output$samWeightResult <- renderPrint(NULL) # Clear this if the plot fails
      shinyjs::html("sam_plot_status_message", paste0("<p style='color: red;'>SAM Prior Plot generation failed: ", e$message, ". Please check inputs and try again.</p>"))
    }, finally = {
      shinyjs::enable("run_sam")
    })
  })


  # ===========================================================================
  # SAM Weight Library
  # ===========================================================================

  # SAM Weight Calculation
  observeEvent(input$run_sam_weight, {
    shinyjs::html("sam_weight_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_sam_weight")
    tryCatch({
      set.seed(123)
      prior_hist <- RBesT::mixbeta(c(1, input$alpha_hist, input$beta_hist))
      weight <- SAMprior::SAM_weight(
        if.prior = prior_hist,
        delta = input$delta,
        method.w = input$sam_method,
        prior.odds = input$prior_odds,
        n = input$n_control,
        r = input$r_control
      )
      output$samWeightResult <- renderPrint({
        print(weight)
      })
      shinyjs::html("sam_weight_status_message", "<p style='color: green;'>SAM Weight calculation ran successfully!</p>")
      updateTabsetPanel(session, "sam_design_tabset_panel", selected = "Results and Conclusion")

      # Dynamic content for Conclusion and Summary tab
      output$samPriorConclusion <- renderUI({
        HTML(
          "<h4>Summary</h4>
          <p>The SAM Prior is a mixture of an informative prior (based on historical data) and a non-informative prior. The calculated mixture weight (W) is <strong>", round(weight, 4), "</strong>. This weight determines the degree of borrowing from the historical data.</p>

          <h4>Statistical Conclusion</h4>
          <p>The weight of <strong>", round(weight, 4), "</strong> suggests that there is a <strong>", round(weight * 100, 2), "%</strong> level of confidence in the historical data's relevance. A higher weight indicates a greater similarity between the historical and concurrent control data, leading to more borrowing and a more informative prior. This approach efficiently utilizes historical information while mitigating the risk of potential conflicts, which could otherwise introduce bias into the study.</p>"
        )
      })

    }, error = function(e){
      output$samWeightResult <- renderPrint({ cat("Error:", conditionMessage(e), '\n') })
      shinyjs::html("sam_weight_status_message", paste0("<p style='color: red;'>SAM Weight calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
    }, finally = {
      shinyjs::enable("run_sam_weight")
    })
  })

  # ===========================================================================
  # Fisher Testing Library
  # ===========================================================================

  # Fisher Power Calculation
  observeEvent(input$run_fisher_power, {
    shinyjs::html("fisher_power_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_fisher_power")
    tryCatch({
      set.seed(input$fp_seed)
      power_result <- BayesianHybridDesign::fisher.power(
        pt = input$fp_pt,
        nt = input$fp_nt,
        pc = input$fp_pc,
        nc = input$fp_nc,
        alpha = input$fp_alpha,
        nsim = input$fp_nsim,
        seed = input$fp_seed
      )
      output$fisherPowerResult <- renderPrint({ cat("Power:", round(power_result, 4), "\n") })
      shinyjs::html("fisher_power_status_message", "<p style='color: green;'>Fisher Power calculation ran successfully!</p>")
      updateTabsetPanel(session, "fisher_power_tabset_panel", selected = "Design Results")

      # Dynamic content for Conclusion and Summary tab
      output$fisherPowerConclusion <- renderUI({
        HTML(
          "<h4>Summary</h4>
          <p>The statistical power of the study was estimated to be <strong>", round(power_result, 4), "</strong>. This value represents the probability of detecting a statistically significant treatment effect given the design parameters (e.g., sample sizes, response rates, alpha level).</p>

          <h4>Statistical Conclusion</h4>
          <p>A power of <strong>", round(power_result, 4), "</strong> indicates a <strong>", round(power_result * 100, 2), "%</strong> chance of correctly rejecting the null hypothesis (i.e., finding a significant effect) if the true effect size is as specified. This value is a crucial measure of the study's ability to avoid a Type II error (a false negative). Typically, a power of 80% or greater is considered acceptable, and your study design either meets or exceeds this standard, suggesting it is well-powered to detect a true difference.</p>"
        )
      })

    }, error = function(e) {
      output$fisherPowerResult <- renderPrint({ cat("Error calculating Fisher Power:", conditionMessage(e), "\n") })
      shinyjs::html("fisher_power_status_message", paste0("<p style='color: red;'>Fisher Power calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
    }, finally = {
      shinyjs::enable("run_fisher_power")
    })
  })

  # Fisher Bound Calculation
  observeEvent(input$run_fisher_bound, {
    shinyjs::html("fisher_bound_status_message", "<p style='color: blue;'>Calculation started. Please wait...</p>")
    shinyjs::disable("run_fisher_bound")
    tryCatch({
      bound_result <- BayesianHybridDesign::fisher.bound(
        pc = input$fb_pc,
        nc = input$fb_nc,
        nt = input$fb_nt,
        alpha = input$fb_alpha
      )
      output$fisherBoundM <- renderPrint({ print(bound_result$M) })
      output$fisherBoundP <- renderPrint({ round(bound_result$p, 4) })
      output$fisherBoundRc <- renderPrint({ bound_result$rc })
      output$fisherBoundNc <- renderPrint({ bound_result$nc })
      output$fisherBoundRt <- renderPrint({ round(bound_result$rt, 4) })
      output$fisherBoundNt <- renderPrint({ bound_result$nt })
      output$fisherBoundDelta <- renderPrint({ round(bound_result$delta, 4) })
      shinyjs::html("fisher_bound_status_message", "<p style='color: green;'>Fisher Bound calculation ran successfully!</p>")
      updateTabsetPanel(session, "fisher_bound_tabset_panel", selected = "Design Results")

      # Dynamic content for Conclusion and Summary tab
      output$fisherBoundConclusion <- renderUI({
        HTML(
          "<h4>Summary</h4>
          <p>The Fisher Bound analysis determined that a minimum of <strong>", bound_result$M, "</strong> responders in the treatment arm are required to achieve statistical significance. This corresponds to a minimum detectable difference of <strong>", round(bound_result$delta, 4), "</strong> and a p-value of <strong>", round(bound_result$p, 4), "</strong> at the boundary.</p>"
        )
      })

      output$fisherBoundStatisticalConclusion <- renderUI({
        HTML(
          "<h4>Statistical Conclusion</h4>
          <p>The number <strong>", bound_result$M, "</strong> is the critical benchmark for this study. If the observed number of responses in the experimental arm is less than this value, the trial will not achieve statistical significance at the specified alpha level. This provides a clear, objective criterion for evaluating the success of the experimental treatment's performance based on the selected study parameters.</p>"
        )
      })

    }, error = function(e) {
      output$fisherBoundM <- renderPrint(NULL)
      output$fisherBoundP <- renderPrint(NULL)
      output$fisherBoundRc <- renderPrint(NULL)
      output$fisherBoundNc <- renderPrint(NULL)
      output$fisherBoundRt <- renderPrint(NULL)
      output$fisherBoundNt <- renderPrint(NULL)
      shinyjs::html("fisher_bound_status_message", paste0("<p style='color: red;'>Fisher Bound calculation failed: ", e$message, ". Please check inputs and try again.</p>"))
    }, finally = {
      shinyjs::enable("run_fisher_bound")
    })
  })
})
