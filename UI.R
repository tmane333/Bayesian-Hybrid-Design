library(shiny)
library(ggplot2)
library(foreach)
library(doParallel)
library(SAMprior)
library(shinyjs)
library(BayesianHybridDesign)
library(future)
library(promises)
library(plotly)

shinyUI(navbarPage(
  id = "main_tabset_panel",
  title = tags$div(
    tags$img(src = "https://i.postimg.cc/59rbFS2r/design1.png", height = "30px"),
    "Bayesian Hybrid Design"
  ),

  tabPanel(
    "Home",
    fluidPage(
      # Title and Tagline under the Logo
      fluidRow(
        column(12,
               div(style = "text-align: center; margin-top: 20px;",
                   tags$img(src = "https://i.postimg.cc/59rbFS2r/design1.png", height = "200px"),
                   h1(strong("Bayesian Hybrid Design")),
                   h4("A Comprehensive R Shiny App for Clinical Study Design and Analysis")
               )
        )
      ),

      # Introduction section
      fluidRow(
        column(10, offset = 1,
               br(),
               h3(strong("Introduction")),
               p("In the pharmaceutical industry, developing a new drug is a high-cost, long-term process. For example, the average cost to develop a cancer drug is approximately $1.2 billion, and the timeline from initial laboratory research to patient use can exceed 10 years."),
               p("To accelerate this process, one promising approach is the hybrid design method. This method allows for the appropriate incorporation of external data into a current study, which enhances go/no-go decision-making for the next phases of development."),
               p("Multiple methods have been proposed within Bayesian hybrid designs. A particularly important method, recently proposed by Lu et al (2025), has demonstrated both strong performance and fast computation. However, without a visual graphical interface, it remains challenging for users to apply this method directly."),
               p("In response, we developed this R Shiny app to provide a user-friendly GUI. It implements the Dynamic Power Prior (DPP) and SAM Prior methods in Bayesian hybrid design. For benchmarking, we have also included the Fisher's Exact approach."),
               p("This application provides a convenient way for users to explore optimal study designs using these hybrid methods. Additionally, the app offers statistical analysis for studies based on Bayesian hybrid design, displaying the posterior distribution of response rates for each arm and implementing comparative statistical inference. The app will be published and available to the public.")
        )
      ),

      hr(),

      # Key Concepts section
      fluidRow(
        column(12,
               h3("Key Concepts", style = "text-align: center;"),
               br(),
               fluidRow(
                 column(4, wellPanel(
                   h4("Dynamic Power Prior (DPP)"),
                   p("A method that allows for the incorporation of historical data into a current study. It dynamically adjusts the weight given to the historical data based on the similarity between the historical and current control groups. This prevents potential bias from conflicting data.")
                 )),
                 column(4, wellPanel(
                   h4("Self-Adapting Mixture (SAM) Prior"),
                   p("This method creates a prior distribution by mixing an informative prior (based on historical data) with a non-informative, or skeptical, prior. The 'self-adapting' feature means it adjusts the borrowing weight dynamically based on conflicts between new and historical data.")
                 )),
                 column(4, wellPanel(
                   h4("Fisher's Exact Method"),
                   p("A statistical test used to analyze the association between two categorical variables, such as treatment response. It is often used as a benchmark for comparison in study design.")
                 ))
               )
        )
      ),

      hr(),

      # Reference section
      fluidRow(
        column(12,
               h4(strong("Reference")),
               p(tags$a(href = "https://onlinelibrary.wiley.com/doi/abs/10.1002/pst.2466",
                        "Lu Z, Toso J, Ayele G, He P. A Bayesian Hybrid Design With Borrowing From Historical Study. Pharm Stat. 2025 Mar-Apr;24(2):e2466. doi: 10.1002/pst.2466. Epub 2024 Dec 27. PMID: 39731333."))
        )
      ),

      # Credit section
      fluidRow(
        column(12,
               br(),
               p(em("This app's design and code were developed with the assistance of Tanvi Mane."))
        )
      )
    )
  ),

  tabPanel(
    "Dynamic Power Prior (DPP)",
    tabsetPanel(
      id = "bha_tabset_panel",
      tabPanel(
        "Library",
        fluidRow(
          column(4,
                 wellPanel(
                   h4("General Terms"),
                   actionLink("term_pt", "pt"), br(), actionLink("term_nt", "nt"), br(),
                   actionLink("term_pc", "pc"), br(), actionLink("term_nc", "nc"), br(),
                   actionLink("term_p_calib", "pc.calib"), br(), actionLink("term_pch", "pch"), br(),
                   actionLink("term_nche", "nche"), br(), actionLink("term_nch", "nch"), br(),
                   actionLink("term_alpha", "alpha"), br(), actionLink("term_tau", "tau"), br(),
                   actionLink("term_a0c", "a0c"), br(), actionLink("term_b0c", "b0c"), br(),
                   actionLink("term_a0t", "a0t"), br(), actionLink("term_b0t", "b0t"), br(),
                   actionLink("term_delta_threshold", "delta_threshold"), br(), actionLink("term_method", "method"), br(),
                   actionLink("term_theta", "theta"), br(), actionLink("term_eta", "eta"), br(),
                   actionLink("term_datamat", "datamat"), br(), actionLink("term_w0", "w0"), br(),
                   actionLink("term_nsim", "nsim"), br(), actionLink("term_seed", "seed")
                 ),
                 # DPP Analysis terms
                 wellPanel(
                   h4("DPP Analysis Terms"),
                   actionLink("dpp_analysis_term_w", "w"), br(),
                   actionLink("dpp_analysis_term_phat", "phat_pt_larger_pc"), br(),
                   actionLink("dpp_analysis_term_apost_c_trial", "apost_c_trial"), br(),
                   actionLink("dpp_analysis_term_bpost_c_trial", "bpost_c_trial"), br(),
                   actionLink("dpp_analysis_term_apost_c_hca", "apost_c_hca"), br(),
                   actionLink("dpp_analysis_term_bpost_c_hca", "bpost_c_hca"), br(),
                   actionLink("dpp_analysis_term_apost_t", "apost_t"), br(),
                   actionLink("dpp_analysis_term_bpost_t", "bpost_t")
                 )
          ),
          column(8,
                 wellPanel(
                   h4("Definition"),
                   uiOutput("definition_output")
                 )
          )
        )
      ),
      tabPanel(
        "Bayesian Hybrid Design",
        shinyjs::useShinyjs(),
        fluidRow(
          column(6,
                 wellPanel(
                   h4("Current Study Design Parameters"),
                   fluidRow(
                     column(6, numericInput("nt", "Experimental Arm Sample Size (nt)", value = 50, min = 1)),
                     column(6, numericInput("pt", "Experimental Arm Response Rate (pt)", value = 0.5, min = 0, max = 1, step = 0.01))
                   ),
                   fluidRow(
                     column(6, numericInput("nc", "Control Arm Sample Size (nc)", value = 50, min = 1)),
                     column(6, numericInput("pc", "Control Arm Response Rate (pc)", value = 0.3, min = 0, max = 1, step = 0.01))
                   ),
                   uiOutput("error_pt"), uiOutput("error_nt"), uiOutput("error_pc"), uiOutput("error_nc")
                 ),
                 wellPanel(
                   h4("Historical Control Study Data"),
                   fluidRow(
                     column(6, numericInput("nch", "Historical Control Sample Size (nch)", value = 100, min = 1)),
                     column(6, numericInput("pch", "Historical Control Number of Responders (pch)", value = 0.3, min = 0, max = 1, step = 0.01))
                   ),
                   uiOutput("error_nch"), uiOutput("error_pch")
                 ),
                 wellPanel(
                   h4("Bayesian Borrowing"),
                   fluidRow(
                     column(6, numericInput("nche", "Maximum Number of Patients from Historical Study (nche)", value = 50, min = 1)),
                     column(6, numericInput("nch_total", "Total number of patients in historical control (nch)", value = 100, min = 1))
                   ),
                   numericInput("delta_threshold", "Delta Threshold (0-1)", value = 0.1, min = 0, max = 1, step = 0.01),
                   uiOutput("error_nche_nch_total"),
                   selectInput("method", "Dynamic Control Method", choices = c("Empirical Bayes", "Bayesian p", "Generalized BC", "JSD")),
                   fluidRow(
                     column(6, uiOutput("theta_ui"))
                   ),
                   fluidRow(
                     column(6, uiOutput("eta_ui"))
                   )
                 )
          ),
          column(6,
                 wellPanel(
                   h4("One-Sided Type I Error"),
                   numericInput("alpha", "Type I Error", value = 0.1, min = 0.001, max = 0.5),
                   uiOutput("error_alpha_tau")
                 ),
                 wellPanel(
                   h4("Response Rate for Calibration (pc.calib)"),
                   numericInput("p_calib", "Default Control Arm Response Rate", value = 0.3, min = 0, max = 1, step = 0.01)
                 ),
                 wellPanel(
                   h4("Hyperpriors in Beta Distribution"),
                   fluidRow(
                     column(6, numericInput("a0c", "Control Arm a0c", value = 0.001, min = 0)),
                     column(6, numericInput("b0c", "Control Arm b0c", value = 0.001, min = 0))
                   ),
                   fluidRow(
                     column(6, numericInput("a0t", "Experimental Arm a0t", value = 0.001, min = 0)),
                     column(6, numericInput("b0t", "Experimental Arm b0t", value = 0.001, min = 0))
                   )
                 ),
                 wellPanel(
                   h4("Simulation"),
                   fluidRow(
                     column(6, numericInput("nsim", "Number of Trials", value = 100000, min = 0, step = 10)),
                     column(6, textInput("seed", "Seed", value = "2000"))
                   )
                 )
          )
        ),
        actionButton("run_dpp", "Study Design Results", class = "btn-primary"),
        br(), br(),
        div(id = "analysis_status_message", style = "text-align: center;")
      ),
      tabPanel(
        "Study Design Results",
        fluidRow(
          column(4, wellPanel(h5("Power"), verbatimTextOutput("dppPower"))),
          column(4, wellPanel(h5("Tau"), verbatimTextOutput("dppTau"))),
          column(4, wellPanel(h5("Statistical critical value in ORR difference (detectable difference)"), verbatimTextOutput("dppDeltaBound")))
        ),
        fluidRow(
          column(6, wellPanel(h5("Posterior mean difference (PMD) between hybrid control and concurrent control"), verbatimTextOutput("dppPcPMD"))),
          column(6, wellPanel(h5("Standard deviation of PMD"), verbatimTextOutput("dppPcSdPMD")))
        ),
        fluidRow(
          column(6,
                 wellPanel(
                   h5("Posterior mean ORR for hybrid control by simulated trials"),
                   plotlyOutput("dpp_mean_hca_hist")
                 )
          ),
          column(6,
                 wellPanel(
                   h5("Posterior mean ORR for concurrent control by simulated trials"),
                   plotlyOutput("dpp_mean_c_hist")
                 )
          )
        ),
        hr(),
        h4("Posterior Mean Difference (PMD) Distribution"),
        plotlyOutput("dpp_plot_pmd")
      ),
      tabPanel(
        "Statistical Analysis (DPP)",
        shinyjs::useShinyjs(),
        fluidRow(
          column(6,
                 wellPanel(
                   h4("Current Study Observed Data"),
                   fluidRow(
                     column(6, numericInput("dpp_analysis_rt", "Number of Responders in Experimental Arm (rt)", value = 30, min = 0)),
                     column(6, numericInput("dpp_analysis_rc", "Number of Responders in Control Arm (rc)", value = 15, min = 0))
                   )
                 ),
                 wellPanel(
                   h4("Hyperpriors in Beta Distribution"),
                   fluidRow(
                     column(6, numericInput("dpp_analysis_a0c", "Control Arm a0c", value = 0.001, min = 0)),
                     column(6, numericInput("dpp_analysis_b0c", "Control Arm b0c", value = 0.001, min = 0))
                   ),
                   fluidRow(
                     column(6, numericInput("dpp_analysis_a0t", "Experimental Arm a0t", value = 0.001, min = 0)),
                     column(6, numericInput("dpp_analysis_b0t", "Experimental Arm b0t", value = 0.001, min = 0))
                   )
                 ),
                 wellPanel(
                   h4("Historical Control Study Data"),
                   fluidRow(
                     column(6, numericInput("dpp_analysis_nch", "Historical Control Sample Size (nch)", value = 87, min = 1)),
                     column(6, numericInput("dpp_analysis_pch", "Historical Control Number of Responders (pch)", value = 0.3, min = 0, max = 1, step = 0.01))
                   )
                 )
          ),
          column(6,
                 wellPanel(
                   h4("Analysis Options"),
                   fluidRow(
                     column(6, numericInput("dpp_analysis_nt", "Experimental Arm Sample Size (nt)", value = 41, min = 1)),
                     column(6,
                            # Non-editable
                            h5("Experimental Arm Response Rate (pt)"),
                            verbatimTextOutput("dpp_analysis_pt_display")
                     )
                   ),
                   fluidRow(
                     column(6, numericInput("dpp_analysis_nc", "Control Arm Sample Size (nc)", value = 44, min = 1)),
                     column(6,
                            # Non-editable
                            h5("Control Arm Response Rate (pc)"),
                            verbatimTextOutput("dpp_analysis_pc_display")
                     )
                   ),
                   wellPanel(
                     h4("Bayesian Borrowing"),
                     fluidRow(
                       column(6, numericInput("dpp_analysis_nche", "Equivalent number of patients borrowed from historical study (nche)", value = 41, min = 1))
                     ),
                     numericInput("dpp_analysis_delta_threshold", "Delta Threshold (0-0.5)", value = 0.1, min = 0, max = 0.5, step = 0.01),
                     selectInput("dpp_analysis_method", "Dynamic Control Method", choices = c("Empirical Bayes", "Bayesian p", "Generalized BC", "JSD")),
                     fluidRow(
                       column(6, uiOutput("dpp_analysis_theta_ui"))
                     ),
                     fluidRow(
                       column(6, uiOutput("dpp_analysis_eta_ui"))
                     )
                   )
                 )
          )
        ),
        actionButton("run_dpp_analysis", "Run DPP Analysis", class = "btn-primary"),
        br(), br(),
        div(id = "dpp_analysis_status_message", style = "text-align: center;")
      ),
      tabPanel(
        "Statistical Analysis Results (DPP)",
        h4("DPP Analysis Results"),
        verbatimTextOutput("dppAnalysisResult"),
        tags$p(style = "color: red; font-size: 1.2em;", tags$strong("Note: For a detailed explanation of the parameters, please check the 'Library' tab.")),
        h4("Posterior Distribution of Response Rate in Experimental Treatment, Concurrent Control, and  Hybrid Control"),
        plotOutput("plotDPP"),

        hr(),

        # Summary
        h3("Summary"),
        fluidRow(
          column(4,
                 wellPanel(
                   h4("1. Concurrent"),
                   verbatimTextOutput("concurrent_summary")
                 )
          ),
          column(4,
                 wellPanel(
                   h4("2. Hybrid Control"),
                   verbatimTextOutput("hybrid_control_summary")
                 )
          ),
          column(4,
                 wellPanel(
                   h4("3. Experimental"),
                   verbatimTextOutput("experimental_summary")
                 )
          )
        ),

        hr(),

        # Statistical Conclusion
        wellPanel(
          h4("Statistical Conclusion"),
          verbatimTextOutput("statistical_conclusion")
        )
      )
    )
  ),

  tabPanel(
    "SAM Prior",
    shinyjs::useShinyjs(),
    tabsetPanel(
      id = "sam_tabset_panel",
      tabPanel(
        "Library",
        fluidRow(
          column(4,
                 wellPanel(
                   h4("SAM Prior Plot Terms"),
                   actionLink("sam_plot_library_term_alpha_hist", "alpha_hist"), br(),
                   actionLink("sam_plot_library_term_beta_hist", "beta_hist"), br(),
                   actionLink("sam_plot_library_term_n_control", "n_control"), br(),
                   actionLink("sam_plot_library_term_p_control", "p_control"), br(),
                   actionLink("sam_plot_library_term_nf_alpha", "nf_alpha"), br(),
                   actionLink("sam_plot_library_term_nf_beta", "nf_beta"), br(),
                   actionLink("sam_plot_library_term_delta", "delta")
                 ),
                 wellPanel(
                   h4("SAM Weight Terms"),
                   actionLink("sam_weight_library_term_alpha_hist_w", "alpha_hist"), br(),
                   actionLink("sam_weight_library_term_beta_hist_w", "beta_hist"), br(),
                   actionLink("sam_weight_library_term_n_control_w", "n_control"), br(),
                   actionLink("sam_weight_library_term_r_control_w", "r_control"), br(),
                   actionLink("sam_weight_library_term_delta_w", "delta"), br(),
                   actionLink("sam_weight_library_term_method_w", "method"), br(),
                   actionLink("sam_weight_library_term_prior_odds", "prior_odds")
                 )
          ),
          column(8,
                 wellPanel(
                   h4("Definition"),
                   uiOutput("sam_plot_definition_output")
                 )
          )
        )
      ),
      tabPanel(
        "SAM Prior Design",
        tabsetPanel(
          id = "sam_design_tabset_panel",
          tabPanel(
            "Inputs",
            fluidRow(
              column(12,
                     wellPanel(
                       h4("1. Used for both SAM Prior Plot and SAM Weight"),
                       fluidRow(
                         column(6, numericInput("alpha_hist", "Alpha (a) for Beta(a, b)", value = 40, min = 0)),
                         column(6, numericInput("beta_hist", "Beta (b) for Beta(a, b)", value = 60, min = 0))
                       ),
                       fluidRow(
                         column(6, numericInput("n_control", "Number of patients (n)", value = 60, min = 1)),
                         column(6, numericInput("delta", "Clinically meaningful difference (delta)", value = 0.15, step = 0.01))
                       )
                     )
              )
            ),
            fluidRow(
              column(6,
                     wellPanel(
                       h4("2. SAM Prior Plot Specific Inputs"),
                       numericInput("p_control", "Simulated Control Rate", value = 0.42, min = 0, max = 1),
                       h4("Non-informative Prior"),
                       fluidRow(
                         column(6, numericInput("nf_alpha", "Alpha for Beta(a, b)", value = 1, min = 0)),
                         column(6, numericInput("nf_beta", "Beta (b) for Beta(a, b)", value = 1, min = 0))
                       )
                     )
              ),
              column(6,
                     wellPanel(
                       h4("3. SAM Weight Specific Inputs"),
                       numericInput("r_control", "Number of Responses (r)", value = 25, min = 0),
                       selectInput("sam_method", "Weight method", choices = c("LRT", "PPR")),
                       numericInput("prior_odds", "Prior Odds H0 vs H1 (PPR only)", value = 1)
                     )
              )
            ),
            actionButton("run_sam", "Generate SAM Prior Plot", class = "btn-primary"),
            actionButton("run_sam_weight", "Calculate SAM Weight", class = "btn-primary"),
            br(), br(),
            div(id = "sam_plot_status_message", style = "text-align: center;"),
            div(id = "sam_weight_status_message", style = "text-align: center;")
          ),
          tabPanel(
            "Results",
            fluidRow(
              column(6,
                     wellPanel(
                       h4("SAM Prior Plot"),
                       plotOutput("samPlot")
                     )
              ),
              column(6,
                     wellPanel(
                       h4("SAM Weight Result"),
                       verbatimTextOutput("samWeightResult")
                     )
              )
            ),
            fluidRow(
              column(12,
                     wellPanel(

                       htmlOutput("samPriorConclusion")
                     )
              )
            )
          )
        )
      )
    )
  ),
  tabPanel(
    "Fisher's Exact Method",
    tabsetPanel(
      id = "fisher_tabset_panel",
      tabPanel(
        "Library",
        fluidRow(
          column(4,
                 wellPanel(
                   h4("Terms"),
                   actionLink("fisher_term_Yc", "Yc (fisher)"), br(), actionLink("fisher_term_nc_fisher", "nc (fisher)"), br(),
                   actionLink("fisher_term_Yt", "Yt (fisher)"), br(), actionLink("fisher_term_nt_fisher", "nt (fisher)"), br(),
                   actionLink("fisher_term_alternative", "alternative (fisher)"), br(), actionLink("fisher_term_pc_bound", "pc (fisher.bound)"), br(),
                   actionLink("fisher_term_nc_bound", "nc (fisher.bound)"), br(), actionLink("fisher_term_nt_bound", "nt (fisher.bound)"), br(),
                   actionLink("fisher_term_alpha_bound", "alpha (fisher.bound)"), br(), actionLink("fisher_term_pt_power", "pt (fisher.power)"), br(),
                   actionLink("fisher_term_nt_power", "nt (fisher.power)"), br(), actionLink("fisher_term_pc_power", "pc (fisher.power)"), br(),
                   actionLink("fisher_term_nc_power", "nc (fisher.power)"), br(), actionLink("fisher_term_alpha_power", "alpha (fisher.power)"), br(),
                   actionLink("fisher_term_nsim_power", "nsim (fisher.power)"), br(), actionLink("fisher_term_seed_power", "seed (fisher.power)")
                 )
          ),
          column(8,
                 wellPanel(
                   h4("Definition"),
                   uiOutput("fisher_definition_output")
                 )
          )
        )
      ),
      tabPanel(
        "Power",
        tabsetPanel(
          id = "fisher_power_tabset_panel",
          tabPanel(
            "Design Inputs",
            fluidRow(
              column(6,
                     wellPanel(
                       h4("Inputs for Fisher Power Calculation"),
                       numericInput("fp_pt", "Probability of success in experimental arm (pt)", value = 0.5, min = 0, max = 1, step = 0.01),
                       numericInput("fp_nt", "Number of subjects in experimental arm (nt)", value = 40, min = 1),
                       numericInput("fp_pc", "Probability of success in control arm (pc)", value = 0.3, min = 0, max = 1, step = 0.01),
                       numericInput("fp_nc", "Number of subjects in control arm (nc)", value = 40, min = 1)
                     )
              ),
              column(6,
                     wellPanel(
                       h4("Additional Inputs"),
                       numericInput("fp_alpha", "One sided type I error rate (alpha)", value = 0.1, min = 0.001, max = 0.5),
                       numericInput("fp_nsim", "Number of replications (nsim)", value = 100000, min = 100, step = 100),
                       numericInput("fp_seed", "Seed for simulations (seed)", value = 2000, min = 1)
                     )
              )
            ),
            actionButton("run_fisher_power", "Calculate Fisher Power", class = "btn-primary"),
            br(), br(),
            div(id = "fisher_power_status_message", style = "text-align: center;")
          ),
          tabPanel(
            "Design Results",
            fluidRow(
              column(6,
                     wellPanel(
                       h4("Fisher Power Result"),
                       verbatimTextOutput("fisherPowerResult")
                     )
              )
            )
          )
        )
      ),
      tabPanel(
        "Bound",
        tabsetPanel(
          id = "fisher_bound_tabset_panel",
          tabPanel(
            "Design Inputs",
            fluidRow(
              column(6,
                     wellPanel(
                       h4("Inputs for Fisher Bound Calculation"),
                       numericInput("fb_pc", "Response rate for control arm (pc)", value = 0.3, min = 0, max = 1, step = 0.01),
                       numericInput("fb_nc", "Number of patients in control arm (nc)", value = 40, min = 1)
                     )
              ),
              column(6,
                     wellPanel(
                       h4("Additional Inputs"),
                       numericInput("fb_nt", "Number of patients in experimental arm (nt)", value = 40, min = 1),
                       numericInput("fb_alpha", "P-value threshold for significance (alpha)", value = 0.1, min = 0.001, max = 0.5)
                     )
              )
            ),
            actionButton("run_fisher_bound", "Calculate Fisher Bound", class = "btn-primary"),
            br(), br(),
            div(id = "fisher_bound_status_message", style = "text-align: center;")
          ),
          tabPanel(
            "Design Results",
            fluidRow(
              column(6,
                     wellPanel(
                       h4("Fisher Bound Result"),
                       wellPanel(h5("Contingency Table (M)"), verbatimTextOutput("fisherBoundM")),
                       wellPanel(h5("P-value at Boundary (p)"), verbatimTextOutput("fisherBoundP")),
                       wellPanel(h5("Responders for Control Arm (rc)"), verbatimTextOutput("fisherBoundRc")),
                       wellPanel(h5("Sample Size in Control Arm (nc)"), verbatimTextOutput("fisherBoundNc")),
                       wellPanel(h5("Response Rate for Experimental Arm (rt)"), verbatimTextOutput("fisherBoundRt")),
                       wellPanel(h5("Sample Size in Experimental Arm (nt)"), verbatimTextOutput("fisherBoundNt")),
                       wellPanel(h5("Minimum Detectable Difference (delta)"), verbatimTextOutput("fisherBoundDelta"))
                     )
              )
            ),
            fluidRow(
              column(12,
                     wellPanel(
                       htmlOutput("fisherBoundConclusion"),
                       htmlOutput("fisherBoundStatisticalConclusion")
                     )
              )
            )
          )
        )
      )
    )
  )
))
