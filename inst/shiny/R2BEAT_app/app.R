# =============================================================================
# R2BEAT - applicazione unica (self-contained)
# Il launcher e i tre workflow (optim / wf1 / wf2) sono incorporati qui.
# Nessuna app esterna viene piu' richiamata: il workflow scelto viene servito
# dalla stessa app tramite la query string  ?app=optim | ?app=wf1 | ?app=wf2.
# =============================================================================

options(shiny.maxRequestSize = 100 * 1024^2)  # 100 MB upload limit

library(shiny)
library(shinyjs)
library(R2BEAT)
library(genalg)
library(DT)
library(ggplot2)
library(shinydashboard)

`%||%` <- function(a, b) if (!is.null(a)) a else b


# ============================ WORKFLOW: OPTIM ============================



# ─── HELPERS ──────────────────────────────────────────────────────────────────


alloc_scatter <- function(result, title) {
  alloc <- result$alloc
  alloc <- alloc[-nrow(alloc), ]
  rng <- c(0, max(alloc$ALLOC, alloc$PROP, na.rm = TRUE))
  fit <- lm(ALLOC ~ PROP, data = alloc)
  r2  <- summary(fit)$r.squared
  coefs <- coef(fit)

  plot(ALLOC ~ PROP, data = alloc,
       main  = title,
       xlab  = "Proportional allocation",
       ylab  = "Optimal allocation",
       xlim  = rng, ylim = rng,
       pch   = 21, bg = "steelblue", col = "black", cex = 1.2)
  abline(a = 0, b = 1, col = "red", lwd = 2, lty = 2)
  abline(fit, col = "darkgreen", lwd = 2)
  legend_x <- rng[2] * 0.05
  legend_y <- rng[2] * 0.95
  text(x = legend_x, y = legend_y,
       labels = sprintf("y = %.2f + %.2f x\nR\u00b2 = %.3f", coefs[1], coefs[2], r2),
       col = "darkgreen", adj = c(0, 1), cex = 0.9)
  text(x = rng[2] * 0.8, y = rng[2] * 0.05,
       labels = "y = x", col = "red", cex = 0.85)
}

round_cv_df <- function(cv_df, digits = 3) {
  cv_cols <- grep("^CV[0-9]+$", names(cv_df))
  cv_df[, cv_cols] <- lapply(cv_df[, cv_cols, drop = FALSE], round, digits = digits)
  cv_df
}

deparse_cv_df <- function(cv_df) {
  n_cv  <- ncol(cv_df) - 1
  doms  <- paste0('"', cv_df$DOM, '"', collapse = ", ")
  lines <- c("as.data.frame(list(", paste0("  DOM = c(", doms, "),"))
  for (j in seq_len(n_cv)) {
    vals  <- paste(cv_df[[paste0("CV", j)]], collapse = ", ")
    comma <- if (j < n_cv) "," else ""
    lines <- c(lines, paste0("  CV", j, " = c(", vals, ")", comma))
  }
  paste(c(lines, "))"), collapse = "\n")
}

# ─── UI ────────────────────────────────────────────────────────────────────────

optim_ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "R2BEAT — CV Optimisation"),

  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("1. Data",           tabName = "tab_data",   icon = icon("database")),
      menuItem("2. CV Constraints", tabName = "tab_cv",     icon = icon("bullseye")),
      menuItem("3. Pareto",         tabName = "tab_pareto", icon = icon("chart-line")),
      menuItem("4. GA Adjustment",  tabName = "tab_ga",     icon = icon("dna")),
      menuItem("5. Final Solution", tabName = "tab_final",  icon = icon("flag-checkered")),
      menuItem("6. R Script",       tabName = "tab_script", icon = icon("file-code"))
    ),
    div(style = "padding:12px 15px; border-top:1px solid #444; margin-top:10px;",
      actionButton("btn_home", " Menu principale", icon = icon("home"),
                   class = "btn-default btn-block",
                   style = "font-weight:bold; margin-bottom:8px;",
                   onclick = "window.location.search='';"),
      actionButton("btn_help", " Help Guide", icon = icon("question-circle"),
                   class = "btn-info btn-block", style = "font-weight:bold;")
    )
  ),

  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .box-body .dataTables_wrapper { overflow-x: auto; }
      pre { max-height: 420px; overflow-y: auto; }
      .cv-grid input[type=number] { width: 72px; }
      .cv-grid td { vertical-align: middle !important; padding: 3px 5px !important; }
      .cv-grid th { text-align: center; padding: 4px 5px; }
    "))),
    tabItems(

      # ── TAB 1: DATA ──────────────────────────────────────────────────────────
      tabItem(tabName = "tab_data",
        fluidRow(
          box(title = "Data Source", width = 4, status = "primary", solidHeader = TRUE,
            radioButtons("data_source", "Choose strata input:",
              choices  = c("Use bundled strata.RData (demo)" = "demo",
                           "Upload .RData file"              = "upload",
                           "Upload delimited file (CSV/TSV)" = "delimited"),
              selected = "demo"),
            conditionalPanel("input.data_source == 'upload'",
              fileInput("upload_rdata", "Upload .RData file",
                        accept = c(".RData", ".rdata", ".rda"))
            ),
            conditionalPanel("input.data_source == 'delimited'",
              fileInput("upload_csv", "Upload CSV/TSV file",
                        accept = c(".csv", ".tsv", ".txt")),
              fluidRow(
                column(6,
                  radioButtons("csv_sep", "Separator:",
                    choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                    selected = ",", inline = FALSE)
                ),
                column(6,
                  checkboxInput("csv_header", "Header row", value = TRUE),
                  checkboxInput("csv_dec_comma", "Decimal comma ( , )", value = FALSE)
                )
              )
            ),
            actionButton("btn_load", "Load Data", class = "btn-primary btn-block",
                         icon = icon("upload"))
          ),
          box(title = "Detected Structure", width = 4, status = "info",
            verbatimTextOutput("txt_structure")
          ),
          box(title = "Target Variable Names (optional)", width = 4, status = "warning",
            p("Rename target variables for display. Leave blank for default V1, V2, …"),
            uiOutput("ui_var_names"),
            actionButton("btn_set_varnames", "Apply Names",
                         class = "btn-warning btn-block", icon = icon("check"))
          )
        ),
        fluidRow(
          box(title = "Strata Data Preview", width = 12, status = "info",
            DTOutput("tbl_strata_head")
          )
        )
      ),

      # ── TAB 2: CV CONSTRAINTS ────────────────────────────────────────────────
      tabItem(tabName = "tab_cv",
        fluidRow(
          box(title = "Initial CV Constraints", width = 5,
              status = "primary", solidHeader = TRUE,
            p("Enter the maximum Coefficient of Variation (CV) for each",
              "domain row and target variable."),
            fluidRow(
              column(7,
                numericInput("cv_equal_val",
                             "Set all cells equal to:", value = 0.05,
                             min = 0.001, max = 1, step = 0.005)
              ),
              column(5,
                br(),
                actionButton("btn_set_equal", "Apply", class = "btn-default btn-block")
              )
            ),
            hr(),
            div(class = "cv-grid", uiOutput("ui_cv_table")),
            br(),
            actionButton("btn_run_init", "Run Initial Allocation",
                         class = "btn-success btn-block", icon = icon("play"))
          ),
          box(title = "Initial Allocation Results", width = 7, status = "info",
            fluidRow(
              column(6,
                h4("Total sample size"),
                verbatimTextOutput("txt_init_total")
              ),
              column(6,
                h4("Cor(Optimal, Proportional)"),
                verbatimTextOutput("txt_init_cor")
              )
            ),
            hr(),
            h4("Allocation by stratum"),
            DTOutput("tbl_init_alloc")
          )
        )
      ),

      # ── TAB 3: PARETO ────────────────────────────────────────────────────────
      tabItem(tabName = "tab_pareto",
        fluidRow(
          box(title = "Pareto Optimisation Settings", width = 4,
              status = "primary", solidHeader = TRUE,
            numericInput("pareto_target",  "Target sample size:",  50000, 100, 1e8),
            numericInput("pareto_tol",     "Tolerance:",               5,   1, 1000),
            numericInput("pareto_maxiter", "Max iterations:",          50,   5,  500),
            numericInput("pareto_maxsame", "Max same iterations:",     20,   5,  200),
            hr(),
            checkboxInput("use_caps", "Enable CV caps", value = FALSE),
            conditionalPanel("input.use_caps == true",
              p(em("Maximum allowable CV per domain / variable:")),
              div(class = "cv-grid", uiOutput("ui_caps_table"))
            ),
            br(),
            actionButton("btn_run_pareto", "Run Pareto",
                         class = "btn-danger btn-block", icon = icon("rocket")),
            br(),
            p(em("Requires Step 2 to be completed first."))
          ),
          box(title = "Pareto Results", width = 8, status = "info",
            fluidRow(
              column(6,
                h4("Total sample size"),
                verbatimTextOutput("txt_pareto_total")
              ),
              column(6,
                h4("Cor(Optimal, Proportional)"),
                verbatimTextOutput("txt_pareto_cor")
              )
            ),
            hr(),
            h4("Pareto CV constraints"),
            DTOutput("tbl_pareto_cvs"),
            hr(),
            fluidRow(
              column(6,
                h4("Convergence: sample size vs iteration"),
                plotOutput("plot_pareto_conv",   height = "280px")
              ),
              column(6,
                h4("Pareto Bisection Search for Feasibility"),
                plotOutput("plot_pareto_bisect", height = "280px")
              )
            ),
            hr(),
            h4("Optimal vs Proportional allocation"),
            plotOutput("plot_pareto_alloc", height = "270px")
          )
        )
      ),

      # ── TAB 4: GA ADJUSTMENT ─────────────────────────────────────────────────
      tabItem(tabName = "tab_ga",
        fluidRow(
          box(title = tags$span("Genetic Algorithm Settings",
                tags$small(tags$span(class = "label label-success",
                                     style = "margin-left:8px;", "optional"))),
              width = 4, status = "primary", solidHeader = TRUE,
            div(class = "alert alert-info", style = "padding:8px; margin-bottom:10px;",
              icon("info-circle"),
              tags$b(" This step is optional."),
              " Skip directly to Step 5 if you only need the Pareto solution."
            ),
            numericInput("ga_iters",      "Iterations:",                 100, 10, 5000),
            numericInput("ga_popsize",    "Population size:",             50, 10,  500),
            numericInput("ga_limit",      "Correlation target (limit):", 0.95, 0.1, 0.999),
            numericInput("ga_min_factor", "CV lower bound (× Pareto CV):", 0.5, 0.1, 1.0),
            numericInput("ga_max_factor", "CV upper bound (× Pareto CV):", 1.5, 1.0, 5.0),
            hr(),
            p(em("The GA searches for CV values that make the optimal",
                 "allocation more proportional while keeping sample size feasible.",
                 "This may take several minutes.")),
            actionButton("btn_run_ga", "Run GA",
                         class = "btn-warning btn-block", icon = icon("dna")),
            br(),
            p(em("Requires Step 3 (Pareto) to be completed first."))
          ),
          box(title = "GA Results", width = 8, status = "info",
            fluidRow(
              column(6,
                h4("Total sample size"),
                verbatimTextOutput("txt_ga_total")
              ),
              column(6,
                h4("Cor(Optimal, Proportional)"),
                verbatimTextOutput("txt_ga_cor")
              )
            ),
            hr(),
            h4("GA CV constraints"),
            DTOutput("tbl_ga_cvs"),
            hr(),
            fluidRow(
              column(6,
                h4("Region of variation for CVs"),
                plotOutput("plot_ga_region", height = "280px")
              ),
              column(6,
                h4("Best and mean evaluation value"),
                plotOutput("plot_ga_fitness", height = "280px")
              )
            )
          )
        )
      ),

      # ── TAB 5: FINAL SOLUTION ────────────────────────────────────────────────
      tabItem(tabName = "tab_final",
        fluidRow(
          box(title = "Final Solution Controls", width = 3,
              status = "primary", solidHeader = TRUE,
            uiOutput("ui_final_controls"),
            hr(),
            downloadButton("btn_download_final_cv",
                           "Download final CVs (CSV)",
                           class = "btn-success btn-block")
          ),
          box(title = "Comparison Summary", width = 9, status = "info",
            h4("Sample size and proportionality by solution"),
            tableOutput("tbl_comparison"),
            hr(),
            h4("Final CV constraints"),
            DTOutput("tbl_final_cvs")
          )
        ),
        fluidRow(
          box(title = "Pareto solution: Optimal vs Proportional",
              width = 6,
            plotOutput("plot_alloc_before", height = "300px")
          ),
          box(uiOutput("ui_final_alloc_title"), width = 6,
            plotOutput("plot_alloc_after", height = "300px")
          )
        ),
        fluidRow(
          box(title = "CV Comparison by Domain", width = 12,
            uiOutput("ui_domain_selector"),
            plotOutput("plot_cv_compare", height = "360px")
          )
        )
      ),

      # ── TAB 6: R SCRIPT ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_script",
        fluidRow(
          box(title = "Generated R Script", width = 12, status = "info",
            p("Reproducible script recording all steps executed in this session."),
            verbatimTextOutput("txt_script_preview")
          )
        ),
        fluidRow(
          box(title = "Download", width = 4, status = "primary", solidHeader = TRUE,
            downloadButton("btn_download_script", "Download R Script",
                           class = "btn-block")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

# ─── SERVER ────────────────────────────────────────────────────────────────────

optim_server <- function(input, output, session) {

  rv <- reactiveValues(
    strata         = NULL,
    n_dom          = NULL,
    n_var          = NULL,
    dom_labels     = NULL,   # "DOM1", "DOM2", ...
    var_names      = NULL,   # "V1", "V2", ... or user-defined labels

    cv_init        = NULL,
    init_result    = NULL,

    cv_pareto      = NULL,
    pareto_out     = NULL,
    pareto_result  = NULL,

    cv_ga          = NULL,
    ga_obj         = NULL,
    ga_result      = NULL,
    ga_vett        = NULL,   # Pareto CV vector used as GA starting point
    ga_vettmin     = NULL,   # lower search bound
    ga_vettmax     = NULL,   # upper search bound

    cv_final       = NULL,
    final_result   = NULL,

    script_lines   = character(0)
  )

  append_script <- function(...) {
    txt <- paste0(...)
    rv$script_lines <- c(rv$script_lines, txt)
  }

  # ── HELP MODAL ────────────────────────────────────────────────────────────────
  observeEvent(input$btn_help, {
    showModal(modalDialog(
      title     = tags$span(icon("question-circle"), " R2BEAT CV Optimisation — Guide"),
      size      = "l",
      easyClose = TRUE,
      footer    = modalButton("Close"),
      tabsetPanel(
        tabPanel(tags$span(icon("database"), " 1. Data"),
          br(),
          tags$h4("Purpose"),
          tags$p("Load a pre-built strata object (output of R2BEAT stratification) that",
                 "contains population means, standard deviations, and sizes for each stratum."),
          tags$h4("Options"),
          tags$ul(
            tags$li(tags$b("Demo dataset:"), " Loads the bundled ", tags$code("strata.RData"),
                    " (24 strata, 3 domain levels, 4 target variables)."),
            tags$li(tags$b("Upload:"), " Any ", tags$code(".RData"), " file containing",
                    " a data frame named ", tags$code("strata"), " with columns ",
                    tags$code("DOM1…DOMk"), ", ", tags$code("M1…Mv"), ", ",
                    tags$code("S1…Sv"), ", ", tags$code("N"), ", ", tags$code("CENS"),
                    ", ", tags$code("COST"), ".")
          ),
          tags$h4("Variable Names"),
          tags$p("Optionally assign meaningful names to the target variables",
                 " (e.g., 'active', 'unemployed'). These are used for plot labels only.")
        ),
        tabPanel(tags$span(icon("bullseye"), " 2. CV Constraints"),
          br(),
          tags$h4("Purpose"),
          tags$p("Define the initial precision constraints — one CV per domain level and",
                 " target variable. The app auto-detects the number of domain columns",
                 " (DOM1, DOM2, …) and variables (M1/S1, M2/S2, …) in the strata object."),
          tags$h4("Workflow"),
          tags$ol(
            tags$li("Enter a starting CV (e.g., 0.05) in the ", tags$b("Set all equal to"),
                    " box and click ", tags$b("Apply"), " to fill every cell at once."),
            tags$li("Adjust individual cells if needed."),
            tags$li("Click ", tags$b("Run Initial Allocation"), " to see the Bethel-style",
                    " sample size and the correlation between optimal and proportional allocation.")
          )
        ),
        tabPanel(tags$span(icon("chart-line"), " 3. Pareto"),
          br(),
          tags$h4("Purpose"),
          tags$p("Find a set of CV constraints that makes the optimal Bethel allocation",
                 " affordable: ", tags$code("pareto()"), " scales the constraints by a",
                 " scalar ", tags$em("t"), " (bisection) until the total sample matches",
                 " the requested ", tags$b("Target sample size"), "."),
          tags$h4("Parameters"),
          tags$ul(
            tags$li(tags$b("Target sample size:"), " Budget-consistent total sample you can afford."),
            tags$li(tags$b("Tolerance:"), " Acceptable deviation from the target (±units)."),
            tags$li(tags$b("Max / Max-same iterations:"), " Stopping criteria for the bisection."),
            tags$li(tags$b("CV caps:"), " Optional per-domain per-variable upper bounds on CVs.",
                    " Prevents any CV from being relaxed beyond a maximum.")
          ),
          tags$h4("Output"),
          tags$p("The Pareto CV table (", tags$code("out$expectedCV"), "), convergence plot,",
                 " and allocation scatter are displayed. The correlation between optimal and",
                 " proportional allocation may still be low — the GA step addresses this.")
        ),
        tabPanel(tags$span(icon("dna"), " 4. GA Adjustment"),
          br(),
          tags$h4("Purpose"),
          tags$p("Run a real-valued genetic algorithm (", tags$code("rbga()"),
                 " from package ", tags$code("genalg"), ") to modify the Pareto CVs so",
                 " that the resulting optimal allocation is more proportional, without",
                 " imposing an exact sample size target."),
          tags$h4("Parameters"),
          tags$ul(
            tags$li(tags$b("Iterations / Population size:"), " Standard GA parameters.",
                    " More iterations → better solution, longer runtime."),
            tags$li(tags$b("Correlation target:"), " The GA's fitness function maximises",
                    " Cor(Optimal, Proportional). Solutions with |r| > target + 0.005 are",
                    " penalised to prevent overshooting."),
            tags$li(tags$b("CV lower/upper bound:"), " The search space for each CV is",
                    " [Pareto_CV × lower, Pareto_CV × upper].")
          ),
          tags$h4("Note"),
          tags$p("The GA ignores the sample size target — the resulting sample size may",
                 " differ from the target. Run Step 5 (Final Pareto) to rescale back.")
        ),
        tabPanel(tags$span(icon("flag-checkered"), " 5. Final Solution"),
          br(),
          tags$h4("Purpose"),
          tags$p("Apply Pareto a second time — now on the GA-adjusted CVs — to rescale",
                 " them back to the target sample size while preserving the improved",
                 " proportionality shape found by the GA."),
          tags$h4("Outputs"),
          tags$ul(
            tags$li(tags$b("Comparison table:"), " Sample size and Cor(Opt, Prop) for",
                    " all four solutions (Initial, Pareto, GA, Final)."),
            tags$li(tags$b("Allocation scatter plots:"), " Pareto (before GA) vs Final (after GA)."),
            tags$li(tags$b("CV bar charts:"), " For each domain, the three sets of CVs",
                    " (Pareto / GA / Final) shown side by side.")
          )
        )
      )
    ))
  })

  # ── STEP 1: LOAD DATA ────────────────────────────────────────────────────────

  observeEvent(input$btn_load, {
    tryCatch({
      if (input$data_source == "demo") {
        env <- new.env()
        load("strata.RData", envir = env)
        if (!exists("strata", envir = env))
          stop("Bundled strata.RData does not contain an object named 'strata'.")
        strata_obj <- env$strata
        data_expr  <- 'load("strata.RData")'

      } else if (input$data_source == "upload") {
        req(input$upload_rdata)
        env <- new.env()
        load(input$upload_rdata$datapath, envir = env)
        if (!exists("strata", envir = env))
          stop("The RData file must contain an object named 'strata'.")
        strata_obj <- env$strata
        data_expr  <- paste0('load("', input$upload_rdata$name, '")')

      } else {   # delimited
        req(input$upload_csv)
        dec <- if (isTRUE(input$csv_dec_comma)) "," else "."
        strata_obj <- read.table(
          input$upload_csv$datapath,
          header    = isTRUE(input$csv_header),
          sep       = input$csv_sep,
          dec       = dec,
          stringsAsFactors = FALSE
        )
        fname     <- input$upload_csv$name
        data_expr <- paste0(
          'strata <- read.table("', fname, '",\n',
          '  header = ', isTRUE(input$csv_header), ', sep = "', input$csv_sep, '",\n',
          '  dec = "', dec, '", stringsAsFactors = FALSE)'
        )
      }

      dom_cols <- sort(grep("^DOM[0-9]+$", names(strata_obj), value = TRUE))
      m_cols   <- sort(grep("^M[0-9]+$",   names(strata_obj), value = TRUE))

      if (length(dom_cols) == 0) stop("No DOMx columns found in the strata object.")
      if (length(m_cols)   == 0) stop("No Mx (mean) columns found in the strata object.")

      rv$strata     <- strata_obj
      rv$n_dom      <- length(dom_cols)
      rv$n_var      <- length(m_cols)
      rv$dom_labels <- paste0("DOM", seq_along(dom_cols))
      rv$var_names  <- paste0("V",   seq_along(m_cols))

      rv$script_lines <- c(
        "# R script generated by R2BEAT CV Optimisation Shiny App",
        "library(R2BEAT)",
        'library("genalg")',
        "",
        "# 1. Load strata data",
        data_expr,
        ""
      )

      showNotification("Strata data loaded successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_structure <- renderText({
    req(rv$strata)
    paste0(
      "Strata rows : ", nrow(rv$strata), "\n",
      "Domain cols : ", rv$n_dom, "  (", paste(rv$dom_labels, collapse = ", "), ")\n",
      "Target vars : ", rv$n_var, "  (", paste(rv$var_names,  collapse = ", "), ")\n",
      "Total N     : ", format(sum(rv$strata$N), big.mark = ","), "\n"
    )
  })

  output$tbl_strata_head <- renderDT({
    req(rv$strata)
    datatable(head(rv$strata, 12), options = list(scrollX = TRUE, pageLength = 6))
  })

  output$ui_var_names <- renderUI({
    req(rv$n_var)
    lapply(seq_len(rv$n_var), function(i) {
      textInput(paste0("varname_", i),
                paste0("Variable ", i, " name:"),
                value = rv$var_names[i], width = "100%")
    })
  })

  observeEvent(input$btn_set_varnames, {
    req(rv$n_var)
    nms <- vapply(seq_len(rv$n_var), function(i) {
      v <- trimws(input[[paste0("varname_", i)]] %||% "")
      if (nchar(v) == 0) paste0("V", i) else v
    }, character(1))
    rv$var_names <- nms
    showNotification("Variable names updated.", type = "message")
  })

  # ── STEP 2: CV CONSTRAINTS ───────────────────────────────────────────────────

  output$ui_cv_table <- renderUI({
    req(rv$n_dom, rv$n_var)
    vn <- isolate(rv$var_names) %||% paste0("V", seq_len(rv$n_var))

    header <- tags$tr(
      tags$th("Domain"),
      lapply(seq_len(rv$n_var), function(j) tags$th(vn[j]))
    )
    rows <- lapply(seq_len(rv$n_dom), function(i) {
      tags$tr(
        tags$td(tags$strong(rv$dom_labels[i])),
        lapply(seq_len(rv$n_var), function(j) {
          tags$td(
            numericInput(paste0("cv_", i, "_", j), NULL,
                         value = 0.05, min = 1e-4, max = 1, step = 0.005, width = "72px")
          )
        })
      )
    })
    tags$table(class = "table table-condensed table-bordered",
               tags$thead(header), tags$tbody(rows))
  })

  observeEvent(input$btn_set_equal, {
    req(rv$n_dom, rv$n_var)
    val <- input$cv_equal_val
    req(!is.na(val))
    for (i in seq_len(rv$n_dom))
      for (j in seq_len(rv$n_var))
        updateNumericInput(session, paste0("cv_", i, "_", j), value = val)
  })

  collect_cv <- function() {
    n_dom <- rv$n_dom; n_var <- rv$n_var
    if (is.null(n_dom) || is.null(n_var)) return(NULL)
    cv_mat <- matrix(0.05, nrow = n_dom, ncol = n_var)
    for (i in seq_len(n_dom))
      for (j in seq_len(n_var)) {
        v <- input[[paste0("cv_", i, "_", j)]]
        if (!is.null(v) && length(v) == 1 && !is.na(v))
          cv_mat[i, j] <- as.numeric(v)
      }
    df <- data.frame(DOM = rv$dom_labels, stringsAsFactors = FALSE)
    for (j in seq_len(n_var)) df[[paste0("CV", j)]] <- cv_mat[, j]
    df
  }

  observeEvent(input$btn_run_init, {
    req(rv$strata)
    tryCatch({
      cv <- collect_cv()
      rv$cv_init <- cv
      withProgress(message = "Running initial allocation…", value = 0.5, {
        rv$init_result <- R2BEAT::beat.1st(rv$strata, cv)
      })
      append_script(paste0(
        "# 2. Initial CV constraints\n",
        "cv_equal <- ", deparse_cv_df(cv), "\n",
        "equal_alloc <- beat.1st(strata, cv_equal)\n",
        "cat('Initial total:', sum(equal_alloc$alloc$ALLOC[-nrow(equal_alloc$alloc)]), '\\n')\n",
        "cor(equal_alloc$alloc$ALLOC[-nrow(equal_alloc$alloc)],\n",
        "    equal_alloc$alloc$PROP[-nrow(equal_alloc$alloc)])\n"
      ))
      showNotification("Initial allocation complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_init_total <- renderText({
    req(rv$init_result)
    alloc <- rv$init_result$alloc
    format(round(sum(alloc$ALLOC[-nrow(alloc)])), big.mark = ",")
  })

  output$txt_init_cor <- renderText({
    req(rv$init_result)
    alloc <- rv$init_result$alloc[-nrow(rv$init_result$alloc), ]
    sprintf("%.4f", cor(alloc$ALLOC, alloc$PROP))
  })

  output$tbl_init_alloc <- renderDT({
    req(rv$init_result)
    datatable(rv$init_result$alloc,
              options = list(scrollX = TRUE, pageLength = 8))
  })

  # ── STEP 3: PARETO ───────────────────────────────────────────────────────────

  output$ui_caps_table <- renderUI({
    req(rv$n_dom, rv$n_var)
    vn <- isolate(rv$var_names) %||% paste0("V", seq_len(rv$n_var))

    header <- tags$tr(
      tags$th("Domain"),
      lapply(seq_len(rv$n_var), function(j) tags$th(vn[j]))
    )
    rows <- lapply(seq_len(rv$n_dom), function(i) {
      tags$tr(
        tags$td(tags$strong(rv$dom_labels[i])),
        lapply(seq_len(rv$n_var), function(j) {
          tags$td(
            numericInput(paste0("cap_", i, "_", j), NULL,
                         value = 0.10, min = 1e-4, max = 1, step = 0.005, width = "72px")
          )
        })
      )
    })
    tags$table(class = "table table-condensed table-bordered",
               style = "font-size:12px;",
               tags$thead(header), tags$tbody(rows))
  })

  collect_caps <- function() {
    if (!isTRUE(input$use_caps)) return(NULL)
    n_dom <- rv$n_dom; n_var <- rv$n_var
    if (is.null(n_dom) || is.null(n_var)) return(NULL)
    dom_vals <- character(0)
    var_vals <- character(0)
    cv_vals  <- numeric(0)
    for (i in seq_len(n_dom))
      for (j in seq_len(n_var)) {
        cap_val <- input[[paste0("cap_", i, "_", j)]]
        if (!is.null(cap_val) && length(cap_val) == 1 && !is.na(cap_val)) {
          dom_vals <- c(dom_vals, rv$dom_labels[i])
          var_vals <- c(var_vals, paste0("V", j))
          cv_vals  <- c(cv_vals,  as.numeric(cap_val))
        }
      }
    if (length(dom_vals) == 0) return(NULL)
    data.frame(DOM = dom_vals, VAR = var_vals, MAX_CV = cv_vals,
               stringsAsFactors = FALSE)
  }

  run_pareto_impl <- function(cv_input) {
    if (is.null(rv$strata)) stop("Strata data not loaded. Complete Step 1 first.")
    caps <- collect_caps()
    args <- list(
      strata        = rv$strata,
      current_cvs   = cv_input,
      target_size   = input$pareto_target,
      tolerance     = input$pareto_tol,
      beat1cv_fun   = beat.1cv_2,
      max_iter      = input$pareto_maxiter,
      max_same_iter = input$pareto_maxsame,
      verbose       = FALSE,
      plot          = FALSE
    )
    if (!is.null(caps)) args$cv_caps <- caps
    do.call(R2BEAT::pareto, args)
  }

  normalize_expectedCV <- function(cv_df) {
    # expectedCV columns come back as V1, V2, … — rename to CV1, CV2, …
    n_cv <- ncol(cv_df) - 1
    colnames(cv_df)[2:ncol(cv_df)] <- paste0("CV", seq_len(n_cv))
    cv_df
  }

  observeEvent(input$btn_run_pareto, {
    req(rv$cv_init)
    tryCatch({
      withProgress(message = "Running Pareto optimisation…", value = 0.5, {
        out <- run_pareto_impl(rv$cv_init)
      })
      rv$pareto_out    <- out
      cv_p             <- normalize_expectedCV(out$expectedCV)
      rv$cv_pareto     <- cv_p
      rv$pareto_result <- R2BEAT::beat.1st(rv$strata, cv_p)

      caps_lines <- ""
      if (isTRUE(input$use_caps)) {
        caps <- collect_caps()
        if (!is.null(caps)) {
          caps_lines <- paste0(
            "caps <- data.frame(\n",
            "  DOM    = c(", paste0('"', caps$DOM, '"', collapse = ", "), "),\n",
            "  VAR    = c(", paste0('"', caps$VAR, '"', collapse = ", "), "),\n",
            "  MAX_CV = c(", paste(caps$MAX_CV, collapse = ", "), ")\n",
            ")\n"
          )
        }
      }
      append_script(paste0(
        "\n# 3. Pareto optimisation\n",
        caps_lines,
        "out_pareto <- pareto(\n",
        "  strata        = strata,\n",
        "  current_cvs   = cv_equal,\n",
        "  target_size   = ", input$pareto_target, ",\n",
        "  tolerance     = ", input$pareto_tol, ",\n",
        "  beat1cv_fun   = beat.1cv_2,\n",
        "  max_iter      = ", input$pareto_maxiter, ",\n",
        "  max_same_iter = ", input$pareto_maxsame, ",\n",
        if (isTRUE(input$use_caps)) "  cv_caps       = caps,\n" else "",
        "  plot          = FALSE\n",
        ")\n",
        "cv_pareto <- out_pareto$expectedCV\n",
        "colnames(cv_pareto)[2:ncol(cv_pareto)] <- paste0('CV', 1:(ncol(cv_pareto)-1))\n",
        "pareto_alloc <- beat.1st(strata, cv_pareto)\n",
        "cat('Pareto total:', sum(pareto_alloc$alloc$ALLOC[-nrow(pareto_alloc$alloc)]), '\\n')\n",
        "cor(pareto_alloc$alloc$ALLOC[-nrow(pareto_alloc$alloc)],\n",
        "    pareto_alloc$alloc$PROP[-nrow(pareto_alloc$alloc)])\n"
      ))
      showNotification("Pareto complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Pareto error:", e$message), type = "error", duration = 8)
    })
  })

  output$tbl_pareto_cvs <- renderDT({
    req(rv$cv_pareto)
    datatable(round_cv_df(rv$cv_pareto), options = list(dom = "t", scrollX = TRUE))
  })

  output$txt_pareto_total <- renderText({
    req(rv$pareto_result)
    alloc <- rv$pareto_result$alloc
    format(round(sum(alloc$ALLOC[-nrow(alloc)])), big.mark = ",")
  })

  output$txt_pareto_cor <- renderText({
    req(rv$pareto_result)
    alloc <- rv$pareto_result$alloc[-nrow(rv$pareto_result$alloc), ]
    sprintf("%.4f", cor(alloc$ALLOC, alloc$PROP))
  })

  output$plot_pareto_conv <- renderPlot({
    req(rv$pareto_out)
    ch <- rv$pareto_out$convergence_history
    req(!is.null(ch$n_vals))
    df <- data.frame(
      iteration = seq_along(ch$n_vals),
      n         = ch$n_vals,
      phase     = ch$phase
    )
    ggplot(df, aes(x = iteration, y = n, colour = phase)) +
      geom_line(aes(group = 1), colour = "grey60") +
      geom_point(size = 2.5) +
      geom_hline(yintercept = input$pareto_target, linetype = "dashed",
                 colour = "red", linewidth = 0.8) +
      geom_hline(yintercept = input$pareto_target + input$pareto_tol,
                 linetype = "dotted", colour = "orange") +
      geom_hline(yintercept = input$pareto_target - input$pareto_tol,
                 linetype = "dotted", colour = "orange") +
      labs(x = "Iteration", y = "Sample size", colour = "Phase") +
      theme_minimal(base_size = 11)
  })

  output$plot_pareto_bisect <- renderPlot({
    req(rv$pareto_out)
    ch     <- rv$pareto_out$convergence_history
    req(!is.null(ch$t_vals), !is.null(ch$n_vals), !is.null(ch$phase))
    target <- input$pareto_target
    tol    <- input$pareto_tol

    df <- data.frame(
      t     = as.numeric(ch$t_vals),
      n     = as.numeric(ch$n_vals),
      phase = as.character(ch$phase),
      stringsAsFactors = FALSE
    )

    # Re-label: first Bracket-up point that goes below target → "Bracket found"
    brk_below <- which(df$phase == "Bracket-up" & df$n < target)
    if (length(brk_below) > 0)
      df$phase[brk_below[1]] <- "Bracket found"

    # Last Bisect → "Converged"
    bisect_rows <- which(df$phase == "Bisect")
    if (length(bisect_rows) > 0)
      df$phase[bisect_rows[length(bisect_rows)]] <- "Converged"

    phase_levels  <- c("Start", "Bracket-up", "Bracket found", "Bisect", "Converged")
    phase_colours <- c("Start"         = "#CC0000",
                       "Bracket-up"    = "#FFAAAA",
                       "Bracket found" = "#7EC8E3",
                       "Bisect"        = "#999999",
                       "Converged"     = "#1F3F7A")
    phase_sizes   <- c("Start"         = 3.5,
                       "Bracket-up"    = 2.5,
                       "Bracket found" = 3.0,
                       "Bisect"        = 2.0,
                       "Converged"     = 5.0)

    df$phase <- factor(df$phase, levels = phase_levels)

    conv_row  <- df[df$phase == "Converged", ]
    n_range   <- max(df$n) - min(df$n)
    t_range   <- max(df$t) - min(df$t)
    ann_nudge <- n_range * 0.07

    p <- ggplot(df, aes(x = t, y = n)) +
      geom_line(aes(group = 1), colour = "grey70", linewidth = 0.5) +
      geom_hline(yintercept = target,       colour = "steelblue",
                 linewidth = 0.9) +
      geom_hline(yintercept = target + tol, colour = "steelblue",
                 linewidth = 0.4, linetype = "dotted") +
      geom_hline(yintercept = target - tol, colour = "steelblue",
                 linewidth = 0.4, linetype = "dotted") +
      geom_point(aes(colour = phase, size = phase)) +
      scale_colour_manual(values = phase_colours, drop = FALSE) +
      scale_size_manual(values  = phase_sizes,    guide = "none") +
      annotate("text",
               x = min(df$t) + t_range * 0.03,
               y = target - n_range * 0.035,
               label = "Feasible", colour = "steelblue",
               size = 3, fontface = "italic", hjust = 0) +
      labs(title  = "Pareto Bisection Search for Feasibility",
           x      = "Scaling Factor (t)",
           y      = "Sample Size (n)",
           colour = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")

    if (nrow(conv_row) > 0) {
      n_bisect <- length(bisect_rows)
      p <- p +
        annotate("text",
                 x     = conv_row$t[1],
                 y     = conv_row$n[1] + ann_nudge,
                 label = sprintf("Converged\nt = %.4f\nn = %s",
                                 conv_row$t[1],
                                 format(round(conv_row$n[1]), big.mark = ",")),
                 colour = "#1F3F7A", size = 2.8, hjust = 0.5,
                 fontface = "italic")
    }
    p
  })

  output$plot_pareto_alloc <- renderPlot({
    req(rv$pareto_result)
    alloc_scatter(rv$pareto_result, "Pareto: Optimal vs Proportional allocation")
  })

  # ── STEP 4: GA ADJUSTMENT ────────────────────────────────────────────────────

  observeEvent(input$btn_run_ga, {
    req(rv$cv_pareto, rv$strata)

    strata_local    <- rv$strata
    cv_pareto_local <- rv$cv_pareto
    n_dom           <- rv$n_dom
    n_var           <- rv$n_var
    limit           <- input$ga_limit

    # CV variation range
    vett <- unlist(lapply(seq_len(n_dom), function(i) {
      as.numeric(cv_pareto_local[i, paste0("CV", seq_len(n_var))])
    }))
    vettmin <- vett * input$ga_min_factor
    vettmax <- vett * input$ga_max_factor

    # Fitness: maximise Cor(Optimal, Proportional)
    fitness_cvs <- function(solution) {
      m       <- matrix(solution, nrow = n_dom, ncol = n_var, byrow = TRUE)
      cv_corr <- as.data.frame(m, stringsAsFactors = FALSE)
      colnames(cv_corr) <- paste0("CV", seq_len(n_var))
      cv_corr <- cbind(DOM = paste0("DOM", seq_len(n_dom)), cv_corr,
                       stringsAsFactors = FALSE)
      tryCatch({
        a     <- R2BEAT::beat.1st(strata_local, cv_corr)
        alloc <- a$alloc[-nrow(a$alloc), ]
        r     <- cor(alloc$ALLOC, alloc$PROP)
        if (!is.finite(r))      return(1)
        if (abs(r) > limit + 0.005) return(0)
        return(-r)
      }, error = function(e) 1)
    }

    tryCatch({
      withProgress(message = "Running Genetic Algorithm…", value = 0.5, {
        ga_obj <- rbga(
          stringMin = vettmin,
          stringMax = vettmax,
          popSize   = input$ga_popsize,
          iters     = input$ga_iters,
          evalFunc  = fitness_cvs,
          verbose   = FALSE
        )
      })
      rv$ga_obj     <- ga_obj
      rv$ga_vett    <- vett
      rv$ga_vettmin <- vettmin
      rv$ga_vettmax <- vettmax

      # Best solution
      evals <- ga_obj$evaluations
      filt  <- evals == min(evals)
      best  <- if (sum(filt) > 1) ga_obj$population[filt, ][1, ] else
                                   ga_obj$population[filt, ]

      m      <- matrix(best, nrow = n_dom, ncol = n_var, byrow = TRUE)
      cv_ga  <- as.data.frame(m, stringsAsFactors = FALSE)
      colnames(cv_ga) <- paste0("CV", seq_len(n_var))
      cv_ga  <- cbind(DOM = paste0("DOM", seq_len(n_dom)), cv_ga,
                      stringsAsFactors = FALSE)

      rv$cv_ga     <- cv_ga
      rv$ga_result <- R2BEAT::beat.1st(rv$strata, cv_ga)

      append_script(paste0(
        "\n# 4. Genetic Algorithm adjustment\n",
        "# Pareto CV vector:\n",
        "vett    <- c(", paste(round(vett, 8), collapse = ", "), ")\n",
        "vettmin <- vett * ", input$ga_min_factor, "\n",
        "vettmax <- vett * ", input$ga_max_factor, "\n",
        "\n",
        "fitness_cvs <- function(solution) {\n",
        "  m       <- matrix(solution, nrow = ", n_dom, ", ncol = ", n_var, ", byrow = TRUE)\n",
        "  cv_corr <- as.data.frame(m)\n",
        "  colnames(cv_corr) <- paste0('CV', 1:", n_var, ")\n",
        "  cv_corr <- cbind(DOM = paste0('DOM', 1:", n_dom, "), cv_corr)\n",
        "  tryCatch({\n",
        "    a <- beat.1st(strata, cv_corr)\n",
        "    al <- a$alloc[-nrow(a$alloc), ]\n",
        "    r <- cor(al$ALLOC, al$PROP)\n",
        "    if (!is.finite(r)) return(1)\n",
        "    if (abs(r) > ", limit, " + 0.005) return(0)\n",
        "    return(-r)\n",
        "  }, error = function(e) 1)\n",
        "}\n",
        "\n",
        "solution_genalg <- rbga(\n",
        "  stringMin = vettmin, stringMax = vettmax,\n",
        "  popSize   = ", input$ga_popsize, ", iters = ", input$ga_iters, ",\n",
        "  evalFunc  = fitness_cvs, verbose = FALSE\n",
        ")\n",
        "\n",
        "# Extract best solution\n",
        "filt <- solution_genalg$evaluations == min(solution_genalg$evaluations)\n",
        "best <- if (sum(filt) > 1) solution_genalg$population[filt,][1,] else\n",
        "        solution_genalg$population[filt,]\n",
        "m    <- matrix(best, nrow = ", n_dom, ", ncol = ", n_var, ", byrow = TRUE)\n",
        "cv_ga <- as.data.frame(m)\n",
        "colnames(cv_ga) <- paste0('CV', 1:", n_var, ")\n",
        "cv_ga <- cbind(DOM = paste0('DOM', 1:", n_dom, "), cv_ga)\n",
        "ga_alloc <- beat.1st(strata, cv_ga)\n",
        "cat('GA total:', sum(ga_alloc$alloc$ALLOC[-nrow(ga_alloc$alloc)]), '\\n')\n",
        "cor(ga_alloc$alloc$ALLOC[-nrow(ga_alloc$alloc)],\n",
        "    ga_alloc$alloc$PROP[-nrow(ga_alloc$alloc)])\n"
      ))
      showNotification("GA complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("GA error:", e$message), type = "error", duration = 8)
    })
  })

  output$tbl_ga_cvs <- renderDT({
    req(rv$cv_ga)
    datatable(round_cv_df(rv$cv_ga), options = list(dom = "t", scrollX = TRUE))
  })

  output$txt_ga_total <- renderText({
    req(rv$ga_result)
    alloc <- rv$ga_result$alloc
    format(round(sum(alloc$ALLOC[-nrow(alloc)])), big.mark = ",")
  })

  output$txt_ga_cor <- renderText({
    req(rv$ga_result)
    alloc <- rv$ga_result$alloc[-nrow(rv$ga_result$alloc), ]
    sprintf("%.4f", cor(alloc$ALLOC, alloc$PROP))
  })

  # "Best and mean evaluation value" — mirrors the article's plot(solution_genalg)
  output$plot_ga_fitness <- renderPlot({
    req(rv$ga_obj)
    best_vals <- rv$ga_obj$best   # per-generation best  (raw = negative correlation)
    mean_vals <- rv$ga_obj$mean   # per-generation mean
    req(length(best_vals) > 0)
    gen <- seq_along(best_vals)
    df <- rbind(
      data.frame(generation = gen, value = best_vals, Series = "best"),
      data.frame(generation = gen, value = mean_vals, Series = "mean")
    )
    y_lim <- c(min(df$value, na.rm = TRUE), max(df$value, na.rm = TRUE))
    ggplot(df, aes(x = generation, y = value, colour = Series)) +
      geom_line(linewidth = 0.8) +
      scale_colour_manual(values = c(best = "black", mean = "blue")) +
      geom_hline(yintercept = -isolate(input$ga_limit),
                 linetype = "dashed", colour = "red", linewidth = 0.6) +
      labs(title = "Best and mean evaluation value",
           x = "generation", y = "evaluation value", colour = NULL) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "top")
  })

  # "Region of variation for CVs" — shows cv_input / cv_min / cv_max per DOM_V
  output$plot_ga_region <- renderPlot({
    req(rv$ga_vett, rv$ga_vettmin, rv$ga_vettmax, rv$n_dom, rv$n_var)
    n_dom  <- rv$n_dom;  n_var <- rv$n_var
    idx    <- seq_len(n_dom * n_var)
    labels <- paste0("DOM", rep(seq_len(n_dom), each = n_var),
                     "_V", rep(seq_len(n_var), times = n_dom))
    df <- rbind(
      data.frame(x = idx, y = rv$ga_vett,    Series = "cv_input"),
      data.frame(x = idx, y = rv$ga_vettmin, Series = "cv_min"),
      data.frame(x = idx, y = rv$ga_vettmax, Series = "cv_max")
    )
    df$Series <- factor(df$Series, levels = c("cv_input", "cv_min", "cv_max"))
    ggplot(df, aes(x = x, y = y, colour = Series)) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.5) +
      scale_x_continuous(breaks = idx, labels = labels) +
      scale_colour_manual(values = c(cv_input = "black",
                                     cv_min   = "red",
                                     cv_max   = "blue")) +
      labs(title = "Region of variation for CVs",
           x = NULL, y = "CV", colour = NULL) +
      theme_minimal(base_size = 11) +
      theme(axis.text.x  = element_text(angle = 90, vjust = 0.5,
                                        hjust = 1, size = 8),
            legend.position = "top")
  })

  # ── STEP 5: FINAL SOLUTION ───────────────────────────────────────────────────

  observeEvent(input$btn_run_final, {
    req(rv$cv_ga)
    tryCatch({
      withProgress(message = "Running final Pareto…", value = 0.5, {
        out2 <- run_pareto_impl(rv$cv_ga)
      })
      # errors_scaled already has CV1, CV2, … column names
      cv_f <- out2$errors_scaled
      rv$cv_final     <- cv_f
      rv$final_result <- R2BEAT::beat.1st(rv$strata, cv_f)

      append_script(paste0(
        "\n# 5. Final Pareto (on GA-adjusted CVs)\n",
        "out_final <- pareto(\n",
        "  strata        = strata,\n",
        "  current_cvs   = cv_ga,\n",
        "  target_size   = ", input$pareto_target, ",\n",
        "  tolerance     = ", input$pareto_tol, ",\n",
        "  beat1cv_fun   = beat.1cv_2,\n",
        "  max_iter      = ", input$pareto_maxiter, ",\n",
        "  max_same_iter = ", input$pareto_maxsame, ",\n",
        "  plot          = FALSE\n",
        ")\n",
        "cv_final     <- out_final$errors_scaled\n",
        "final_alloc  <- beat.1st(strata, cv_final)\n",
        "cat('Final total:', sum(final_alloc$alloc$ALLOC[-nrow(final_alloc$alloc)]), '\\n')\n",
        "cor(final_alloc$alloc$ALLOC[-nrow(final_alloc$alloc)],\n",
        "    final_alloc$alloc$PROP[-nrow(final_alloc$alloc)])\n"
      ))
      showNotification("Final Pareto complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Final Pareto error:", e$message), type = "error", duration = 8)
    })
  })

  # Helpers: fall back to Pareto solution when GA / Final Pareto were skipped
  final_cv_display     <- reactive({ rv$cv_final     %||% rv$cv_pareto     })
  final_result_display <- reactive({ rv$final_result %||% rv$pareto_result })

  output$ui_final_controls <- renderUI({
    if (is.null(rv$cv_ga)) {
      tagList(
        div(class = "alert alert-info", style = "padding:8px;",
          icon("info-circle"),
          " GA step was skipped — the",
          tags$b(" Pareto solution"), " is used as the final result below."
        ),
        p(em("Run Step 4 (GA Adjustment) if you want to improve",
             " the proportionality of the allocation."))
      )
    } else if (is.null(rv$cv_final)) {
      tagList(
        p("GA adjustment complete. Apply Pareto to rescale the GA CVs",
          "back to the target sample size."),
        actionButton("btn_run_final", "Run Final Pareto",
                     class = "btn-danger btn-block", icon = icon("flag-checkered"))
      )
    } else {
      tagList(
        div(class = "alert alert-success", style = "padding:8px;",
          icon("check-circle"), tags$b(" Final solution ready.")
        ),
        actionButton("btn_run_final", "Re-run Final Pareto",
                     class = "btn-default btn-block", icon = icon("flag-checkered"))
      )
    }
  })

  output$ui_final_alloc_title <- renderUI({
    if (is.null(rv$cv_ga)) {
      tags$div(class = "box-header with-border",
        tags$h3(class = "box-title",
                "Pareto solution: Optimal vs Proportional (GA skipped)"))
    } else {
      tags$div(class = "box-header with-border",
        tags$h3(class = "box-title",
                "Final solution: Optimal vs Proportional (after GA)"))
    }
  })

  output$tbl_comparison <- renderTable({
    req(rv$init_result)
    solutions <- list(Initial = rv$init_result,
                      Pareto  = rv$pareto_result)
    if (!is.null(rv$ga_result))    solutions[["GA"]]    <- rv$ga_result
    if (!is.null(rv$final_result)) solutions[["Final"]] <- rv$final_result

    rows <- lapply(names(solutions), function(nm) {
      res    <- solutions[[nm]]
      alloc  <- res$alloc
      allsub <- alloc[-nrow(alloc), ]
      n <- sum(alloc$ALLOC[-nrow(alloc)])
      r <- cor(allsub$ALLOC, allsub$PROP)
      data.frame(Solution = nm,
                 `Sample size`    = format(round(n), big.mark = ","),
                 `Cor(Opt, Prop)` = sprintf("%.4f", r),
                 check.names = FALSE, stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
  }, striped = TRUE, hover = TRUE, spacing = "s")

  output$tbl_final_cvs <- renderDT({
    cv <- final_cv_display()
    req(cv)
    datatable(round_cv_df(cv), options = list(dom = "t", scrollX = TRUE))
  })

  output$btn_download_final_cv <- downloadHandler(
    filename = function() {
      if (is.null(rv$cv_final)) "pareto_cvs.csv" else "final_cvs.csv"
    },
    content = function(file) {
      cv <- final_cv_display()
      req(cv)
      write.csv(cv, file, row.names = FALSE)
    }
  )

  output$plot_alloc_before <- renderPlot({
    req(rv$pareto_result)
    alloc_scatter(rv$pareto_result, "Pareto solution")
  })

  output$plot_alloc_after <- renderPlot({
    res <- final_result_display()
    req(res)
    title <- if (is.null(rv$cv_ga)) "Pareto solution (GA skipped)" else "Final solution (after GA)"
    alloc_scatter(res, title)
  })

  output$ui_domain_selector <- renderUI({
    req(rv$cv_pareto)
    selectInput("sel_domain", "Domain:",
                choices = rv$cv_pareto$DOM, selected = rv$cv_pareto$DOM[1])
  })

  output$plot_cv_compare <- renderPlot({
    req(rv$cv_pareto, input$sel_domain)
    d    <- input$sel_domain
    n_var <- rv$n_var
    vn    <- rv$var_names %||% paste0("V", seq_len(n_var))
    vcols <- paste0("CV", seq_len(n_var))

    cvsets <- list(Pareto = rv$cv_pareto)
    if (!is.null(rv$cv_ga))    cvsets[["GA"]]    <- rv$cv_ga
    if (!is.null(rv$cv_final)) cvsets[["Final"]] <- rv$cv_final

    rows <- lapply(names(cvsets), function(nm) {
      cv_df <- cvsets[[nm]]
      row   <- cv_df[cv_df$DOM == d, vcols, drop = FALSE]
      if (nrow(row) == 0) return(NULL)
      data.frame(Solution = nm, Variable = vn,
                 CV = as.numeric(row[1, ]),
                 stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)

    df_plot          <- do.call(rbind, rows)
    df_plot$Solution <- factor(df_plot$Solution, levels = names(cvsets))
    df_plot$Variable <- factor(df_plot$Variable, levels = vn)

    ggplot(df_plot, aes(x = Variable, y = CV, fill = Solution)) +
      geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = c(Pareto = "#e41a1c", GA = "#377eb8",
                                   Final  = "#4daf4a")) +
      labs(title = paste("CV by variable — Domain:", d),
           x = NULL, y = "Coefficient of Variation") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top")
  })

  # ── STEP 6: R SCRIPT ─────────────────────────────────────────────────────────

  output$txt_script_preview <- renderText({
    if (length(rv$script_lines) == 0)
      return("# No steps recorded yet. Complete the workflow to generate a script.")
    paste(rv$script_lines, collapse = "\n")
  })

  output$btn_download_script <- downloadHandler(
    filename = function() "r2beat_cv_optimization.R",
    content  = function(file) writeLines(rv$script_lines, file)
  )

} # end server


# ============================ WORKFLOW: WF1 ============================



# ─── HELPERS ──────────────────────────────────────────────────────────────────


wf1_weight_boxplot <- function(df, grp_var, title_txt) {
  if (!is.null(grp_var) && nchar(grp_var) > 0 && grp_var %in% names(df)) {
    df[[grp_var]] <- factor(df[[grp_var]])
    ggplot(df, aes(x = .data[[grp_var]], y = weight)) +
      geom_boxplot(fill = "orange", colour = "grey30", outlier.size = 0.6) +
      labs(title = title_txt, x = NULL, y = "Weight") +
      coord_flip() +
      theme_minimal(base_size = 11)
  } else {
    ggplot(df, aes(x = "", y = weight)) +
      geom_boxplot(fill = "orange", colour = "grey30") +
      labs(title = title_txt, x = NULL, y = "Weight") +
      theme_minimal(base_size = 11)
  }
}

# ─── UI ────────────────────────────────────────────────────────────────────────

wf1_ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "R2BEAT - Workflow 1"),

  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("1. Frame Data",    tabName = "tab_data",   icon = icon("database")),
      menuItem("2. Parameters",    tabName = "tab_params", icon = icon("sliders-h")),
      menuItem("3. Sensitivity",   tabName = "tab_sens",   icon = icon("chart-line")),
      menuItem("4. Prepare Input", tabName = "tab_prep",   icon = icon("cogs")),
      menuItem("5. Allocation",    tabName = "tab_alloc",  icon = icon("bullseye")),
      menuItem("6. Select Sample", tabName = "tab_samp",   icon = icon("users")),
      menuItem("7. R Script",      tabName = "tab_script", icon = icon("file-code"))
    ),
    div(style = "padding:12px 15px; border-top:1px solid #444; margin-top:10px;",
      actionButton("btn_home", " Menu principale", icon = icon("home"),
                   class = "btn-default btn-block",
                   style = "font-weight:bold; margin-bottom:8px;",
                   onclick = "window.location.search='';"),
      actionButton("btn_help", " Help Guide", icon = icon("question-circle"),
                   class = "btn-info btn-block", style = "font-weight:bold;")
    )
  ),

  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .box-body .dataTables_wrapper { overflow-x: auto; }
      pre { max-height: 380px; overflow-y: auto; }
      .cv-grid input[type=number] { width: 72px; }
      .cv-grid td { vertical-align: middle !important; padding: 3px 5px !important; }
      .cv-grid th { text-align: center; padding: 4px 5px; }
    "))),
    tabItems(

      # ── TAB 1: FRAME DATA ─────────────────────────────────────────────────
      tabItem(tabName = "tab_data",
        fluidRow(
          box(title = "Data Source", width = 4, status = "primary", solidHeader = TRUE,
            radioButtons("data_source", "Choose population frame:",
              choices  = c("Use bundled pop.RData (demo)" = "demo",
                           "Upload .RData file"           = "upload",
                           "Upload delimited file (CSV/TSV)" = "delimited"),
              selected = "demo"),
            conditionalPanel("input.data_source == 'upload'",
              fileInput("upload_rdata", "Upload .RData file",
                        accept = c(".RData", ".rdata", ".rda"))
            ),
            conditionalPanel("input.data_source == 'delimited'",
              fileInput("upload_csv", "Upload CSV/TSV file",
                        accept = c(".csv", ".tsv", ".txt")),
              fluidRow(
                column(6,
                  radioButtons("csv_sep", "Separator:",
                    choices  = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                    selected = ",", inline = FALSE)
                ),
                column(6,
                  checkboxInput("csv_header",    "Header row",        value = TRUE),
                  checkboxInput("csv_dec_comma", "Decimal comma ( , )", value = FALSE)
                )
              )
            ),
            div(class = "alert alert-warning", style = "padding:8px; font-size:12px;",
              icon("exclamation-triangle"),
              " The demo dataset has 2.2 M rows — loading may take ~10 s."
            ),
            actionButton("btn_load", "Load Frame", class = "btn-primary btn-block",
                         icon = icon("upload"))
          ),
          box(title = "Frame Summary", width = 4, status = "info",
            verbatimTextOutput("txt_frame_summary")
          ),
          box(title = "Structure", width = 4, status = "info",
            verbatimTextOutput("txt_frame_str")
          )
        ),
        fluidRow(
          box(title = "Frame Preview (first 100 rows)", width = 12, status = "info",
            DTOutput("tbl_frame_head")
          )
        )
      ),

      # ── TAB 2: PARAMETERS ────────────────────────────────────────────────
      tabItem(tabName = "tab_params",
        fluidRow(
          box(title = "Variable Assignment & Numeric Parameters",
              width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("ui_id_psu"),
            uiOutput("ui_id_ssu"),
            uiOutput("ui_strata_var"),
            uiOutput("ui_target_vars"),
            uiOutput("ui_deff_var"),
            uiOutput("ui_domain_vars"),
            hr(),
            numericInput("minimum",     "Minimum SSU per PSU:",   25, 2, 500),
            numericInput("delta",       "Delta (intraclass corr.):", 1, 0, 10),
            numericInput("deff_sugg",   "Suggested DEFF:",       1.5, 1, 20),
            hr(),
            actionButton("btn_set_params", "Confirm Parameters",
                         class = "btn-success btn-block", icon = icon("check"))
          ),
          box(title = "CV Constraints", width = 8, status = "primary", solidHeader = TRUE,
            p("One row per domain level: ", tags$b("DOM1"), " = national,",
              tags$b(" DOM2"), " = Level-1 domain variable,",
              tags$b(" DOM3"), " = Level-2, …",
              " One column per selected target variable.",
              " Set ", tags$b("Domain levels = N\u2009+\u20091"),
              " where N is the number of domain variables selected on the left."),
            fluidRow(
              column(3,
                numericInput("n_dom", "Domain levels:", 2, 1, 6)
              ),
              column(4,
                numericInput("cv_equal_val", "Set all equal to:", 0.05,
                             0.001, 1, 0.005)
              ),
              column(3, br(),
                actionButton("btn_set_equal", "Apply",
                             class = "btn-default btn-block")
              )
            ),
            uiOutput("ui_ndom_warning"),
            fluidRow(
              column(5,
                fileInput("upload_cv_file", "Load CV values from CSV:",
                          accept = c(".csv", ".txt"),
                          buttonLabel = "Browse…",
                          placeholder = "DOM, CV1, CV2, … format")
              ),
              column(7,
                div(style = "margin-top:28px;",
                  div(class = "alert alert-info", style = "padding:6px 10px; font-size:12px; margin-bottom:0;",
                    icon("info-circle"),
                    " CSV must contain a ", tags$b("DOM"), " column and one or more ",
                    tags$b("CV1, CV2, …"), " columns. The number of domain levels",
                    " is updated automatically."
                  )
                )
              )
            ),
            hr(),
            div(class = "cv-grid", uiOutput("ui_cv_table")),
            hr(),
            h5("Current CV dataframe:"),
            verbatimTextOutput("txt_cv_preview")
          )
        )
      ),

      # ── TAB 3: SENSITIVITY ANALYSIS ──────────────────────────────────────
      tabItem(tabName = "tab_sens",
        fluidRow(
          box(title = tags$span("Sensitivity Analysis",
                tags$small(tags$span(class = "label label-success",
                                     style = "margin-left:8px;", "optional"))),
              width = 4, status = "primary", solidHeader = TRUE,
            div(class = "alert alert-info", style = "padding:8px;",
              icon("info-circle"),
              " This step is ", tags$b("optional."),
              " It helps choose a good minimum SSU per PSU value."
            ),
            numericInput("sens_min", "Min SSU to test:", 30,  2, 500),
            numericInput("sens_max", "Max SSU to test:", 80, 10, 1000),
            hr(),
            div(class = "alert alert-warning", style = "padding:8px; font-size:12px;",
              icon("clock"),
              " May take several minutes for large frames."
            ),
            actionButton("btn_run_sens", "Run Sensitivity",
                         class = "btn-warning btn-block", icon = icon("chart-line")),
            br(),
            p(em("Requires Step 2 to be confirmed first."))
          ),
          box(title = "Sensitivity Results", width = 8, status = "info",
            h4("PSU count and total SSU vs minimum SSU per PSU"),
            plotOutput("plot_sens", height = "320px"),
            hr(),
            h4("Sensitivity table"),
            DTOutput("tbl_sens")
          )
        )
      ),

      # ── TAB 4: PREPARE INPUT ─────────────────────────────────────────────
      tabItem(tabName = "tab_prep",
        fluidRow(
          box(title = "Prepare Allocation Input", width = 3,
              status = "primary", solidHeader = TRUE,
            p("Aggregates the frame into strata-level statistics and computes",
              "design-effect parameters needed by ", tags$code("beat.2st()"), "."),
            div(class = "alert alert-warning", style = "padding:8px; font-size:12px;",
              icon("clock"),
              " May take 1–2 minutes for large frames."
            ),
            actionButton("btn_run_prep", "Run Preparation",
                         class = "btn-danger btn-block", icon = icon("cogs")),
            br(),
            p(em("Requires Step 2 to be confirmed first."))
          ),
          box(title = "Execution Log", width = 9, status = "warning",
            verbatimTextOutput("txt_prep_log")
          )
        ),
        fluidRow(
          box(title = "Strata (head)", width = 6, status = "info",
            DTOutput("tbl_inp_strata")
          ),
          box(title = "Design File (des_file)", width = 6, status = "info",
            DTOutput("tbl_inp_des")
          )
        ),
        fluidRow(
          box(title = "PSU File (head)", width = 6, status = "info",
            DTOutput("tbl_inp_psu")
          ),
          box(title = "Design Parameters", width = 6, status = "info",
            h5("DEFF"),
            verbatimTextOutput("txt_inp_deff"),
            h5("Effective stratum size (effst — head)"),
            verbatimTextOutput("txt_inp_effst"),
            h5("Intraclass correlation (rho — head)"),
            verbatimTextOutput("txt_inp_rho")
          )
        )
      ),

      # ── TAB 5: ALLOCATION ────────────────────────────────────────────────
      tabItem(tabName = "tab_alloc",
        fluidRow(
          box(title = "Allocation Settings", width = 3,
              status = "primary", solidHeader = TRUE,
            numericInput("alloc_minPSU",      "Min PSU per stratum (minPSUstrat):", 2, 1, 50),
            numericInput("alloc_minnumstrat", "Min strata (minnumstrat):",          2, 1, 20),
            br(),
            actionButton("btn_run_alloc", "Run Allocation (beat.2st)",
                         class = "btn-danger btn-block", icon = icon("bullseye")),
            br(),
            p(em("Requires Step 4 to be completed first."))
          ),
          box(title = "Allocation Results", width = 9, status = "info",
            fluidRow(
              column(4,
                h4("Total PSUs (1st stage)"),
                verbatimTextOutput("txt_alloc_total")
              ),
              column(4,
                h4("Total SSUs (2nd stage)"),
                verbatimTextOutput("txt_alloc_ssu")
              )
            ),
            hr(),
            h4("Allocation by stratum"),
            DTOutput("tbl_alloc")
          )
        ),
        fluidRow(
          box(title = "Expected CVs", width = 6, status = "info",
            DTOutput("tbl_expected")
          ),
          box(title = "Sensitivity to CV constraints", width = 6, status = "info",
            p(em("Values > 1 indicate that a constraint is active (binding).")),
            DTOutput("tbl_sensitivity")
          )
        )
      ),

      # ── TAB 6: SELECT SAMPLE ─────────────────────────────────────────────
      tabItem(tabName = "tab_samp",
        fluidRow(
          box(title = "Sampling Settings", width = 3,
              status = "primary", solidHeader = TRUE,
            numericInput("seed", "Random seed:", 1234, 1, 99999),
            hr(),
            h4("Step 1 — Select PSUs"),
            selectInput("psu_type", "Allocation type:",
              choices  = c("ALLOC", "PROP", "EQUAL"), selected = "ALLOC"),
            checkboxInput("pps", "PPS sampling", value = TRUE),
            actionButton("btn_select_psu", "Select PSUs",
                         class = "btn-warning btn-block", icon = icon("map-marker")),
            hr(),
            h4("Step 2 — Select SSUs"),
            uiOutput("ui_ssu_code"),
            actionButton("btn_select_ssu", "Select SSUs",
                         class = "btn-success btn-block", icon = icon("users")),
            hr(),
            downloadButton("btn_download_samp", "Download Sample (CSV)",
                           class = "btn-block"),
            br(),
            p(em("Requires Step 5 to be completed first."))
          ),
          box(title = "Sample Results", width = 9, status = "info",
            h4("PSU statistics"),
            verbatimTextOutput("txt_psu_stats"),
            hr(),
            fluidRow(
              column(6,
                h4("Sample size (SSU)"),
                verbatimTextOutput("txt_samp_size")
              ),
              column(6,
                h4("Reconstructed population total (sum of weights)"),
                verbatimTextOutput("txt_pop_n")
              )
            ),
            hr(),
            h4("Sample preview"),
            DTOutput("tbl_samp_head")
          )
        ),
        fluidRow(
          box(title = "Weight distribution — total sample", width = 6,
            plotOutput("plot_wt_total",    height = "260px")
          ),
          box(title = uiOutput("ui_wt_domain_title"), width = 6,
            plotOutput("plot_wt_domain",   height = "260px")
          )
        ),
        fluidRow(
          box(title = "Weight distribution — by province", width = 6,
            plotOutput("plot_wt_province", height = "260px")
          ),
          box(title = uiOutput("ui_wt_stratum_title"), width = 6,
            plotOutput("plot_wt_stratum",  height = "260px")
          )
        )
      ),

      # ── TAB 7: R SCRIPT ──────────────────────────────────────────────────
      tabItem(tabName = "tab_script",
        fluidRow(
          box(title = "Generated R Script", width = 12, status = "info",
            p("Reproducible script recording all steps executed in this session."),
            verbatimTextOutput("txt_script")
          )
        ),
        fluidRow(
          box(title = "Download", width = 4, status = "primary", solidHeader = TRUE,
            downloadButton("btn_dl_script", "Download R Script", class = "btn-block")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

# ─── SERVER ────────────────────────────────────────────────────────────────────

wf1_server <- function(input, output, session) {

  rv <- reactiveValues(
    frame        = NULL,
    params_set   = FALSE,
    cv           = NULL,
    domain_vars  = NULL,   # confirmed domain variable(s), set by btn_set_params
    target_vars  = NULL,   # confirmed target variables (snapshot at btn_set_params)
    id_psu       = NULL,   # confirmed PSU column
    id_ssu       = NULL,   # confirmed SSU column
    strata_var   = NULL,   # confirmed stratum column
    deff_var     = NULL,   # confirmed DEFF variable
    minimum      = NULL,   # confirmed minimum SSU per PSU
    delta        = NULL,   # confirmed delta
    deff_sugg    = NULL,   # confirmed suggested DEFF
    pending_cv   = NULL,
    sens_result  = NULL,
    inp          = NULL,
    prep_log     = NULL,
    alloc        = NULL,
    sample_1st   = NULL,
    samp         = NULL,
    script_lines = character(0)
  )

  append_script <- function(...) {
    rv$script_lines <- c(rv$script_lines, paste0(...))
  }

  # ── HELP MODAL ──────────────────────────────────────────────────────────────
  observeEvent(input$btn_help, {
    showModal(modalDialog(
      title = tags$span(icon("question-circle"), " R2BEAT Workflow 1 — Guide"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      tabsetPanel(
        tabPanel(tags$span(icon("database"), " 1. Frame"),
          br(),
          tags$h4("Purpose"),
          tags$p("Load a row-level population frame (one row per individual / SSU) that",
                 " contains the PSU identifier, a stratum column, estimation domain,",
                 " and the target survey variables."),
          tags$h4("Demo dataset"),
          tags$p(tags$code("pop.RData"), " holds a synthetic population of 2.2 M individuals",
                 " with 13 columns: region (domain), province, municipality (PSU),",
                 " id_ind (SSU), stratum, and four target variables:"),
          tags$code("income_hh, active, inactive, unemployed"),
          tags$h4("Upload options"),
          tags$p("Any ", tags$code(".RData"), " file (first object is used) or a delimited",
                 " text file (CSV, TSV, etc.).")
        ),
        tabPanel(tags$span(icon("sliders-h"), " 2. Parameters"),
          br(),
          tags$h4("Variable assignment"),
          tags$ul(
            tags$li(tags$b("PSU identifier:"), " Column identifying Primary Sampling Units (e.g., municipality)."),
            tags$li(tags$b("SSU identifier:"), " Column for Secondary Sampling Units (e.g., id_ind)."),
            tags$li(tags$b("Stratum column:"), " Pre-defined stratification variable."),
            tags$li(tags$b("Target variables:"), " Columns to be estimated (drive CV constraints)."),
            tags$li(tags$b("DEFF variable:"), " Variable used to compute the design effect."),
            tags$li(tags$b("Domain variable:"), " Column defining estimation domains (e.g., region).")
          ),
          tags$h4("Numeric parameters"),
          tags$ul(
            tags$li(tags$b("Minimum SSU per PSU:"), " Minimum within-PSU sample size (also used in sensitivity analysis)."),
            tags$li(tags$b("Delta:"), " Intra-cluster correlation coefficient for SSU within PSU."),
            tags$li(tags$b("Suggested DEFF:"), " Design effect estimate used in initial allocation.")
          ),
          tags$h4("CV constraints"),
          tags$p("One row per domain level:",
                 " DOM1 = national, DOM2 = Level-1 domain variable, DOM3 = Level-2, …",
                 " One column per target variable. Values are the maximum acceptable",
                 " Coefficient of Variation (e.g., 0.05 = 5 %)."),
          tags$p(tags$b("Rule: Domain levels = N\u2009+\u20091"),
                 " where N is the number of domain variables selected.",
                 " The extended ", tags$code("prepareInputToAllocation1_v2()"),
                 " creates DOM1 (national constant) + one DOM column per domain variable,",
                 " all of which are used by ", tags$code("beat.2st()"), ".")
        ),
        tabPanel(tags$span(icon("chart-line"), " 3. Sensitivity"),
          br(),
          tags$h4("Purpose (optional)"),
          tags$p("Calls ", tags$code("sensitivity_min_SSU()"), " to explore how the",
                 " number of PSUs and total SSUs change as the minimum-SSU-per-PSU",
                 " threshold varies. This helps choose a practical minimum that balances",
                 " field costs (PSU count) and precision (total SSU)."),
          tags$h4("Plot"),
          tags$p("The ", tags$code("plot_sens()"), " output shows PSU count and SSU total",
                 " over the tested range of minimum values. Look for the 'elbow' where",
                 " adding more SSUs per PSU yields little gain in total sample size.")
        ),
        tabPanel(tags$span(icon("cogs"), " 4. Prepare Input"),
          br(),
          tags$h4("Purpose"),
          tags$p("Calls ", tags$code("prepareInputToAllocation1()"), " which aggregates the",
                 " row-level frame into the strata-level summary objects required by",
                 tags$code("beat.2st()"), "."),
          tags$h4("Output objects"),
          tags$ul(
            tags$li(tags$b("strata:"), " Mean (M1…Mv), SD (S1…Sv), population size (N) per stratum."),
            tags$li(tags$b("des_file:"), " Design file — stratum MOS, delta, minimum SSU."),
            tags$li(tags$b("psu_file:"), " PSU-level measure of size."),
            tags$li(tags$b("deff / effst / rho:"), " Design-effect components passed to beat.2st.")
          )
        ),
        tabPanel(tags$span(icon("bullseye"), " 5. Allocation"),
          br(),
          tags$h4("Purpose"),
          tags$p("Calls ", tags$code("beat.2st()"), " — the two-stage Bethel optimal",
                 " allocation that minimises total PSUs subject to all CV constraints."),
          tags$h4("Key parameters"),
          tags$ul(
            tags$li(tags$b("minPSUstrat:"), " Minimum PSUs per stratum."),
            tags$li(tags$b("minnumstrat:"), " Minimum number of strata in each domain.")
          ),
          tags$h4("Results"),
          tags$ul(
            tags$li(tags$b("Allocation table:"), " ALLOC (optimal), PROP, EQUAL per stratum."),
            tags$li(tags$b("Expected CVs:"), " Precision achieved per domain / variable."),
            tags$li(tags$b("Sensitivity:"), " 1 = constraint is active (binding).")
          )
        ),
        tabPanel(tags$span(icon("users"), " 6. Sample"),
          br(),
          tags$h4("Step 1 — Select PSUs"),
          tags$p(tags$code("select_PSU()"), " draws primary sampling units according to",
                 " the chosen allocation type (ALLOC, PROP, EQUAL) with or without PPS."),
          tags$h4("Step 2 — Select SSUs"),
          tags$p(tags$code("select_SSU()"), " draws secondary units within the selected",
                 " PSUs and attaches first-stage, second-stage, and total inclusion",
                 " probabilities plus the sampling weight."),
          tags$h4("Weight distribution plots"),
          tags$p("Four boxplots check weight variability: overall sample, by domain,",
                 " by province, and by stratum. Large spread may signal stratification",
                 " or allocation issues.")
        )
      )
    ))
  })

  # ── STEP 1: LOAD FRAME ──────────────────────────────────────────────────────
  observeEvent(input$btn_load, {
    tryCatch({
      if (input$data_source == "demo") {
        withProgress(message = "Loading pop.RData…", value = 0.5, {
          env <- new.env()
          load("pop.RData", envir = env)
          if (!exists("pop", envir = env))
            stop("pop.RData must contain an object named 'pop'.")
          rv$frame <- env$pop
        })
        data_expr <- 'load("pop.RData")\nsamp_frame <- pop'

      } else if (input$data_source == "upload") {
        req(input$upload_rdata)
        env <- new.env()
        load(input$upload_rdata$datapath, envir = env)
        objs <- ls(envir = env)
        if (length(objs) == 0) stop("No objects found in the RData file.")
        rv$frame  <- env[[objs[1]]]
        data_expr <- paste0('load("', input$upload_rdata$name, '")\n',
                            'samp_frame <- ', objs[1])

      } else {
        req(input$upload_csv)
        dec <- if (isTRUE(input$csv_dec_comma)) "," else "."
        rv$frame <- read.table(
          input$upload_csv$datapath,
          header           = isTRUE(input$csv_header),
          sep              = input$csv_sep,
          dec              = dec,
          stringsAsFactors = FALSE
        )
        data_expr <- paste0(
          'samp_frame <- read.table("', input$upload_csv$name, '",\n',
          '  header = ', isTRUE(input$csv_header),
          ', sep = "', input$csv_sep,
          '", dec = "', dec, '", stringsAsFactors = FALSE)'
        )
      }

      rv$script_lines <- c(
        "# R2BEAT workflow 1 — generated by Shiny app",
        "library(R2BEAT)",
        "",
        "# 1. Load frame",
        data_expr,
        ""
      )
      showNotification("Frame loaded successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_frame_summary <- renderPrint({
    req(rv$frame)
    df <- rv$frame
    cat("Rows :", format(nrow(df), big.mark = ","), "\n")
    cat("Cols :", ncol(df), "\n\n")
    for (v in names(df)) {
      x <- df[[v]]
      if (is.numeric(x))
        cat(sprintf("  %-22s numeric  [%s, %s]\n", v,
                    format(min(x, na.rm = TRUE), big.mark = ","),
                    format(max(x, na.rm = TRUE), big.mark = ",")))
      else
        cat(sprintf("  %-22s factor   (%d levels)\n", v, length(unique(x))))
    }
  })

  output$txt_frame_str <- renderPrint({
    req(rv$frame)
    str(rv$frame, give.attr = FALSE)
  })

  output$tbl_frame_head <- renderDT({
    req(rv$frame)
    datatable(head(rv$frame, 100), options = list(scrollX = TRUE, pageLength = 8))
  })

  # ── STEP 2: PARAMETERS ──────────────────────────────────────────────────────
  col_choices <- reactive({ req(rv$frame); names(rv$frame) })

  is_demo <- reactive({ identical(input$data_source, "demo") })
  output$ui_id_psu      <- renderUI({ selectInput("id_psu",     "PSU identifier column:", col_choices(),
                                                   selected = if (is_demo() && "municipality" %in% col_choices()) "municipality") })
  output$ui_id_ssu      <- renderUI({ selectInput("id_ssu",     "SSU identifier column:", col_choices(),
                                                   selected = if (is_demo() && "id_ind" %in% col_choices()) "id_ind") })
  output$ui_strata_var  <- renderUI({ selectInput("strata_var", "Stratum column:",        col_choices(),
                                                   selected = if (is_demo() && "stratum" %in% col_choices()) "stratum") })
  output$ui_deff_var    <- renderUI({ selectInput("deff_var",   "DEFF variable:", col_choices()) })

  # ── Multi-level domain variable selector ─────────────────────────────────────
  output$ui_domain_vars <- renderUI({
    cols <- col_choices()
    tagList(
      tags$label("Domain variable(s):"),
      numericInput("n_domain_vars", "Number of domain variables:", 1, min = 1, max = 5),
      uiOutput("ui_domain_rows"),
      tags$small(class = "text-muted",
        "Level 1 = primary domain (e.g. Region), Level 2 = finer (e.g. Province), …",
        tags$br(),
        tags$b("Note:"), " all selected variables are used in computation.",
        " The finest level (last) drives strata statistics;",
        " coarser levels are looked up from the frame via nested merge.",
        tags$br(),
        tags$b("CV table tip:"), " set ", tags$b("Domain levels = N\u2009+\u20091"),
        " where N = number of domain variables above",
        " (DOM1 = national + one row per variable)."
      )
    )
  })

  output$ui_domain_rows <- renderUI({
    req(input$n_domain_vars)
    n  <- as.integer(input$n_domain_vars)
    cols <- col_choices()
    rows <- lapply(seq_len(n), function(i) {
      fluidRow(style = "margin-bottom:0px;",
        column(3,
          div(style = "margin-top:6px; text-align:center;",
            tags$span(class = "badge badge-primary",
                      style = "font-size:12px; padding:4px 7px;",
                      paste0("L", i))
          )
        ),
        column(9,
          selectInput(paste0("domain_col_", i),
                      label    = NULL,
                      choices  = cols,
                      selected = if (is_demo() && i == 1 && "region" %in% cols) "region"
                                 else if (length(cols) >= i) cols[i] else cols[1],
                      width    = "100%")
        )
      )
    })
    tagList(rows)
  })

  # Helper: returns ordered character vector of ALL selected domain columns
  collect_domain_vars <- function() {
    n <- isolate(input$n_domain_vars)   # isolate avoids reactive dependency
    n <- if (is.null(n) || is.na(n)) 1L else as.integer(n)
    vars <- character(n)
    for (i in seq_len(n)) {
      v <- input[[paste0("domain_col_", i)]]
      vars[i] <- if (!is.null(v) && nchar(v) > 0) v else ""
    }
    vars <- vars[nchar(vars) > 0]
    if (length(vars) == 0) NULL else vars
  }

  # Helper: returns ONLY the first (primary) domain variable as a scalar string.
  # R2BEAT functions (sensitivity_min_SSU, prepareInputToAllocation1, beat.2st)
  # accept a single domain_var; higher-level domains are tracked in the CV table.
  primary_domain_var <- function() {
    dv <- collect_domain_vars()
    if (is.null(dv)) NULL else dv[1]
  }

  output$ui_target_vars <- renderUI({
    checkboxGroupInput("target_vars", "Target variables (select ≥ 1):",
                       choices = col_choices(),
                       selected = if (is_demo())
                         intersect(c("active", "inactive", "unemployed", "income_hh"),
                                   col_choices()),
                       inline = FALSE)
  })

  output$ui_cv_table <- renderUI({
    req(input$n_dom, input$target_vars, length(input$target_vars) > 0)
    n_dom  <- as.integer(input$n_dom)
    n_var  <- length(input$target_vars)
    vnames <- input$target_vars

    header <- tags$tr(
      tags$th("Domain"),
      lapply(seq_len(n_var), function(j) tags$th(vnames[j]))
    )
    rows <- lapply(seq_len(n_dom), function(i) {
      tags$tr(
        tags$td(tags$strong(paste0("DOM", i))),
        lapply(seq_len(n_var), function(j) {
          tags$td(numericInput(paste0("cv_", i, "_", j), NULL,
                               value = 0.05, min = 1e-4, max = 1,
                               step = 0.005, width = "72px"))
        })
      )
    })
    tags$table(class = "table table-condensed table-bordered",
               tags$thead(header), tags$tbody(rows))
  })

  observeEvent(input$btn_set_equal, {
    req(input$n_dom, input$target_vars)
    val <- input$cv_equal_val
    req(!is.na(val))
    for (i in seq_len(as.integer(input$n_dom)))
      for (j in seq_len(length(input$target_vars)))
        updateNumericInput(session, paste0("cv_", i, "_", j), value = val)
  })

  # ── Inline Domain levels consistency indicator ───────────────────────────────
  output$ui_ndom_warning <- renderUI({
    n_dom   <- input$n_dom
    n_dvars <- input$n_domain_vars
    if (is.null(n_dom) || is.na(n_dom)) return(NULL)
    n_dom   <- as.integer(n_dom)
    n_dvars <- if (is.null(n_dvars) || is.na(n_dvars)) 1L else as.integer(n_dvars)
    expected <- n_dvars + 1L

    if (n_dom == expected) {
      div(class = "alert alert-success",
          style = "padding:6px 10px; font-size:12px; margin-bottom:6px;",
          icon("check-circle"),
          tags$b(" Correct:"), " ", n_dvars,
          " domain variable(s) \u2192 ", expected,
          " domain row(s) (DOM1 national + ", n_dvars, " sub-national).")
    } else {
      div(class = "alert alert-warning",
          style = "padding:8px 12px; font-size:12px; margin-bottom:6px;",
          icon("exclamation-triangle"),
          tags$b(" Mismatch:"), " ", n_dvars,
          " domain variable(s) require ",
          tags$b(paste0(expected, " domain rows")),
          " (DOM1 national + one per variable),",
          " but Domain levels is set to ", tags$b(n_dom), ".",
          tags$br(),
          "Adjust Domain levels to ", tags$b(expected),
          " before clicking Confirm Parameters.")
    }
  })

  observeEvent(input$upload_cv_file, {
    req(input$upload_cv_file)
    tryCatch({
      cv_data <- read.csv(input$upload_cv_file$datapath, stringsAsFactors = FALSE)

      # Validate structure
      if (!"DOM" %in% names(cv_data))
        stop("The CSV must contain a 'DOM' column.")
      cv_cols <- grep("^CV[0-9]+$", names(cv_data))
      if (length(cv_cols) == 0)
        stop("The CSV must contain at least one CV column named CV1, CV2, …")

      n_rows <- nrow(cv_data)
      n_cvs  <- length(cv_cols)

      # Persist for use after UI re-render
      rv$pending_cv <- cv_data

      # Update number of domain rows (triggers ui_cv_table to re-render)
      updateNumericInput(session, "n_dom", value = n_rows)

      # After a short delay (to allow the CV grid inputs to be (re-)created),
      # push the loaded values into every cv_{i}_{j} input.
      shinyjs::delay(700, {
        cv_d <- isolate(rv$pending_cv)
        if (!is.null(cv_d)) {
          cv_c <- grep("^CV[0-9]+$", names(cv_d))
          for (i in seq_len(nrow(cv_d)))
            for (j in seq_along(cv_c)) {
              val <- cv_d[i, cv_c[j]]
              if (!is.na(val))
                updateNumericInput(session, paste0("cv_", i, "_", j), value = val)
            }
          rv$pending_cv <- NULL
        }
      })

      showNotification(
        paste0("CV file loaded: ", n_rows, " domain row(s), ", n_cvs, " CV column(s)."),
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("CV file error:", e$message), type = "error", duration = 8)
    })
  })

  collect_cv <- function() {
    n_dom <- as.integer(input$n_dom)
    n_var <- length(input$target_vars)
    if (is.null(n_dom) || n_var == 0) return(NULL)
    cv_mat <- matrix(0.05, nrow = n_dom, ncol = n_var)
    for (i in seq_len(n_dom))
      for (j in seq_len(n_var)) {
        v <- input[[paste0("cv_", i, "_", j)]]
        if (!is.null(v) && length(v) == 1 && !is.na(v)) cv_mat[i, j] <- as.numeric(v)
      }
    df <- data.frame(DOM = paste0("DOM", seq_len(n_dom)), stringsAsFactors = FALSE)
    for (j in seq_len(n_var)) df[[paste0("CV", j)]] <- cv_mat[, j]
    df
  }

  output$txt_cv_preview <- renderPrint({
    cv <- collect_cv()
    if (is.null(cv)) return(cat("(set target variables first)"))
    print(cv)
  })

  observeEvent(input$btn_set_params, {
    req(rv$frame, input$target_vars, length(input$target_vars) > 0)
    cv  <- collect_cv()
    dvs <- collect_domain_vars()
    if (is.null(cv)) {
      showNotification("Select at least one target variable.", type = "warning")
      return()
    }
    if (is.null(dvs)) {
      showNotification("Select at least one domain variable.", type = "warning")
      return()
    }
    # Non-blocking check: n_dom should equal length(dvs) + 1
    expected_ndom <- length(dvs) + 1L
    actual_ndom   <- nrow(cv)
    if (actual_ndom != expected_ndom) {
      showNotification(
        tagList(
          icon("exclamation-triangle"),
          tags$b(" Domain levels mismatch: "),
          actual_ndom, " CV row(s) vs. ", expected_ndom,
          " expected (", length(dvs), " domain variable(s) + 1 national).",
          tags$br(),
          "Parameters confirmed anyway — adjust Domain levels and re-confirm if needed."
        ),
        type = "warning", duration = 12
      )
    }

    rv$cv          <- cv
    rv$domain_vars <- dvs            # snapshot — safe to use from any tab
    rv$target_vars <- input$target_vars
    rv$id_psu      <- input$id_psu
    rv$id_ssu      <- input$id_ssu
    rv$strata_var  <- input$strata_var
    rv$deff_var    <- input$deff_var
    rv$minimum     <- input$minimum
    rv$delta       <- input$delta
    rv$deff_sugg   <- input$deff_sugg
    rv$params_set  <- TRUE

    cv_r <- paste0(
      "cv <- as.data.frame(list(\n",
      "  DOM = c(", paste(shQuote(cv$DOM), collapse = ", "), "),\n",
      paste(sapply(seq_len(ncol(cv) - 1), function(j) {
        paste0("  CV", j, " = c(", paste(cv[[paste0("CV", j)]], collapse = ", "), ")")
      }), collapse = ",\n"),
      "\n))"
    )
    append_script(paste0(
      "# 2. Parameters\n",
      "id_PSU      <- ", shQuote(input$id_psu),      "\n",
      "id_SSU      <- ", shQuote(input$id_ssu),      "\n",
      "strata_var  <- ", shQuote(input$strata_var),  "\n",
      "target_vars <- c(", paste(shQuote(input$target_vars), collapse = ", "), ")\n",
      "deff_var    <- ", shQuote(input$deff_var),    "\n",
      "domain_var  <- c(", paste(shQuote(collect_domain_vars()), collapse = ", "), ")\n",
      "minimum     <- ", input$minimum,              "\n",
      "delta       <- ", input$delta,                "\n",
      "deff_sugg   <- ", input$deff_sugg,            "\n",
      cv_r, "\n"
    ))
    showNotification("Parameters confirmed.", type = "message")
  })

  # ── STEP 3: SENSITIVITY ─────────────────────────────────────────────────────
  observeEvent(input$btn_run_sens, {
    req(rv$params_set, rv$frame, rv$cv, rv$domain_vars, rv$target_vars,
        rv$id_psu, rv$id_ssu, rv$strata_var, rv$deff_var)
    tryCatch({
      frame <- rv$frame
      withProgress(message = "Running sensitivity analysis…", value = 0.5, {
        sens <- sensitivity_min_SSU(
          samp_frame  = frame,
          errors      = rv$cv,
          id_PSU      = rv$id_psu,
          id_SSU      = rv$id_ssu,
          strata_var  = rv$strata_var,
          target_vars = rv$target_vars,
          deff_var    = rv$deff_var,
          domain_var  = rv$domain_vars,
          delta       = rv$delta,
          deff_sugg   = rv$deff_sugg,
          min         = input$sens_min,
          max         = input$sens_max,
          plot        = FALSE
        )
      })
      rv$sens_result <- list(result  = sens,
                             min_val = input$sens_min,
                             max_val = input$sens_max)
      append_script(paste0(
        "\n# 3. Sensitivity analysis (optional)\n",
        "sens_min_SSU <- sensitivity_min_SSU(\n",
        "  samp_frame=samp_frame, errors=cv,\n",
        "  id_PSU=", shQuote(rv$id_psu), ", id_SSU=", shQuote(rv$id_ssu), ",\n",
        "  strata_var=", shQuote(rv$strata_var), ",\n",
        "  target_vars=c(", paste(shQuote(rv$target_vars), collapse = ", "), "),\n",
        "  deff_var=", shQuote(rv$deff_var), ",\n",
        "  domain_var=c(", paste(shQuote(rv$domain_vars), collapse = ", "), "),\n",
        "  delta=", rv$delta, ", deff_sugg=", rv$deff_sugg, ",\n",
        "  min=", input$sens_min, ", max=", input$sens_max, ", plot=TRUE)\n"
      ))
      showNotification("Sensitivity analysis complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Sensitivity error:", e$message), type = "error", duration = 8)
    })
  })

  output$plot_sens <- renderPlot({
    req(rv$sens_result)
    plot_sens(rv$sens_result$result,
              rv$sens_result$min_val,
              rv$sens_result$max_val)
  })

  output$tbl_sens <- renderDT({
    req(rv$sens_result)
    psu_v <- rv$sens_result$result[[1]]
    ssu_v <- rv$sens_result$result[[2]]
    n     <- length(psu_v)
    df    <- data.frame(
      `Min SSU` = seq(rv$sens_result$min_val, length.out = n),
      `PSU`     = psu_v,
      `SSU`     = ssu_v,
      check.names = FALSE
    )
    datatable(df, options = list(dom = "t", scrollX = TRUE, pageLength = 20))
  })

  # ── STEP 4: PREPARE INPUT ───────────────────────────────────────────────────
  observeEvent(input$btn_run_prep, {
    req(rv$params_set, rv$frame, rv$cv, rv$domain_vars, rv$target_vars,
        rv$id_psu, rv$id_ssu, rv$strata_var, rv$deff_var)
    tryCatch({
      frame <- rv$frame
      withProgress(message = "Running prepareInputToAllocation1()…", value = 0.5, {
        log_txt <- capture.output({
          inp <- prepareInputToAllocation1(
            samp_frame  = frame,
            id_PSU      = rv$id_psu,
            id_SSU      = rv$id_ssu,
            strata_var  = rv$strata_var,
            target_vars = rv$target_vars,
            deff_var    = rv$deff_var,
            domain_var  = rv$domain_vars,
            minimum     = rv$minimum,
            delta       = rv$delta,
            deff_sugg   = rv$deff_sugg
          )
        })
      })
      rv$inp      <- inp
      rv$prep_log <- log_txt
      append_script(paste0(
        "\n# 4. Prepare allocation input\n",
        "inp <- prepareInputToAllocation1(\n",
        "  samp_frame=samp_frame,\n",
        "  id_PSU=", shQuote(rv$id_psu), ", id_SSU=", shQuote(rv$id_ssu), ",\n",
        "  strata_var=", shQuote(rv$strata_var), ",\n",
        "  target_vars=c(", paste(shQuote(rv$target_vars), collapse = ", "), "),\n",
        "  deff_var=", shQuote(rv$deff_var),
        ", domain_var=c(", paste(shQuote(rv$domain_vars), collapse = ", "), "),\n",
        "  minimum=", rv$minimum, ", delta=", rv$delta,
        ", deff_sugg=", rv$deff_sugg, "\n)\n"
      ))
      showNotification("Preparation complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Preparation error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_prep_log   <- renderText({ req(rv$prep_log); paste(rv$prep_log, collapse = "\n") })
  output$tbl_inp_strata <- renderDT({ req(rv$inp); datatable(head(rv$inp$strata, 30), options = list(scrollX = TRUE, pageLength = 8)) })
  output$tbl_inp_des    <- renderDT({ req(rv$inp); datatable(rv$inp$des_file,        options = list(scrollX = TRUE, pageLength = 8)) })
  output$tbl_inp_psu    <- renderDT({ req(rv$inp); datatable(head(rv$inp$psu_file, 30), options = list(scrollX = TRUE, pageLength = 8)) })
  output$txt_inp_deff   <- renderPrint({ req(rv$inp); print(rv$inp$deff)  })
  output$txt_inp_effst  <- renderPrint({ req(rv$inp); print(head(rv$inp$effst, 10)) })
  output$txt_inp_rho    <- renderPrint({ req(rv$inp); print(head(rv$inp$rho,   10)) })

  # ── STEP 5: ALLOCATION ──────────────────────────────────────────────────────
  observeEvent(input$btn_run_alloc, {
    req(rv$inp, rv$cv)
    tryCatch({
      inp <- rv$inp
      inp$des_file$MINIMUM <- input$minimum
      withProgress(message = "Running beat.2st()…", value = 0.5, {
        alloc <- beat.2st(
          stratif     = inp$strata,
          errors      = rv$cv,
          des_file    = inp$des_file,
          psu_file    = inp$psu_file,
          rho         = inp$rho,
          deft_start  = NULL,
          effst       = inp$effst,
          minPSUstrat = input$alloc_minPSU,
          minnumstrat = input$alloc_minnumstrat
        )
      })
      rv$alloc <- alloc
      append_script(paste0(
        "\n# 5. Two-stage allocation\n",
        "inp$des_file$MINIMUM <- ", input$minimum, "\n",
        "alloc <- beat.2st(\n",
        "  stratif=inp$strata, errors=cv,\n",
        "  des_file=inp$des_file, psu_file=inp$psu_file,\n",
        "  rho=inp$rho, deft_start=NULL, effst=inp$effst,\n",
        "  minPSUstrat=", input$alloc_minPSU,
        ", minnumstrat=", input$alloc_minnumstrat, "\n)\n"
      ))
      showNotification("Allocation complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Allocation error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_alloc_total <- renderText({
    req(rv$alloc)
    iters <- rv$alloc$iterations
    psu_col <- "PSU Total"
    if (!is.null(iters) && is.data.frame(iters) && nrow(iters) > 0 && psu_col %in% names(iters))
      format(iters[[psu_col]][nrow(iters)], big.mark = ",")
    else "N/A"
  })

  output$txt_alloc_ssu <- renderText({
    req(rv$alloc)
    iters <- rv$alloc$iterations
    if (!is.null(iters) && is.data.frame(iters) && nrow(iters) > 0 && "SSU" %in% names(iters))
      format(iters$SSU[nrow(iters)], big.mark = ",")
    else "N/A"
  })

  output$tbl_alloc <- renderDT({
    req(rv$alloc)
    datatable(rv$alloc$alloc, options = list(scrollX = TRUE, pageLength = 10))
  })
  output$tbl_expected <- renderDT({
    req(rv$alloc)
    datatable(rv$alloc$expected, options = list(dom = "t", scrollX = TRUE))
  })
  output$tbl_sensitivity <- renderDT({
    req(rv$alloc)
    datatable(rv$alloc$sensitivity, options = list(dom = "t", scrollX = TRUE))
  })

  # ── STEP 6: SELECT SAMPLE ───────────────────────────────────────────────────
  observeEvent(input$btn_select_psu, {
    req(rv$alloc)
    tryCatch({
      set.seed(input$seed)
      withProgress(message = "Selecting PSUs…", value = 0.5, {
        s1 <- select_PSU(rv$alloc, type = input$psu_type,
                         pps = isTRUE(input$pps), plot = FALSE)
      })
      rv$sample_1st <- s1
      append_script(paste0(
        "\n# 6a. Select PSUs\n",
        "set.seed(", input$seed, ")\n",
        "sample_1st <- select_PSU(alloc, type=", shQuote(input$psu_type),
        ", pps=", isTRUE(input$pps), ", plot=FALSE)\n"
      ))
      showNotification("PSU selection complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("PSU selection error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_psu_stats <- renderPrint({
    req(rv$sample_1st); print(rv$sample_1st$PSU_stats)
  })

  # SSU code per select_SSU (puo' differire dall'id_ssu usato per le statistiche).
  output$ui_ssu_code <- renderUI({
    req(rv$frame)
    selectInput("ssu_select_code", "SSU selection variable (frame):",
                choices  = names(rv$frame),
                # Demo: la selezione fisica delle SSU usa id_ind (vedi esempio R2BEAT
                # select_SSU); le statistiche/allocazione usano id_hh (passo 2).
                selected = if (identical(input$data_source, "demo") &&
                               "id_ind" %in% names(rv$frame)) "id_ind"
                           else rv$id_ssu %||% names(rv$frame)[1])
  })

  observeEvent(input$btn_select_ssu, {
    req(rv$sample_1st, rv$frame, input$ssu_select_code)
    tryCatch({
      withProgress(message = "Selecting SSUs…", value = 0.5, {
        samp <- select_SSU(
          df          = rv$frame,
          PSU_code    = input$id_psu,
          SSU_code    = input$ssu_select_code,
          PSU_sampled = rv$sample_1st$sample_PSU
        )
      })
      rv$samp <- samp
      append_script(paste0(
        "\n# 6b. Select SSUs\n",
        "samp <- select_SSU(\n",
        "  df=samp_frame, PSU_code=", shQuote(input$id_psu),
        ", SSU_code=", shQuote(input$ssu_select_code), ",\n",
        "  PSU_sampled=sample_1st$sample_PSU\n)\n",
        "cat('Sample size (SSU):', nrow(samp), '\\n')\n",
        "cat('Pop. total (sum weights):', round(sum(samp$weight)), '\\n')\n"
      ))
      showNotification("SSU selection complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("SSU selection error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_samp_size <- renderText({ req(rv$samp); format(nrow(rv$samp), big.mark = ",") })
  output$txt_pop_n     <- renderText({ req(rv$samp); format(round(sum(rv$samp$weight)), big.mark = ",") })

  output$tbl_samp_head <- renderDT({
    req(rv$samp)
    datatable(head(rv$samp, 50), options = list(scrollX = TRUE, pageLength = 8))
  })

  output$btn_download_samp <- downloadHandler(
    filename = function() "r2beat_sample.csv",
    content  = function(file) { req(rv$samp); write.csv(rv$samp, file, row.names = FALSE) }
  )

  # Weight distribution plots
  output$ui_wt_domain_title  <- renderUI({
    dv <- collect_domain_vars()
    paste("Weight distribution — by", if (!is.null(dv)) dv[1] else "domain")
  })
  output$ui_wt_stratum_title <- renderUI({ req(input$strata_var); paste("Weight distribution — by", input$strata_var) })

  output$plot_wt_total    <- renderPlot({ req(rv$samp); wf1_weight_boxplot(rv$samp, NULL, "Weight distribution — total sample") })
  output$plot_wt_domain   <- renderPlot({
    req(rv$samp)
    dv <- collect_domain_vars()
    d1 <- if (!is.null(dv)) dv[1] else NULL
    wf1_weight_boxplot(rv$samp, d1, paste("Weights by", d1 %||% "domain"))
  })
  output$plot_wt_province <- renderPlot({
    req(rv$samp)
    dv      <- collect_domain_vars()
    # Level 2 variable if available, else fall back to "province" column, then level 1
    prov_col <- if (!is.null(dv) && length(dv) >= 2) dv[2] else
                  if ("province" %in% names(rv$samp)) "province" else
                    if (!is.null(dv)) dv[1] else NULL
    wf1_weight_boxplot(rv$samp, prov_col, paste("Weights by", prov_col %||% "province"))
  })
  output$plot_wt_stratum  <- renderPlot({ req(rv$samp); wf1_weight_boxplot(rv$samp, input$strata_var, paste("Weights by", input$strata_var))             })

  # ── STEP 7: R SCRIPT ────────────────────────────────────────────────────────
  output$txt_script <- renderText({
    if (length(rv$script_lines) == 0)
      return("# No steps recorded yet. Complete the workflow to generate a script.")
    paste(rv$script_lines, collapse = "\n")
  })

  output$btn_dl_script <- downloadHandler(
    filename = function() "r2beat_workflow.R",
    content  = function(file) writeLines(rv$script_lines, file)
  )

} # end server


# ============================ WORKFLOW: WF2 ============================



# ─── HELPERS ──────────────────────────────────────────────────────────────────


wf2_weight_boxplot <- function(df, grp_var, title_txt) {
  if (!is.null(grp_var) && nchar(grp_var) > 0 && grp_var %in% names(df)) {
    df[[grp_var]] <- factor(df[[grp_var]])
    ggplot(df, aes(x = .data[[grp_var]], y = weight)) +
      geom_boxplot(fill = "steelblue", colour = "grey30", outlier.size = 0.6) +
      labs(title = title_txt, x = NULL, y = "Weight") +
      coord_flip() + theme_minimal(base_size = 11)
  } else {
    ggplot(df, aes(x = "", y = weight)) +
      geom_boxplot(fill = "steelblue", colour = "grey30") +
      labs(title = title_txt, x = NULL, y = "Weight") +
      theme_minimal(base_size = 11)
  }
}

# ─── UI ───────────────────────────────────────────────────────────────────────

wf2_ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "R2BEAT - Workflow 2"),

  dashboardSidebar(
    sidebarMenu(id = "tabs",
      menuItem("1. Load Data",     tabName = "tab_data",   icon = icon("database")),
      menuItem("2. Parameters",    tabName = "tab_params", icon = icon("sliders-h")),
      menuItem("3. Prepare Input", tabName = "tab_prep",   icon = icon("cogs")),
      menuItem("4. Allocation",    tabName = "tab_alloc",  icon = icon("bullseye")),
      menuItem("5. Select Sample", tabName = "tab_samp",   icon = icon("users")),
      menuItem("6. R Script",      tabName = "tab_script", icon = icon("file-code"))
    ),
    div(style = "padding:12px 15px; border-top:1px solid #444; margin-top:10px;",
      actionButton("btn_home", " Menu principale", icon = icon("home"),
                   class = "btn-default btn-block",
                   style = "font-weight:bold; margin-bottom:8px;",
                   onclick = "window.location.search='';"),
      actionButton("btn_help", " Help Guide", icon = icon("question-circle"),
                   class = "btn-info btn-block", style = "font-weight:bold;")
    )
  ),

  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .box-body .dataTables_wrapper { overflow-x: auto; }
      pre { max-height: 380px; overflow-y: auto; }
      .cv-grid input[type=number] { width: 72px; }
      .cv-grid td { vertical-align: middle !important; padding: 3px 5px !important; }
      .cv-grid th { text-align: center; padding: 4px 5px; }
    "))),
    tabItems(

      # ── TAB 1: LOAD DATA ──────────────────────────────────────────────────
      tabItem(tabName = "tab_data",
        fluidRow(
          box(title = "ReGenesees Objects", width = 4, status = "primary", solidHeader = TRUE,
            p("Load the ReGenesees design and calibrated objects from a previous survey round."),
            radioButtons("rg_source", "Source:",
              choices  = c("Use bundled demo (R2BEAT_ReGenesees)" = "demo",
                           "Upload .RData file"                   = "upload"),
              selected = "demo"),
            conditionalPanel("input.rg_source == 'upload'",
              fileInput("upload_rg", "Upload .RData file (must contain 'des' and 'cal'):",
                        accept = c(".RData", ".rdata", ".rda")),
              div(class = "alert alert-info", style = "padding:6px 10px; font-size:12px;",
                icon("info-circle"),
                " The file must contain a ReGenesees design object named ",
                tags$b("des"), " and a calibrated object named ", tags$b("cal"), "."
              )
            ),
            hr(),
            p("Population frame (needed for SSU selection):"),
            radioButtons("frame_source", "Source:",
              choices  = c("Use bundled pop.RData (same as demo)" = "demo",
                           "Upload .RData file"                   = "upload_rdata",
                           "Upload delimited file (CSV/TSV)"      = "delimited"),
              selected = "demo"),
            conditionalPanel("input.frame_source == 'upload_rdata'",
              fileInput("upload_frame_rdata", "Upload .RData file:",
                        accept = c(".RData", ".rdata", ".rda"))
            ),
            conditionalPanel("input.frame_source == 'delimited'",
              fileInput("upload_frame_csv", "Upload CSV/TSV:",
                        accept = c(".csv", ".tsv", ".txt")),
              fluidRow(
                column(6,
                  radioButtons("csv_sep", "Separator:",
                    choices = c(Comma=",", Semicolon=";", Tab="\t"),
                    selected = ",", inline = FALSE)
                ),
                column(6,
                  checkboxInput("csv_header",    "Header row",         value = TRUE),
                  checkboxInput("csv_dec_comma", "Decimal comma ( , )", value = FALSE)
                )
              )
            ),
            actionButton("btn_load", "Load Data", class = "btn-primary btn-block",
                         icon = icon("upload"))
          ),
          box(title = "ReGenesees — Design Object Variables", width = 4, status = "info",
            verbatimTextOutput("txt_des_str")
          ),
          box(title = "Checks & Frame Summary", width = 4, status = "info",
            h5("Lonely strata in design object:"),
            verbatimTextOutput("txt_lonely"),
            hr(),
            h5("Population frame:"),
            verbatimTextOutput("txt_frame_summary")
          )
        ),
        fluidRow(
          box(title = "Frame Preview (first 100 rows)", width = 12, status = "info",
            DTOutput("tbl_frame_head")
          )
        )
      ),

      # ── TAB 2: PARAMETERS ─────────────────────────────────────────────────
      tabItem(tabName = "tab_params",
        fluidRow(
          box(title = "Variable Assignment & Numeric Parameters",
              width = 4, status = "primary", solidHeader = TRUE,
            uiOutput("ui_id_psu"),
            uiOutput("ui_id_ssu"),
            uiOutput("ui_strata_var"),
            uiOutput("ui_target_vars"),
            uiOutput("ui_deff_var"),
            uiOutput("ui_domain_vars"),
            hr(),
            numericInput("minimum", "Minimum SSU per PSU:", 25, 2, 500),
            numericInput("delta",   "Delta (intraclass corr.):", 1, 0, 10),
            hr(),
            actionButton("btn_set_params", "Confirm Parameters",
                         class = "btn-success btn-block", icon = icon("check"))
          ),
          box(title = "CV Constraints", width = 8, status = "primary", solidHeader = TRUE,
            p("One row per domain level: ", tags$b("DOM1"), " = national,",
              tags$b(" DOM2"), " = Level-1 domain variable,",
              tags$b(" DOM3"), " = Level-2, …",
              " One column per selected target variable.",
              " Set ", tags$b("Domain levels = N\u2009+\u20091"),
              " where N is the number of domain variables selected on the left."),
            fluidRow(
              column(3,
                numericInput("n_dom", "Domain levels:", 2, 1, 10)
              ),
              column(4,
                numericInput("cv_equal_val", "Set all equal to:", 0.05, 0.001, 1, 0.005)
              ),
              column(3, br(),
                actionButton("btn_set_equal", "Apply", class = "btn-default btn-block")
              )
            ),
            uiOutput("ui_ndom_warning"),
            fluidRow(
              column(5,
                fileInput("upload_cv_file", "Load CV values from CSV:",
                          accept = c(".csv", ".txt"),
                          buttonLabel = "Browse…",
                          placeholder = "DOM, CV1, CV2, … format")
              ),
              column(7,
                div(style = "margin-top:28px;",
                  div(class = "alert alert-info",
                      style = "padding:6px 10px; font-size:12px; margin-bottom:0;",
                    icon("info-circle"),
                    " CSV must contain a ", tags$b("DOM"), " column and one or more ",
                    tags$b("CV1, CV2, …"), " columns."
                  )
                )
              )
            ),
            hr(),
            div(class = "cv-grid", uiOutput("ui_cv_table")),
            hr(),
            h5("Current CV dataframe:"),
            verbatimTextOutput("txt_cv_preview")
          )
        )
      ),

      # ── TAB 3: PREPARE INPUT ──────────────────────────────────────────────
      tabItem(tabName = "tab_prep",
        fluidRow(
          box(title = "Prepare Allocation Input", width = 3,
              status = "primary", solidHeader = TRUE,
            p("Calls ", tags$code("prepareInputToAllocation2()"),
              " using the ReGenesees objects to compute strata statistics,",
              " design effects, and intraclass correlations."),
            p("Then runs ", tags$code("check_input()"),
              " to reconcile population counts between the strata and PSU files."),
            div(class = "alert alert-warning", style = "padding:8px; font-size:12px;",
              icon("clock-o"),
              " May take 1–2 minutes."
            ),
            actionButton("btn_run_prep", "Run Preparation",
                         class = "btn-danger btn-block", icon = icon("cogs")),
            br(),
            p(em("Requires Step 2 to be confirmed first."))
          ),
          box(title = "Execution Log", width = 9, status = "warning",
            verbatimTextOutput("txt_prep_log")
          )
        ),
        fluidRow(
          box(title = "Strata (head)", width = 6, status = "info",
            DTOutput("tbl_inp_strata")
          ),
          box(title = "Design File (des_file)", width = 6, status = "info",
            DTOutput("tbl_inp_des")
          )
        ),
        fluidRow(
          box(title = "PSU File (head)", width = 6, status = "info",
            DTOutput("tbl_inp_psu")
          ),
          box(title = "Design Parameters", width = 6, status = "info",
            h5("DEFF (head)"),
            verbatimTextOutput("txt_inp_deff"),
            h5("Effective stratum size (effst — head)"),
            verbatimTextOutput("txt_inp_effst"),
            h5("Intraclass correlation (rho — head)"),
            verbatimTextOutput("txt_inp_rho")
          )
        )
      ),

      # ── TAB 4: ALLOCATION ─────────────────────────────────────────────────
      tabItem(tabName = "tab_alloc",
        fluidRow(
          box(title = "Allocation Settings", width = 3,
              status = "primary", solidHeader = TRUE,
            numericInput("alloc_minimum",     "Minimum SSU per PSU (MINIMUM):", 25, 2, 500),
            numericInput("alloc_minPSU",      "Min PSU per stratum (minPSUstrat):", 2, 1, 50),
            numericInput("alloc_minnumstrat", "Min strata (minnumstrat):", 2, 1, 20),
            br(),
            actionButton("btn_run_alloc", "Run Allocation (beat.2st)",
                         class = "btn-danger btn-block", icon = icon("bullseye")),
            br(),
            p(em("Requires Step 3 to be completed first."))
          ),
          box(title = "Allocation Results", width = 9, status = "info",
            fluidRow(
              column(4,
                h4("Total PSUs (1st stage)"),
                verbatimTextOutput("txt_alloc_total")
              ),
              column(4,
                h4("Total SSUs (2nd stage)"),
                verbatimTextOutput("txt_alloc_ssu")
              )
            ),
            hr(),
            h4("Allocation by stratum"),
            DTOutput("tbl_alloc")
          )
        ),
        fluidRow(
          box(title = "Expected CVs", width = 6, status = "info",
            DTOutput("tbl_expected")
          ),
          box(title = "Sensitivity to CV constraints", width = 6, status = "info",
            p(em("Values > 1 indicate that a constraint is active (binding).")),
            DTOutput("tbl_sensitivity")
          )
        )
      ),

      # ── TAB 5: SELECT SAMPLE ──────────────────────────────────────────────
      tabItem(tabName = "tab_samp",
        fluidRow(
          box(title = "Sampling Settings", width = 3,
              status = "primary", solidHeader = TRUE,
            numericInput("seed", "Random seed:", 1234, 1, 99999),
            hr(),
            h4("Step 1 — Select PSUs"),
            selectInput("psu_type", "Allocation type:",
              choices = c("ALLOC", "PROP", "EQUAL"), selected = "ALLOC"),
            checkboxInput("pps", "PPS sampling", value = TRUE),
            actionButton("btn_select_psu", "Select PSUs",
                         class = "btn-warning btn-block", icon = icon("map-marker")),
            hr(),
            h4("Step 2 — Select SSUs"),
            uiOutput("ui_ssu_code"),
            actionButton("btn_select_ssu", "Select SSUs",
                         class = "btn-success btn-block", icon = icon("users")),
            hr(),
            downloadButton("btn_download_samp", "Download Sample (CSV)",
                           class = "btn-block"),
            br(),
            p(em("Requires Step 4 to be completed first."))
          ),
          box(title = "Sample Results", width = 9, status = "info",
            h4("PSU statistics"),
            verbatimTextOutput("txt_psu_stats"),
            hr(),
            fluidRow(
              column(6,
                h4("Sample size (SSU)"),
                verbatimTextOutput("txt_samp_size")
              ),
              column(6,
                h4("Reconstructed population total (sum of weights)"),
                verbatimTextOutput("txt_pop_n")
              )
            ),
            hr(),
            h4("Sample preview"),
            DTOutput("tbl_samp_head")
          )
        ),
        fluidRow(
          box(title = "Weight distribution — total sample", width = 6,
            plotOutput("plot_wt_total",  height = "260px")
          ),
          box(title = uiOutput("ui_wt_domain_title"), width = 6,
            plotOutput("plot_wt_domain", height = "260px")
          )
        ),
        fluidRow(
          box(title = "Weight distribution — by province", width = 6,
            plotOutput("plot_wt_province", height = "260px")
          ),
          box(title = uiOutput("ui_wt_stratum_title"), width = 6,
            plotOutput("plot_wt_stratum", height = "260px")
          )
        )
      ),

      # ── TAB 6: R SCRIPT ───────────────────────────────────────────────────
      tabItem(tabName = "tab_script",
        fluidRow(
          box(title = "Generated R Script", width = 12, status = "info",
            p("Reproducible script recording all steps executed in this session."),
            verbatimTextOutput("txt_script")
          )
        ),
        fluidRow(
          box(title = "Download", width = 4, status = "primary", solidHeader = TRUE,
            downloadButton("btn_dl_script", "Download R Script", class = "btn-block")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
)

# ─── SERVER ───────────────────────────────────────────────────────────────────

wf2_server <- function(input, output, session) {

  rv <- reactiveValues(
    des          = NULL,   # ReGenesees design object
    cal          = NULL,   # ReGenesees calibrated object
    frame        = NULL,   # population frame (for SSU selection)
    params_set   = FALSE,
    cv           = NULL,
    domain_vars  = NULL,
    target_vars  = NULL,
    id_psu       = NULL,
    id_ssu       = NULL,
    strata_var   = NULL,
    deff_var     = NULL,
    minimum      = NULL,
    delta        = NULL,
    pending_cv   = NULL,
    inp          = NULL,
    prep_log     = NULL,
    alloc        = NULL,
    sample_1st   = NULL,
    samp         = NULL,
    script_lines = character(0)
  )

  append_script <- function(...) {
    rv$script_lines <- c(rv$script_lines, paste0(...))
  }

  # ── HELP MODAL ──────────────────────────────────────────────────────────────
  observeEvent(input$btn_help, {
    showModal(modalDialog(
      title = tags$span(icon("question-circle"), " R2BEAT Sample Workflow — Guide"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      tabsetPanel(
        tabPanel(tags$span(icon("database"), " 1. Load Data"),
          br(),
          tags$h4("ReGenesees objects"),
          tags$p("Load the ", tags$b("des"), " (survey design) and ", tags$b("cal"),
                 " (calibrated) objects produced by the ", tags$code("ReGenesees"),
                 " package from a previous survey round."),
          tags$p("The demo uses the ", tags$code("R2BEAT_ReGenesees"), " dataset bundled",
                 " with R2BEAT, which contains a synthetic survey design and calibration."),
          tags$h4("Population frame"),
          tags$p("The full population frame is needed only for the SSU selection step (Step 5).",
                 " It must contain the PSU code, the SSU code, and the stratum variable.")
        ),
        tabPanel(tags$span(icon("sliders-h"), " 2. Parameters"),
          br(),
          tags$h4("Variable assignment"),
          tags$ul(
            tags$li(tags$b("PSU identifier:"), " Column identifying Primary Sampling Units."),
            tags$li(tags$b("SSU identifier:"), " Column for Secondary Sampling Units used in",
                    " strata statistics (e.g. households)."),
            tags$li(tags$b("Stratum column:"), " Pre-defined stratification variable."),
            tags$li(tags$b("Target variables:"), " Columns to be estimated (drive CV constraints)."),
            tags$li(tags$b("DEFF grouping variable:"), " A ", tags$b("categorical"), " column",
                    " used to group the design-effect computation (typically the same as the",
                    " stratum column, e.g. ", tags$code("stratum"), ")."),
            tags$li(tags$b("Domain variable(s):"), " Column(s) defining estimation domains.")
          ),
          tags$h4("CV constraints"),
          tags$p("One row per domain level (DOM1 = national, DOM2 = Level-1 domain, …),",
                 " one column per target variable."),
          tags$p(tags$b("Rule: Domain levels = N\u2009+\u20091"),
                 " where N = number of domain variables selected.")
        ),
        tabPanel(tags$span(icon("cogs"), " 3. Prepare Input"),
          br(),
          tags$h4("prepareInputToAllocation2()"),
          tags$p("Uses the ReGenesees design and calibrated objects to compute:"),
          tags$ul(
            tags$li(tags$b("strata:"), " Population counts (N), means (M1…Mv), SDs (S1…Sv),",
                    " domain assignments (DOM1…DOMk) per stratum."),
            tags$li(tags$b("deff:"), " Design effects per stratum / variable."),
            tags$li(tags$b("effst:"), " Estimator effects."),
            tags$li(tags$b("rho:"), " Intraclass correlations (AR and NAR)."),
            tags$li(tags$b("psu_file:"), " PSU-level measure of size."),
            tags$li(tags$b("des_file:"), " Stratum-level design parameters.")
          ),
          tags$h4("check_input()"),
          tags$p("Reconciles population counts between the strata and PSU files.",
                 " The corrected strata file replaces the original before allocation.")
        ),
        tabPanel(tags$span(icon("bullseye"), " 4. Allocation"),
          br(),
          tags$h4("beat.2st()"),
          tags$p("Two-stage Bethel optimal allocation minimising total PSUs",
                 " subject to CV constraints."),
          tags$ul(
            tags$li(tags$b("MINIMUM:"), " Minimum SSUs per PSU (can differ from Step 2 value)."),
            tags$li(tags$b("minPSUstrat:"), " Minimum PSUs per stratum."),
            tags$li(tags$b("minnumstrat:"), " Minimum number of strata in each domain.")
          )
        ),
        tabPanel(tags$span(icon("users"), " 5. Select Sample"),
          br(),
          tags$h4("Step 1 — select_PSU()"),
          tags$p("PPS first-stage selection of PSUs according to allocation type",
                 " (ALLOC, PROP, or EQUAL)."),
          tags$h4("Step 2 — select_SSU()"),
          tags$p("Second-stage selection of SSUs within selected PSUs.",
                 " The ", tags$b("SSU selection variable"),
                 " may differ from the SSU identifier used in Step 2",
                 " (e.g., households for statistics vs. individuals for selection)."),
          tags$h4("Weight distribution plots"),
          tags$p("Four boxplots check weight variability: overall sample, by domain,",
                 " by province, and by stratum.")
        )
      )
    ))
  })

  # ── STEP 1: LOAD DATA ───────────────────────────────────────────────────────
  observeEvent(input$btn_load, {
    tryCatch({
      # Load ReGenesees objects
      if (input$rg_source == "demo") {
        withProgress(message = "Loading R2BEAT_ReGenesees…", value = 0.3, {
          env_rg <- new.env()
          load("R2BEAT_ReGenesees.RData", envir = env_rg)
          if (!exists("des", envir = env_rg) || !exists("cal", envir = env_rg))
            stop("Demo data must contain 'des' and 'cal' objects.")
          rv$des <- env_rg$des
          rv$cal <- env_rg$cal
        })
        rg_expr <- 'load("R2BEAT_ReGenesees.RData")\nRGdes <- des\nRGcal <- cal'
      } else {
        req(input$upload_rg)
        env_rg <- new.env()
        load(input$upload_rg$datapath, envir = env_rg)
        if (!exists("des", envir = env_rg) || !exists("cal", envir = env_rg))
          stop("Uploaded file must contain objects named 'des' and 'cal'.")
        rv$des <- env_rg$des
        rv$cal <- env_rg$cal
        rg_expr <- paste0('load("', input$upload_rg$name, '")\nRGdes <- des\nRGcal <- cal')
      }

      # Load population frame
      if (input$frame_source == "demo") {
        withProgress(message = "Loading pop.RData…", value = 0.6, {
          env_fr <- new.env()
          load("pop.RData", envir = env_fr)
          if (!exists("pop", envir = env_fr))
            stop("pop.RData must contain an object named 'pop'.")
          rv$frame <- env_fr$pop
        })
        frame_expr <- 'load("pop.RData")\nsamp_frame <- pop'
      } else if (input$frame_source == "upload_rdata") {
        req(input$upload_frame_rdata)
        env_fr <- new.env()
        load(input$upload_frame_rdata$datapath, envir = env_fr)
        objs <- ls(envir = env_fr)
        if (length(objs) == 0) stop("No objects found in the frame RData file.")
        rv$frame  <- env_fr[[objs[1]]]
        frame_expr <- paste0('load("', input$upload_frame_rdata$name, '")\nsamp_frame <- ',
                             objs[1])
      } else {
        req(input$upload_frame_csv)
        dec <- if (isTRUE(input$csv_dec_comma)) "," else "."
        rv$frame <- read.table(input$upload_frame_csv$datapath,
                               header = isTRUE(input$csv_header),
                               sep = input$csv_sep, dec = dec,
                               stringsAsFactors = FALSE)
        frame_expr <- paste0('samp_frame <- read.table("', input$upload_frame_csv$name,
                             '", header=', isTRUE(input$csv_header),
                             ', sep="', input$csv_sep, '", dec="', dec, '",',
                             ' stringsAsFactors=FALSE)')
      }

      rv$script_lines <- c(
        "# R2BEAT sample workflow — generated by Shiny app",
        "library(R2BEAT)",
        "",
        "# 1. Load data",
        rg_expr, frame_expr, ""
      )
      showNotification("Data loaded successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Load error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_des_str <- renderPrint({
    req(rv$des)
    cat("Class:", class(rv$des), "\n\n")
    cat("Variables in des$variables:\n")
    str(rv$des$variables, give.attr = FALSE)
  })

  output$txt_lonely <- renderPrint({
    req(rv$des)
    tryCatch({
      if (requireNamespace("ReGenesees", quietly = TRUE)) {
        ls_res <- ReGenesees::find.lon.strata(rv$des)
        if (is.null(ls_res) || nrow(as.data.frame(ls_res)) == 0)
          cat("No lonely strata found.\n")
        else
          print(ls_res)
      } else {
        cat("ReGenesees package not available — cannot check lonely strata.\n")
      }
    }, error = function(e) cat("Check skipped:", e$message, "\n"))
  })

  output$txt_frame_summary <- renderPrint({
    req(rv$frame)
    df <- rv$frame
    cat("Rows:", format(nrow(df), big.mark = ","), "\n")
    cat("Cols:", ncol(df), "\n")
    cat("Columns:", paste(names(df), collapse = ", "), "\n")
  })

  output$tbl_frame_head <- renderDT({
    req(rv$frame)
    datatable(head(rv$frame, 100), options = list(scrollX = TRUE, pageLength = 8))
  })

  # ── STEP 2: PARAMETERS ──────────────────────────────────────────────────────
  # Column choices come from des$variables (survey design variables)
  col_choices <- reactive({
    req(rv$des)
    names(rv$des$variables)
  })

  is_demo <- reactive({ identical(input$rg_source, "demo") })
  output$ui_id_psu     <- renderUI({ selectInput("id_psu",     "PSU identifier column:", col_choices(),
                                                  selected = if (is_demo() && "municipality" %in% col_choices()) "municipality") })
  output$ui_id_ssu     <- renderUI({ selectInput("id_ssu",     "SSU identifier column:", col_choices(),
                                                  selected = if (is_demo() && "id_ind" %in% col_choices()) "id_ind") })
  output$ui_strata_var <- renderUI({ selectInput("strata_var", "Stratum column:",        col_choices(),
                                                  selected = if (is_demo() && "stratum" %in% col_choices()) "stratum") })
  output$ui_deff_var   <- renderUI({
    cols <- col_choices()
    # Default deff_var to the stratum column (categorical) — same as the vignette (deff_var = "stratum")
    sel  <- if (!is.null(input$strata_var) && input$strata_var %in% cols) input$strata_var else cols[1]
    tagList(
      selectInput("deff_var", "DEFF grouping variable:", cols, selected = sel),
      tags$small(class = "text-muted",
        icon("info-circle"), " Must be a ", tags$b("categorical"), " column (typically the stratum",
        " or domain variable). Used to compute design effects by group."
      )
    )
  })

  # ── Multi-level domain variable selector ─────────────────────────────────────
  output$ui_domain_vars <- renderUI({
    cols <- col_choices()
    tagList(
      tags$label("Domain variable(s):"),
      numericInput("n_domain_vars", "Number of domain variables:", 1, min = 1, max = 5),
      uiOutput("ui_domain_rows"),
      tags$small(class = "text-muted",
        "Level 1 = primary domain (e.g. Region), Level 2 = finer (e.g. Province), …",
        tags$br(),
        tags$b("Note:"), " all selected variables are used in computation.",
        " The finest level (last) drives strata statistics;",
        " coarser levels are attached via the nested merge in ",
        tags$code("input_to_beat.2st_1()"), ".",
        tags$br(),
        tags$b("CV table tip:"), " set ", tags$b("Domain levels = N\u2009+\u20091"),
        " where N = number of domain variables above."
      )
    )
  })

  output$ui_domain_rows <- renderUI({
    req(input$n_domain_vars)
    n    <- as.integer(input$n_domain_vars)
    cols <- col_choices()
    rows <- lapply(seq_len(n), function(i) {
      fluidRow(style = "margin-bottom:0px;",
        column(3,
          div(style = "margin-top:6px; text-align:center;",
            tags$span(class = "badge badge-primary",
                      style = "font-size:12px; padding:4px 7px;", paste0("L", i))
          )
        ),
        column(9,
          selectInput(paste0("domain_col_", i), label = NULL,
                      choices  = cols,
                      selected = if (is_demo() && i == 1 && "region" %in% cols) "region"
                                 else if (length(cols) >= i) cols[i] else cols[1],
                      width    = "100%")
        )
      )
    })
    tagList(rows)
  })

  collect_domain_vars <- function() {
    n <- isolate(input$n_domain_vars)
    n <- if (is.null(n) || is.na(n)) 1L else as.integer(n)
    vars <- character(n)
    for (i in seq_len(n)) {
      v <- input[[paste0("domain_col_", i)]]
      vars[i] <- if (!is.null(v) && nchar(v) > 0) v else ""
    }
    vars <- vars[nchar(vars) > 0]
    if (length(vars) == 0) NULL else vars
  }

  output$ui_target_vars <- renderUI({
    checkboxGroupInput("target_vars", "Target variables (select ≥ 1):",
                       choices = col_choices(),
                       selected = if (is_demo())
                         intersect(c("active", "inactive", "unemployed", "income_hh"),
                                   col_choices()),
                       inline = FALSE)
  })

  # ── CV table ─────────────────────────────────────────────────────────────────
  output$ui_cv_table <- renderUI({
    req(input$n_dom, input$target_vars, length(input$target_vars) > 0)
    n_dom  <- as.integer(input$n_dom)
    n_var  <- length(input$target_vars)
    vnames <- input$target_vars
    header <- tags$tr(
      tags$th("Domain"),
      lapply(seq_len(n_var), function(j) tags$th(vnames[j]))
    )
    rows <- lapply(seq_len(n_dom), function(i) {
      tags$tr(
        tags$td(tags$strong(paste0("DOM", i))),
        lapply(seq_len(n_var), function(j) {
          tags$td(numericInput(paste0("cv_", i, "_", j), NULL,
                               value = 0.05, min = 1e-4, max = 1,
                               step = 0.005, width = "72px"))
        })
      )
    })
    tags$table(class = "table table-condensed table-bordered",
               tags$thead(header), tags$tbody(rows))
  })

  observeEvent(input$btn_set_equal, {
    req(input$n_dom, input$target_vars)
    val <- input$cv_equal_val; req(!is.na(val))
    for (i in seq_len(as.integer(input$n_dom)))
      for (j in seq_len(length(input$target_vars)))
        updateNumericInput(session, paste0("cv_", i, "_", j), value = val)
  })

  # ── Domain levels consistency indicator ──────────────────────────────────────
  output$ui_ndom_warning <- renderUI({
    n_dom   <- input$n_dom
    n_dvars <- input$n_domain_vars
    if (is.null(n_dom) || is.na(n_dom)) return(NULL)
    n_dom   <- as.integer(n_dom)
    n_dvars <- if (is.null(n_dvars) || is.na(n_dvars)) 1L else as.integer(n_dvars)
    expected <- n_dvars + 1L
    if (n_dom == expected) {
      div(class = "alert alert-success",
          style = "padding:6px 10px; font-size:12px; margin-bottom:6px;",
          icon("check-circle"),
          tags$b(" Correct:"), " ", n_dvars,
          " domain variable(s) \u2192 ", expected,
          " domain row(s) (DOM1 national + ", n_dvars, " sub-national).")
    } else {
      div(class = "alert alert-warning",
          style = "padding:8px 12px; font-size:12px; margin-bottom:6px;",
          icon("exclamation-triangle"),
          tags$b(" Mismatch:"), " ", n_dvars,
          " domain variable(s) require ",
          tags$b(paste0(expected, " domain rows")),
          " (DOM1 national + one per variable),",
          " but Domain levels is set to ", tags$b(n_dom), ".",
          tags$br(),
          "Adjust Domain levels to ", tags$b(expected),
          " before clicking Confirm Parameters.")
    }
  })

  observeEvent(input$upload_cv_file, {
    req(input$upload_cv_file)
    tryCatch({
      cv_data <- read.csv(input$upload_cv_file$datapath, stringsAsFactors = FALSE)
      if (!"DOM" %in% names(cv_data))
        stop("The CSV must contain a 'DOM' column.")
      cv_cols <- grep("^CV[0-9]+$", names(cv_data))
      if (length(cv_cols) == 0)
        stop("The CSV must contain at least one CV column named CV1, CV2, …")
      n_rows <- nrow(cv_data); n_cvs <- length(cv_cols)
      rv$pending_cv <- cv_data
      updateNumericInput(session, "n_dom", value = n_rows)
      shinyjs::delay(700, {
        cv_d <- isolate(rv$pending_cv)
        if (!is.null(cv_d)) {
          cv_c <- grep("^CV[0-9]+$", names(cv_d))
          for (i in seq_len(nrow(cv_d)))
            for (j in seq_along(cv_c)) {
              val <- cv_d[i, cv_c[j]]
              if (!is.na(val))
                updateNumericInput(session, paste0("cv_", i, "_", j), value = val)
            }
          rv$pending_cv <- NULL
        }
      })
      showNotification(
        paste0("CV file loaded: ", n_rows, " domain row(s), ", n_cvs, " CV column(s)."),
        type = "message")
    }, error = function(e) {
      showNotification(paste("CV file error:", e$message), type = "error", duration = 8)
    })
  })

  collect_cv <- function() {
    n_dom <- as.integer(input$n_dom)
    n_var <- length(input$target_vars)
    if (is.null(n_dom) || n_var == 0) return(NULL)
    cv_mat <- matrix(0.05, nrow = n_dom, ncol = n_var)
    for (i in seq_len(n_dom))
      for (j in seq_len(n_var)) {
        v <- input[[paste0("cv_", i, "_", j)]]
        if (!is.null(v) && length(v) == 1 && !is.na(v)) cv_mat[i, j] <- as.numeric(v)
      }
    df <- data.frame(DOM = paste0("DOM", seq_len(n_dom)), stringsAsFactors = FALSE)
    for (j in seq_len(n_var)) df[[paste0("CV", j)]] <- cv_mat[, j]
    df
  }

  output$txt_cv_preview <- renderPrint({
    cv <- collect_cv()
    if (is.null(cv)) return(cat("(set target variables first)"))
    print(cv)
  })

  observeEvent(input$btn_set_params, {
    req(rv$des, input$target_vars, length(input$target_vars) > 0)
    cv  <- collect_cv()
    dvs <- collect_domain_vars()
    if (is.null(cv))  { showNotification("Select at least one target variable.", type = "warning"); return() }
    if (is.null(dvs)) { showNotification("Select at least one domain variable.",  type = "warning"); return() }

    # Non-blocking mismatch check
    expected_ndom <- length(dvs) + 1L
    if (nrow(cv) != expected_ndom)
      showNotification(
        tagList(icon("exclamation-triangle"),
                tags$b(" Domain levels mismatch: "),
                nrow(cv), " CV row(s) vs. ", expected_ndom,
                " expected. Adjust Domain levels and re-confirm if needed."),
        type = "warning", duration = 12)

    rv$cv          <- cv
    rv$domain_vars <- dvs
    rv$target_vars <- input$target_vars
    rv$id_psu      <- input$id_psu
    rv$id_ssu      <- input$id_ssu
    rv$strata_var  <- input$strata_var
    rv$deff_var    <- input$deff_var
    rv$minimum     <- input$minimum
    rv$delta       <- input$delta
    rv$params_set  <- TRUE

    cv_r <- paste0(
      "cv <- as.data.frame(list(\n",
      "  DOM = c(", paste(shQuote(cv$DOM), collapse = ", "), "),\n",
      paste(sapply(seq_len(ncol(cv) - 1), function(j)
        paste0("  CV", j, " = c(", paste(cv[[paste0("CV", j)]], collapse = ", "), ")")),
        collapse = ",\n"),
      "\n))"
    )
    append_script(paste0(
      "# 2. Parameters\n",
      "id_PSU      <- ", shQuote(input$id_psu),     "\n",
      "id_SSU      <- ", shQuote(input$id_ssu),     "\n",
      "strata_var  <- ", shQuote(input$strata_var), "\n",
      "target_vars <- c(", paste(shQuote(input$target_vars), collapse = ", "), ")\n",
      "deff_var    <- ", shQuote(input$deff_var),   "\n",
      "domain_var  <- c(", paste(shQuote(collect_domain_vars()), collapse = ", "), ")\n",
      "minimum     <- ", input$minimum,             "\n",
      "delta       <- ", input$delta,               "\n",
      cv_r, "\n"
    ))
    updateNumericInput(session, "alloc_minimum", value = input$minimum)
    showNotification("Parameters confirmed.", type = "message")
  })

  # ── STEP 3: PREPARE INPUT ───────────────────────────────────────────────────
  observeEvent(input$btn_run_prep, {
    req(rv$params_set, rv$des, rv$cal, rv$cv, rv$domain_vars, rv$target_vars,
        rv$id_psu, rv$id_ssu, rv$strata_var, rv$deff_var)
    tryCatch({
      withProgress(message = "Running prepareInputToAllocation2()…", value = 0.4, {
        log_txt <- capture.output({
          inp <- prepareInputToAllocation2(
            samp_frame  = rv$frame %||% {
              env_fr <- new.env(); load("pop.RData", envir = env_fr); env_fr$pop },
            RGdes       = rv$des,
            RGcal       = rv$cal,
            id_PSU      = rv$id_psu,
            id_SSU      = rv$id_ssu,
            strata_var  = rv$strata_var,
            target_vars = rv$target_vars,
            deff_var    = rv$deff_var,
            domain_var  = rv$domain_vars,
            delta       = rv$delta,
            minimum     = rv$minimum
          )
        })
      })

      # check_input: reconcile population counts
      withProgress(message = "Running check_input()…", value = 0.8, {
        log_txt2 <- capture.output({
          newstrata <- check_input(
            strata           = inp$strata,
            des              = inp$des_file,
            strata_var_strata = "STRATUM",
            strata_var_des   = "STRATUM"
          )
        })
        inp$strata <- newstrata
      })

      rv$inp      <- inp
      rv$prep_log <- c(log_txt,
                       "",
                       "--- check_input() ---",
                       if (length(log_txt2) > 0) log_txt2 else "No discrepancies found.")

      append_script(paste0(
        "\n# 3. Prepare allocation input\n",
        "inp <- prepareInputToAllocation2(\n",
        "  samp_frame=samp_frame, RGdes=RGdes, RGcal=RGcal,\n",
        "  id_PSU=", shQuote(rv$id_psu), ", id_SSU=", shQuote(rv$id_ssu), ",\n",
        "  strata_var=", shQuote(rv$strata_var), ",\n",
        "  target_vars=c(", paste(shQuote(rv$target_vars), collapse = ", "), "),\n",
        "  deff_var=", shQuote(rv$deff_var),
        ", domain_var=c(", paste(shQuote(rv$domain_vars), collapse = ", "), "),\n",
        "  delta=", rv$delta, ", minimum=", rv$minimum, "\n)\n",
        "inp$strata <- check_input(strata=inp$strata, des=inp$des_file,\n",
        "  strata_var_strata='STRATUM', strata_var_des='STRATUM')\n"
      ))
      showNotification("Preparation complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Preparation error:", e$message), type = "error", duration = 10)
    })
  })

  output$txt_prep_log   <- renderText({ req(rv$prep_log); paste(rv$prep_log, collapse = "\n") })
  output$tbl_inp_strata <- renderDT({ req(rv$inp); datatable(head(rv$inp$strata, 30), options = list(scrollX = TRUE, pageLength = 8)) })
  output$tbl_inp_des    <- renderDT({ req(rv$inp); datatable(rv$inp$des_file,        options = list(scrollX = TRUE, pageLength = 8)) })
  output$tbl_inp_psu    <- renderDT({ req(rv$inp); datatable(head(rv$inp$psu_file, 30), options = list(scrollX = TRUE, pageLength = 8)) })
  output$txt_inp_deff   <- renderPrint({ req(rv$inp); print(head(rv$inp$deff,  10)) })
  output$txt_inp_effst  <- renderPrint({ req(rv$inp); print(head(rv$inp$effst, 10)) })
  output$txt_inp_rho    <- renderPrint({ req(rv$inp); print(head(rv$inp$rho,   10)) })

  # ── STEP 4: ALLOCATION ──────────────────────────────────────────────────────
  observeEvent(input$btn_run_alloc, {
    req(rv$inp, rv$cv)
    tryCatch({
      inp <- rv$inp
      inp$des_file$MINIMUM <- input$alloc_minimum
      withProgress(message = "Running beat.2st()…", value = 0.5, {
        alloc <- beat.2st(
          stratif     = inp$strata,
          errors      = rv$cv,
          des_file    = inp$des_file,
          psu_file    = inp$psu_file,
          rho         = inp$rho,
          deft_start  = NULL,
          effst       = inp$effst,
          minPSUstrat = input$alloc_minPSU,
          minnumstrat = input$alloc_minnumstrat
        )
      })
      rv$alloc <- alloc
      append_script(paste0(
        "\n# 4. Two-stage allocation\n",
        "inp$des_file$MINIMUM <- ", input$alloc_minimum, "\n",
        "alloc <- beat.2st(\n",
        "  stratif=inp$strata, errors=cv,\n",
        "  des_file=inp$des_file, psu_file=inp$psu_file,\n",
        "  rho=inp$rho, deft_start=NULL, effst=inp$effst,\n",
        "  minPSUstrat=", input$alloc_minPSU,
        ", minnumstrat=", input$alloc_minnumstrat, "\n)\n"
      ))
      showNotification("Allocation complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Allocation error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_alloc_total <- renderText({
    req(rv$alloc)
    iters <- rv$alloc$iterations
    psu_col <- "PSU Total"
    if (!is.null(iters) && is.data.frame(iters) && nrow(iters) > 0 && psu_col %in% names(iters))
      format(iters[[psu_col]][nrow(iters)], big.mark = ",")
    else "N/A"
  })
  output$txt_alloc_ssu <- renderText({
    req(rv$alloc)
    iters <- rv$alloc$iterations
    if (!is.null(iters) && is.data.frame(iters) && nrow(iters) > 0 && "SSU" %in% names(iters))
      format(iters$SSU[nrow(iters)], big.mark = ",")
    else "N/A"
  })
  output$tbl_alloc     <- renderDT({ req(rv$alloc); datatable(rv$alloc$alloc,        options = list(scrollX = TRUE, pageLength = 10)) })
  output$tbl_expected  <- renderDT({ req(rv$alloc); datatable(rv$alloc$expected,     options = list(dom = "t", scrollX = TRUE)) })
  output$tbl_sensitivity <- renderDT({ req(rv$alloc); datatable(rv$alloc$sensitivity, options = list(dom = "t", scrollX = TRUE)) })

  # ── STEP 5: SELECT SAMPLE ───────────────────────────────────────────────────
  # SSU code for select_SSU (may differ from id_ssu used in prepareInput)
  output$ui_ssu_code <- renderUI({
    req(rv$frame)
    selectInput("ssu_select_code", "SSU selection variable (frame):",
                choices  = names(rv$frame),
                # Demo: la selezione fisica delle SSU usa id_ind (vedi esempio R2BEAT
                # select_SSU); le statistiche/allocazione usano id_hh (passo 2).
                selected = if (identical(input$frame_source, "demo") &&
                               "id_ind" %in% names(rv$frame)) "id_ind"
                           else rv$id_ssu %||% names(rv$frame)[1])
  })

  observeEvent(input$btn_select_psu, {
    req(rv$alloc)
    tryCatch({
      set.seed(input$seed)
      withProgress(message = "Selecting PSUs…", value = 0.5, {
        s1 <- select_PSU(rv$alloc, type = input$psu_type,
                         pps = isTRUE(input$pps), plot = FALSE)
      })
      rv$sample_1st <- s1
      append_script(paste0(
        "\n# 5a. Select PSUs\n",
        "set.seed(", input$seed, ")\n",
        "sample_1st <- select_PSU(alloc, type=", shQuote(input$psu_type),
        ", pps=", isTRUE(input$pps), ", plot=FALSE)\n"
      ))
      showNotification("PSU selection complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("PSU selection error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_psu_stats <- renderPrint({ req(rv$sample_1st); print(rv$sample_1st$PSU_stats) })

  observeEvent(input$btn_select_ssu, {
    req(rv$sample_1st, rv$frame, rv$id_psu, input$ssu_select_code)
    tryCatch({
      withProgress(message = "Selecting SSUs…", value = 0.5, {
        samp <- select_SSU(
          df          = rv$frame,
          PSU_code    = rv$id_psu,
          SSU_code    = input$ssu_select_code,
          PSU_sampled = rv$sample_1st$sample_PSU
        )
      })
      rv$samp <- samp
      append_script(paste0(
        "\n# 5b. Select SSUs\n",
        "samp <- select_SSU(\n",
        "  df=samp_frame,\n",
        "  PSU_code=", shQuote(rv$id_psu), ",\n",
        "  SSU_code=", shQuote(input$ssu_select_code), ",\n",
        "  PSU_sampled=sample_1st$sample_PSU\n)\n"
      ))
      showNotification("SSU selection complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("SSU selection error:", e$message), type = "error", duration = 8)
    })
  })

  output$txt_samp_size <- renderText({
    req(rv$samp)
    paste(format(nrow(rv$samp), big.mark = ","), "SSUs selected")
  })
  output$txt_pop_n <- renderText({
    req(rv$samp)
    paste(format(round(sum(rv$samp$weight), 0), big.mark = ","))
  })
  output$tbl_samp_head <- renderDT({
    req(rv$samp)
    datatable(head(rv$samp, 100), options = list(scrollX = TRUE, pageLength = 8))
  })
  output$btn_download_samp <- downloadHandler(
    filename = function() paste0("sample_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(rv$samp); write.csv(rv$samp, file, row.names = FALSE)
    }
  )

  # Weight distribution plots
  domain_var_primary <- reactive({
    dv <- rv$domain_vars
    if (is.null(dv) || length(dv) == 0) NULL else dv[1]
  })

  output$ui_wt_domain_title  <- renderUI({
    dv <- domain_var_primary()
    paste("Weight distribution — by", if (is.null(dv)) "domain" else dv)
  })
  output$ui_wt_stratum_title <- renderUI({
    req(rv$strata_var)
    paste("Weight distribution — by", rv$strata_var)
  })

  output$plot_wt_total    <- renderPlot({ req(rv$samp); wf2_weight_boxplot(rv$samp, NULL,                  "Weights — total sample") })
  output$plot_wt_domain   <- renderPlot({ req(rv$samp); wf2_weight_boxplot(rv$samp, domain_var_primary(),  paste("Weights by", domain_var_primary() %||% "domain")) })
  output$plot_wt_province <- renderPlot({ req(rv$samp); wf2_weight_boxplot(rv$samp, "province",            "Weights by province") })
  output$plot_wt_stratum  <- renderPlot({ req(rv$samp, rv$strata_var); wf2_weight_boxplot(rv$samp, rv$strata_var, paste("Weights by", rv$strata_var)) })

  # ── STEP 6: R SCRIPT ────────────────────────────────────────────────────────
  output$txt_script <- renderText({
    paste(rv$script_lines, collapse = "\n")
  })
  output$btn_dl_script <- downloadHandler(
    filename = function() paste0("R2BEAT_sample_workflow_", Sys.Date(), ".R"),
    content  = function(file) writeLines(rv$script_lines, file)
  )
}


# ============================ LAUNCHER ============================

# =============================================================================
# R2BEAT — Launcher (landing page)
# Apre il workflow scelto nella stessa app via query string (?app=optim|wf1|wf2).
# =============================================================================


# ─── Deployed sub-app URLs ────────────────────────────────────────────────────

app_urls <- c(
  "1" = "?app=optim",
  "2" = "?app=wf1",
  "3" = "?app=wf2"
)

# ─── CSS ──────────────────────────────────────────────────────────────────────

launcher_css <- HTML("
  body {
    background: #f4f6f9;
    font-family: 'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif;
  }
  .launcher-wrapper {
    max-width: 760px;
    margin: 55px auto 40px auto;
    padding: 0 20px;
  }
  .launcher-header { text-align: center; margin-bottom: 34px; }
  .launcher-header h1 {
    font-size: 2.7em; font-weight: 700; color: #2c3e50;
    margin-bottom: 4px; letter-spacing: -0.5px;
  }
  .launcher-header p { font-size: 1.08em; color: #666; }

  .option-card {
    display: block; width: 100%; background: #fff;
    border: 2px solid #dde3ea; border-radius: 10px;
    padding: 20px 22px; margin-bottom: 14px; cursor: pointer;
    transition: border-color .15s, box-shadow .15s, background .15s;
    position: relative; user-select: none;
  }
  .option-card:hover {
    border-color: #3c8dbc;
    box-shadow: 0 3px 12px rgba(60,141,188,.15);
  }
  .option-card.selected {
    border-color: #3c8dbc; background: #eaf4fb;
    box-shadow: 0 3px 14px rgba(60,141,188,.22);
  }
  .card-number {
    display: inline-block; width: 32px; height: 32px; line-height: 32px;
    text-align: center; border-radius: 50%; background: #3c8dbc;
    color: #fff; font-weight: 700; font-size: .95em;
    margin-right: 12px; vertical-align: middle;
  }
  .option-card.selected .card-number { background: #2171a0; }
  .card-title { font-size: 1.08em; font-weight: 700; color: #2c3e50; vertical-align: middle; }
  .card-desc  { margin-top: 9px; margin-left: 44px; font-size: .95em; color: #555; line-height: 1.5; }
  .card-badge {
    position: absolute; top: 13px; right: 16px;
    font-size: .76em; padding: 2px 9px; border-radius: 20px; font-weight: 600;
  }
  .badge-green  { background: #dff0d8; color: #3c763d; }
  .badge-blue   { background: #d9edf7; color: #31708f; }
  .badge-orange { background: #fcf8e3; color: #8a6d3b; }

  #btn_launch {
    margin-top: 26px; width: 100%; padding: 13px;
    font-size: 1.12em; font-weight: 700; border-radius: 8px;
  }
  .footer-note { text-align: center; color: #aaa; font-size: .83em; margin-top: 28px; }
")

# ─── UI ───────────────────────────────────────────────────────────────────────

launcher_ui <- fluidPage(
  useShinyjs(),
  tags$head(tags$style(launcher_css)),

  div(class = "launcher-wrapper",
    div(class = "launcher-header",
      h1("R2BEAT"),
      p("Select the workflow you want to perform, then click ", tags$b("Launch"), ".")
    ),

    div(id = "card_1", class = "option-card",
        onclick = "Shiny.setInputValue('card_clicked','1',{priority:'event'})",
      span(class = "card-badge badge-green", "Optimisation"),
      span(class = "card-number", "1"),
      span(class = "card-title", "Optimise precision constraints"),
      div(class = "card-desc",
        "Find the optimal set of CV constraints for a given set of target variables",
        " and domain levels, minimising total sample size while meeting user-defined requirements."
      )
    ),

    div(id = "card_2", class = "option-card",
        onclick = "Shiny.setInputValue('card_clicked','2',{priority:'event'})",
      span(class = "card-badge badge-blue", "Workflow 1"),
      span(class = "card-number", "2"),
      span(class = "card-title",
           "Two-stage design & selection \u2014 from the population frame"),
      div(class = "card-desc",
        "Design a two-stage stratified sample starting from a complete population frame.",
        " Covers sensitivity analysis, optimal two-stage allocation (beat.2st),",
        " and PSU / SSU selection."
      )
    ),

    div(id = "card_3", class = "option-card",
        onclick = "Shiny.setInputValue('card_clicked','3',{priority:'event'})",
      span(class = "card-badge badge-orange", "Workflow 2"),
      span(class = "card-number", "3"),
      span(class = "card-title",
           "Two-stage design & selection \u2014 from a previous survey round"),
      div(class = "card-desc",
        "Design a two-stage stratified sample using ReGenesees design and calibration",
        " objects from a previous survey round.",
        " Covers allocation, PSU and SSU selection with updated weights."
      )
    ),

    disabled(
      actionButton("btn_launch",
                   label = tagList(icon("rocket"), " Launch"),
                   class = "btn-primary btn-lg")
    ),

    div(class = "footer-note",
      "R2BEAT \u2014 R to Bethel Allocation for Estimation in Two-stage sampling"
    )
  ),

  tags$script(HTML("
    function highlightCard(val) {
      document.querySelectorAll('.option-card').forEach(function(c) {
        c.classList.remove('selected');
      });
      var el = document.getElementById('card_' + val);
      if (el) el.classList.add('selected');
    }
  "))
)

# ─── Server ───────────────────────────────────────────────────────────────────

launcher_server <- function(input, output, session) {
  current_choice <- reactiveVal(NULL)

  observeEvent(input$card_clicked, {
    current_choice(input$card_clicked)
    enable("btn_launch")
    runjs(paste0("highlightCard('", input$card_clicked, "');"))
  })

  observeEvent(input$btn_launch, {
    req(current_choice())
    url <- app_urls[[current_choice()]]
    # Naviga nella stessa scheda; dal workflow si torna qui col pulsante "Menu principale".
    runjs(paste0("window.location.href='", url, "';"))
  })
}


# ============================ ROUTER ============================
# ui e' una funzione della richiesta: legge ?app=... e restituisce la UI giusta.
ui <- function(req) {
  q <- parseQueryString(req$QUERY_STRING)
  switch(q$app %||% "",
    optim = optim_ui,
    wf1   = wf1_ui,
    wf2   = wf2_ui,
    launcher_ui)
}

server <- function(input, output, session) {
  q   <- parseQueryString(isolate(session$clientData$url_search))
  app <- q$app %||% ""
  if (app == "optim")    optim_server(input, output, session)
  else if (app == "wf1") wf1_server(input, output, session)
  else if (app == "wf2") wf2_server(input, output, session)
  else                   launcher_server(input, output, session)
}

shinyApp(ui, server)
