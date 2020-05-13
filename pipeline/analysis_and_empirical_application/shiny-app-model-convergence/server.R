library(shiny)
library(shinyjs)
library(sleepsimReval)
library(dplyr)
# Define server logic for random distribution app ----
server <- function(input, output) {
  # Observe event: when user clicks converge / not converged
  observeEvent(input$models_converged, {
    # Process the current dataset
    history[[f]] <<- list(
      "iteration_id" = gsub("\\.rds", "", f),
      "converged" = 1,
      "n" = n,
      "n_t" = n_t,
      "zeta" = zeta,
      "Q" = Q,
      "emiss_1_mu_1" = input$conv_em1_mu1,
      "emiss_1_mu_2" = input$conv_em1_mu2,
      "emiss_1_mu_3" = input$conv_em1_mu3,
      "emiss_2_mu_1" = input$conv_em2_mu1,
      "emiss_2_mu_2" = input$conv_em2_mu2,
      "emiss_2_mu_3" = input$conv_em2_mu3,
      "emiss_3_mu_1" = input$conv_em3_mu1,
      "emiss_3_mu_2" = input$conv_em3_mu2,
      "emiss_3_mu_3" = input$conv_em3_mu3,
      "emiss_1_varmu_1" = input$conv_em1_varmu1,
      "emiss_1_varmu_2" = input$conv_em1_varmu2,
      "emiss_1_varmu_3" = input$conv_em1_varmu3,
      "emiss_2_varmu_1" = input$conv_em2_varmu1,
      "emiss_2_varmu_2" = input$conv_em2_varmu2,
      "emiss_2_varmu_3" = input$conv_em2_varmu3,
      "emiss_3_varmu_1" = input$conv_em3_varmu1,
      "emiss_3_varmu_2" = input$conv_em3_varmu2,
      "emiss_3_varmu_3" = input$conv_em3_varmu3,
      "gamma_int_bar_S1toS2" = input$gamma_int_bar12,
      "gamma_int_bar_S1toS3" = input$gamma_int_bar13,
      "gamma_int_bar_S2toS2" = input$gamma_int_bar22,
      "gamma_int_bar_S2toS3" = input$gamma_int_bar23,
      "gamma_int_bar_S3toS2" = input$gamma_int_bar32,
      "gamma_int_bar_S3toS3" = input$gamma_int_bar33
    )
    # Save
    saveRDS(history, "history.rds")
    # Load new datasets
    needs_processing <<- needs_processing[-1]
    print(length(needs_processing))
    # File name
    f <<- needs_processing[1]
    # Load
    mods <<- read_models(f)
    # Get model details
    scen_cur <- scen %>%
      filter(iteration_id == gsub("\\.rds", "", f))
    # Parameter settings
    n <- scen_cur$n
    n_t <- scen_cur$n_t
    zeta <- scen_cur$zeta
    Q <- scen_cur$Q
    # Set html output
    output$model_id <- renderUI({HTML(paste0("<p>Iteration ID: ", scen_cur$iteration_id, "</p>"))})
    output$model_N <- renderUI({HTML(paste0("<p>N: ", n, "</p>"))})
    output$model_NT <- renderUI({HTML(paste0("<p>N_t: ", n_t, "</p>"))})
    output$model_zeta <- renderUI({HTML(paste0("<p>Zeta: ", zeta, "</p>"))})
    output$model_Q <- renderUI({HTML(paste0("<p>Q: ", Q, "</p>"))})
    # Reset selectinput
    reset("conv_em1_mu1")
    reset("conv_em1_mu2")
    reset("conv_em1_mu3")
    reset("conv_em2_mu1")
    reset("conv_em2_mu2")
    reset("conv_em2_mu3")
    reset("conv_em3_mu1")
    reset("conv_em3_mu2")
    reset("conv_em3_mu3")
    # Varmu
    reset("conv_em1_varmu1")
    reset("conv_em1_varmu2")
    reset("conv_em1_varmu3")
    reset("conv_em2_varmu1")
    reset("conv_em2_varmu2")
    reset("conv_em2_varmu3")
    reset("conv_em3_varmu1")
    reset("conv_em3_varmu2")
    reset("conv_em3_varmu3")
    # gamm int bar
    reset("gamma_int_bar12")
    reset("gamma_int_bar13")
    reset("gamma_int_bar22")
    reset("gamma_int_bar23")
    reset("gamma_int_bar32")
    reset("gamma_int_bar33")
  })
  # When not converged
  observeEvent(input$models_not_converged, {
    # Process the current dataset
    history[[f]] <<- list(
      "iteration_id" = gsub("\\.rds", "", f),
      "converged" = 0,
      "n" = n,
      "n_t" = n_t,
      "zeta" = zeta,
      "Q" = Q,
      "emiss_1_mu_1" = input$conv_em1_mu1,
      "emiss_1_mu_2" = input$conv_em1_mu2,
      "emiss_1_mu_3" = input$conv_em1_mu3,
      "emiss_2_mu_1" = input$conv_em2_mu1,
      "emiss_2_mu_2" = input$conv_em2_mu2,
      "emiss_2_mu_3" = input$conv_em2_mu3,
      "emiss_3_mu_1" = input$conv_em3_mu1,
      "emiss_3_mu_2" = input$conv_em3_mu2,
      "emiss_3_mu_3" = input$conv_em3_mu3,
      "emiss_1_varmu_1" = input$conv_em1_varmu1,
      "emiss_1_varmu_2" = input$conv_em1_varmu2,
      "emiss_1_varmu_3" = input$conv_em1_varmu3,
      "emiss_2_varmu_1" = input$conv_em2_varmu1,
      "emiss_2_varmu_2" = input$conv_em2_varmu2,
      "emiss_2_varmu_3" = input$conv_em2_varmu3,
      "emiss_3_varmu_1" = input$conv_em3_varmu1,
      "emiss_3_varmu_2" = input$conv_em3_varmu2,
      "emiss_3_varmu_3" = input$conv_em3_varmu3,
      "gamma_int_bar_S1toS2" = input$gamma_int_bar12,
      "gamma_int_bar_S1toS3" = input$gamma_int_bar13,
      "gamma_int_bar_S2toS2" = input$gamma_int_bar22,
      "gamma_int_bar_S2toS3" = input$gamma_int_bar23,
      "gamma_int_bar_S3toS2" = input$gamma_int_bar32,
      "gamma_int_bar_S3toS3" = input$gamma_int_bar33
    )
    # Save
    saveRDS(history, "history.rds")
    # Load new datasets
    needs_processing <<- needs_processing[-1]
    print(length(needs_processing))
    # File name
    f <<- needs_processing[1]
    # Load
    mods <<- read_models(f)
    # Get model details
    scen_cur <- scen %>%
      filter(iteration_id == gsub("\\.rds", "", f))
    # Parameter settings
    n <- scen_cur$n
    n_t <- scen_cur$n_t
    zeta <- scen_cur$zeta
    Q <- scen_cur$Q
    # Set html output
    output$model_id <- renderUI({HTML(paste0("<p>Iteration ID: ", scen_cur$iteration_id, "</p>"))})
    output$model_N <- renderUI({HTML(paste0("<p>N: ", n, "</p>"))})
    output$model_NT <- renderUI({HTML(paste0("<p>N_t: ", n_t, "</p>"))})
    output$model_zeta <- renderUI({HTML(paste0("<p>Zeta: ", zeta, "</p>"))})
    output$model_Q <- renderUI({HTML(paste0("<p>Q: ", Q, "</p>"))})
    # Reset selectinput
    reset("conv_em1_mu1")
    reset("conv_em1_mu2")
    reset("conv_em1_mu3")
    reset("conv_em2_mu1")
    reset("conv_em2_mu2")
    reset("conv_em2_mu3")
    reset("conv_em3_mu1")
    reset("conv_em3_mu2")
    reset("conv_em3_mu3")
    # Varmu
    reset("conv_em1_varmu1")
    reset("conv_em1_varmu2")
    reset("conv_em1_varmu3")
    reset("conv_em2_varmu1")
    reset("conv_em2_varmu2")
    reset("conv_em2_varmu3")
    reset("conv_em3_varmu1")
    reset("conv_em3_varmu2")
    reset("conv_em3_varmu3")
    # gamm int bar
    reset("gamma_int_bar12")
    reset("gamma_int_bar13")
    reset("gamma_int_bar22")
    reset("gamma_int_bar23")
    reset("gamma_int_bar32")
    reset("gamma_int_bar33")
  })
  # Create plots
  output$trace_plot_mu <- renderPlot({
    param <- input$parameter_input
    # If param is var1, var2, var3
    if(param %in% c("var1", "var2", "var3")) {
      param_in <- "emiss_mu_bar"
      var <- as.numeric(gsub("var", "", param))
    } else {
      param_in <- "gamma_int_bar"
      var <- 1
    }
    sleepsimReval::tpp(mods[[1]], mods[[2]],
                       param = param_in,
                       var=var)
  })
  output$trace_plot_varmu <- renderPlot({
    param <- input$parameter_input
    # If param is var1, var2, var3
    if(param %in% c("var1", "var2", "var3")) {
      param_in <- "emiss_varmu_bar"
      var <- as.numeric(gsub("var", "", param))
    } else {
      param_in <- "gamma_int_bar"
      var <- 1
      return(NULL)
    }
    sleepsimReval::tpp(mods[[1]], mods[[2]],
                       param = param_in,
                       var=var)
  })
  # GR statistic
  output$gelman_rubin <- renderTable({
    param <- input$parameter_input
    # If param is var1, var2, var3
    if(param %in% c("var1", "var2", "var3")) {
      param_in <- "emiss_mu_bar"
      var <- as.numeric(gsub("var", "", param))
      m <- sleepsimReval::compute_grs(mods[[1]], mods[[2]],
                                      param = param_in,
                                      var=var)
      m2 <- sleepsimReval::compute_grs(mods[[1]], mods[[2]],
                                       param = "emiss_varmu_bar",
                                       var=var)
      # Bind
      od <- as.data.frame(rbind(m$psrf, m2$psrf))
      od$var <- c(row.names(m$psrf), row.names(m2$psrf))
      od
    } else {
      param_in <- "gamma_int_bar"
      var <- 1
      m <- sleepsimReval::compute_grs(mods[[1]], mods[[2]],
                                      param = param_in,
                                      var=var)
      od <- as.data.frame(m$psrf)
      od$var <- row.names(m$psrf)
      od
    }
  })
  # Density plot
  output$dens_plot <- renderPlot({
    param <- input$parameter_input
    # If param is var1, var2, var3
    if(param %in% c("var1", "var2", "var3")) {
      param_in <- "emiss_mu_bar"
      var <- as.numeric(gsub("var", "", param))
    } else {
      param_in <- "gamma_int_bar"
      var <- 1
    }
    sleepsimReval::dens_plot(mods[[1]], mods[[2]],
                       param = param_in,
                       var=var)    
  })
  output$dens2_plot <- renderPlot({
    param <- input$parameter_input
    # If param is var1, var2, var3
    if(param %in% c("var1", "var2", "var3")) {
      param_in <- "emiss_varmu_bar"
      var <- as.numeric(gsub("var", "", param))
      sleepsimReval::dens_plot(mods[[1]], mods[[2]],
                               param = param_in,
                               var=var) 
    } else {
      NULL
    }
  })
}