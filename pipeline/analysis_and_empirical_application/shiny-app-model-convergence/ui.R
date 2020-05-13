library(shiny)
library(shinyjs)
# UI logic goes here
ui <- fluidPage(
    useShinyjs(),
    # App title ----
    titlePanel("Model evaluation"),
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
        # Sidebar panel for inputs ----
        sidebarPanel(
            h3("Model details"),
            htmlOutput("model_id"),
            fluidRow(
                column(6,
                       htmlOutput("model_N"),
                       htmlOutput("model_NT")),
                column(6,
                       htmlOutput("model_zeta"),
                       htmlOutput("model_Q")),
            ),
            # Show the ID of the current model
            hr(),
            h3("Parameter selection"),
            # Select parameter
            selectInput("parameter_input", "Select parameter", c(
                "Emiss. dist. 1" = "var1",
                "Emiss. dist. 2" = "var2",
                "Emiss. dist. 3" = "var3",
                "Gamma int bar" = "gamma"
            )),
            hr(),
            h3("Select parameters that did not converge"),
            fluidRow(
                column(
                    3,
                    h5("Emission variable 1"),
                    selectInput("conv_em1_mu1", "State 1 (mean)", c("Yes", "No")),
                    selectInput("conv_em1_mu2", "State 2 (mean)", c("Yes", "No")),
                    selectInput("conv_em1_mu3", "State 3 (mean)", c("Yes", "No")),
                    selectInput("conv_em1_varmu1", "State 1 (var)", c("Yes", "No")),
                    selectInput("conv_em1_varmu2", "State 2 (var)", c("Yes", "No")),
                    selectInput("conv_em1_varmu3", "State 3 (var)", c("Yes", "No"))
                ),
                column(
                    3,
                    h5("emission variable 2"),
                    selectInput("conv_em2_mu1", "State 1 (mean)", c("Yes", "No")),
                    selectInput("conv_em2_mu2", "State 2 (mean)", c("Yes", "No")),
                    selectInput("conv_em2_mu3", "State 3 (mean)", c("Yes", "No")),
                    selectInput("conv_em2_varmu1", "State 1 (var)", c("Yes", "No")),
                    selectInput("conv_em2_varmu2", "State 2 (var)", c("Yes", "No")),
                    selectInput("conv_em2_varmu3", "State 3 (var)", c("Yes", "No"))
                ),
                column(
                    3,
                    h5("emission variable 3"),
                    selectInput("conv_em3_mu1", "State 1 (mean)", c("Yes", "No")),
                    selectInput("conv_em3_mu2", "State 2 (mean)", c("Yes", "No")),
                    selectInput("conv_em3_mu3", "State 3 (mean)", c("Yes", "No")),
                    selectInput("conv_em3_varmu1", "State 1 (var)", c("Yes", "No")),
                    selectInput("conv_em3_varmu2", "State 2 (var)", c("Yes", "No")),
                    selectInput("conv_em3_varmu3", "State 3 (var)", c("Yes", "No"))
                ),
                column(
                    3,
                    h5("Gamma intercepts"),
                    selectInput("gamma_int_bar12", "S1toS2", c("Yes", "No")),
                    selectInput("gamma_int_bar13", "S1toS3", c("Yes", "No")),
                    selectInput("gamma_int_bar22", "S2toS2", c("Yes", "No")),
                    selectInput("gamma_int_bar23", "S2toS3", c("Yes", "No")),
                    selectInput("gamma_int_bar32", "S3toS2", c("Yes", "No")),
                    selectInput("gamma_int_bar33", "S3toS3", c("Yes", "No"))
                )
            ),
            fluidRow(
                column(
                    6,
                    h4("Model convergence"),
                    # Button to click (converged/not converged)
                    actionButton("models_converged", "Converged"),
                    actionButton("models_not_converged", "Not converged")
                )
            ),
        ),
        # Main panel for displaying outputs ----
        mainPanel(
            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(type = "tabs",
                        tabPanel("Trace plots", plotOutput("trace_plot_mu"),
                                 plotOutput("trace_plot_varmu")),
                        tabPanel("Gelman-Rubin", tableOutput("gelman_rubin")),
                        tabPanel("Density plots", list(
                            plotOutput("dens_plot"),
                            plotOutput("dens2_plot")
                        )
                        )
            )
            
        )
    )
)