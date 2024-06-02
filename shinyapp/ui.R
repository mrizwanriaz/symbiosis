library(pacman)
p_load(shiny,
       ggplot2,
       gplots,
       dplyr,
       formattable,
       tidyverse,
       magrittr,
       shinyjs)

################# Define UI for the Shiny app ################
ui <- fluidPage(
  useShinyjs(),
  
  tags$head(# Note the wrapping of the string in HTML()
    tags$style(
      HTML(
        "
                    .well{
                    background-color:#fcf9f0;
                    border-top:0px;
                    border-radius:0px;
                    border-color:#f0dfaf;
                    }
                      .nav-tabs>li>a {
                        color:black;
                      }
                      .nav-tabs>li.active>a, .nav-tabs>li.active>a:focus, .nav-tabs>li.active>a:hover {
                      font-weight:600;
                    color:black;
                    background-color:#fcf9f0;
                      }
                    .nav-tabs>li.active {
                        border-botom:0px;
                      }
      "
      )
    )),
  
  titlePanel(
    fluidRow(
      h2("Systems Genetics of Legume-Rhizobia Partnership Quality"),
      tags$h4("Dual-RNAseq Analysis of Nodule Transcriptome", tags$i("(Medicago truncatula & Sinorhizobium meliloti)")),
      tags$a(
        href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212235",
        target = "_BLANK",
        "[GEO: GSE212235]",
        style = "display:inline;padding-left:10px;font-size:16px;color:#bb5566;"
      ),
      tags$a(
        href = "",
        target = "",
        "[GitHub]",
        style = "display:inline;padding-left:10px;font-size:16px;color:#bb5566;"
      ),
      tags$a(
        href = "",
        target = "",
        "[DOI]",
        style = "display:inline;padding-left:10px;font-size:16px;color:#bb5566;"
      )
      ,
      style = "padding-left:10px;background-color:#dddddd;padding-bottom:10px;"
    ),
    windowTitle = "Symbiosis"
  ),
  #tags$br(),
  
  tabsetPanel(
    type = "tabs",
    
    #single gene panel
    tabPanel("Single Gene",
             sidebarLayout(
               sidebarPanel(
                 radioButtons(
                   "sel_specie",
                   "Select Specie:",
                   c(
                     "Sinorhizobium meliloti" = "smel",
                     "Medicago truncatula" = "med"
                   )
                 ),
                 
                 
                 selectizeInput(
                   inputId = "gene_select",
                   label = "Select a Gene:",
                   multiple = FALSE,
                   choices = character(0)
                 )
                 ,
                 fluidRow(column(
                   12,
                   actionButton(
                     inputId = "gene_slc_btn",
                     label = "Plot",
                     icon("chart-simple"),
                     style = "color: #fff; background-color: #4A9759; border-color: ##4A9759;float:right;"
                   )
                 )),
                 width = 3
               ),
               
               mainPanel(fluidRow(column(
                 12, div(uiOutput(outputId = "gene_feature_info"), style = "margin-top:20px;")
               )),
               fluidRow(column(
                 12, div(plotOutput(outputId = "expression_plot"), style = "margin-top:60px;")
               )))
             )),
    
    #multiple genes panel
    tabPanel(
      "Multiple Genes",
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel(
          radioButtons(
            "sel_multi_specie",
            "Select Specie:",
            c(
              "Sinorhizobium meliloti" = "smel",
              "Medicago truncatula" = "med"
            )
          ),
          textAreaInput(
            inputId = "gene_list",
            label = "Enter Gene IDs (one per line):",
            value = "",
            rows = 5
          ),
          fluidRow(column(
            12,
            
            actionButton(
              inputId = "clear_txt",
              label = "Clear",
              icon("broom"),
              style = "color: #fff; background-color: #8E9590; border-color: #8E9590;"
            ),
            actionButton(
              inputId = "load_exmpl",
              label = "Load Example IDs",
              icon("pencil"),
              style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
            ),
            
            actionButton(
              inputId = "plot_btn",
              label = "Plot",
              icon("chart-simple"),
              style = "color: #fff;background-color: #4A9759; border-color: #4A9759;display:inline-block;float:right;"
            )
          )),
          width = 3
        ),
        
        mainPanel(
          fluidRow(
            tags$i(id = "highlow2",
              "Strains are sorted from ",
              tags$b("high", style = "color:#B155A5;"),
              "to ",
              tags$b("low", style = "color:#419460;"),
              " quality partners"
            ),
            style = "float:right;font-size:10px;padding:10px;", 
          ),
          fluidRow(column(
            12, plotOutput(outputId = "heatmap_plot"), style = "min-height:420px;"
          )),
          fluidRow(column(
            12,
            downloadButton(outputId = 'download_multi', "Export Table as CSV", style =
                             "margin-top:40px;margin-bottom:10px;color:#fff;background-color:#997700;float:right;"),
            div(dataTableOutput('multi_dto'), style =
                  "margin-top:10px;"),
            style = "padding:5px"
          )),
        )
      )
    ),
    
    #module-wise
    tabPanel(
      "WGCNA Modules",
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel(
          radioButtons(
            "sel_mod_specie",
            "Select Specie:",
            c(
              "Sinorhizobium meliloti" = "smel",
              "Medicago truncatula" = "med"
            )
          ),
          selectInput(
            inputId = "mod_select",
            label = "Select a Module:",
            multiple = FALSE,
            choices = paste("M", c(1:40), sep = "")
          )
          ,
          fluidRow(id = "mm_cutoff_row", column(
            12,
            sliderInput(
              "mm_cutoff",
              "Module Membership:",
              min = -1,
              max = 1,
              value = c(-1, 1),
              step = 0.2
            )
          ), style = "display:none;"),
          fluidRow(id = "gs_cutoff_row", column(
            12,
            sliderInput(
              "gs_cutoff",
              "Gene Significance:",
              min = -1,
              max = 1,
              value = c(-1, 1),
              step = 0.2
            )
          ), style = "display:none;"),
          fluidRow(column(
            6,
            actionButton(
              inputId = "mod_filt_btn",
              label = "Show/hide parameters",
              icon("gear"),
              style = "color: #fff; background-color: grey; border-color: grey; float:left;"
            )
          ),
          column(
            6,
            actionButton(
              inputId = "mod_slc_btn",
              label = "Plot",
              icon("chart-simple"),
              style = "color: #fff; background-color: #4A9759; border-color: ##4A9759; float:right;"
            )
          )),
          width = 3
        ),
        
        mainPanel(
          fluidRow(
            tags$i(id = "highlow",
              "Strains are sorted from ",
              tags$b("high", style = "color:#B155A5;"),
              "to ",
              tags$b("low", style = "color:#419460;"),
              " quality partners"
            ),
            style = "float:right;font-size:10px;padding:10px;"
          ),
          
          fluidRow(column(
            12, plotOutput(outputId = "mod_heatmap_plot"), style = "min-height:420px;"
          )),
          fluidRow(column(
            12,
            downloadButton(outputId = 'download_mod', "Export Table as CSV", style =
                             "margin-top:40px;margin-bottom:10px;color:#fff;background-color:#997700;float:right;"),
            div(dataTableOutput('mod_dto'), style =
                  "margin-top:10px;"),
            style = "padding:5px"
          )),
        )
      )
    ),
  )
)
