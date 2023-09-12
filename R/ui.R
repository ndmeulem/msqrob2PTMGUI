#' Shiny app server object
#' @import shiny shinyhelper shinyWidgets DT plotly shinybusy BiocManager
# create the shiny application user interface

library(shiny)
library(shinyhelper)
library(DT)
library(plotly)
library(shinybusy)
library(shinyWidgets)


ui <- function() {
    shiny::addResourcePath("msqrob2PTMGUI", system.file("helpfiles", package="msqrob2PTMGUI"))
    shinyjs::useShinyjs()


    (navbarPage(

    # Application title
    title = "msqrob2PTM",

    # Add the tabs you want present in the app
    #First tab: upload necessary data files
    tabPanel(title = "Data input",

         add_busy_spinner(spin = "flower",
                          position = "full-page"),

         #Example data button
         fluidRow(
             column(width=4,
                    helper(shinyWidgets::actionBttn(inputId = "example_data",
                                 label = "Read in example data",
                                 color = "royal",
                                 style="bordered",
                                 size = "sm"),
                    type = "markdown", content = "example_data"))
         ),

         #Text output when example data has been read in
         fluidRow(
             column(width=4,
                    textOutput("read_in_example_data"))
         ),

         br(),

         # Show two input data columns

         # Input intensity file
         fluidRow(
             column(width = 4,
                helper(fileInput(inputId = "data", label = "Upload intensity file"),
                    type = "markdown", content = "peptidesfile")),

         #Input protein file
              column(width = 4,
                helper(fileInput(inputId = "proteindata", label = "Upload non-enriched intensity file"),
                       type = "markdown", content = "proteinfile")),

         # Input metadata file
            column(width = 4,
                helper(fileInput(inputId = "metadata", label = "Upload metadata file"),
                       type = "markdown", content = "metadatafile"))
         ),

        #Now add another row for input reading options
        fluidRow(
            #Options for the peptides intensity file
            column(width = 4,
                   tags$div(tags$h4("Options")),
                   numericInput(inputId = "skip", label = "Number of lines to skip", value=0, min=0),
                   radioButtons(inputId = "separator", label = "File separator",
                                choices = c("comma" = ",",
                                            "space" = " ",
                                            "semicolon" = ";",
                                            "tab" = "\t")),
                   textInput(inputId = "intensityIdentifier", label = "Intensity columns identifier"),
                   numericInput(inputId = "proteinColumn", label = "Protein column", value=1, min=1),
                   numericInput(inputId = "sequenceColumn", label = "Sequence column", value=2, min=1),
                   numericInput(inputId = "modificationsColumn", label = "Modifications column", value=3, min=1),
                   textInput(inputId = "modificationSplit", label = "Modifications separator", value = "|")
                   ),
            #Options for the protein intensity file
            column(width = 4,
                   tags$div(tags$h4("Options")),
                   numericInput(inputId = "proteinskip", label = "Number of lines to skip", value=0, min=0),
                   radioButtons(inputId = "proteinseparator", label = "File separator",
                                choices = c("comma" = ",",
                                            "space" = " ",
                                            "semicolon" = ";",
                                            "tab" = "\t")),
                   textInput(inputId = "proteinintensityIdentifier", label = "Intensity columns identifier"),
                   numericInput(inputId = "proteinproteinColumn", label = "Protein column", value=1, min=1)
            ),
            #Option for the metadata file
            column(width = 4,
                   tags$div(tags$h4("Options")),
                   checkboxInput("SDRF", label = tags$strong("Metadata in SDRF format"),
                                 value = F),
                   radioButtons(inputId = "separatorMetadata", label = "file separator",
                                choices = c("comma" = ",",
                                            "space" = " ",
                                            "semicolon" = ";",
                                            "tab" = "\t"))
                   ),
                  actionButton(inputId="generateMetaData",
                               label="Generate metadata (sdrf) file",
                               style="background-color: transparent; color: #5f9c68; border: 2px solid #5f9c68"),
                  htmlOutput("downloadButtonDownloadAnnot"),
        ),
        #Add action button for reading in the data
        fluidRow(
            column(width=4,offset=8,
                   shinyWidgets::actionBttn(inputId = "go",
                                label = "Read in your own data",
                                style="bordered",
                                color="royal",
                                size="sm"))
        ),
        #Text output for when data has been read in
        fluidRow(
            column(width=4,offset=8,
                   textOutput("read_in_data"))
        )),

    #Second tab: preprocessing: options + visualisations
    tabPanel(title = "Preprocessing",

         add_busy_spinner(spin = "fading-circle"),
         #Add row for preprocessing options
         fluidRow(
             column(width = 3,
                helper(tags$div(tags$h4("Preprocessing options")),
                       type = "markdown", content = "preprocessingoptions"),

                    radioButtons(inputId = "logTransform",
                      label = "Logarithmic transformation",
                      choiceNames = c("none","log2", "log10", "natural"),
                      choiceValues = c("none", 2, 10, exp(1))),

                    radioButtons(inputId = "summarisation",
                      label = "Summarise to PTMs",
                      choices = c("yes", "no"))
                ),
             column(width = 3,
                    numericInput(inputId = "nnonzero",
                                 label = "Minimum number of nonzero columns", value=2, min=0),
                    radioButtons(inputId = "normalisationMethodGlobal",
                                 label = "Normalisation",
                                 choiceNames = c("none","mean centering (DA)", "median centering (DA)", "against protein (DU)"),
                                 choiceValues = c("none", "center.mean", "center.median", "robust")
                                 )
                ),
             column(width = 3,
                    tags$div(tags$h4("Stats before preprocessing")),
                    tableOutput("statstable")),
             column(width = 3,
                    tags$div(tags$h4("Stats after preprocessing")),
                    tableOutput("statstableafter")
                    )
            ),
         #Add row for preprocessing button
         fluidRow(
             column(width=4,offset=4,
                    actionButton(inputId = "preprocess",
                                 label = "Preprocess!",
                                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
         ),
         #Add row for before plots
         fluidRow(
             column(width = 6,
                    tags$div(tags$h4("Density plot before preprocessing")),
                    plotlyOutput("densityBefore")
                    ),
             column(width = 6,
                    tags$div(tags$h4("Density plot after preprocessing")),
                    plotlyOutput("densityAfter")
                    )
            ),
         #Add row for after plots
         fluidRow(
             column(width = 6,
                    tags$div(tags$h4("Boxplot before preprocessing")),
                    plotOutput("boxplotBefore")
             ),
             column(width = 6,
                    tags$div(tags$h4("Boxplot after preprocessing")),
                    plotOutput("boxplotAfter")
             )
         )
         ),

    #Third tab: actual data visualisation
    tabPanel(title = "Data exploration",

        add_busy_spinner(spin = "fading-circle"),

        #Add row for protein dropdown menu, normalisation options and data table
        fluidRow(
            #Protein dropdown menu and normalisation options
            column(width = 3,
                tags$div(tags$h4("Protein dropdown menu")),
                selectInput(inputId = "protein", label = "select protein", choices = ""),
                helper(radioButtons(inputId = "normalisationMethod",
                                    label = "Usage options",
                             choiceNames = c("absolute abundance",
                                             "relative abundance (center mean)",
                                             "relative abundance (center median)"),
                             choiceValues = c("none","center.mean", "center.median")),
                       type = "markdown", content = "normalisationfile")
                ),
            #Data table
            column(width = 9,
                   tags$div(tags$h4("Corresponding protein data table")),
                   DTOutput("proteinDataTable"))
            ),
        fluidRow(column(width = 3, offset=2,
                        actionButton(inputId = "deselect",
                                     label = "Deselect all features",
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
        #Add row for the plots themselves
        fluidRow(
            column(width = 6,
                   tags$div(tags$h4("Lineplot of log-normalised features")),
                   plotlyOutput("lineplot")
                   ),
            column(width = 6,
                   tags$div(tags$h4("Boxplot of selected features")),
                   plotlyOutput("boxplot")
                   )
            ),
        #Add row for the arrange x-axis dropdown menu
        fluidRow(
            column(width = 3,
                   tags$div(tags$h4("Arrange x-axis based on")),
                   selectInput(inputId = "x_axis", label = "select variable", choices = ""))
        )
        ),

    #Fourth tab: model building with formula
    tabPanel(title = "Model",

        add_busy_spinner(spin = "fading-circle"),

        #Add row to build formula and visualise colData
        fluidRow(
          column(width = 5,
                 helper(tags$div(tags$h4("Build model formula")),
                        type = "markdown", content = "build_model_formula"),
                 #list available variables
                 tags$h5("Following variables can be selected to build the model: "),
                  tags$div(tags$em(textOutput("available_modelvariables"))),
                 br(),
                 textInput("designformula", label = "Design formula",
                           placeholder = "~ var1 + var2*var3"),
                 radioButtons(inputId = "robustRegression",
                              label = "Robust Regression",
                              choiceNames = c("yes", "no"),
                              choiceValues = c(TRUE, FALSE)),
                 actionButton("fitModel", "Fit Model",
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                 br(),
                 textOutput("model_fitted")

        ),
        column(width = 7,
               tags$div(tags$h4("Design variables")),
               dataTableOutput("designVariables")
               )

    ),

      #Add row to visualise design matrix
        fluidRow(
          column(width = 11,
          tags$div(tags$h4("Visualise design")),
          uiOutput("designmatrix", width = "100%")
        ))

    ),

    #Fifth tab: inference tab: make contrasts, plot volcano and table with
    #significant ptms/pepforms (option)
    tabPanel(title = "Inference",

        add_busy_spinner(spin = "fading-circle"),

        #Add row to specify contrasts, choose significance level and choose to
        #display ptms and/or pepforms
        fluidRow(
          #Specify contrasts
          column(width = 4,
                 helper(tags$div(tags$h4("Specify contrast")),
                        type = "markdown", content = "specify_contrast"),
                 #list available variables
                 tags$p("Following parameters can be used in contrasts for hypothesis tests: "),
                 tags$p(textOutput("available_parameters")),
                 ),
          column(width = 4,
                 textInput("contrast", label = "Null hypothesis",
                           placeholder = "parameter1 = 0"),
                 actionButton("testcontrast", "Test contrasts",
                              style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
          #Choose options
          column(width = 4,
                 numericInput("significancelevel", "Choose significance level",
                              value = 0.05, min = 0, max = 1, step = 0.01),
                 checkboxInput("onlysignificant", label = tags$strong("Only significant features in table"),
                               value = T),
                 )
          # column(width = 4,
          #        radioButtons("whichfeatures", "Choose which features to display",
          #                     choiceNames = c("PTMs",
          #                                     "Peptidoforms"),
          #                     choiceValues = c("ptmRel","peptideLogNorm")))

        ),

        #Add row for volcanoplot and significance table
        fluidRow(
          #volcanoplot
          column(width = 6,
                 tags$div(tags$h4("Volcano plot")),
                 plotlyOutput("volcano")
                 ),
          #Results table
          helper(column(width = 6,
                 tags$div(tags$h4("Results table")),
                 DTOutput("significanceTable")),
                 type = "markdown", content = "results_table")
        ),
        fluidRow(
          #Results table for selected peptidoforms
          column(width = 6, offset = 6,
                 DTOutput("significanceTablePepforms"))
        )

             ),

    #Sixth panel: detail plots of selected features
    tabPanel("Detail plot",

      add_busy_spinner(spin = "fading-circle"),
      #LinePlot


      fluidRow(
        column(width= 12,
          tags$div(tags$h4("Select a feature in the results table of the inference tab.
                           Please only select one feature at a time.")),
          plotlyOutput("lineplot_feature", height = "780"),
          ),
        ),

      #Add row for the arrange x-axis dropdown menu
      fluidRow(
        column(width = 3,
               tags$div(tags$h4("Arrange x-axis based on")),
               selectInput(inputId = "x_axis1", label = "select 1st variable", choices = "")),
        column(width = 3,
               tags$div(tags$h4("Arrange x-axis based on")),
               selectInput(inputId = "x_axis2", label = "select 2nd variable", choices = "")),
        column(width = 3,
               tags$div(tags$h4("Arrange x-axis based on")),
               selectInput(inputId = "x_axis3", label = "select 3rd variable", choices = ""))
      )
    ),
    tabPanel("Report",
        helper(fluidRow(
          column(width = 6,
            numericInput("maxPlot", "Number of significant features for which you want to have detail plots",
                       value = 10, min = 1, max = NA, step = 1))
          ), type = "markdown", content = "report"),
        fluidRow(
          column(width = 6,
                 actionButton("generate", "Generate Report"),
                 htmlOutput("downloadButtonDownloadReport"))
        )

      )
    )
)
}
