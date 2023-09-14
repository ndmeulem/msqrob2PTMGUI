#' Shiny app server function
#'
#' @param input provided by shiny
#' @param output provided by shiny
#' @param session provided by shiny
#' @import tidyverse shiny shinyhelper shinymeta rmarkdown knitr msqrob2 QFeatures limma plotly ggplot2 DT wesanderson BiocManager utils ExploreModelMatrix forcats MsCoreUtils

library(shiny)
library(shinyhelper)
library(shinymeta)
library(QFeatures)
library(tidyverse)
library(DT)
library(ggplot2)
library(plotly)
library(limma)
library(msqrob2)
library(dplyr)

server <- (function(input, output, session) {
    options(shiny.maxRequestSize=400*1024^2)
    #Activate the helper files
    shiny::addResourcePath("msqrob2PTMGUI", system.file("helpfiles", package="msqrob2PTMGUI"))
    shinyhelper::observe_helpers(help_dir = system.file("helpfiles", package="msqrob2PTMGUI"),
                    withMathJax = TRUE)

    #---------------------------
    #Data input
    #---------------------------

    #add variables to work with
    variables <- reactiveValues(pe = NULL)

    #When the user clicks read data or example data, it will trigger a number of events:
    #get ecols, read in files, get coldata, do filtering steps

    observeEvent(input$example_data, {

          intensity_file = paste0(system.file("example_data", package="msqrob2PTMGUI"),
                                  "/example_intensitiesfile.csv")

          metadata_file = paste0(system.file("example_data", package="msqrob2PTMGUI"),
                                  "/example_metadatafile.csv")
          df <- read.delim(intensity_file,
                           sep = ",", header = TRUE)
          df$peptidoform <- paste(df[,input$sequenceColumn],
                                  df[,input$modificationsColumn], sep = "_")

          #Get intensity columns
          ecols <- grep("2021",
                        names(df))

          #Read in peptides intensity file into QFeatures object
          pe <- readQFeatures(df,
                              ecol = ecols,
                              name = "peptideRaw",
                              fnames = "peptidoform")

          #Read in metadatafile
          metadataFile <- read.delim(metadata_file,
                                     sep = ",")

          variables$pe <- pe
          variables$metadataFile <- metadataFile
          variables$ecols <- ecols
          variables$protein_included <- FALSE
          variables$input_example <- TRUE

          output$read_in_example_data <- renderText({
            "Data read in complete. You may now continue to the preprocessing tab."
          })
    })


    #If necessary, generate metadatafile
     observeEvent(input$generateMetaData, {
      req(input$data$name)
      req(input$intensityIdentifier)

      runs = grep(input$intensityIdentifier,
                  strsplit(readLines(input$data$datapath, 1), split = input$separator)[1][[1]],
                  value = T)
      pattern_ <- paste0("\"", input$intensityIdentifier, "|\"")
      runs = gsub(pattern_, "", runs)

      metadataFileToDownload = data.frame("comment[data.file]" = runs,
                                          "comment[technical.replicate]" = "to fill out",
                                          "factor.value[toFillOut]" = "to fill out")
      variables$metadataFileToDownload = metadataFileToDownload
      #browser()

    output$DownloadAnnot <- downloadHandler(
      filename = function() {
        paste0(gsub(" |:","-",Sys.time()),"_sdrf.tsv")
      },
      content = function(file) {
        write.table(variables$metadataFileToDownload, file, col.names = TRUE,
                    row.names = FALSE)
      }
    )

    output$downloadButtonDownloadAnnot<- renderUI({
      if(!is.null(metadataFileToDownload)) {
        downloadButton("DownloadAnnot", "Download Metadata File")}
    })
     })

    observeEvent(input$go, {
        #make sure files are actually uploaded
        req(input$data$name, input$metadata$name, input$intensityIdentifier)
        variables$protein_included <- FALSE
        variables$input_example <- FALSE

        df <- read.delim(input$data$datapath, skip = input$skip,
                         sep = input$separator, header = T)
        df$peptidoform <- paste(df[,input$sequenceColumn],
                                df[,input$modificationsColumn], sep = "_")
        #Get intensity columns
        ecols <- grep(input$intensityIdentifier,
                      names(df))

        #Read in peptides intensity file into QFeatures object
        pe <- readQFeatures(df,
                            ecol = ecols,
                            name = "peptideRaw",
                            fnames = "peptidoform")


        #Read in metadatafile
        ##if SDRF
        if (input$SDRF == T){
          metadataFile <- read.delim(input$metadata$datapath,
                                     sep = "\t")
          colnames_to_keep <- c(grep("data.file", colnames(metadataFile), value=T),
                                grep("technical.replicate", colnames(metadataFile), value=T),
                                grep("factor.value", colnames(metadataFile), value=T))
          metadataFile <- metadataFile %>% select(all_of(colnames_to_keep))
          colnames(metadataFile) <- c("datafile", "technicalreplicate",
                                      gsub(".*factor\\.value\\.(.*)\\.$", "\\1",
                                           grep("factor.value",
                                           colnames(metadataFile), value=T)))
        } else {
        metadataFile <- read.delim(input$metadata$datapath,
                                   sep = input$separatorMetadata)}

        #Read in proteinfile
        if(is.character(input$proteindata$datapath)){

          df_prot <- read.delim(input$proteindata$datapath, skip = input$proteinskip,
                                sep = input$proteinseparator, header = T)
          #Get intensity columns
          ecols_prot <- grep(input$proteinintensityIdentifier,
                        names(df_prot))

          #Read in peptides intensity file into QFeatures object
          prot <- readQFeatures(df_prot,
                              ecol = ecols_prot,
                              name = "proteinRaw")

          variables$prot <- prot
          variables$protein_included <- TRUE
        }
        variables$pe <- pe
        variables$metadataFile <- metadataFile

        output$read_in_data <- renderText({
          "Data read in complete. You may now continue to the preprocessing tab."
        })

    })

    observeEvent({input$go
                  input$example_data}, {

        #Make coldata for the pe object based on metadatafile
        idcols <- colnames(variables$metadataFile)[2:length(colnames(variables$metadataFile))]
        for (col in idcols){
            colData(variables$pe)[as.character(col)] <- as.factor(variables$metadataFile[[as.character(col)]])
        }

        #Update selectInput for arranging of x-axis
        updateSelectInput(session, "x_axis",
                          label = "select variable",
                          choices = variables$idcols,
                          selected = variables$idcols[1])

        #Filtering steps: already calculate nNonZero, zero -> NA is necessary
        rowData(variables$pe[["peptideRaw"]])$nNonZero <- rowSums(assay(variables$pe[["peptideRaw"]]) > 0, na.rm = T)
        pe <- zeroIsNA(variables$pe, i = "peptideRaw")
        if(is.character(input$proteindata$datapath)){
          rowData(variables$prot[["proteinRaw"]])$nNonZero <- rowSums(assay(variables$prot[["proteinRaw"]]) > 0, na.rm = T)
          prot <- zeroIsNA(variables$prot, i = "proteinRaw")
        }

        #Calculate some quick stats about the data
        features = paste(rowData(pe[["peptideRaw"]])[[input$sequenceColumn]],
                         rowData(pe[["peptideRaw"]])[[input$modificationsColumn]], sep="_")
        stats = tibble(stats = c("number of unique proteins",
                                 "number of unique peptides",
                                 "number of unique peptidoforms"))
        stats$values = c(length(unique(rowData(pe[["peptideRaw"]])[[input$proteinColumn]])),
                         length(unique(rowData(pe[["peptideRaw"]])[[input$sequenceColumn]])),
                         length(unique(features)))

        output$statstable <- renderTable(stats)

        #Plot the before preprocessing plots
        #for ggplot: pivot the intensity assay to long format
        assay_pe_long <- assay(pe[["peptideRaw"]]) %>%
          as_tibble() %>%
          pivot_longer(cols = everything(), names_to = "sample", values_to = "intensity")
        #set the colour palette
        pal <- wesanderson::wes_palette("Darjeeling1",
                                        length(unique(assay_pe_long$sample)),
                                        type = "continuous")
        #Density plot
        output$densityBefore <- renderPlotly({

          p1 <- ggplot(assay_pe_long, aes(x=intensity, color=sample)) +
            geom_density(show.legend = F) +
            scale_colour_manual(values = pal)
          p1 <- ggplotly(p1) %>% layout(showlegend=F) %>%
            config(toImageButtonOptions = list(
              format = "png", filename = "densityplot_before", width = 1920, height = 1080
            ),
            modeBarButtonsToRemove = list("toImage"))
          p1
        })
        #Boxplot
        output$boxplotBefore <- renderPlot({
          p1 = ggplot(assay_pe_long, aes(x = sample, y=intensity, color=sample)) +
            geom_boxplot(show.legend = F) +
            scale_color_manual(values=pal)+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =10))
          p1
        })

        variables$pe <- pe
        variables$idcols <- idcols
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    #---------------
    #Preprocessing
    #--------------
    #User will click on preprocessing button, then do preprocessing as indicated by user
    observeEvent(input$preprocess, {
      #Each time the preprocessing changes, this will be triggered,
      #so data needs to start from scratch every time
      pe2 <- variables$pe
      if(is.character(input$proteindata$datapath)){
      prot <- variables$prot}

      #logtransformation
      if (input$logTransform != "none"){
        pe2 <- logTransform(pe2, base = as.numeric(input$logTransform),
                           i = "peptideRaw", name = "peptideLog")
        if(is.character(input$proteindata$datapath)){
          prot <- logTransform(prot, base = as.numeric(input$logTransform),
                              i = "proteinRaw", name = "proteinLog")
        }
      }
      else if (input$logTransform == "none"){
        #This probably means the data are already logtransformed, or it is not necessary
        #So I will manually add the peptideRaw as peptideLog
        pe2 <- addAssay(pe2, pe2[["peptideRaw"]], "peptideLog")
        if(is.character(input$proteindata$datapath)){
          prot <- addAssay(prot, prot[["proteinRaw"]], "proteinLog")
        }
      }
      colData(pe2[["peptideLog"]]) <- colData(pe2)

      #filter on number of nonzero columns
      pe2 <- pe2[rowData(pe2[["peptideLog"]])$nNonZero>=input$nnonzero,,]
      if(is.character(input$proteindata$datapath)){
        prot <- prot[rowData(prot[["proteinLog"]])$nNonZero>=input$nnonzero,,]
      }

      protein_column <- colnames(rowData(pe2[["peptideLog"]]))[input$proteinColumn]
      #Normalisation
      #Either it is a robust normalisation, a centering or no normalisation
      if (input$normalisationMethodGlobal == "robust"){
        #If normalisation against protein is chosen, I do need to normalise against
        #technical variation as well, I typically use median centering for that
        pe2 <- normalize(pe2,
                         i = "peptideLog",
                         name = "peptideNorm",
                         method = "center.median")
        if(is.character(input$proteindata$datapath)){
          prot <- normalize(prot,
                           i = "proteinLog",
                           name = "proteinNorm",
                           method = "center.median")
          prot <- aggregateFeatures(prot,
                                    i = "proteinNorm",
                                    fcol = colnames(rowData(prot[["proteinLog"]]))[input$proteinproteinColumn],
                                    na.rm = T,
                                    name = "proteinRobust",
                                    fun = MsCoreUtils::robustSummary)
          pe2 <- addAssay(pe2, prot[["proteinRobust"]], "proteinRobust")

        }else{pe2 <- aggregateFeatures(pe2,
                                       i = "peptideNorm",
                                       fcol = protein_column,
                                       na.rm = TRUE,
                                       name = "proteinRobust",
                                       fun = MsCoreUtils::robustSummary)
        }


        #Do not change pe2 yet, gives problems further down the line
        #Normalisation step: subtract protein value from corresponding pepform value
        if(is.character(input$proteindata$datapath)){
          #Only take pepforms that have a parent protein
          pepWithProtein <- rowData(pe2[["peptideNorm"]])[[input$proteinColumn]] %in%
                                  rownames(prot[["proteinRobust"]])
          pePepWithProtein  <- pe2[["peptideNorm"]][pepWithProtein,]
          pe2 <- addAssay(pe2,pePepWithProtein,"peptideLogNorm")
        }else{pe2 <- addAssay(pe2,pe2[["peptideNorm"]],"peptideLogNorm")}

        assay(pe2[["peptideLogNorm"]]) <- assay(pe2[["peptideLogNorm"]]) -
          assay(pe2[["proteinRobust"]])[rowData(pe2[["peptideLogNorm"]])[[input$proteinColumn]],
                                        colnames(assay(pe2[["peptideLogNorm"]]))]
        #Filter out peptidoforms that now have only NA or 0 intensities
        filtering <- rowSums(is.na(assay(pe2[["peptideLogNorm"]])) |
                               assay(pe2[["peptideLogNorm"]])==0) != ncol(assay(pe2[["peptideLogNorm"]]))
        pe2[["peptideLogNorm"]] <- pe2[["peptideLogNorm"]][filtering,]
      }

      else if (input$normalisationMethodGlobal != "none"){
        pe2 <- normalize(pe2,
                        i = "peptideLog",
                        name = "peptideLogNorm",
                        method = input$normalisationMethodGlobal)
      }
      else if (input$normalisationMethodGlobal == "none"){
        #I do continue with the name peptideLognorm, so I will have to add this
        #in this case as well, it is then just unnormalised
        pe2 <- addAssay(pe2, pe2[["peptideLog"]], "peptideLogNorm")
      }
      #PTM summarisation
      if (input$summarisation == "yes"){
          #ptm is the form of protein + modification + location (in the protein)
          #so add ptm column
          rowData(pe2[["peptideLogNorm"]])$ptm <- mapply(function(x, y) {
          vector2 <- strsplit(y, paste0("\\",input$modificationSplit))[[1]]
          result <- lapply(vector2, function(i) paste(x, i, sep = "__"))
          final_result <- paste(unlist(result), collapse = ", ")
        }, rowData(pe2[["peptideLogNorm"]])[[input$proteinColumn]],
        rowData(pe2[["peptideLogNorm"]])[[input$modificationsColumn]])
          #Get all unique PTMs
          ptms <- unique(unlist(strsplit(rowData(pe2[["peptideLogNorm"]])$ptm, split = ", ", fixed = TRUE)))

          #For each ptm do
          ptm_x_assay <- sapply(unique(ptms), function(i){

            #ptm_sub <- filterFeatures(pe2,~grepl(ptm,pattern=i,fixed = T))[["peptideLogNorm"]]
            ptm_sub <- pe2[["peptideLogNorm"]][grepl(rowData(pe2[["peptideLogNorm"]])$ptm, pattern = i, fixed = T),]
            
            #Get intensity values of those peptidoforms
            z <- assay(ptm_sub)
            z <- filter(as.data.frame(z), rowSums(is.na(z) | z==0) != ncol(z))
            #And summarise them to 1 row of intensity values: 1 value per sample for that ptm
            ptm_y <- try(MsCoreUtils::robustSummary(as.matrix(z)), silent = T)
            if (is(ptm_y, "try-error"))
            {ptm_y <- matrix(rep(NA,ncol(ptm_sub)), nrow = 1,
                             dimnames = list(c(i), c(colnames(ptm_sub))))}

            ptm_y
          })
          #Then we get the intensity assay on ptm level
          ptm_x_assay <- t(ptm_x_assay)
          
          #Filter out ptms that have too many missing values
          filtering <- (rowSums(ptm_x_assay != 0, na.rm=TRUE)) > 0 &
            (rowSds(ptm_x_assay, na.rm=TRUE) > 1e-4) &
            (rowSums(is.na(ptm_x_assay)) < (ncol(ptm_x_assay)-1))
          ptm_x_assay <- ptm_x_assay[filtering,]
          colnames(ptm_x_assay) <- colnames(assay(pe2[["peptideLogNorm"]]))
          #Add to QFeatures object
          rowdata <- data.frame(ptm = rownames(ptm_x_assay))
          rowdata[[protein_column]] <- sapply(str_split(rowdata$ptm, pattern="_"),
                                                    function(x) x[1])
          ptm_y_assay <- SummarizedExperiment(assays=as.matrix(ptm_x_assay),
                                              rowData=rowdata, colData=colData(pe2)[colnames(ptm_x_assay),])
          pe2 <- addAssay(pe2, ptm_y_assay, "ptmRel")
      }
      #Update selectInput for data table as some of proteins may have been filtered out
      updateSelectInput(session, "protein",
                        label = "Protein",
                        choices = rowData(pe2[["peptideLogNorm"]])[[input$proteinColumn]],
                        selected = rowData(pe2[["peptideLogNorm"]])[[input$proteinColumn]][1])

      #Update selectInput for arranging of x-axis
      updateSelectInput(session, "x_axis",
                        label = "select variable",
                        choices = variables$idcols,
                        selected = variables$idcols[1])

      #Calculate some quick stats about the data after preprocessing
      features = paste(rowData(pe2[["peptideLogNorm"]])[[input$sequenceColumn]],
                       rowData(pe2[["peptideLogNorm"]])[[input$modificationsColumn]], sep="_")
      stats = tibble(stats = c("number of unique proteins",
                               "number of unique peptides",
                               "number of unique peptidoforms"))
      stats$values = c(length(unique(rowData(pe2[["peptideLogNorm"]])[[input$proteinColumn]])),
                       length(unique(rowData(pe2[["peptideLogNorm"]])[[input$sequenceColumn]])),
                       length(unique(features)))
      if (input$summarisation == "yes"){
        stats = rbind(stats, (tibble(stats = c("number of unique ptms"),
                                     values = c(length(unique(rowData(pe2[["ptmRel"]])$ptm))))))
      }

      output$statstableafter <- renderTable(stats)

      #Plot the after preprocessing plots
      #for ggplot: pivot the intensity assay to long format
      assay_pe_long <- assay(pe2[["peptideLogNorm"]]) %>%
        as_tibble() %>%
        pivot_longer(cols = everything(), names_to = "sample", values_to = "intensity")
      #set the colour palette
      pal <- wesanderson::wes_palette("Darjeeling1",
                                      length(unique(assay_pe_long$sample)),
                                      type = "continuous")

      #Density plot
      output$densityAfter <- renderPlotly({

        p1 <- ggplot(assay_pe_long, aes(x=intensity, color=sample)) +
          geom_density(show.legend = F) +
          scale_colour_manual(values = pal)
        p1 <- ggplotly(p1) %>% layout(showlegend=F) %>%
          config(toImageButtonOptions = list(
            format = "png", filename = "densityplot_after", width = 1920, height = 1080
          ),
          modeBarButtonsToRemove = list("toImage"))
        p1
      })

      #Boxplot
      output$boxplotAfter <- renderPlot({
        p1 = ggplot(assay_pe_long, aes(x = sample, y=intensity, color=sample)) +
          geom_boxplot(show.legend = F) +
          scale_color_manual(values=pal)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size =10))
        p1
      })

      rowData(pe2[["peptideLogNorm"]])$pepform <- paste(
                      rowData(pe2[["peptideLogNorm"]])[,input$sequenceColumn],
                      rowData(pe2[["peptideLogNorm"]])[,input$modificationsColumn])
      variables$pe2 <- pe2
      variables$protein_column <- protein_column

    })

    #----------------------------
    #Data exploration
    #----------------------------

    #If user clicks on radiobuttons for normalisation or selects different protein,
    #or user changes preprocessing settings, or chooses to arrange x-axis differently on plots,
    #variables$proteindf needs to update to the correct version
    observeEvent({input$protein
                 input$normalisationMethod
                 input$preprocess
                 input$x_axis}, {

        req(input$x_axis)

        #Get data for particular protein
        proteinpe <- variables$pe2[["peptideLogNorm"]][grepl(input$protein,
                                rowData(variables$pe2[["peptideLogNorm"]])[[input$proteinColumn]],
                                fixed = T)]
        #Normalisation + get dataset instead of QFeatures object
        #get all metadata + intensity values together in df
        if( input$normalisationMethod != "none"){
            proteinpe <- normalize(proteinpe,
                            method = input$normalisationMethod)
            df <- longFormat(proteinpe)
        } else if (input$normalisationMethod == "none"){
            df <- longFormat(proteinpe)
        }
        colH <- colData(proteinpe)
        for (col in variables$idcols){
            df[as.character(col)] = as.factor(colH[df$colname, col])
        }
        #add id column
        df <- df %>% as_tibble() %>% tidyr::unite("id", all_of(variables$idcols),sep="_", remove = F)
        #add biorepeat column
        df$biorepeat <- df$id %>% as.factor %>% as.double
        #add features column
        df$features <-
            rep(paste(rowData(proteinpe)[,input$sequenceColumn],
                      rowData(proteinpe)[,input$modificationsColumn],sep="_"),
                      length(unique(df$biorepeat)))
        #arrange according to input$x_axis
        df <- as.data.frame(df)
        df <- df[order(df[,input$x_axis]),]
        df$id = factor(df$id, unique(df$id))
        #Save dataset
        variables$proteindf <- df
        }, ignoreInit = TRUE)

    #Visualisation
    #Plot data table: wide format so that users can easily see the features
    output$proteinDataTable <- DT::renderDataTable(server = FALSE, {
        #Transform dataset into wide format
        proteindf_wide <- variables$proteindf %>% tibble::as_tibble() %>%
            tidyr::pivot_wider(id_cols = c("features"), names_from = "id", values_from = "value")
        variables$proteindf_wide = proteindf_wide
        DT::datatable(proteindf_wide,
            extensions = "Buttons",
            options = list(
            paging = TRUE,
            pageLength = 4,
            searching = TRUE,
            dom = "Bfrtip",
            buttons = list("copy", list(
              extend = "collection",
              buttons = c("csv", "excel"),
              text = "Download",
              exportOptions = list(
                modifier = list(page="all")
              ))
            ))
            ) %>% DT::formatStyle(names(proteindf_wide), lineHeight="80%")
    })

    #If user clicks deselect button -> all selected rows are deselected
    proxy = dataTableProxy('proteinDataTable')
    observeEvent(input$deselect, {
      proxy %>% selectRows(NULL)
    })


    #delay input$proteinDataTable_rows_selected so that it does not
    #recalculate every time a row is clicked
    rows_selected <- reactive(input$proteinDataTable_rows_selected)
    rows_selected_d <- debounce(rows_selected, 1000)

    #Plot lineplot
    #Use plotly to make the plot interactive
    output$lineplot <- renderPlotly({
      #set colour palette
      pal <- wesanderson::wes_palette("Darjeeling1",
                                      variables$proteindf %>% pull(rowname) %>%
                                        unique %>% length,
                                      type = "continuous")

      #base lineplot (ggplot)
      p1 <- ggplot(
        data = variables$proteindf %>% as.data.frame,
        aes(x = as.factor(id), y = value, group = rowname,colour=rowname)) +
        geom_line(show.legend = F)  +
        scale_colour_manual(values = pal) +
        ggtitle("log2-normalised data") +
        xlab("sample id") +
        ylab("Intensity (log2)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


      #Plot base lineplot already when nothing is selected
      if (is.null(rows_selected_d())){
        #turn ggplot into plotly object
        p1 <- ggplotly(p1) %>% layout(showlegend=F) %>%
          config(toImageButtonOptions = list(
            format = "png", filename = "lineplot", width = 1920, height = 1080
          ),
          modeBarButtonsToRemove = list("toImage"))
        print(p1)
        }

      #If a row is selected, highlight that row
        else if (!is.null(rows_selected_d())){
          #Get dataset with selected rows
          features_selected <- variables$proteindf_wide[rows_selected_d(),]$features
          filter_proteindf <- variables$proteindf %>% as_tibble %>%
              filter(features %in% features_selected)

          hlp1 <- p1 +
              geom_line(
                  aes(x=as.factor(id),
                      y=value),
                  size=2.3,
                  color="palevioletred4",
                  show.legend=FALSE,
                  data = as_tibble(filter_proteindf)) +
              geom_line(
                  aes(x=as.factor(id),
                      y=value),
                  size=0.7,
                  color="grey85",
                  show.legend=FALSE,
                  data = as_tibble(filter_proteindf))

          hlp1 <- ggplotly(hlp1) %>% layout(showlegend=F) %>%
            config(toImageButtonOptions = list(
              format = "png", filename = "lineplot", width = 1920, height = 1080
            ),
            modeBarButtonsToRemove = list("toImage"))
          print(hlp1)
        }
    })

    #Plot boxplot
    output$boxplot <- renderPlotly({
      #set colour palette
      pal <- wesanderson::wes_palette("Darjeeling1",
                                      length(unique(variables$proteindf$id)),
                                      type = "continuous")
      #base boxplot
      boxplot <- variables$proteindf %>% as_tibble() %>%
        ggplot(aes(x=as.factor(id),y=value,col=as.factor(id))) +
        geom_boxplot() +
        #geom_jitter(aes(col=features),size=0.7) +
        xlab("sample id") +
        ylab("Intensity (log2)") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.title = element_blank(), legend.position = "None") +
        scale_colour_manual(values = pal)

      #If nothing is selected, already print base boxplot
      if(is.null(rows_selected_d())){
        boxplot <- ggplotly(boxplot) %>% layout(showlegend=F) %>%
          config(toImageButtonOptions = list(
            format = "png", filename = "boxplot", width = 1920, height = 1080
          ),
          modeBarButtonsToRemove = list("toImage"))
        print(boxplot)
      }
      #If something is selected, highlight selected dots
       else if (!is.null(rows_selected_d())){
         #Get dataset with selected rows
            features_selected <- variables$proteindf_wide[rows_selected_d(),]$features
            filter_proteindf <- variables$proteindf %>% as_tibble %>%
                filter(features %in% features_selected)

            boxplot_select <- boxplot +
              geom_point(data=filter_proteindf,
                         aes(x=as.factor(id),y=value), color = "red", size = 1.2)

            boxplot_select <- ggplotly(boxplot_select) %>% layout(showlegend=F) %>%
              config(toImageButtonOptions = list(
                format = "png", filename = "boxplot", width = 1920, height = 1080
              ),
              modeBarButtonsToRemove = list("toImage"))
            print(boxplot_select)
            }
    })

    #--------------------------
    #Modeling tab
    #--------------------------

    ##output possible modeling variables
    output$available_modelvariables <- renderText({
      colnames(variables$metadataFile)},
      sep = ", ")

    ##output metadata
    output$designVariables <- renderDataTable(variables$metadataFile,
                                              options = list(
                                                pageLength = 5
                                              ))

    ##output designmatrix
    design <- reactive(
    #If the formula contains a random effect, remove it in order to use VisualizeDesign
    if (any(grepl("\\|",attr(terms(as.formula(input$designformula)), "term.labels")))){
      design <- ExploreModelMatrix::VisualizeDesign(colData(variables$pe),update(as.formula(input$designformula), as.formula(paste("~. -",paste0("(",attr(terms(as.formula(input$designformula)), "term.labels")[grepl("\\|", attr(terms(as.formula(input$designformula)), "term.labels"))], ")")))))
    } else{
      design <- ExploreModelMatrix::VisualizeDesign(variables$metadataFile,input$designformula)}
    )

    output$designmatrix <- renderUI({
      lapply(1:length(design()[[2]]),
             function(j){
               renderPlot(design()[[2]][[j]])
             })
    })

    ##model the data when user clicks button
    observeEvent(input$fitModel,{
      pe2 <- variables$pe2
      #check whether peptideLogNorm actually exists
      if (!"peptideLogNorm" %in% names(pe2)){
        output$model_fitted <- renderText({"Please go through the preprocessing step first"})
      }
      req(pe2[["peptideLogNorm"]])

      pe2 <- msqrob2::msqrob(object = pe2, i = "peptideLogNorm",
                             formula = stats::as.formula(input$designformula),
                             robust=input$robustRegression, overwrite = T)
      output$model_fitted <- renderText({"Model fitting complete"})
      variables$pe2 <- pe2
      #if we have a ptm assay, also do modelling for ptm
      req(pe2[["ptmRel"]])
      pe2 <- msqrob2::msqrob(object = pe2, i = "ptmRel",
                             formula = stats::as.formula(input$designformula),
                             robust=input$robustRegression, overwrite = T)

      output$model_fitted <- renderText({"Model fitting complete"})
      variables$pe2 <- pe2
    })


    #------------------------
    #Inference
    #------------------------

    #Inference tab
    ##output possible parameter names
    output$available_parameters <- renderText({
      colnames(design()$designmatrix)
      }, sep = ", ")

    observeEvent(input$testcontrast, {
      pe2 <- variables$pe2
      contrasts <- strsplit(input$contrast, ", ")[[1]] %>% as.vector()
      L <- makeContrast(contrasts, parameterNames = colnames(design()$designmatrix))
      pe2 <- hypothesisTest(object = pe2, i = "peptideLogNorm", contrast = L, overwrite = TRUE)

      # If we have a ptm assay, also do inference for ptm
      if (input$summarisation == "yes") {
        pe2 <- hypothesisTest(object = pe2, i = "ptmRel", contrast = L, overwrite = TRUE)
      }

      variables$pe2 <- pe2
      variables$contrasts <- contrasts

      # Return a value to trigger the event for the second observeEvent
      pe2
    })

    clicked_feature <- reactive(input$significanceTable_rows_selected)

    observeEvent({input$testcontrast
                 input$significancelevel
                 input$onlysignificant}, {
      req(variables$pe2[["peptideLogNorm"]])

      pe2 <- variables$pe2
      contrasts <- variables$contrasts
      contrasts <- gsub("\\s*=\\s*0", "", contrasts)

      #Render significance table
      filter_alpha <- ifelse(input$onlysignificant, input$significancelevel, 1)
      # if (input$onlysignificant) {
      #   filter_alpha <- input$significancelevel
      # } else {
      #   filter_alpha <- 1
      # }
      result_assay <- ifelse(input$summarisation == "yes", "ptmRel", "peptideLogNorm")
      # if (input$summarisation == "yes"){
      #   result_assay <- "ptmRel"
      # } else {result_assay <- "peptideLogNorm"}

      sign_table <- do.call(rbind, lapply(contrasts, function(contrast) {
        rowData(pe2[[result_assay]])[[contrast]] %>%
          na.exclude() %>%
          filter(adjPval <= filter_alpha) %>%
          arrange(pval) %>%
          mutate(Contrast = contrast)
      }))

      output$significanceTable <- DT::renderDataTable(server = FALSE, {
        DT::datatable(sign_table,
                      extensions = "Buttons",
                      options = list(
                        paging = TRUE,
                        pageLength = 5,
                        searching = TRUE,
                        dom = "Bfrtip",
                        buttons = list("copy", list(
                          extend = "collection",
                          buttons = c("csv", "excel"),
                          text = "Download",
                          exportOptions = list(
                            modifier = list(page = "all")
                          ))
                        ))
        ) %>% DT::formatStyle(names(sign_table), lineHeight = "80%") %>%
              DT::formatRound(1:6,
                              digits = 3)
      })
      variables$sign_table <- sign_table

      annotation_dots <- ifelse(input$summarisation == "yes", "ptm", "pepform")

      #Plot volcano plot

      # volcanos <- list()
      # volcanos[[contrast]] <- p1
      # combined_volcanos <- subplot(plotlist = volcanos, nrows = length(volcanos))
      # for (contrast in contrasts){

        p1 <- rowData(pe2[[result_assay]])[[contrasts[1]]] %>%
          ggplot(aes(x = logFC, y = -log10(pval),
                     color = adjPval <= input$significancelevel,
                     annotation = rowData(pe2[[result_assay]])[, annotation_dots])) +
          geom_point(cex = 2.5) +
          scale_color_manual(values = alpha(c("black", "red"), 0.5)) +
          theme_minimal() +
          ylab("-log10(pvalue)") +
          ggtitle(contrasts[1])
        p1_ly <- ggplotly(p1) %>%
          layout(showlegend = F) %>%
          config(toImageButtonOptions = list(
            format = "png", filename = "lineplot", width = 1920, height = 1580
          ),
          modeBarButtonsToRemove = list("toImage"))


      output$volcano <- renderPlotly({p1_ly})
      variables$volcanoplot <- p1
      variables$contrasts_short <- contrasts
      variables$result_assay <- result_assay
    })

    selected_feature <- reactiveVal()

    observeEvent(input$significanceTable_rows_selected, {
      #Update selectInput for arranging of x-axis
      updateSelectInput(session, "x_axis1",
                        label = "select 1st variable",
                        choices = variables$idcols,
                        selected = variables$idcols[1])
      if(length(variables$idcols) > 1){
        updateSelectInput(session, "x_axis2",
                          label = "select 2nd variable",
                          choices = variables$idcols,
                          selected = variables$idcols[2])
      }
      if(length(variables$idcols) > 2){
        updateSelectInput(session, "x_axis3",
                          label = "select 3rd variable",
                          choices = variables$idcols,
                          selected = variables$idcols[3])
      }
      selected_rows <- input$significanceTable_rows_selected
      if (!is.null(selected_rows) && length(selected_rows) > 0) {
        selected_features = variables$sign_table[selected_rows,] %>% rownames()
        selected_feature(selected_features)
      } else {
        selected_feature(NULL)
      }
      p1 <- variables$volcanoplot
      pe2 <- variables$pe2
      contrasts <- variables$contrasts_short
      result_assay <- variables$result_assay
      annotation_dots <- ifelse(input$summarisation == "yes", "ptm", "pepform")

      p1 <- p1 +
        geom_point(data = rowData(pe2[[result_assay]])[[contrasts[1]]][selected_feature(),],
                   aes(x = logFC, y = -log10(pval),
                       annotation = rowData(pe2[[result_assay]])[selected_feature(), annotation_dots]),
                   color = ifelse(!is.null(selected_feature()), "chartreuse",
                              ifelse(rowData(pe2[[result_assay]])[[contrasts[1]]][selected_rows,]$adjPval <= input$significancelevel,
                                      "red", "black")),
                   cex = 2.5)

      p1_ly <- ggplotly(p1) %>%
        layout(showlegend = F) %>%
        config(toImageButtonOptions = list(
          format = "png", filename = "lineplot", width = 1920, height = 1580
        ),
        modeBarButtonsToRemove = list("toImage"))
      output$volcano <- renderPlotly({p1_ly})
      #print significance table corresponding peptidoforms in case of ptms
      if(input$summarisation == "yes"){
        pepform_table <- rowData(pe2[["peptideLogNorm"]])[grepl(selected_feature(),
                              rowData(pe2[["peptideLogNorm"]])$ptm, fixed = T),] %>%
                        as.data.frame() %>% select(starts_with(contrasts[1])) %>%
                        arrange(!!sym(paste(contrasts[1], "adjPval", sep = ".")))

        output$significanceTablePepforms <- DT::renderDataTable(server = FALSE, {
          DT::datatable(pepform_table,
                        extensions = "Buttons",
                        options = list(
                          paging = TRUE,
                          pageLength = 3,
                          searching = TRUE,
                          dom = "Bfrtip",
                          buttons = list("copy", list(
                            extend = "collection",
                            buttons = c("csv", "excel"),
                            text = "Download",
                            exportOptions = list(
                              modifier = list(page = "all")
                            ))
                          ))
          ) %>% DT::formatStyle(names(pepform_table), lineHeight = "80%") %>%
            DT::formatRound(1:6,
                            digits = 3)
        })
        variables$pepform_table <- pepform_table
      }


    })


    #------------------------
    #Detail Plots
    #------------------------

    clicked_feature <- reactive(input$significanceTable_rows_selected)

    output$lineplot_feature <- renderPlotly({
      #add error handling for zero or two or more features selected

      #Get dataset with selected rows
      pe2 <- variables$pe2
      sign_table <- variables$sign_table
      feature_selected <- sign_table[clicked_feature(),] %>% rownames()
      result_assay <- variables$result_assay


      if (result_assay == "peptideLogNorm"){
        pep <- feature_selected
        prot <- rowData(pe2[[result_assay]])[pep,input$proteinColumn]
      } else if (result_assay == "ptmRel"){
        prot <- strsplit(feature_selected, "__")[[1]][1]
        pep <- rowData(pe2[["peptideLogNorm"]]) %>% as.data.frame() %>%
          filter(grepl(feature_selected, ptm, fixed=T)) %>% rownames
      }
        #get protein value by aggregating subset
        pe2_sub_prot <- pe2[["peptideLog"]][rowData(pe2[["peptideLog"]])[[variables$protein_column]] == prot,]
        pe2_sub_prot <- aggregateFeatures(pe2_sub_prot,
                                 i = "peptideLog",
                                 fcol = variables$protein_column,
                                 na.rm = TRUE,
                                 name = "proteinLog",
                                 fun = MsCoreUtils::robustSummary)

        protein <- assay(pe2_sub_prot) %>% as.data.frame() %>% pivot_longer(
          cols = colnames(assay(pe2_sub_prot)), values_to = "LogIntensities",
          names_to = "run"
        ) %>% mutate(type = "protein", sequence = prot)
        #add biorepeat column
        #df$biorepeat <- df$id %>% as.factor %>% as.double

        #peptidoform abundance
        pepform <- assay(pe2[["peptideLog"]])[pep,] %>% as.data.frame()
        if (length(pep)>1){
          pepform <- pepform %>% rownames_to_column("sequence") %>%
            pivot_longer(cols = colnames(pepform), names_to = "run",
                         values_to = "LogIntensities") %>%
            mutate(type = "pepform -log")
        } else{
          colnames(pepform) <- "LogIntensities"
          pepform <- pepform %>% rownames_to_column("run") %>%
            mutate(sequence = pep, type = "pepform -log")
        }
        #peptidoform abundance - normalised
        pepform_norm <- assay(pe2[["peptideLogNorm"]])[pep,] %>% as.data.frame()
        if (length(pep)>1){
          pepform_norm <- pepform_norm %>% rownames_to_column("sequence") %>%
            pivot_longer(cols = colnames(pepform_norm), names_to = "run",
                         values_to = "LogIntensities")
        } else{
          colnames(pepform_norm) <- "LogIntensities"
          pepform_norm <- pepform_norm %>% rownames_to_column("run") %>%
            mutate(sequence = pep)
        }
        pepform_norm <- pepform_norm %>% mutate(type = "pepform - norm",
                            sequence = paste(sequence, "norm"))
        #ptm abundance - normalised
        if(result_assay == "ptmRel"){
          ptm_df <- assay(pe2[["ptmRel"]])[feature_selected,] %>% as.data.frame()
          colnames(ptm_df) <- "LogIntensities"
          ptm_df <- ptm_df %>%
            rownames_to_column("run") %>% mutate(type = "ptm",
                                                 sequence = feature_selected)
        }

        if(result_assay == "ptmRel"){
          df_plot <- rbind(protein, pepform, pepform_norm, ptm_df)
          sign_pepforms <- variables$pepform_table %>% filter(
            !!sym(paste(gsub("\\s*=\\s*0", "", variables$contrasts[1]), "adjPval",
                        sep = ".")) <= input$significancelevel
          ) %>% rownames %>% paste("norm")
          df_plot[df_plot$sequence==sign_pepforms,]$type = "pepform - norm - significant"
        } else {df_plot <- rbind(protein, pepform, pepform_norm)}

        colH <- colData(pe2)

        for (col in variables$idcols){
          df_plot[as.character(col)] = as.factor(colH[df_plot$run, col])
        }
        #add id column
        df_plot <- df_plot %>% as_tibble() %>%
          tidyr::unite("id", all_of(variables$idcols),sep="_", remove = F)

        if (result_assay == "ptmRel"){
        colour_values <- c("pepform - norm" = "gray72",
                           "pepform - norm - significant" = "pink",
                           "protein" = "dodgerblue2", "ptm" = "seagreen3",
                           "pepform -log" = "grey44")}
        else{
          colour_values <- c("pepform - norm" = "seagreen3",
                             "protein" = "dodgerblue2",
                             "pepform -log" = "palevioletred2")}

        if(length(variables$idcols) == 1 ){
          df_plot <- df_plot %>% arrange(!!sym(input$x_axis1))}
        else if(length(variables$idcols) == 2 ){
          df_plot <- df_plot %>% arrange(!!sym(input$x_axis1), !!sym(input$x_axis2))
        }
        else if(length(variables$idcols) > 2){
          df_plot <- df_plot %>%
            arrange(!!sym(input$x_axis1), !!sym(input$x_axis2), !!sym(input$x_axis3))
        }
        df_plot$id <- forcats::fct_inorder(df_plot$id)

        p1 <-  df_plot %>% ggplot() +
          geom_line(aes(x = id, y = LogIntensities, group = sequence,
                        color = type), size = 1) +
          geom_point(aes(x = id, y = LogIntensities , group = sequence,
                         color = type), size = 2.5) +
          scale_colour_manual(values = colour_values)  +
          labs(title = feature_selected, x = "BioReplicate", y = "LogIntensity") +
          theme_bw() +
           theme(axis.text.x = element_text(angle = 60, hjust=1))


        p1 <- ggplotly(p1, width = 1520, height= 780)  %>%
          config(toImageButtonOptions = list(
            format = "png", filename = paste("detail_plot", feature_selected),
            width = 1920, height = 1080
          ),
          modeBarButtonsToRemove = list("toImage"))
        return(p1)
    })

    getDataPath <- function(datapath){
      if(Sys.info()['sysname']=="Windows"){
        datapath <- gsub("\\","/",datapath, fixed=TRUE)
      }
      return(datapath)
    }

    #REPORT
    ### use metaReactive to write options to file
    intensityDatapath <- metaReactive({ifelse(variables$input_example == T,
                                              ..(paste0(system.file("example_data", package="msqrob2PTMGUI"),
                                                          "/example_intensitiesfile.csv")),
                                              getDataPath(..(input$data$datapath)))})
    annotationDatapath <- metaReactive({ifelse(variables$input_example == T,
                                               ..(paste0(system.file("example_data", package="msqrob2PTMGUI"),
                                                        "/example_metadatafile.csv")),
                                               getDataPath(..(input$metadata$datapath)))})
    SDRF_file <- metaReactive({..(input$SDRF)}, varname = "SDRF_file")
    proteinDatapath <- metaReactive({getDataPath(..(input$proteindata$datapath))})
    protein_included <- metaReactive({..(variables$protein_included)}, varname = "protein_included")
    input_example <- metaReactive({..(variables$input_example)}, varname = "input_example")
    intensityIdentifier <- metaReactive({..(input$intensityIdentifier)}, varname = "intensityIdentifier")
    skip <- metaReactive({..(input$skip)}, varname = "skip")
    separator <- metaReactive({..(input$separator)}, varname = "separator")
    proteinskip <- metaReactive({..(input$proteinskip)}, varname = "proteinskip")
    proteinseparator <- metaReactive({..(input$proteinseparator)}, varname = "proteinseparator")
    proteinintensityIdentifier <- metaReactive({..(input$proteinintensityIdentifier)}, varname = "proteinintensityIdentifier")
    annotationSep <- metaReactive({..(input$separatorMetadata)}, varname = "annotationSep")
    sequenceColumn <- metaReactive({..(input$sequenceColumn)}, varname = "sequenceColumn")
    modificationsColumn <- metaReactive({..(input$modificationsColumn)}, varname = "modificationsColumn")
    modificationSplit <- metaReactive({..(input$modificationSplit)}, varname = "modificationSplit")
    proteinColumn <- metaReactive({..(input$proteinColumn)}, varname = "proteinColumn")
    proteinproteinColumn <- metaReactive({..(input$proteinproteinColumn)}, varname = "proteinproteinColumn")
    logTrans <- metaReactive({..(input$logTransform)}, varname = "logTrans")
    minObsFeat <- metaReactive({..(input$nnonzero)}, varname ="minObsFeat")
    normMethod <- metaReactive({..(input$normalisationMethodGlobal)}, varname = "normMethod")
    summarisation <- metaReactive({..(input$summarisation)}, varname ="summarisation")
    formula <- metaReactive({..(input$designformula)}, varname = "formula")
    robustRegression <- metaReactive({..(input$robustRegression)}, varname = "robustRegression")
    contrast <- metaReactive({..(input$contrast)}, varname = "contrast")
    sigLevel <- metaReactive({..(input$significancelevel)}, varname = "sigLevel")
    selectedFeature <- metaReactive({..(input$significanceTable_rows_selected)}, varname = "selectedFeature")
    maxPlot <- metaReactive({..(input$maxPlot)}, varname = "maxPlot")
    x_axis1 <- metaReactive({..(input$x_axis1)}, varname = "x_axis1")
    x_axis2 <- metaReactive({..(input$x_axis2)}, varname ="x_axis2")
    x_axis3 <- metaReactive({..(input$x_axis3)}, varname = "x_axis3")
    
    #report <- reactiveValues(filepath = NULL) #This creates a short-term storage location for a filepath
    variables$filepath = NULL
    
    observeEvent(input$generate, {
     
      # filename = function() {
      #   paste0(input$project_name,"-report-", gsub(" |:","-",Sys.time()),".zip")
      # }
      # content = function(file) {
      print(intensityDatapath())
      print(annotationDatapath())
         file.copy(from = intensityDatapath(), to = "intensityFile.txt", overwrite = TRUE)
         file.copy(from = annotationDatapath(), to = "annotationFile.csv", overwrite = TRUE)
         if(protein_included()==T){
           file.copy(from = proteinDatapath(), to = "intensity_nonenriched.txt", overwrite = TRUE)
         }
        
        if(protein_included()==T){
        inputfiles <- expandChain(
          quote({
            intensityFile <- "intensityFile.txt"
            annotationFile <- "annotationFile.csv"
            proteinFile <- "intensity_nonenriched.txt"
          }))}
        else if (protein_included()==F){
          inputfiles <- expandChain(
            quote({
              intensityFile <- "intensityFile.txt"
              annotationFile <- "annotationFile.csv"
            }))}

        inputparameters <- expandChain(
          invisible(protein_included()),
          invisible(input_example()),
          invisible(intensityIdentifier()),
          invisible(skip()),
          invisible(separator()),
          invisible(proteinintensityIdentifier()),
          invisible(proteinskip()),
          invisible(proteinseparator()),
          invisible(annotationSep()),
          invisible(SDRF_file()),
          invisible(sequenceColumn()),
          invisible(modificationsColumn()),
          invisible(modificationSplit()),
          invisible(proteinColumn()),
          invisible(proteinproteinColumn())
        )
        preprocessing <- expandChain(
          invisible(logTrans()),
          invisible(minObsFeat()),
          invisible(normMethod())
        )
        summarization <- expandChain(invisible(summarisation()))
        model <- expandChain(
          invisible(formula()),
          invisible(robustRegression())
        )
        inference <- expandChain(
          invisible(contrast()),
          invisible(sigLevel()))
        report <- expandChain(
          invisible(selectedFeature()),
          invisible(maxPlot()),
          invisible(x_axis1()),
          invisible(x_axis2()),
          invisible(x_axis3())
        )
        
        tmp_file <- getDataPath(paste0(tempfile(), ".zip"))#Creating the temp where the .pdf is going to be stored
        if(Sys.info()['sysname']=="Windows"){
          tmp_file <- getDataPath(paste0(tempfile(tmpdir = system.file(package="msqrob2PTMGUI")), ".zip"))
        }
        
        if(protein_included()==T){
          buildRmdBundle(
            system.file("data/Report.Rmd",package="msqrob2PTMGUI"),
            tmp_file,
            list(
              inputfiles = inputfiles,
              inputparameters = inputparameters,
              preprocessing = preprocessing,
              summarization = summarization,
              model = model,
              inference = inference,
              report = report
            ),
            render=TRUE,
            include_files = c("intensityFile.txt",'annotationFile.csv',"intensity_nonenriched.txt")
          )}
        else if (protein_included()==F){
          buildRmdBundle(
            system.file("data/Report.Rmd",package="msqrob2PTMGUI"),
            tmp_file,
            list(
              inputfiles = inputfiles,
              inputparameters = inputparameters,
              preprocessing = preprocessing,
              summarization = summarization,
              model = model,
              inference = inference,
              report = report
            ),
            render=TRUE,
            include_files = c("intensityFile.txt",'annotationFile.csv')
          )}
        
        variables$filepath <- tmp_file #Assigning in the temp file where the .pdf is located to the reactive file created above
        print(variables$filepath)
        
      })
    
    # # Hide download button until report is generated
    # output$reportbuilt <- reactive({
    #   return(!is.null(report$filepath))
    # })
    # outputOptions(output, 'reportbuilt', suspendWhenHidden= FALSE)
    # 
    #Download report  
    output$DownloadReport <- downloadHandler(
      
      # This function returns a string which tells the client
      # browser what name to use when saving the file.
      filename = function() {
        paste0(input$project_name,"-report-", gsub(" |:","-",Sys.time()),".zip")
      },
      
      # This function should write data to a file given to it by
      # the argument 'file'.
      content = function(file) {
        
        file.copy(getDataPath(variables$filepath), file, overwrite = T)
        
      }
    )
   #outputOptions(output, 'download', suspendWhenHidden= FALSE)
    
    output$downloadButtonDownloadReport <- renderUI({
      if(!is.null(variables$filepath)) {
        downloadButton("DownloadReport", "Download Report")}
    })


})
