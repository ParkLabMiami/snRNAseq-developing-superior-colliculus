
# packages ----------------------------------------------------------------

library(shiny)
library(shinythemes)
# library(ggplot2)
library(plotly)
library(hdf5r)
library(RColorBrewer)
library(patchwork)
library(ggbeeswarm)
library(dplyr)
library(reshape2)
library(shinymanager)
# library(ggrepel)

# Load
# load(file = './data/appdata.RData') # This is faster by ~10%
source(file = 'helper.R')



# Password protection -----------------------------------------------------

# To enable password protection, I followed the advice here: https://stackoverflow.com/questions/28987622/starting-shiny-app-after-password-input


inactivity <- "function idleTimer() {
var t = setTimeout(logout, 120000);
window.onmousemove = resetTimer; // catches mouse movements
window.onmousedown = resetTimer; // catches mouse movements
window.onclick = resetTimer;     // catches mouse clicks
window.onscroll = resetTimer;    // catches scrolling
window.onkeypress = resetTimer;  //catches keyboard actions

function logout() {
window.close();  //close the window
}

function resetTimer() {
clearTimeout(t);
t = setTimeout(logout, 300000);  // time is in milliseconds (1000 is 1 second)
}
}
idleTimer();"


IDs = c('kpark','jsc228','aca136','fxf201','rxm1687','klevay','mtapia2','ncm110')
allowed.users = c(paste0(IDs, '@miami.edu'),
                  paste0(IDs, '@med.miami.edu'))
password <- 'axonregen'

credentials <- data.frame(
  user = allowed.users,
  password = rep(password, length(allowed.users))
)





# UI ----------------------------------------------------------------------

ui <- secure_app(ui = fluidPage(
  
  head_auth = tags$script(inactivity),
  
#   # jQuery chunk for ensuring DOM is fully loaded before user can interact
#   tags$script(HTML(
#     '$(document).ready(function () {
#     $.getJSON("https://ipapi.co/json/", function (client) {
#         Shiny.onInputChange("client", client);
#     });
# });'
#   )),
  
  # Set layout colors/theme
  theme = shinytheme("yeti"),
  
  # App title
  titlePanel(
    fluidRow(
      column(
        width = 9, 
        div(
          h2(
            window_title, 
            style="display:inline;"
          ),
          tags$a(
            h3(
              title_link_text, 
              style="display:inline;"
            ),
            href = title_link_url, 
            target = "_blank"
          ),
          style='padding-left:10px;'
        )
      ),
      column(
        width = 3, 
        div(
          tags$a(
            h5(
              "Browser App", 
              style = "display:inline;color:black;vertical-align:middle;"
            ),
            # href = "https://github.com/JaeLeeLab/sci_scRNAseq_portal", 
            target = "_blank"
          ),
          tags$a(
            img(
              height = 24, 
              width = 24, 
              src = "GitHub-Mark-64px.png", 
              style="display:inline;vertical-align:middle;"
            ),
            # href = "https://github.com/JaeLeeLab/sci_scRNAseq_portal", 
            target="_blank"
          )
        ),
        align="right"
      )
    ),
    
    # This is to display on the headers of web browser.
    windowTitle = window_title
  ),
  
  # Create panel of tabs for different queries
  tabsetPanel(
    
    # About study panel
    tabPanel(
      title = 'About',
      br(),
      fluidRow(
        column(
          width = 12,
          br(),
          HTML(
            text = "<p style='font-size:15px'>This website accompanies the manuscript:<br><br>Ana C. Ayupe, James S. Choi, Felipe Beckedorff, Robyn McCartan, Konstantin Levay, Kevin K. Park. Single-Nucleus RNA Sequencing of Developing and Mature Superior Colliculus Identifies Neuronal Diversity and Candidate Mediators of Circuit Assembly. <i>bioRxiv</i> [INSERT DATE OF PUB]; <a href='https://parklabmiami.com'>DOI: [INSERT DOI]</a></p>"
          ),
          HTML(
            text = '<br><p><b>Click the "Gene expression" tab to start exploring the data presented in the paper.</b></p>'
          ),
          HTML(
            text = '<br><h4>Summary:</h4>The superior colliculus (SC) is a sensorimotor structure in the midbrain that integrates inputs from multiple sensory modalities to initiate motor commands. It undergoes well-characterized steps of circuit assembly during development, rendering the mouse SC a popular model to study establishment and refinement of neural connectivity. Here we performed single nucleus RNA-sequencing analysis of the mouse SC isolated at various developmental time points. Our study provides a transcriptomic landscape of the cell-types that comprise the SC across murine development with particular emphasis on neuronal heterogeneity. We used these data to identify Pax7 as a marker for an anatomically homogeneous population of GABAergic neurons. Lastly, we report a repertoire of genes differentially expressed across the different postnatal ages, many of which are known for regulating axon guidance and synapse formation. Our data provide a valuable resource for interrogating the mechanisms of circuit development, and identifying markers for manipulating specific neuronal populations and circuits in the SC.'
          ),
          HTML(
            text = '<h4>Methods:</h4><p>To profile gene expression of SC during development and maturation, we analyzed the transcriptomes of single nuclei across four time points: embryonic day 19 (E19), postnatal day 4 (P4), P8 and P21. We selected these time points because together they encompass discrete yet overlapping developmental events including axonal outgrowth, axonal targeting, topographic mapping, synaptogenesis, oligodendrocyte differentiation, myelination, and synaptic refinement and maturation. For each age, we micro-dissected the SC territory from several animals, which were pooled. After detergent-based digestion and mechanical dissociation followed by nuclei isolation using a sucrose density gradient centrifugation, the single nucleus suspension was processed and sequenced using the 10X Genomics droplet-based snRNA-seq platform. Downstream clustering and gene expression analyses were done using Seurat and the Bioconductor suite of bioinformatics tools.'
          ),
          style='padding-left:10px;'
        ),
      ),
      hr(),
      fluidRow(
        h3(
          "Data availability", 
          style="display:inline;"
        ),
        br(),
        br(),
        column(
          width = 12,
          HTML(
            text = "<p>Code used to analyze the snRNA-seq data presented in <i>Single-Nucleus RNA Sequencing of Developing and Mature Superior Colliculus Identifies Neuronal Diversity and Candidate Mediators of Circuit Assembly</i> are available on <a href='https://parklabmiami.com'>Github.</a> Direct download of SeuratObjects can be found in the repo.</p><p>Raw sequencing data are available from the SRA (Sequence Read Archive) database under study <a href='https://www.ncbi.nlm.nih.gov/sra'>[INSERT SRA ACCESSION]</a>. Gene-count matrices are available from the Gene Expression Omnibus under accession <a href='https://www.ncbi.nlm.nih.gov/geo/'>[INSERT GEO ACCESSION]</a>. Relevant sample-level and cell-level metadata are available under the GEO accession metadata.</p>"
          )
        ),
        hr(),
        style='padding-left:10px;'
      ),
      style='padding-left:20px; padding-right:20px'
    ),
    # Tab to query features 
    tabPanel(
      title = 'Gene expression',
      # shinyjs::useShinyjs(),
      sidebarLayout(
        fluid = TRUE,
        # `choices = NULL` so that choices are set in server-side
        sidebarPanel(
          tags$head(tags$script(
          'var dimension = [0, 0];
          $(document).on("shiny:connected", function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange("dimension", dimension);
          });
          $(window).resize(function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange("dimension", dimension);
          });')),
          tags$style(type = 'text/css', ".selectize-input { font-size: 16px; line-height: 22px;} .selectize-dropdown { font-size: 16px; line-height: 22px; }"),
          width = 2,
          selectizeInput(
            inputId = "selected_dataset",
            choices = NULL,
            label = 'Select dataset:'
          ),
          selectizeInput(
            inputId = "selected_feature", 
            choices = NULL, 
            options = list(placeholder = 'Select a gene'),
            label = 'Select gene:'
          ),
          selectizeInput(
            inputId = "selected_groupby",
            choices = NULL,
            label = 'Group cells by:'
          ),
          checkboxInput(
            inputId = 'draw_labels',
            label = 'Overlay group labels',
            value = TRUE
          ),
          HTML(
            text = "
            <p style='font-size:12px'><b>How to use:</b>
            <br>
            To get started, select a dataset from the \"Select dataset\" drop-down menu. Select genes of interest from the \"Select Gene\" drop-down menu. Alternatively, start typing your gene of interest for matching items.
            </p>
            <br>
            <p style='font-size:12px'>
            <b>Drop-down legend:</b>
            Cells can be grouped and colored by celltype, subtype, time-point of collection, and other metadata from the \"Group cells by\" menu.
            <ul style=\"padding-left:10px\">
            <li><span style='font-family:courier'>SuperiorColliculus</span>: all cells from study</li>
            <li><span style='font-family:courier'>Neurons</span>: neuronal cells only</li>
            <li><span style='font-family:courier'>time</span>: development time-point of tissue collection</li>
            <li><span style='font-family:courier'>celltype</span>: celltype class</li>
            <li><span style='font-family:courier'>subtype</span>: neuronal subtypes and further glial subtypes not reported in manuscript</li>
            <li><span style='font-family:courier'>integrated_snn_res.0.8</span>: original SuperiorColliculus cluster output by Seurat</li>
            <li><span style='font-family:courier'>SingleR_Transseq</span>: SingleR-based prediction of Trans-seq subtypes in SuperiorColliculus dataset</li>
            <li><span style='font-family:courier'>SingleR_Vectorseq</span>: SingleR-based prediction of Vector-seq subtypes in SuperiorColliculus dataset</li>
            </ul>
            </p>"
          )
        ),
        mainPanel(
          fluid = TRUE, 
          width = 10,
          fluidRow(
            width = 12,
            column(
              width = 12,
              plotOutput(
                outputId = 'jointplot',
                width = '100%',
                height = '500px'
              )
            )
          ),
          fluidRow(
            width = 12,
            plotOutput(
              outputId = 'splitfeatureplot',
              width = 'auto',
              height = '350px'
            )
          ),
          fluidRow(
            # column(
            #   width = 6,
            #   # plotOutput(
            #   #   outputId = 'featuredotplot',
            #   #   height = '550px'
            #   # )
            # ),
            column(
              width = 12,
              plotOutput(
                outputId = 'featureviolinplot',
                height = '350px'
              )
            ),
          )
        )
      )
    ),
    tabPanel(
      title = 'Multiple genes',
      sidebarLayout(
        fluid = TRUE,
        sidebarPanel = sidebarPanel(
          width = 2,
          tags$head(tags$script(
          'var dimension = [0, 0];
          $(document).on("shiny:connected", function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange("dimension", dimension);
          });
          $(window).resize(function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange("dimension", dimension);
          });')),
          selectizeInput(
            inputId = "selected_dataset_multiple",
            choices = NULL,
            label = 'Select dataset:'
          ),
          selectizeInput(
            inputId = "selected_feature_multiple", 
            choices = NULL, 
            options = list(placeholder = 'Select a gene'),
            multiple = TRUE,
            label = 'Select gene:'
          ),
          selectizeInput(
            inputId = "selected_groupby_multiple",
            choices = NULL,
            label = 'Group cells by:'
          ),
          actionButton(
            inputId = 'submit_multiple',
            label = 'Submit'
          ),
          HTML(
            text = "
            <p style='font-size:12px'><b>How to use:</b>
            <br>
            Type multiple genes of interest into the \"Select Gene\" box. Click \"Submit\" to update dot plots (note: dot plots may reload without clicking \"Submit\" but will not reflect new query list).
            </p>
            <br>
            <p style='font-size:12px'>
            <b>Drop-down legend:</b>
            Cells can be grouped and colored by celltype, subtype, time-point of collection, and other metadata from the \"Group cells by\" menu.
            <ul style=\"padding-left:10px\">
            <li><span style='font-family:courier'>SuperiorColliculus</span>: all cells from study</li>
            <li><span style='font-family:courier'>Neurons</span>: neuronal cells only</li>
            <li><span style='font-family:courier'>time</span>: development time-point of tissue collection</li>
            <li><span style='font-family:courier'>celltype</span>: celltype class</li>
            <li><span style='font-family:courier'>subtype</span>: neuronal subtypes and further glial subtypes not reported in manuscript</li>
            <li><span style='font-family:courier'>integrated_snn_res.0.8</span>: original SuperiorColliculus cluster output by Seurat</li>
            <li><span style='font-family:courier'>SingleR_Transseq</span>: SingleR-based prediction of Trans-seq subtypes in SuperiorColliculus dataset</li>
            <li><span style='font-family:courier'>SingleR_Vectorseq</span>: SingleR-based prediction of Vector-seq subtypes in SuperiorColliculus dataset</li>
            </ul>
            </p>"
          )
        ),
        mainPanel = mainPanel(
          fluid = TRUE, 
          width = 10,
          fluidRow(
            column(
              width = 12,
              plotOutput(
                outputId = 'multiplefeaturedotplot',
                width = '100%',
                height = '500px'
              )
            )
          ),
          fluidRow(
            column(
              width = 12,
              plotOutput(
                outputId = 'featuresplitdotplot',
                width = '100%',
                height = '600px'
              )
            )
          )
        )
      )
    )
  )
))


# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  
  #Logging
  observeEvent(
    eventExpr = {
      input$client 
    },
    handlerExpr = {
      logging::loginfo("New client with ip: %s", input$client$ip)
    },
    ignoreNULL = TRUE, 
    ignoreInit = FALSE
  )
  
  result_auth <- secure_server(check_credentials = check_credentials(credentials))
  output$res_auth <- renderPrint({
    reactiveValuesToList(result_auth)
  })
  
  
  # Initialize dataset
  updateSelectizeInput(
    session = session, 
    label = 'Select dataset:',
    inputId = "selected_dataset",
    choices = names(dataset_dict),
    selected = 'Neurons'
  )
  
  # Initialize cell groupby
  updateSelectizeInput(
    session = session,
    inputId = "selected_groupby", 
    label = "Group cells by:",
    choices = categoricalVars,
    selected = 'subtype'
  )
  
  # Initialize selected feature
  updateSelectizeInput(
    session = session,
    inputId = "selected_feature", 
    label = "Select Gene:",
    choices = all_features,
    server = TRUE,
    selected = 'Etv1'
  )
  
  # Initialize dataset multiple
  updateSelectizeInput(
    session = session, 
    label = 'Select dataset:',
    inputId = "selected_dataset_multiple",
    choices = names(dataset_dict),
    selected = 'Neurons'
  )
  
  # Initialize cell groupby multiple
  updateSelectizeInput(
    session = session,
    inputId = "selected_groupby_multiple", 
    label = "Group cells by:",
    choices = categoricalVars,
    selected = 'subtype'
  )
  
  # Initialize selected feature multiple
  updateSelectizeInput(
    session = session,
    inputId = "selected_feature_multiple", 
    label = "Select Gene:",
    choices = all_features,
    server = TRUE,
    selected = c('Etv1','Slc17a6')
  )
  
  # Upon change to selected_dataset, take subset of expression matrix
  dataset_value <- eventReactive(
    eventExpr = {
      input$selected_dataset
    },
    valueExpr = {
      dataset_value <- dataset_dict[input$selected_dataset]
      return(dataset_value)
    },
    ignoreInit = FALSE,
    ignoreNULL = TRUE
  )
  
  # Joined featureplot + dimplot
  observeEvent(
    eventExpr = {
      c(input$selected_dataset, input$selected_groupby, input$selected_feature, input$draw_labels)
    },
    handlerExpr = {
      tmp_group <- req(input$selected_groupby)
      tmp_feature <- req(input$selected_feature)
      logging::loginfo("loaded dataset %s with dimplot cells group by %s.", 
                       dataset_value(), tmp_group)
      gc(verbose = FALSE)
      p1 <- drawFeaturePlot(
        dataset = dataset_value(),
        feature = tmp_feature,
        reduction = 'UMAP'
      )
      p2 <- drawDimPlot(
        dataset = dataset_value(),
        groupby = tmp_group,
        reduction = 'UMAP',
        draw_labels = input$draw_labels
      )
      output$jointplot <- renderPlot(
        expr = {((p1 | p2) + plot_layout(guides = 'collect'))}
      )
    }
  )
  
  
  # Split feature plot - separate plot per injury time-point
  observeEvent(
    eventExpr = {
      c(input$selected_dataset, input$selected_feature)
    },
    handlerExpr = {
      tmp_feature <- req(input$selected_feature)
      output$splitfeatureplot <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          drawSplitFeaturePlot(
            dataset = dataset_value(),
            feature = tmp_feature
          )
        }
      )
    }
  )
  
  # Feature violin plot with dot rows split by `input$selected_groupby`
  observeEvent(
    eventExpr = {
      c(input$selected_dataset, input$selected_feature, input$selected_groupby)
    },
    handlerExpr = {
      tmp_feature <- req(input$selected_feature)
      tmp_groupby <- req(input$selected_groupby)
      output$featureviolinplot <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          drawFeatureViolinPlot(
            dataset = dataset_value(),
            feature = tmp_feature,
            groupby = tmp_groupby
          )
        }
      )
    }
  )
  
  # multi feature dot plot
  observeEvent(
    eventExpr = input$submit_multiple,
    handlerExpr = {
      output$multiplefeaturedotplot <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          tmp_feature <- req(input$selected_feature_multiple)
          tmp_groupby <- req(input$selected_groupby_multiple)
          drawDotPlot(
            dataset = dataset_value(),
            feature = tmp_feature,
            groupby = tmp_groupby
          )
        }
      )
    }
  )
  # output$multiplefeaturedotplot <- renderPlot(
  #   expr = {
  #         gc(verbose = FALSE)
  #         drawDotPlot(
  #           dataset_value = dataset_value(),
  #           feature = plot_feats(),
  #           groupby = group(),
  #         )
  #       }
  # )
  # 
  # observeEvent(
  #   eventExpr = {
  #     c(input$submit_multiple)
  #   },
  #   handlerExpr = {
  #     tmp_feature <- req(input$selected_feature_multiple)
  #     tmp_groupby <- req(input$selected_groupby_multiple)
  #     logging::loginfo(msg = 'dataset_value(): %s', dataset_value())
  #     logging::loginfo(msg = 'dataset_value(): %s', tmp_feature[2])
  #     logging::loginfo(msg = 'group: %s', tmp_groupby)
  #     output$multiplefeaturedotplot <- renderPlot(
  #       expr = {
  #         gc(verbose = FALSE)
  #         drawDotPlot(
  #           dataset = dataset_value(),
  #           feature = plot_feats(),
  #           groupby = tmp_groupby,
  #         )
  #       }
  #     )
  #   }
  # )
  
  # multi feature split dot plot
  observeEvent(
    eventExpr = {
      c(input$submit_multiple,
        input$selected_dataset_multiple,
        input$selected_groupby_multiple)
    },
    handlerExpr = {
      tmp_feature <- req(input$selected_feature_multiple)
      tmp_groupby <- req(input$selected_groupby_multiple)
      tmp_ncol <- switch(
        EXPR = LETTERS[which.max(input$dimension[1] > c(1008, 640, 0))],
        A = 2,
        B = 2,
        C = 1
      )
      output$featuresplitdotplot <- renderPlot(
        expr = {
          gc(verbose = FALSE)
          drawSplitDotPlot(
            dataset = dataset_value(),
            feature = tmp_feature,
            groupby = tmp_groupby,
            ncol = tmp_ncol
          )
        }
      )
    }
  )
}

shinyApp(ui = ui, server = server)

# To deploy, run the following two lines:
# Current working diretory should contain the project directory.
# setwd('..')
# rsconnect::deployApp(appDir = 'superior-colliculus-snRNAseq/', appName = 'superior-colliculus-snRNAseq', account = 'parklabmiami')
# rsconnect::accounts()
# rsconnect::accountInfo()


# To allow bioconductor packages to be sourced. Run the following in your 
# current session before pushing.
# library(BiocManager)
# options(repos = BiocManager::repositories())
