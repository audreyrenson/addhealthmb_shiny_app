
# Setup -------------------------------------------------------------------

library(shiny)
source("helpers.R")

# Actual app --------------------------------------------------------------

ui <- fluidPage(

    titlePanel("Results: Taxa-biomarker associations"),

    sidebarLayout(
        sidebarPanel(
            width = 2,
            sliderInput(inputId = "lfc_cutoff",
                        label = "Log-fold-change cutoff (absolute value)",
                        min = 0,
                        max = 2.5,
                        value = 0.8, 
                        step=0.1),
            sliderInput(inputId = "sval_cutoff",
                       label = "-log10(s-value) cutoff",
                       min = 0, 
                       max = 10,
                       value = 2, 
                       step = 0.5),
            uiOutput("selectTaxa"),
            checkboxInput(inputId = "label_taxa_yn",
                          label = "Label taxa?", 
                          value = TRUE),
            selectInput(inputId = "select_biomarkers", 
                        label = "Focal biomarkers",
                        choices = c("All", "DBS", "Hallmark", "Immunological"), 
                        selected = "All"),
            radioButtons(inputId = "sort_taxa_by", 
                         label = "Sort by", 
                         choices = c("appearance","alphabetic"), 
                         selected = "appearance"),
            checkboxInput(inputId = "collapse_biomarkers",
                          label = "Collapse biomarker names if not labeled", value=TRUE),
            radioButtons(inputId = "descriptive_plot_type", 
                         label = "Descriptive plot type", 
                         choices = c("density", "histogram")),
            checkboxInput(inputId = "nonzero_only", 
                          label = "Show only nonzero abundances in descriptive plot")
        ),

        mainPanel(
            tabsetPanel(
                
                tabPanel(
                    title = "Manhattan plot",
                    "Contents",
                    plotOutput("manhattan_plot")
                ),
                
                tabPanel(
                    title = "Volcano plot",
                    
                    plotOutput("volcano_plot")
                ),
                
                tabPanel(
                    title = "Taxon Information",
                    h2(textOutput("taxon_heading")),
                    fluidRow(
                        column(
                            width=5, 
                            plotOutput("distribution_plot")
                        ),
                        column(
                            width=3, 
                            plotOutput("prevalence_piechart")
                        ),
                        column(
                            h4("Taxonomic information"),
                            width=2,
                            tableOutput("taxonomy_table")
                        )
                    )
                )
            )
        )
    )
)

# Define server logic
server <- function(input, output) {

    getFilteredResultsData <- reactive({
        out_data <- full_raw_data %>% 
            make_key_variables() %>%
            label_by_lfc_sval(
                lfc_cutoff = input$lfc_cutoff,
                sval_cutoff = input$sval_cutoff
            ) %>%
            mutate(biomarker_group = make_biomarker_groups(biomarker_label)) %>%
            filter_biomarker_groups(input$select_biomarkers) %>%
            collapse_biomarker_labels(
                collapse = input$collapse_biomarkers
            )
        if(!input$label_taxa_yn) {
            out_data <- out_data %>% 
                mutate(taxon_label = NA)
        }
        out_data
    })
    
    
    output$selectTaxa <- renderUI ({
        taxa_list <- get_sorted_list_of_labelled_taxa(
            getFilteredResultsData(), 
            sort_by = input$sort_taxa_by
        )
        
        selectInput(inputId = "labelled_taxa",
                    label = "Focal taxon", 
                    choices = c("All", taxa_list),
                    multiple = FALSE,
                    selected = "All")
    })
    
    output$volcano_plot <- renderPlot({
        
        volcano_data <- getFilteredResultsData() %>%
            drop_taxa_labels_except(taxa_names = input$labelled_taxa) 
        
        figure_panel_A <- volcano_data %>% 
            make_volcano_plot(
                biomrkr_source = "gene",
                legend_position = c(1,1)
            )
        figure_panel_B <- volcano_data %>%
            make_volcano_plot(
                biomrkr_source = "dbs", 
                legend_position = c(0.7, 1)
            )
        
        figure_panel_A + ggtitle("Gene expression biomarkers") | 
            figure_panel_B + ggtitle("Dried blood spot biomarkers")
        
    }, height=600, width=1200, res=120)
    
    output$manhattan_plot <- renderPlot({
        
        manhattan_data <- getFilteredResultsData()
        
        manhattan_plot <- manhattan_data %>%
            drop_taxa_labels_except(taxa_names = input$labelled_taxa) %>% 
            make_manhattan_plot()
        
        labeled_tree_plot <- tree %>% 
            format_and_collapse_tree_phyla(taxonomy) %>%
            make_tree_plot() %>%
            add_labels_to_tree_plot()
        
        (manhattan_plot / labeled_tree_plot) + plot_layout(heights = c(4,1), design = TREE_SQUISHER)
    }, height=600, width=1200, res=120)
    
    output$taxon_heading <- renderText({
        paste0("Info about ", input$labelled_taxa)
    })
    
    output$distribution_plot <- renderPlot({
        if(input$descriptive_plot_type == "histogram") {
            plot_data <- if(input$nonzero_only) {
                histogram_nonzero_data 
            } else {
                histogram_data   
            }
            distribution_plot <- make_histogram_plot(
                histogram_data = plot_data,
                taxon = input$labelled_taxa
            )
        } else {
            plot_data <- if(input$nonzero_only) {
                density_nonzero_data 
            } else {
                density_data   
            }
            distribution_plot <- make_density_plot(
                density_data = plot_data,
                taxon = input$labelled_taxa
            )
        }
        
        distribution_plot
    }, height=600, width=600, res=120)
    
    output$prevalence_piechart <- renderPlot({
        getFilteredResultsData() %>%
            make_pie_data(taxon = input$labelled_taxa) %>%
            make_pie_chart()
    }, height=300, width=300, res=120)
    
    output$taxonomy_table <- renderTable({
        get_taxonomy_table(glommed_taxonomies, taxon = input$labelled_taxa)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
