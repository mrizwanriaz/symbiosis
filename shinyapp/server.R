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

high_to_low_colors <-
  rev(
    c(
      "#419460",
      "#419462",
      "#419467",
      "#419469",
      "#5F966B",
      "#65986A",
      "#65986A",
      "#77996B",
      "#A59D6A",
      "#B59D6A",
      "#BA9D69",
      "#BE9E69",
      "#EC9E69",
      "#EEAD6B",
      "#DE9B78",
      "#CE8084",
      "#B15696",
      "#B15697",
      "#B15698",
      "#B155A5"
    )
  )

# Define the color palette for the heatmap and barplots
my_color <-
  colorpanel(50, "#0077bb", "white", "#EE7733") 
bar_color <- high_to_low_colors
#c(rep("grey", 20))


############# DATA ################
strain_sort <- readRDS("high_to_low.RDS")

# Sinorhizobium meliloti ############
sm_expr <- readRDS("smel_expr.RDS")
sm_data <- readRDS("smel_data2.RDS")
rhizobia_sym_genes <-
  readRDS("rhizobia_symgenes.RDS") #example gene list of rhizobia

#gene list for dropdown menu
dd_genes <-
  c(sm_expr$Expr_Id, sm_expr$LocusTag) %>% unique() %>% sort()

# Medicago truncatula ############
md_data <-
  readRDS("med_data.RDS") #%>% filter(BioType == "mRNA") #only mRNAs
md_expr <-
  readRDS("med_expr.RDS") %>% filter(LocusTag %in% md_data$LOCUS_TAG)
medicago_sym_genes <-
  readRDS("medicago_symgenes.RDS") #example gene list of medicago



################## Define server logic for the Shiny app ##############3
server <- function(session, input, output) {
  #hide download buttons
  shinyjs::hide("download_multi")
  shinyjs::hide("download_mod")
  
  #hide high/low legends
  shinyjs::hide("highlow")
  shinyjs::hide("highlow2")
  
  ######## SINGLE GENE ##############
  
  # Update gene dropdown menu
  observeEvent(input$sel_specie, {
    switch (
      input$sel_specie,
      "smel" = updateSelectizeInput(
        session,
        "gene_select",
        choices = c("", dd_genes),
        server = TRUE,
        selected = "SM_RS27005"
      ),
      "med" = updateSelectizeInput(
        session,
        "gene_select",
        choices = c("", md_expr$LocusTag),
        server = TRUE,
        selected = "MtrunA17_Chr1g0155701"
      )
    )
    
  })
  
  #update module list
  observeEvent(input$sel_mod_specie, {
    output$mod_dto <- NULL
    output$mod_heatmap_plot <- NULL
    shinyjs::hide("download_mod")
    shinyjs::hide("highlow")
    
    switch (
      input$sel_mod_specie,
      "smel" = updateSelectInput(
        session,
        "mod_select",
        choices = paste("r-M", c(1:40), sep = ""),
        selected = "r-M1"
      ),
      "med" = updateSelectInput(
        session,
        "mod_select",
        choices = paste("p-M", c(1:65), sep = ""),
        selected = "p-M1"
      )
    )
  })
  
  
  # Show gene feature info
  output$gene_feature_info <- renderUI({
    req(input$gene_slc_btn)
    
    if (!is.null(input$gene_select) && input$gene_select != "") {
      get_gene_feature_info(input$gene_select, input$sel_specie) %>% paste(sep = "<br/>") %>% HTML()
      
    }
  })
  
  # Plot expression barplot
  observeEvent(input$gene_slc_btn, {
    output$expression_plot <- renderPlot({
      req(input$gene_slc_btn)
      
      if (!is.null(input$gene_select) && input$gene_select != "") {
        plot_gene_expression(input$gene_select, input$sel_specie)
      }
    })
  })
  
  
  # Helper function to get gene feature info
  get_gene_feature_info <- function(gene, specie) {
    if (specie == "smel") {
      gene_record <-
        sm_data %>% filter(locus_tag %in% c(gene) |
                             proteinId %in% c(gene))
      
      
      if (!is.null(gene_record)) {
        # Constructing feature_info string
        feature_info <- paste(
          bold("Locus Tag:"),
          gene_record[, "locus_tag"],
          "<br/>",
          bold("Refseq Protein:"),
          paste(
            "<a href='https://www.ncbi.nlm.nih.gov/protein/",
            gene_record[, "proteinId"],
            "' target=_BLANK>",
            gene_record[, "proteinId"],
            "</a>",
            sep = ""
          ),
          "<br/>",
          bold("Name:"),
          gene_record["name"],
          "<br/>",
          bold("Location:"),
          ifelse(
            gene_record["seq_type"] != "",
            paste(
              gene_record["seq_type"],
              " (",
              gene_record["strand"],
              ") [",
              gene_record["start"],
              "..",
              gene_record["end"],
              "]",
              sep = ""
            ),
            ""
          ),
          "<br/>",
          bold("Feature length:"),
          gene_record["feature_interval_length"],
          "<br/>",
          bold("Product length:"),
          gene_record["product_length"],
          "<br/><br/>",
          bold("WGCNA Module:"),
          gene_record["WGCNA_module"],
          "<br/>",
          bold("Module Membership:"),
          formatValue(as.numeric(gene_record["WGCNA_MM"]), 0.6),
          "(p =",
          formatC(as.numeric(gene_record["WGCNA_p.MM"]), digits = 2),
          ") <br/>",
          bold("Gene Significance (shoot biomass):"),
          formatValue(as.numeric(gene_record["WGCNA_GS.shoot"]), 0.6),
          "(p =",
          formatC(as.numeric(gene_record["WGCNA_p.GS.shoot"]), digits = 2),
          ")",
          "<br/><br/><i style=float:right;>Strains are sorted from <b style=color:#B155A5; title='high plant shoot biomass';>high</b> to <b style=color:#419460; title='low plant shoot biomass';>low</b> quality partners</i>"
        )
        
      } else {
        feature_info <- "No feature information available."
      }
    }
    else{
      #medicago
      gene_record <- md_data %>% filter(LOCUS_TAG %in% c(gene))
      
      
      if (!is.null(gene_record)) {
        # Constructing feature_info string
        feature_info <- paste(
          bold("Locus Tag:"),
          paste(
            "<a href='https://www.ncbi.nlm.nih.gov/gene/?term=",
            gene_record[, "LOCUS_TAG"],
            "' target=_BLANK>",
            gene_record[, "LOCUS_TAG"],
            "</a>",
            sep = ""
          ),
          "<br/>",
          bold("ACRONYM:"),
          gene_record["SYMBOL"],
          "<br/>",
          bold("Name:"),
          gene_record["PRODUCT"],
          "<br/><br/>",
          bold("WGCNA Module:"),
          gene_record["WGCNA_module"],
          "<br/>",
          bold("Module Membership:"),
          formatValue(as.numeric(gene_record["WGCNA_MM"]), 0.6),
          "(p =",
          formatC(as.numeric(gene_record["WGCNA_p.MM"]), digits = 2),
          ") <br/>",
          bold("Gene Significance (shoot biomass):"),
          formatValue(as.numeric(gene_record["WGCNA_GS.shoot"]), 0.6),
          "(p =",
          formatC(as.numeric(gene_record["WGCNA_p.GS.shoot"]), digits = 2),
          ")",
          "<br/><br/><i style=float:right;>Strains are sorted from <b style=color:#B155A5; title='high plant shoot biomass';>high</b> to <b style=color:#419460; title='low plant shoot biomass';>low</b> quality partners</i>"
        )
        
      } else {
        feature_info <- "No feature information available."
      }
    }
    
    return(feature_info)
  }
  
  
  # Helper function to plot gene expression
  plot_gene_expression <- function(gene, specie) {
    if (specie == "smel") {
      gene_expr <- sm_expr %>%
        subset(Expr_Id %in% c(gene) | LocusTag %in% c(gene)) %>%
        select(all_of(strain_sort)) %>% t() %>% data.frame() #aveLogCPM smel
    } else{
      gene_expr <- md_expr %>%
        subset(LocusTag %in% c(gene)) %>%
        select(all_of(strain_sort)) %>% t() %>% data.frame() #aveLogCPM med
    }
    
    
    colnames(gene_expr) <- "expr"
    
    
    ggplot(data = gene_expr, aes(
      x = rownames(gene_expr),
      y = expr,
      fill = bar_color
    )) +
      xlab('Strains') +
      ylab("log2(Average CPM)") +
      geom_bar(stat = "identity",
               color = "black",
               width = 0.8) +
      scale_fill_identity() +
      #geom_text(aes(label = formattable(expr,format="f",digits=2),vjust = 0))+
      scale_x_discrete(limits = rownames(gene_expr)) +
      guides(fill = "none") +
      theme_classic(
        base_size = 18,
        base_line_size = 0.2,
        base_family = "aerial"
      ) +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.2,
          hjust = 1,
          size = 14
        ),
        axis.title = element_text(
          size = 16,
          face = "bold",
          family = "sans"
        ),
        plot.title = element_text(hjust = 0.5,),
        line = element_blank()
      )
  }
  
  
  ############# MULTIPLE GENES #########################
  
  # Plot heatmap of multiple genes
  observeEvent(input$plot_btn, {
    if (!is.null(input$gene_list) && input$gene_list != "") {
      gene_list <- strsplit(input$gene_list, "\n")[[1]]
      
      if (input$sel_multi_specie == "smel") {
        #rhizobia
        #expression heatmap
        output$heatmap_plot <- renderPlot({
          sm_expr %>% subset(Expr_Id %in% gene_list |
                               LocusTag %in% gene_list) %>% plot_heatmap()
          
        })
        
        #data table
        thedata <-
          sm_data %>% filter(locus_tag %in% c(gene_list) |
                               proteinId %in% c(gene_list)) %>% select(c(1:3, 8, 13:17)) %>% reactive()
        
      }
      else{
        #medicago
        
        output$heatmap_plot <- renderPlot({
          #expression
          md_expr %>% subset(LocusTag %in% gene_list) %>% plot_heatmap()
          
        })
        
        thedata <-
          md_data %>% filter(LOCUS_TAG %in% gene_list) %>% select(c(1, 4, 8:12)) %>% reactive()
        
        
      }
      
      output$multi_dto <- renderDataTable({
        thedata()
      })
      
      #show download button
      shinyjs::show("download_multi")
      
      #show high/low legend
      shinyjs::show("highlow2")
      
      
      #download data
      output$download_multi <- downloadHandler(
        filename = function() {
          "genes_data.csv"
        },
        content = function(fname) {
          write.csv(thedata(), fname)
        }
      )
      
      
    }
    else{
      #show modal dialog
      show_message("Please provide a list of genes.")
    }
  })
  
  #clear text field
  observeEvent(input$clear_txt, {
    updateTextAreaInput(session, "gene_list", value = "")
    output$multi_dto <- NULL
    output$heatmap_plot <- NULL
    shinyjs::hide("download_multi")
    #hide high/low legend
    shinyjs::hide("highlow2")
    
  })
  
  
  #load example Ids in the text field
  observeEvent(input$load_exmpl, {
    switch(
      input$sel_multi_specie,
      "smel" = updateTextAreaInput(session, "gene_list", value = gsub(
        ", ", "\n", toString(rhizobia_sym_genes)
      )),
      "med" = updateTextAreaInput(session, "gene_list", value = gsub(
        ", ", "\n", toString(medicago_sym_genes)
      ))
    )
  })
  
  #when user changes specie
  observeEvent(input$sel_multi_specie, {
    updateTextAreaInput(session, "gene_list", value = "")
    output$multi_dto <- NULL
    output$heatmap_plot <- NULL
    shinyjs::hide("download_multi")
    #hide high/low legend
    shinyjs::hide("highlow2")
  })
  
  
  
  
  ############# MODULE GENES #########################
  
  # Plot heatmap of module genes
  observeEvent(input$mod_slc_btn, {
    #req(input$mm_cutoff)
    
    #expression heatmap
    if (!is.null(input$mod_select) && input$mod_select != "") {
      if (input$sel_mod_specie == "smel") {
        #rhizobia
        mod <- gsub("r-M", "", input$mod_select)
        
        mod_data <- sm_data %>%
          subset(WGCNA_module == mod) %>% filter(
            between(
              as.numeric(WGCNA_MM),
              input$mm_cutoff[1],
              input$mm_cutoff[2]
            ),
            between(
              as.numeric(WGCNA_GS.shoot),
              input$gs_cutoff[1],
              input$gs_cutoff[2]
            )
          )
        
        output$mod_heatmap_plot <- renderPlot({
          #expression
          sm_expr %>% subset(Expr_Id %in% mod_data$locus_tag |
                               LocusTag %in% mod_data$locus_tag) %>% plot_heatmap()
          
        })
        
        thedata <- mod_data[, c(1:3, 8,  13:17)] %>% reactive()
        
      }
      else{
        #medicago
        mod <- gsub("p-M", "", input$mod_select)
        
        mod_data <-
          md_data %>%
          subset(WGCNA_module == mod) %>% filter(
            between(
              as.numeric(WGCNA_MM),
              input$mm_cutoff[1],
              input$mm_cutoff[2]
            ),
            between(
              as.numeric(WGCNA_GS.shoot),
              input$gs_cutoff[1],
              input$gs_cutoff[2]
            )
          )
        
        output$mod_heatmap_plot <- renderPlot({
          #expression
          md_expr %>% subset(LocusTag %in% mod_data$LOCUS_TAG) %>% plot_heatmap()
          
        })
        
        thedata <- mod_data[, c(1, 4, 8:12)] %>% reactive()
        
      }
      output$mod_dto <- renderDataTable({
        thedata()
      })
      
      
      
      #show download button
      shinyjs::show("download_mod")
      
      
      #show high low legend
      shinyjs::show("highlow")
      
      #download data
      output$download_mod <- downloadHandler(
        filename = function() {
          paste(input$sel_mod_specie,
                "_",
                "Module",
                mod,
                "_data.csv",
                sep = "")
        },
        content = function(fname) {
          write.csv(thedata(), fname)
        }
      )
    }
  })
  
  
  # show/hide parameters for module genes
  observeEvent(input$mod_filt_btn, {
    shinyjs::toggle("mm_cutoff_row", anim = TRUE, animType = "slide")
    shinyjs::toggle("gs_cutoff_row", anim = TRUE, animType = "slide")
    
    
  })
  
  
  
  ## Helper function to plot heatmap
  plot_heatmap <- function(gene_expr) {
    gene_expr %<>% data.frame()
    rownames(gene_expr) <- gene_expr[, 1]
    
    gene_expr.filt <-
      gene_expr[rowSums(gene_expr > 0) >= 4, strain_sort] #plot expr for genes with >0 in atleast 4 samples
    
    if (nrow(gene_expr.filt) > 1) {
      heatmap.2(
        as.matrix(gene_expr.filt),
        dendrogram = "row",
        scale = "row",
        Colv = FALSE,
        Rowv = TRUE,
        key = TRUE,
        col = my_color,
        colCol = high_to_low_colors,
        density.info = "none",
        key.title = NA,
        trace = "none",
        lwid = c(2, 10),
        lhei = c(2, 9),
        cexCol = 1.3,
        srtCol = 55,
        #labRow = FALSE,
        margins = c(10, 12)
      )
    }
    else{
      show_message("No expressed gene found.")
      return(0)
      
    }
  }
  
  #helper functions
  show_message <- function(msg) {
    showModal(modalDialog(
      msg,
      easyClose = TRUE,
      title = NULL,
      footer = NULL
    ))
  }
  
  bold <- function(text) {
    paste("<b>", text, "</b>")
  }
  
  italic <- function(text) {
    paste("<i>", text, "</i>")
  }
  
  formatValue <- function(value, threshold) {
    v <- as.character(formatC(value, digits = 3))
    ifelse(abs(value) > threshold, italic(v), v)
  }
  
}
