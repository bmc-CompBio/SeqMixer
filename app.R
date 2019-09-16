




library(shiny)
library(shinydashboard)
library(gdata)



#################################################################################################################################

ui <- dashboardPage(
      
      skin = "black",

      dashboardHeader(title = "SeqMixer",
                      titleWidth = 250,

                      tags$li(class = "dropdown",
                              tags$style(".main-header {max-height: 80px}"),
                              tags$style(".skin-black .main-header .logo {
                                         height: 80px;
                                         font-size: 32px;
                                         font-weight: bold;
                                         color: white;
                                         background-color: #1e282c;
                                         padding-top: 15px;
                                         }"),
                            tags$style(".skin-black .main-header .logo:hover {background-color: #1e282c;"),
                            tags$style(".skin-black .sidebar-toggle {
                                       height: 80px;
                                       font-size: 20px;
                                       padding: 25px;
                            }"),
                            tags$style(".navbar {min-height:80px}")
                            )
                      ),
      dashboardSidebar(width = 250,
            tags$style(".left-side, .main-sidebar {padding-top: 100px}"),
            
            column(width = 12, align = "center", offset = 0,
                   helpText("SeqMixer reshuffles DNA sequences with synonymous substitutions")
            ),
            
            column(width = 12, br()),
            
            # fileInput(inputId = "my_data", 
            #           label = "Upload Sequence File (txt)", 
            #           accept = c("txt")),
            
            #column(width = 12, br()),
            
            column(width = 12, align = "left", offset = 0.5,
                   h4("Rules")),
            
            column(width = 12, align = "left", offset = 0,
                   HTML("<font color=\"#6E6E6E\">
                        <ul>
                        <li>Paste plain text</li>
                        <li>Use only capital GATC letters</li>
                        <li>Sequence has to be in frame</li>
                        <ul/>
                        </font>"
                   )
            ),
            
            column(width = 12, br()),
            
            column(width = 12, align = "left", offset = 0.5,
                   h4("Summary")),
            
            column(width = 12, br()),
            
            column(6, 
                  valueBoxOutput("total", width = 32)),
            
            column(6, 
                  valueBoxOutput("changed", width = 32)),
            

            
            column(width = 12, br())
            
 
      ),
      dashboardBody(
            
            tags$head(tags$style(HTML("pre { white-space: pre-wrap; word-break: keep-all; }"))),
            
            shiny::tags$head(shiny::tags$style(shiny::HTML(
                  "#text { font-family: Courier New; height: 500px; overflow: auto; }"))),


            
            fluidRow(
                  box(title = "Input", solidHeader = T, status = "primary", width = 6, height = 562, 
                      
                          textAreaInput("my_data", label = "Enter Text", height = 462,
                                        value = "ATGGAGGAGGAGCCTGTTGCCATTGAAAA")
                          
                      #span(textOutput(outputId = "textout1"), style = "font-family: Courier New;")
                      
  
                      ),
                  box(title = "Output", solidHeader = T, status = "primary", width = 6, #height = 500,
                      
                      
                      div(id = "text",
                          htmlOutput(outputId = "textout2")
                          )
                  )
                  
            ) #,
            # fluidRow(
            #       # box(title = "Summary", solidHeader = T, status = "primary", width = 6, #height = 500, 
            #       #     
            #       #     span(textOutput(outputId = "textout3"), style = "font-family: Courier New;")
            #       # )
            #       
            #       #infoBoxOutput("summary", width = 6)
            # )
      )
)



#################################################################################################################################



server <- function(input, output) {
      
      
      my_seq_input <- reactive({
            
            # inFile <- input$my_data
            # 
            # if (is.null(inFile)){
            #       return("ATGGAGGAGGAGCCTGTTGCCATTGAAAA")
            # }
            # 
            # my_seq_input <- readChar(inFile$datapath, file.info(inFile$datapath)$size)
            # 
            my_seq_input <- input$my_data
            
            gsub("[^GATC]+", "", my_seq_input)

      })
      
      
      my_mixer <- function(){
            
            
            codon_table <- read.xls("CodonUsageGenScriptDrosphila.xlsx", sheet = 1, header = TRUE, row.names = NULL, stringsAsFactors = FALSE, na.strings = "") 
            
            
            my_seq_input <- my_seq_input()
            

            
            my_seq_split <- sapply(seq(from=1, to=nchar(my_seq_input), by=3), function(i) substr(my_seq_input, i, i+2))
            
            
            for(i in seq_along(my_seq_split)){
                  codon_table_sub <- codon_table[codon_table$AminoAcid == codon_table[codon_table$Triplet == my_seq_split[i],]$AminoAcid,]
                  
                  
                  if(nrow(codon_table_sub) > 1){
                        codon_table_sub <- codon_table_sub[codon_table_sub$Triplet != my_seq_split[i],]  
                        
                        my_seq_split[i] <- codon_table_sub$Triplet[which.max(codon_table_sub$Fraction)]
                  }
                  
            }
            
            my_seq_out <- paste(my_seq_split, collapse = "")
            

            # my_seq_out_col <- paste0(ifelse(strsplit(my_seq_input,"")[[1]] != strsplit(my_seq_out, "")[[1]],  "<font color=\"#FF0000\">", "<font color=\"#000000\">"),
            #                      strsplit(my_seq_out, "")[[1]])
            
            my_seq_out_col <- paste0(ifelse(strsplit(my_seq_input,"")[[1]] != strsplit(my_seq_out, "")[[1]],  "<font color=\"#FF0000\">", "<font color=\"#000000\">"),
                                     strsplit(my_seq_out, "")[[1]],#
                                     ifelse(strsplit(my_seq_input,"")[[1]] != strsplit(my_seq_out, "")[[1]],  "</font>", "</font>"))
            

            my_seq_out_col <- paste(my_seq_out_col, collapse = "")

            #my_seq_out_split <- sapply(seq(from=1, to=nchar(my_seq_out_col), by=230), function(i) substr(my_seq_out_col, i, i+229))
            my_seq_out_split <- sapply(seq(from=1, to=nchar(my_seq_out_col), by=300), function(i) substr(my_seq_out_col, i, i+299))
            
            #my_seq_out_split <- sapply(seq(from=1, to=nchar(my_seq_out), by=10), function(i) substr(my_seq_out, i, i+9))
            
            paste(my_seq_out_split, collapse = "\n")

   
      }
      
      my_unmixed <- function(){
            
            my_seq_input <- my_seq_input()
            
            my_seq_out2 <- sapply(seq(from=1, to=nchar(my_seq_input), by=10), function(i) substr(my_seq_input, i, i+9))
            
            paste(my_seq_out2, collapse = "\n")
            
      }
      
      
      my_total <- function(){
            
            my_seq_out <- my_unmixed()
            my_seq_out <- gsub("[^GATC]+", "", my_seq_out)
            
            nchar(my_seq_out)

      }
      
      
      my_changed <- function(){
            
            my_seq_out <- my_unmixed()
            my_seq_out <- gsub("[^GATC]+", "", my_seq_out)
            
            
            my_seq_out2 <- my_mixer()
            my_seq_out2 <- gsub("[^GATC]+", "", my_seq_out2)
            
            sum(strsplit(my_seq_out,"")[[1]] != strsplit(my_seq_out2, "")[[1]])
      }
      
      
      output$textout1 <- renderText(my_unmixed())
      
      output$textout2 <- renderText(my_mixer())
      
      
      #output$textout3 <- renderText(my_summary())
      
      output$total <- renderValueBox({
            valueBox(value = tags$p("Total", style = "font-size: 18px;"), 
                     subtitle = tags$p(my_total(), style = "font-size: 18px;"), 
                     color = "blue"
            )
      })
      
      output$changed <- renderValueBox({
            valueBox(value = tags$p("Changed", style = "font-size: 18px;"), 
                     subtitle = tags$p(my_changed(), style = "font-size: 18px;"), 
                     color = "red"
            )
      })
      
      
      
}







#################################################################################################################################


shinyApp(ui, server)





