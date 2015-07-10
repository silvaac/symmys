library(shiny)
library(ggplot2)
library(reshape2)

dataset <- diamonds
source("utility.R")

server <- function(input,output) {
  xx<-reactive({
      switch(input$model,
      "bm" = BM(T=input$nYears,
                  sigma=as.numeric(input$sigma),
                  mu   =as.numeric(input$mu),
                  N    =input$sampleSize,
                  step =input$timeStep),
      "vg" = VG(T=input$nYears,
                sigma=as.numeric(input$sigma),
                mu   =as.numeric(input$mu),
                N    =input$sampleSize,
                step =input$timeStep,
                kappa=as.numeric(input$kappa))
      )
  })
  
  output$plot <- renderPlot({
    yy = xx()$value
    yy$period = xx()$period
    zz = melt(yy,id.vars = "period")
    p1<- ggplot(zz, aes(x=period, y=value, colour=variable))
    p1<- p1+geom_line() #+theme_bw()
    p1<- p1+theme(legend.position="none")+geom_hline(aes(yintercept=0))
    print(p1)
    #qplot(xx()$period,yy,geom="line")
  })
  
}

ui<-
  fluidPage(
  
  plotOutput('plot'),
  
  hr(),
  fluidRow(
    column(4,
           selectInput("model", "Model",
                      choices = c("BM" = "bm",
                                  "VG" = "vg"),
                      selected = "bm"),
           sliderInput('sampleSize', '# of Simulations', 
                       min=1, max=100, value=1, 
                       step=1, round=0)
    ),
    column(4,
           sliderInput('timeStep', 'Time unit (1/252)', 
                       min=1, max=100, value=1, 
                       step=1, round=0),
           sliderInput('nYears', '# of 252 time units', 
                       min=1, max=10, value=2, 
                       step=1, round=0)
    ),
    column(3, 
      # This outputs the dynamic UI component
      textInput("mu", "Deterministic Drift (\\mu)",value="0"),
      conditionalPanel(
        condition = "input.model == \"bm\" || input.model == \"vg\"",
        textInput("sigma", "Std dev (\\sigma)",value="0.2")
      ),
      conditionalPanel(
        condition = "input.model == \"vg\"",
        textInput("kappa", "Jump(\\kappa)",value="1")
      )
    )
  )
)

shinyApp(ui = ui, server = server)
