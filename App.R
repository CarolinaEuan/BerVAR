library(shiny)
library(here)
library(shinythemes)
here()
source("SourceBerVAR.R")

ui = fluidPage(theme = shinytheme("cerulean"),
               column(3,img(src="1.png",width=300,height=100)),
               column(3,offset=6,a(href="https://es.kaust.edu.sa/Pages/Home.aspx",
                                   img(src="envstat.png",width=100,height=100))),
               withMathJax(),
               fluidRow(
                 column(10,offset=.1,h1("Ber-VAR model")),
                 column(6,offset = 0.1,p('This app shows the effect of the 
                                          parameters for the latent processes \\((Y_{t},Z_{t})\\) 
                                          on the covariance structure of the observed process \\(X_{t}\\).'),
               a(href="https://es.kaust.edu.sa/Pages/CarolinaEuan.aspx",h5("PhD, Carolina Euan.")))
               ),
               hr(),
               sidebarPanel(    selectInput(inputId="LP",label="Select a process",
                                c( 'Y' = "y", 'Z'  = "z")),
                                helpText("Parameters:"), 
                                conditionalPanel(condition = "input.LP=='y' ", 
                                            sliderInput("PY1", withMathJax('Marginal probability $$(p_1^Y,p_2^Y,
                                                                           p_3^Y)$$'),
                                            min = 0, max = 1, value = .5,step = .05
                                ),
                                sliderInput("PY2", " ",
                                            min = 0, max = 1, value = .5,step = .05
                                ),
                                sliderInput("PY3", " ",
                                            min = 0, max = 1, value = .5,step = .05
                                ),
                                sliderInput("RY12", withMathJax('Correlation $$(\\rho_{12}^Y,\\rho_{13}^Y,
                                                                           \\rho_{23}^Y)$$'),
                                            min = -1, max = 1, value = 0,step = .1
                                ),
                                sliderInput("RY13", " ",
                                            min = -1, max = 1, value = 0,step = .1
                                ),
                                sliderInput("RY23", " ",
                                            min = -1, max = 1, value = 0,step = .1
                                )),
                                conditionalPanel(condition = "input.LP=='z' ", 
                                sliderInput("PZ1", withMathJax('Marginal probability $$(p_1^Z,p_2^Z,
                                                                           p_3^Z)$$'),
                                            min = 0, max = 1, value = .5,step = .05
                                                 ),
                                sliderInput("PZ2", " ",
                                            min = 0, max = 1, value = .5,step = .05
                                                 ),
                                sliderInput("PZ3", " ",
                                            min = 0, max = 1, value = .5,step = .05
                                                 ),
                                sliderInput("RZ12", withMathJax('Correlation $$(\\rho_{12}^Z,\\rho_{13}^Z,
                                                                           \\rho_{23}^Z)$$'),
                                            min = -1, max = 1, value = 0,step = .1
                                                 ),
                                sliderInput("RZ13", " ",
                                            min = -1, max = 1, value = 0,step = .1
                                                 ),
                                sliderInput("RZ23", " ",
                                            min = -1, max = 1, value = 0,step = .1
                                                 )),
                                actionButton("Run", "Run Setting"),
                                h4(htmlOutput("AR2Par"))),
               #output
               mainPanel(
                 tabsetPanel(
                  tabPanel("Correlation Plot",
                           plotOutput(outputId = "CorPlot",width = '100%',height = "600px")         
                           ),
                  tabPanel("Covariance Plot",
                           plotOutput(outputId = "CovPlot",width = '100%',height = "600px")         
                  )
                 )
               )
               
)

server<-function(input, output, session){
  pY<-reactive(c(input$PY1,input$PY2,input$PY3))
  pZ<-reactive(c(input$PZ1,input$PZ2,input$PZ3))
  Omega<-reactive(matrix(as.numeric(unlist(as.binary((1:2^3)-1,n=3,littleEndian=TRUE))),ncol=3,byrow = TRUE))
  observe({
    Bound1 <- Check.pho(pY(),Omega())
    updateSliderInput(session, "RY12", value = sum(Bound1[1,])/2,
                      min = max(-1,ceiling(Bound1[1,1]*100)/100), max = min(1,floor(Bound1[1,2]*100)/100), step = .1)
    updateSliderInput(session, "RY13", value = sum(Bound1[2,])/2,
                      min = max(-1,ceiling(Bound1[2,1]*100)/100), max = min(1,floor(Bound1[2,2]*100)/100), step = .1)
    updateSliderInput(session, "RY23", value = sum(Bound1[3,])/2,
                      min = max(-1,ceiling(Bound1[3,1]*100)/100), max = min(1,floor(Bound1[3,2]*100)/100), step = .1)
  })
  observe({
    Bound2 <- Check.pho(pZ(),Omega())
    updateSliderInput(session, "RZ12", value = sum(Bound2[1,])/2,
                      min = max(-1,ceiling(Bound2[1,1]*100)/100), max = min(1,floor(Bound2[1,2]*100)/100), step = .1)
    updateSliderInput(session, "RZ13", value = sum(Bound2[2,])/2,
                      min = max(-1,ceiling(Bound2[2,1]*100)/100), max = min(1,floor(Bound2[2,2]*100)/100), step = .1)
    updateSliderInput(session, "RZ23", value = sum(Bound2[3,])/2,
                      min = max(-1,ceiling(Bound2[3,1]*100)/100), max = min(1,floor(Bound2[3,2]*100)/100), step = .1)
  })
  observeEvent(input$Run,
               {
  Cov<-CovMat(pY(),pZ(),c(input$RY12,input$RY13,input$RY23),c(input$RZ12,input$RZ13,input$RZ23))
  output$CovPlot<-renderPlot({
    indexij<-cbind(rep(1:3,3),rep(1:3,each=3))
    par(mfrow=c(3,3))
    for(i in 1:9)plot(Cov[,i],type="b",pch=20,ylim=c(min(Cov),max(Cov)),xlab="Lag",
                      ylab = "Cov",main=mymain(indexij[i,1],indexij[i,2],FALSE),cex.main=2.5)
  })
  output$CorPlot<-renderPlot({
    indexij<-cbind(rep(1:3,3),rep(1:3,each=3))
    ind.equal<-which(indexij[,1]==indexij[,2])
    gamma0<-Cov[1,ind.equal]
    par(mfrow=c(3,3))
    for(i in 1:9)plot(Cov[,i]/sqrt(prod(gamma0[indexij[i,]])),type="b",pch=20,ylim=c(-1,1),xlab="Lag",
                      ylab="Corr",main=mymain(indexij[i,1],indexij[i,2]),cex.main=2.5)
  })
  })
}

shinyApp(ui=ui,server=server)


