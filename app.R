library(shiny)
library(shinythemes)
library(shinycustomloader)
library(ggplot2)
library(ggpubr)
library(clam)

calSampleApprox <- function(x,t1,t2,r){
    n <- length(x)
    funs <- lapply(x,approxfun)
    y_list <- lapply(1:n,function(j)funs[[j]](seq(t1,t2,r)))
    y_mat <- do.call(cbind,y_list)
    y_mat[which(is.na(y_mat))] <- 0
    return(y_mat)
}

makeSpdf <- function(input){
    enddate <- input$startdate - input$span + 1
    if(input$dgp == 1){
        dgpfun <- exp((1:input$span) * input$exprate)
    }else if(input$dgp == 2){
        dgpfun <- (input$A *
                    sin(2 *
                        pi *
                        input$f *
                        (1:input$span) +
                        input$phase)) + (input$A + 1e-5)
    }else if(input$dgp == 3){
        lmid <- input$logmid * input$span
        dgpfun <- input$logmax / ( 1 + exp(-input$lograte * ((1:input$span) - lmid)) )
    }else{
        dgpfun <- rep(1,input$span)
    }
    dgp_df <- data.frame(
                        Process = dgpfun,
                        YBP = input$startdate:enddate
                    )
    simdates <- sample(
                    x=input$startdate:enddate,
                    size=input$ndates,
                    replace=T,
                    prob=dgpfun)
    event_counts <- as.data.frame(table(simdates))
    names(event_counts) <- c("YBP","Count")
    event_counts$YBP <- as.numeric(levels(event_counts$YBP))
    event_counts$Count <- as.numeric(event_counts$Count)
    event_counts <- merge(data.frame(YBP=input$startdate:enddate),
                        event_counts,
                        by=1,
                        all=T
                        )
    event_counts[which(is.na(event_counts[,2])),2] <- 0
    event_counts$Time <- rev(1:input$span)
    simc14 <- t(sapply(simdates,calBP.14C))
    c14post <- lapply(1:input$ndates,
                        function(x){
                            capture.output(d <- calibrate(simc14[x,1],
                                                        simc14[x,2],
                                                        graph=F)$calib)
                            return(d)
                        })
    sample_date_range <- range(
                            unlist(
                                lapply(
                                    c14post,
                                    function(x)range(x[,1]))
                                )
                            )
    c14_matrix <- calSampleApprox(
                                c14post,
                                sample_date_range[1],
                                sample_date_range[2],
                                r=1
                                )
    spdf <- data.frame(
                    YBP = seq(
                            sample_date_range[1],
                            sample_date_range[2]
                            ),
                    SumDensity = rowSums(c14_matrix),
                    Time = rev(1:(diff(sample_date_range) + 1))
                    )
    return(
        list(
            dgp_df=dgp_df,
            event_counts=event_counts,
            spdf=spdf,
            sample_date_range=sample_date_range
            )
    )
}


source("./content.R")

# Define UI ----
ui <- fluidPage(
    fluidRow(
        style="background-color: #303030",
        pagetitle,
        shorthead,
        column(12,
            fluidRow(
                style="background-color: #FFFFFF",
                plotspacetitle,
                withLoader(
                    plotOutput(
                        outputId = "spdfPlot"
                        ),
                    type="image",
                    loader="spdf_loader.gif"
                ),
                hr()
                ),
            fluidRow(
                style="background-color: #AEABAD",
                column(11,
                    paramstitle
                ),
                column(1,
                    actionButton("draw","Draw")
                )
                ),
            fluidRow(
                style="background-color: #AEABAD",
                column(3,
                    numericInput(
                        "startdate",
                        h3("Start Date (cal. BP)"),
                        value = floor(runif(
                                    n = 1,
                                    min = 2000,
                                    max = 48000
                                    ))
                        ),
                    numericInput(
                        "span",
                        h3("Span"),
                        value = 1000
                        ),
                    numericInput(
                        "ndates",
                        h3("N. Dates"),
                        value=500
                        )
                    ),
                column(3,
                    selectInput(
                        "dgp",
                        h3("True Process"),
                        choices = list(
                                        "Exponential" = 1,
                                        "Sinusoidal" = 2,
                                        "Logistic" = 3),
                        selected = 1
                        ),
                    conditionalPanel(
                        style = "border: 1px solid #d9d9d9; border-radius: 5px; background-color: #EEE9EC",
                        condition = "input.dgp == 1",
                        numericInput(
                            "exprate",
                            h4("Rate"),
                            value = 0.004
                            )
                        ),
                    conditionalPanel(
                        style = "border: 1px solid #d9d9d9; border-radius: 5px; background-color: #EEE9EC",
                        condition = "input.dgp == 2",
                        numericInput(
                            "A",
                            h4("Amplitude"),
                            value = 1
                            ),
                        numericInput(
                            "f",
                            h4("Common Frequency"),
                            value = 0.01
                            ),
                        numericInput(
                            "phase",
                            h4("Phase"),
                            value = 0
                            )
                        ),
                    conditionalPanel(
                        style = "border: 1px solid #d9d9d9; border-radius: 5px; background-color: #EEE9EC",
                        condition = "input.dgp == 3",
                        numericInput(
                            "lograte",
                            h4("Rate"),
                            value = 0.01
                            ),
                        numericInput(
                            "logmax",
                            h4("Max."),
                            value = 1
                            ),
                        sliderInput(
                            "logmid",
                            h4("Midpoint"),
                            min=0,
                            max=1,
                            value = 0.5
                            )
                        )
                    ),
                column(6,
                    style = "background-color: #636162",
                    explainertext1,
                    explainertext2,
                    contact
                    )
                )
            )
        ),
    theme = shinytheme("cosmo")
)

# Define server logic ----
server <- function(input, output) {

    userinputs <- reactive({
        validate(
            need(
                input$startdate <= 48000 & (input$startdate - input$span) >= 500,
                "Error: Select a start date and span that produce an interval between 48,000 cal. BP and 500 cal. BP")
        )
    })

    spdfdata <- eventReactive(input$draw,{
        makeSpdf(input)
    })

    output$spdfPlot <- renderPlot({
        if(input$draw == 0) {
            instruct_img <- readPNG("./www/instructions.png")
            plot(as.raster(instruct_img))
        }else{
            userinputs()
            spdfdat <- spdfdata()
            sample_date_range <- spdfdat$sample_date_range
            ##plots
            p1 <- ggplot(data=spdfdat$dgp_df
                            ) +
                        geom_path(
                                mapping=aes(y=Process,x=YBP),
                                alpha = 0.8
                                ) +
                        xlim(c(sample_date_range[2],sample_date_range[1])
                            ) +
                        labs(
                            y="Process",
                            x="Year BP"
                            ) +
                        theme_minimal() +
                        theme(
                            plot.margin=unit(c(0,1,0,1),"cm"),
                            text = element_text(size=16))
            p2 <- ggplot(data=spdfdat$event_counts) +
                        geom_col(
                                mapping=aes(y=Count,x=YBP),
                                position="identity",
                                alpha=0.8,
                                colour=NA
                                ) +
                        xlim(c(sample_date_range[2],sample_date_range[1])
                            ) +
                        labs(
                            y="Count",
                            x="Year BP"
                            ) +
                        theme_minimal() +
                        theme(
                            plot.margin=unit(c(0,1,0,1),"cm"),
                            text = element_text(size=16),
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            p3 <- ggplot(data=spdfdat$spdf) +
                        geom_area(
                                mapping=aes(y=SumDensity,x=YBP),
                                alpha=0.7,
                                fill="steelblue",
                                colour="black"
                                ) +
                        labs(
                            y="Summed Density",
                            x="Year BP"
                            ) +
                        xlim(c(sample_date_range[2],sample_date_range[1])
                            ) +
                        theme_minimal() +
                        theme(
                            plot.margin=unit(c(0,1,0,1),"cm"),
                            text = element_text(size=16),
                            axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.ticks.x=element_blank())
            fig <- ggarrange(
                            p3,p2,p1,
                            ncol=1,
                            nrow=3,
                            align="v"
                            )
            fig
        }
    })

}

# Run the app ----
shinyApp(ui = ui, server = server)
