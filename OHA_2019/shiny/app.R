library(tidyverse)

attr_data <- read.csv("data/attribution_data.csv")
sts = attr_data %>%
  group_by(ST) %>% count(Source) %>% spread(Source, n, fill=0) %>%
  ungroup()

top20 <- sts %>% mutate(ST = fct_lump(factor(ST), n=20, w=Human)) %>% gather(Source, Count, -ST) %>%
  group_by(ST, Source) %>% summarise(Count = sum(Count)) %>% group_by(Source) %>% mutate(Count = Count/sum(Count)) %>% ungroup() %>%
  spread(Source, Count) %>% mutate(ST = fct_reorder(ST, Human, .fun = identity, .desc=TRUE),
                                   ST = fct_relevel(ST, "Other", after = 23)) 

top20_hum <- top20 %>% select(ST, Human)
top20_anim <- top20 %>% select(-Human) %>% gather(Source, Count, -ST) %>%
  mutate(Source = fct_relevel(Source, "Other", after=2))

library(shiny)

ui <- fluidPage(
  fluidRow(plotOutput("sourceplot", height="160px")),
  fluidRow(column(3,sliderInput("s1", NULL, 0, 100, 25, step=1, ticks=FALSE, post='%'), align='center'),
           column(3,sliderInput("s2", NULL, 0, 100, 25, step=1, ticks=FALSE, post='%'), align='center'),
           column(3,sliderInput("s3", NULL, 0, 100, 25, step=1, ticks=FALSE, post='%'), align='center'),
           column(3,sliderInput("s4", NULL, 0, 100, 25, step=1, ticks=FALSE, post='%'), align='center')),
  fluidRow(column(8, offset=2, plotOutput("humanplot", height="330px")))
)

server <- function(input, output, session) {
  state <- reactiveValues(balance = rep(0.25, 4))
  observeEvent(input$s1, {
    # normalise the others
    state$balance[1] <- input$s1 / 100
    state$balance[-1] <- state$balance[-1] / sum(state$balance[-1]) * (1-state$balance[1])
  })
  observeEvent(input$s2, {
    # normalise the others
    state$balance[2] <- input$s2 / 100
    state$balance[-2] <- state$balance[-2] / sum(state$balance[-2]) * (1-state$balance[2])
  })
  observeEvent(input$s3, {
    # normalise the others
    state$balance[3] <- input$s3 / 100
    state$balance[-3] <- state$balance[-3] / sum(state$balance[-3]) * (1-state$balance[3])
  })
  observeEvent(input$s4, {
    # normalise the others
    state$balance[4] <- input$s4 / 100
    state$balance[-4] <- state$balance[-4] / sum(state$balance[-4]) * (1-state$balance[4])
  })
  
  observeEvent(state$balance, {
    updateSliderInput(session, "s1", value = state$balance[1] * 100)
    updateSliderInput(session, "s2", value = state$balance[2] * 100)
    updateSliderInput(session, "s3", value = state$balance[3] * 100)
    updateSliderInput(session, "s4", value = state$balance[4] * 100)
  })

  output$sourceplot = renderPlot({
    ggplot(top20_anim, aes(x=ST, y=Count, fill=Source)) +
      geom_col() +
      facet_wrap(~Source, nrow=1) +
      guides(fill='none') +
      theme_bw(base_size=14) +
      coord_cartesian(ylim = c(0,0.20)) +
      scale_y_continuous(expand=c(0,0)) +
      scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank())
  })
  output$humanplot = renderPlot({
    balance <- tibble(Source = levels(top20_anim$Source), Balance = state$balance)
    attrib_data <- top20_anim %>% left_join(balance, by="Source") %>% mutate(Total = Count * Balance)
    ggplot(attrib_data, aes(x=ST)) + geom_col(aes(y=Total, fill=Source)) +
      geom_col(data=top20_hum, aes(y=Human), fill=NA, col='black') +
      theme_bw(base_size=16) +
      guides(fill='none') +
      coord_cartesian(ylim = c(0,0.25)) +
      scale_y_continuous(name = "Attributed human cases", expand=c(0,0), labels=scales::percent_format(accuracy = 1)) +
      scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
      theme(axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 90, vjust=0.5))
  })
}

shinyApp(ui, server)
