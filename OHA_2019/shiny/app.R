library(tidyverse)
library(patchwork)


attr_data <- read.csv("attribution_data.csv")
sts = attr_data %>% filter(Source != "Human" | Year >= 2008) %>%
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
  tags$head(
    tags$style(HTML("
      span.irs-min,span.irs-max,span.irs-single {
        visibility: hidden !important;
      }
    "))),
  fluidRow(plotOutput("sourceplot", height="160px")),
  fluidRow(column(3,sliderInput("s1", NULL, 0, 100, 25, step=1, ticks=FALSE), align='center'),
           column(3,sliderInput("s2", NULL, 0, 100, 25, step=1, ticks=FALSE), align='center'),
           column(3,sliderInput("s3", NULL, 0, 100, 25, step=1, ticks=FALSE), align='center'),
           column(3,sliderInput("s4", NULL, 0, 100, 25, step=1, ticks=FALSE), align='center')),
  fluidRow(plotOutput("humanplot", height="330px"))
)

server <- function(input, output, session) {
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
    bal <- c(input$s1, input$s2, input$s3, input$s4)
    balance <- tibble(Source = levels(top20_anim$Source), Balance = bal/sum(bal)) %>%
      mutate(Source = fct_inorder(Source))
    attrib_data <- top20_anim %>% left_join(balance, by="Source") %>% mutate(Total = Count * Balance)
    g1 = ggplot(balance, aes(x=Source)) + geom_col(aes(y=Balance, fill=Source)) +
      theme_bw(base_size=16) +
      guides(fill='none') +
      scale_y_continuous(name = "Overall attribution", expand=c(0,0), labels=scales::percent_format(accuracy = 1), limits=c(0,1)) +
      scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
      theme(axis.title.x = element_blank(),
            axis.ticks = element_blank())
    g2 = ggplot(attrib_data, aes(x=ST)) + geom_col(aes(y=Total, fill=Source)) +
      geom_col(data=top20_hum, aes(y=Human), fill=NA, col='black') +
      theme_bw(base_size=16) +
      guides(fill='none') +
      coord_cartesian(ylim = c(0,0.25)) +
      scale_y_continuous(name = "Attributed human cases", expand=c(0,0), labels=scales::percent_format(accuracy = 1)) +
      scale_fill_manual(name=NULL, values = c(Poultry="brown", Ruminants="steelblue2", Other="plum4", Water="green4")) +
      theme(axis.title.x = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = 90, vjust=0.5))
    g1 + g2 + plot_layout(nrow = 1, widths = c(1,3))
  })
}

shinyApp(ui, server)
