library(gtools)
library(gganimate)
library(transformr)
library(gifski)

set.seed(3)
alpha <- top50$Count
random <- as.data.frame(t(rdirichlet(20, alpha)))
#random <- bind_cols(random, V21=random[,1])

uncertainty <- bind_cols(top50, random) %>%
  select(-Count, -Proportion) %>% gather(Iteration, Proportion, V1:V20) %>%
  extract(Iteration, into="Iteration", regex="([0-9]+)", convert=TRUE)

# plot
anim <- ggplot(uncertainty, aes(x=ST, y=Proportion, fill=Source)) + geom_col() +
  scale_y_continuous("Percent", labels = scales::percent, expand=c(0,0)) +
  scale_fill_manual(values=c(Poultry = "steelblue", Prior = "grey50")) +
  guides(fill='none') +
  coord_cartesian(ylim = c(0,0.15)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=0.5)) +
  transition_states(Iteration, transition_length = 1, state_length = 0.1) +
  ease_aes('cubic-in-out')

animate(anim, nframes=220, fps=20, width=800, height=480, renderer=gifski_renderer(file="anims/dirichlet.gif"))
