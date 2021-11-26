library(tidyverse)


make_density_plot <- function(density_data, taxon) {
  density_data %>% 
    filter(tax_id == taxon) %>%
    mutate(y = y/max(y)) %>%
    ggplot(aes(x=x)) + 
    geom_ribbon(aes(ymin=y, ymax=1), fill="grey90") + 
    geom_ribbon(aes(ymin=0, ymax=y), fill="grey70", color="black")  +
    scale_x_continuous(labels = scales::percent_format(scale = 100))
}

make_histogram_plot <- function(histogram_data, taxon) {
  hist_data %>%
    filter(tax_id == taxon) %>%
    rename(right_break = breaks) %>%
    mutate(left_break =right_break - min(right_break),
           widths = right_break - left_break) %>%
    ggplot() + 
    geom_rect(aes(xmin=left_break, xmax=right_break, ymin=0, ymax=counts),
              fill="grey60", color="grey30")
}

histogram_data <- read_csv("histograms.csv")[,-1]
histogram_nonzero_data <- read_csv("histograms_nonzero.csv")[,-1]
density_data <- read_csv("densities.csv")[,-1]
density_nonzero_data <- read_csv("densities_nonzero.csv")[,-1]

make_density_plot(density_nonzero_data, taxon = "s__muciniphila1")
