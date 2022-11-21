library(tidyverse)

# read in our data
min_dist <- read_csv(here::here("data/distance.csv"))

md <- as.matrix(min_dist)

md[md < 5] <- 0
md[md > 0] <- 1

# ok, now sample rows randomly and then pick the best genes
h <- 100
 
g <- function(h, prob_h) {
  sum(colSums(md[sample(nrow(md), size=h, prob=prob_h),,drop=FALSE]) == 0)
}

g_max <- function(h, prob_h, samples=1000) {
  cat("doing h=", h, "\n")
  replicate(n=samples, g(h, prob_h)) |>
    max()
}

h <- seq(1, nrow(md), 21)
g_h <- map(h, g_max, prob_h = ncol(md) - rowSums(md),
           samples=10000) # favour big h
plot(h, g_h, type='l')
