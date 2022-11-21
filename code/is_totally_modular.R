library(tidyverse)
library(ROI)
library(ROI.plugin.glpk)
library(furrr)
library(glue)

min_dist <- read_csv(here::here("data/distance.csv"))

dist_to_nearest <- 5 # allow this much distance before we count as being 'unique'
humans_required <- seq(0, 651, by=31) # the number of humans that we require (we'll find the maximum number of genes we can use)
#seq(0, 651-11, by=21)+11 # the number of humans that we require (we'll find the maximum number of genes we can use)
#humans_required <- expand_grid(s=seq(0,651-21,by=21), x=1:20) |>
#  mutate(h = s + x) |> pull(h)

humans_required <- seq(0, 50)

plan(multisession, workers = 4) # use up to 7 cores for the future_map_dfr() call. NOTE: The outside loop is sequential

is_tm <- function(dist_to_hearest, humans_required) {
  # construct the graph incidence matrix we'll need. We do this by first pivoting long
  # then filtering out the edges we want
  edges <- min_dist %>% tibble::rowid_to_column("human") %>%
    pivot_longer(-human, names_to="gene", values_to="edge", names_prefix="V") %>%
    mutate(gene = as.numeric(gene)) %>%
    filter(edge > dist_to_nearest)

  num_humans <- nrow(min_dist)
  num_genes <- ncol(min_dist)

  # identify all the edges and replace them with uniques?
  num_edges <- nrow(edges)

  incidence <- matrix(0, nrow=num_humans+num_genes, ncol=num_edges)
  incidence[cbind(edges$human, 1:nrow(edges))] <- 1
  incidence[cbind(edges$gene + num_humans, 1:nrow(edges))] <- 1

  A <- t(incidence)

  is_tm_k <- function(k) {
    constr <- matrix(c(rep(1, num_humans), rep(0, num_genes)), nrow=1)
    A_comb <- rbind(A, constr)
    return(data.frame(k=k, tm = is_totally_unimodular(A_comb)))
  }

  tm <- future_map_dfr(humans_required, is_tm_k, .progress = TRUE)
  
  tm
}

everything <- map_dfr(dist_to_nearest, is_tm, humans_required = humans_required)

ggplot(everything) + geom_line(aes(x=humans, y=genes, col=dist_to_nearest))
