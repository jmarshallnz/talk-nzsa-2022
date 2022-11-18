library(tidyverse)
library(seqinr)

Dat <- read.csv(here::here("raw/SNP_data_Full_JM.csv"),header=TRUE)
SourceDat <- read.csv(here::here("raw/LabID_Source_Species_560.csv"))  #Animal isolate numbers

Source <- SourceDat$LabID
DatLong <- Dat[,c(3,4,5,7)]
DatLong$Isolate <- factor(DatLong$Isolate)
DatWide <- spread(DatLong[,-3], Gene, SNPs)
df <- as.data.frame(nchar(DatWide[1,-1]))
for(i in 2:1211){
  df1 <- as.data.frame(nchar(DatWide[i,-1]))
  df <- cbind(df, df1)
}
x <- which(is.na(as.data.frame(apply(df,1,function(x) {max(x)-min(x)}))))
length(which(is.na(as.data.frame(apply(df,1,function(x) {max(x)-min(x)})))))  # 41 genes have NAs
names(DatWide)[x+1]
#DatWideBU <- DatWide
for(i in x){
  a <- DatWide[,i+1]
  b <- max(nchar(a), na.rm=TRUE)
  d <- c2s(rep("-", times=b))
  e <- which(is.na(a))
  DatWide[e,i+1] <- d
}
HumanIsolatesWide <- DatWide[!(DatWide$Isolate %in% Source),]
AnimalIsolatesWide <- DatWide[DatWide$Isolate %in% Source,]

library(stringdist)
min_dist <- matrix(0, nrow(HumanIsolatesWide), 1343)
for (j in 1:1343) {
  # Find distance between each human and source
  dist <- stringdistmatrix(HumanIsolatesWide[,j+1], AnimalIsolatesWide[,j+1],
                           method="hamming", useBytes=TRUE)
  min_dist[,j] <- apply(dist, 1, min) # check row vs column
  cat("Done column", j, "\n")
}
write_csv(as.data.frame(min_dist), "data/distance.csv")

#plot_this <- list()
#for (j in 0:5) {

# We want to find the maximum number of genes we can use if we want to attribute at least K humans (K = 1..651)
# This turns out to be a famous problem. Essentially we want to be able to delete rows and/or columns while maximising
# the amount of data we have left. In general, this is an NP-hard problem known as the Maximal biclique in a bipartite graph.
# This optimises the 'area' (i.e. nrows*ncol left after deleting rows/columns).
# But, we want to be able to vary the number of humans from 1 to 651 and find the number of columns to use then.
# If instead of focusing on the rows/columns to keep, we focus on the rows/columns to remove, then the problem reduces to
# a vertex covering problem. If we're only interested in minimising the sum of (num removed columns + num removed rows) then
# that is polynomial-solvable. Turns out we can reframe it (like lots of problems) as a linear programming problem. Treat
# the humans and genes as nodes in a bipartite (two-coloured) graph, and link them with an edge whenever a human has a unique
# allele for a particular gene. We then want a vertex-covering: Pick the smallest number of vertices such that when we remove
# them, the graph is made completely disconnected (all edges are gone). Once we've done that, we have no edges, and thus
# no unique alleles left.
# This is equivalent to the linear program.
# Let v be a binary vector of length num_humans + num_genes (one per vertex). If v=0 we keep it. If v=1 we drop it.
# Then:
# minimise sum(v) for vertices V
# such that v_i + v_j >= 1 if (i,j) is an Edge (if a human is connected to a gene - has a unique allele)
# and such that v_i is binary.

library(tidyverse)
library(ROI)
library(ROI.plugin.glpk)
library(furrr)
library(glue)

min_dist <- read_csv(here::here("data/distance.csv"))

dist_to_nearest <- 0:5 # allow this much distance before we count as being 'unique'
humans_required <- seq(0, 651, by=21) # the number of humans that we require (we'll find the maximum number of genes we can use)
#seq(0, 651-11, by=21)+11 # the number of humans that we require (we'll find the maximum number of genes we can use)
humans_required <- expand_grid(s=seq(0,651-21,by=21), x=1:20) |>
  mutate(h = s + x) |> pull(h)

plan(multisession, workers = 4) # use up to 7 cores for the future_map_dfr() call. NOTE: The outside loop is sequential

compute_max_genes_by_distance <- function(dist_to_nearest, humans_required) {
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

  nm <- c(paste0("H",1:num_humans),paste0("G",1:num_genes))

  objective <- L_objective(L = c(rep(0, num_humans), rep(1, num_genes)), names=nm)

  constraint_base <- L_constraint(L = t(incidence), dir = rep('>=', num_edges), rhs = rep(1, num_edges), names=nm)

  # This is repeated for k = 1..num_humans
  find_max_genes <- function(k) {
    # we want to allow removal of up to k humans
    constraint_rows <- L_constraint(L = matrix(c(rep(1, num_humans), rep(0, num_genes)), nrow=1), dir = "<=", rhs=num_humans - k)

    lp <- OP(objective, c(constraint_base, constraint_rows), types = rep("B", num_humans+num_genes), maximum=FALSE)

    cat("Solving linear programming problem for k=", k, "...")

    foo = system.time({
      (lp_sol <- ROI_solve(lp, solver="glpk"))
    })
    cat("took", foo[3], "s\n")

    out <- solution(lp_sol) %>% enframe %>% filter(value == 1)

    humans <- out %>% filter(str_starts(name, "H")) %>%
      mutate(row = as.numeric(str_sub(name, 2))) %>%
      pull(row)

    genes <- out %>% filter(str_starts(name, "G")) %>%
      mutate(col = as.numeric(str_sub(name, 2))) %>%
      pull(col)

    data.frame(humans = k,
               genes = num_genes - length(genes))
  }

  distances <- future_map_dfr(humans_required, find_max_genes, .progress = TRUE) %>%
    mutate(dist_to_nearest = dist_to_nearest)
  
  distances
}

everything <- map_dfr(dist_to_nearest, compute_max_genes_by_distance, humans_required = humans_required)

ggplot(everything) + geom_line(aes(x=humans, y=genes, col=dist_to_nearest))

write_csv(everything, here::here("data/max_genes_by_distance.csv"))

read_csv(here::here('data/max_genes_by_distance.csv')) |>
  ggplot() +
  geom_line(aes(x=humans, y=genes, col=as_factor(dist_to_nearest)))
