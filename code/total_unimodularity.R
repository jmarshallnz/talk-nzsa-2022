# playing with totally unimodular matrices

h = 3
g = 6

library(igraph)
library(lintools)

incidence_matrix <- function(g) {
  A <- matrix(0, nrow=length(V(g)), ncol=length(E(g)))
  edge <- ends(g, E(g))
  A[cbind(edge[,1], 1:nrow(edge))] <- 1
  A[cbind(edge[,2], 1:nrow(edge))] <- 1
  A
}

graph <- sample_bipartite(h, g, p=0.2)
A <- graph |>
  incidence_matrix() |> t() |>
  rbind(c(rep(-1, h), rep(0, g)))
A |> is_totally_unimodular()

# sample a bunch of square submatrices until we hit det
subsample <- function(A, n=1000) {
  for (i in 1:n) {
    size <- sample(2:min(nrow(A), ncol(A)), 1)
    rows <- sample(nrow(A), size)
    cols <- sample(ncol(A), size)
    if (!det(A[rows, cols]) %in% c(-1, 0, 1)) {
      return(list(rows = rows, cols = cols))
    }
  }
  return(NULL)
}

a <- subsample(A); a

det(A[a$rows, a$cols])
A
# check for total unimod

# This matrix will do it: 3rd row is row of -1's
A = matrix(1, 3, 3);
Ab = A - diag(diag(A))
det(Ab)

