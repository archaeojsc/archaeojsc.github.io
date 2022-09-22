require(tidyverse)
require(tidygraph)
require(igraph)
require(ggplot2)
require(ggraph)
# require(visNetwork)
# require(lsa)
# require(WGCNA)
# library(readr)
# require(bipartite)
require(ade4)
# require(tnet)



# Import data from file ---------------------------------------------------

dat <- read_csv("Catalog_SiteA.csv")

# Create bipartite graph from unique provenience & artifact code pairs ----

G_assemblages_bpg <-
  graph_from_data_frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                        directed = FALSE)

V(G_assemblages_bpg)$type <- bipartite_mapping(G_assemblages_bpg)$type

plot(
  G_assemblages_bpg,
  layout = layout.bipartite,
  vertex.size = 5,
  vertex.label = NA,
  vertext.color = V(G_assemblages_bpg)$type
)

# Create incidence matrix from bipartite graph ----------------------------

assemblages_inc <- as_incidence_matrix(G_assemblages_bpg)

# Project graph of proveniences -------------------------------------------

assemblages_projection <-
  bipartite_projection(G_assemblages_bpg, multiplicity = TRUE)

prov_g <- assemblages_projection$proj1

artifact_g <- assemblages_projection$proj2

## Project Overlap by counts ----------------------------------------------
prov_adj_overlap <-
  assemblages_inc %*% t(assemblages_inc)

## Project by Jaccard index similarity ------------------------------------
prov_adj_jacc <-
  as.matrix(1 - dist.binary(
    assemblages_inc,
    method = 1,
    upper = TRUE,
    diag = FALSE
  ) ^ 2)

## Find density of non-zero similarity scores -----------------------------
plot(density(c(prov_adj_jacc[prov_adj_jacc > 0])))

## Find threshold for strongly similar nodes ------------------------------
prov_sim_thresh <- quantile(c(prov_adj_jacc[prov_adj_jacc > 0]), 0.95)

## Binarize by threshold value --------------------------------------------
prov_adj_jacc_thresh <-
  ifelse(prov_adj_jacc < prov_sim_thresh, 0, prov_adj_jacc)

# Create provenience graph ------------------------------------------------
G_prov <- graph_from_adjacency_matrix(
  prov_adj_jacc,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

plot(
  G_prov,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold
)

## Find and remove isolated nodes -----------------------------------------
G_prov_isolates <- which(degree(G_prov) == 0)

G_prov_conn <- delete.vertices(G_prov, G_prov_isolates)

plot(
  G_prov_conn,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout_nicely
)


# Project graph of artifacts ----------------------------------------------


## Project by overlap -----------------------------------------------------

artifact_adj_overlap <- t(assemblages_inc) %*% assemblages_inc

## Project by Jaccard index similarity ------------------------------------

artifact_adj_jacc <- # by Jaccard index
  as.matrix(1 - dist.binary(
    t(assemblages_inc),
    method = 1,
    upper = TRUE,
    diag = FALSE
  ) ^ 2)

## Find density of non-zero similarity scores -----------------------------

plot(density(c(artifact_adj_jacc[artifact_adj_jacc > 0])))

## Find threshold for strongly similar nodes ------------------------------
art_sim_thresh <- quantile(c(artifact_adj_jacc[artifact_adj_jacc > 0]), 0.95)

## Binarize by threshold value --------------------------------------------
artifact_adj_jacc_thresh <-
  ifelse(artifact_adj_jacc < art_sim_thresh, 0, artifact_adj_jacc)

# Create artifact graph ---------------------------------------------------

G_artifact <-
  graph_from_adjacency_matrix(
    artifact_adj_jacc_thresh,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

plot(
  G_artifact,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold
)

## Find and remove isolated nodes -----------------------------------------
G_artifact_isolates <- which(degree(G_artifact) == 0)

G_artifact_conn <- delete.vertices(G_artifact, G_artifact_isolates)

plot(
  G_artifact_conn,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold
)
