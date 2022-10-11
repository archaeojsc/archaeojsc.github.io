require(tidyverse)
require(igraph)
require(ggraph)

# Import data from file ---------------------------------------------------

dat <- read_csv("Catalog_SiteA.csv",
                col_select = c(LEVEL_ID, CODE))

# Create un-weighted bipartite graph --------------------------------------

g_assemblages_bpg <-
  graph_from_data_frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                        directed = FALSE)

V(g_assemblages_bpg)$type <-
  bipartite_mapping(g_assemblages_bpg)$type

g_assemblages_bpg_layout <-
  g_assemblages_bpg %>% layout_as_bipartite()

g_assemblages_bpg %>%
  ggraph(layout = g_assemblages_bpg_layout) +
  geom_edge_link(color = "gray", alpha = 0.25) +
  geom_node_point(aes(color = type), size = 2) +
  scale_color_manual(
    values = c("green", "blue"),
    name = "Node Type",
    breaks = c(FALSE, TRUE),
    labels = c("Provenience", "Artifact")
  )


# Create incidence matrix from bipartite graph ----------------------------

g_assemblages_bpg_inc <- as_incidence_matrix(g_assemblages_bpg)

# Project one-mode graphs igraph ------------------------------------------

g_assemblages_proj <-
  bipartite_projection(g_assemblages_bpg, multiplicity = TRUE)

g_assemblages_proj_prov <- g_assemblages_proj$proj1

g_assemblages_proj_artifact <- g_assemblages_proj$proj2

g_assemblages_proj_prov %>%
  ggraph(layout = "fr") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "green", size = 2) +
  ggtitle("Network of Proveniences")

g_assemblages_proj_artifact %>%
  ggraph(layout = "fr") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "blue", size = 2) +
  ggtitle("Network of Artifact Types")


# Project one-mode graphs, Szymkiewicz-Simpson -----------------------------

# adapted from `rbiosUtils` function `overlapCoefficient`

overlap_coef <- function(x, y = NULL) {
  if (is.null(y))
    y <- x
  
  t_mat <- t(x) %*% y
  
  x_count <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  y_count <- apply(y, 2, function(yy)
    sum(yy != 0))
  
  t_mat_pmin <- outer(x_count, y_count, FUN = pmin)
  
  res <- t_mat / t_mat_pmin
  
  if (is.null(y)) {
    diag(res) <- 1L
  }
  
  dimnames(res) <- list(colnames(x), colnames(y))
  
  return(res)
}


## Project provenience -----------------------------------------------------

prov_adj_ssoc <- overlap_coef(t(g_assemblages_bpg_inc))
diag(prov_adj_ssoc) <- 0

prov_ssoc_vals <- prov_adj_ssoc[lower.tri(prov_adj_ssoc)]

ggplot(data = data.frame(x = c(prov_ssoc_vals)), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Szymkiewicz-Simpson Similarity for Proveniences")

g_assemblages_proj_prov_oc <-
  graph_from_adjacency_matrix(prov_adj_ssoc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_proj_prov_oc %>% 
#   ggraph(layout = "kk") + 
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "green", size = 2) +
#   ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(
  g_assemblages_proj_prov_oc)), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Szymkiewicz-Simpson Degree for Proveniences")

ggplot(
  data = data.frame(x = E(g_assemblages_proj_prov_oc)$weight), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Szymkiewicz-Simpson Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

artifact_adj_ssoc <- overlap_coef(g_assemblages_bpg_inc)
diag(artifact_adj_ssoc) <- 0

artifact_ssoc_vals <-
  artifact_adj_ssoc[lower.tri(artifact_adj_ssoc)]

ggplot(data = data.frame(x = c(artifact_ssoc_vals)), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Szymkiewicz-Simpson Similarity for Artifacts")

g_assemblages_proj_artifact_oc <-
  graph_from_adjacency_matrix(
    artifact_adj_ssoc,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

# g_assemblages_proj_artifact_oc %>% 
#   ggraph(layout = "kk") +
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "blue", size = 2) +
#   ggtitle("Network of Artifacts")

ggplot(data = data.frame(x = degree(
  g_assemblages_proj_artifact_oc)), aes(x = x)) + 
  geom_density(color = "blue") + 
  ggtitle("Distribution of Szymkiewicz-Simpson Degree for Artifacts")

ggplot(
  data = data.frame(x = E(g_assemblages_proj_artifact_oc)$weight), aes(x = x)) + 
  geom_density(color = "blue") + 
  ggtitle("Distribution of Szymkiewicz-Simpson Edge Weight for Artifacts")


# Project one-mode graphs, Sorenson-Dice -----------------------------------

soren_dice_sim <- function(x, y = NULL) {
  if (is.null(y))
    y <- x
  
  t_mat <- t(x) %*% y
  
  x_count <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  y_count <- apply(y, 2, function(yy)
    sum(yy != 0))
  
  t_mat_sum <- outer(x_count, y_count, FUN = "+")
  
  res <- (2 * t_mat) / t_mat_sum
  
  if (is.null(y)) {
    diag(res) <- 1L
  }
  
  dimnames(res) <- list(colnames(x), colnames(y))
  return(res)
}

## Project provenience -----------------------------------------------------

prov_adj_sdsim <- soren_dice_sim(t(g_assemblages_bpg_inc))
diag(prov_adj_sdsim) <- 0

prov_sdsim_val <-
  prov_adj_sdsim[lower.tri(prov_adj_sdsim)]

ggplot(data = data.frame(x = prov_sdsim_val), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Sorenson-Dice Similarity for Provenience")

g_assemblages_proj_prov_sdsim <-
  graph_from_adjacency_matrix(prov_adj_sdsim,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_proj_prov_sdsim %>% 
#   ggraph(layout = "kk") + 
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "green", size = 2) +
#   ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(
  g_assemblages_proj_prov_sdsim)), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Sorenson-Dice Degree for Proveniences")

ggplot(
  data = data.frame(x = E(g_assemblages_proj_prov_sdsim)$weight), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Sorenson-Dice Edge Weight for Proveniences")


## Project artifact types --------------------------------------------------

artifact_adj_sdsim <- soren_dice_sim(g_assemblages_bpg_inc)
diag(artifact_adj_sdsim) <- 0

artifact_sdsim_val <-
  artifact_adj_sdsim[lower.tri(artifact_adj_sdsim)]

ggplot(data = data.frame(x = artifact_sdsim_val), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Sorenson-Dice Similarity for Artifacts")


# Project one-mode graphs, Jaccard ----------------------------------------

jaccard_sim <- function(x, y = NULL) {
  if (is.null(y))
    y <- x
  
  t_mat <- t(x) %*% y
  
  x_count <- apply(x, 2, function(xx)
    sum(xx != 0))
  
  y_count <- apply(y, 2, function(yy)
    sum(yy != 0))
  
  t_mat_sum <- outer(x_count, y_count, FUN = "+")
  
  res <- t_mat / (t_mat_sum - t_mat)
  
  if (is.null(y)) {
    diag(res) <- 1L
  }
  
  dimnames(res) <- list(colnames(x), colnames(y))
  
  return(res)
}

## Project provenience -----------------------------------------------------

prov_adj_jacc <- jaccard_sim(t(g_assemblages_bpg_inc))
diag(prov_adj_jacc) <- 0

prov_jacc_val <-
  prov_adj_jacc[lower.tri(prov_adj_jacc)]

ggplot(data = data.frame(x = prov_jacc_val), aes(x = x)) +
  geom_density(color = "green") +
  ggtitle("Distribution of Jaccard Similarity for Provenience")

g_assemblages_proj_prov_jacc <-
  graph_from_adjacency_matrix(prov_adj_jacc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_proj_prov_jacc %>% 
#   ggraph(layout = "kk") + 
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "green", size = 2) +
#   ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(
  g_assemblages_proj_prov_jacc)), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Jaccard Degree for Proveniences")

ggplot(
  data = data.frame(x = E(g_assemblages_proj_prov_jacc)$weight), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Jaccard Edge Weight for Proveniences")



## Project artifact types --------------------------------------------------

artifact_adj_jacc <- jaccard_sim(g_assemblages_bpg_inc)
diag(artifact_adj_jacc) <- 0

artifact_jacc_val <-
  artifact_adj_jacc[lower.tri(artifact_adj_jacc)]

ggplot(data = data.frame(x = artifact_jacc_val), aes(x = x)) +
  geom_density(color = "blue") +
  ggtitle("Distribution of Jaccard Similarity for Artifacts")

g_assemblages_proj_artifact_jacc <-
  graph_from_adjacency_matrix(artifact_adj_jacc,
                              mode = "undirected",
                              weighted = TRUE,
                              diag = FALSE)

# g_assemblages_proj_prov_jacc %>% 
#   ggraph(layout = "kk") + 
#   geom_edge_link(color = "gray", aes(alpha = weight)) +
#   geom_node_point(color = "green", size = 2) +
#   ggtitle("Network of Proveniences")

ggplot(data = data.frame(x = degree(
  g_assemblages_proj_artifact_jacc)), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Jaccard Degree for Proveniences")

ggplot(
  data = data.frame(x = E(g_assemblages_proj_artifact_jacc)$weight), aes(x = x)) + 
  geom_density(color = "green") + 
  ggtitle("Distribution of Jaccard Edge Weight for Proveniences")

# TOM Adjacency -----------------------------------------------------------

tom_adjacency_matrix <- function(adj_mat) {
  l_mat <- adj_mat %*% t(adj_mat)
  
  k_row <- rowSums(adj_mat)
  k_col <- colSums(adj_mat)
  k_min <- outer(k_row, k_col, FUN = pmin)
  
  tom_adj <- (l_mat + adj_mat) / (k_min + 1 - adj_mat)
  
  diag(tom_adj) <- 1
  
  dimnames(tom_adj) <- list(colnames(adj_mat), colnames(adj_mat))
  
  return(tom_adj)
}

## Provenience Overlap Matrix -----------------------------------------------------



## Project artifact types --------------------------------------------------

