require(tidyverse)
require(igraph)
require(ggraph)
require(ade4) # for jaccard similarity index
require(lsa) #for cosine similarity index

# Import data from file ---------------------------------------------------

dat <- read_csv("Catalog_SiteA.csv",
                col_select = c(LEVEL_ID, CODE))

# Create un-weighted bipartite graph --------------------------------------

g_assemblages_bpg <-
  graph_from_data_frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                        directed = FALSE)

V(g_assemblages_bpg)$type <-
  bipartite_mapping(g_assemblages_bpg)$type

g_assemblages_bpg %>%
  ggraph(layout = "bipartite") +
  geom_edge_link(color = "gray", alpha = 0.25) +
  geom_node_point(aes(color = type), size = 2) +
  scale_color_manual(
    values = c("green", "blue"),
    name = "Node Type",
    breaks = c(FALSE, TRUE),
    labels = c("Provenience", "Artifact")
  )


# Create incidence matrix from bipartite graph ----------------------------

g_assemblages_inc <- as_incidence_matrix(g_assemblages_bpg)

# Project one-mode graphs igraph ------------------------------------------

g_assemblages_proj <-
  bipartite_projection(g_assemblages_bpg, multiplicity = TRUE)

g_assemblage_prov <- g_assemblages_proj$proj1

g_assemblage_artifact <- g_assemblages_proj$proj2

g_assemblage_prov %>%
  ggraph(layout = "fr") +
  geom_edge_link(color = "gray", aes(alpha = weight)) +
  geom_node_point(color = "green", size = 2) +
  ggtitle("Network of Proveniences")

g_assemblage_artifact %>%
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

prov_adj_ssoc <- overlap_coef(t(g_assemblages_inc))

artifact_adj_ssoc <- overlap_coef(g_assemblages_inc)


# Project one-mode graphs, Sorenson-Dice -----------------------------------

soren_dice_sim <- function(x, y = NULL) {
  if (is.null(y))
    y <- x
  t_mat <- (t(x) %*% y)
  x_count <- apply(x, 2, function(xx)
    sum(xx != 0))
  y_count <- apply(y, 2, function(yy)
    sum(yy != 0))
  res <- (2 * t_mat) / (x_count + y_count)
  if (is.null(y)) {
    diag(res) <- 1L
  }
  dimnames(res) <- list(colnames(x), colnames(y))
  return(res)
}

prov_adj_sdsim <- soren_dice_sim(t(g_assemblages_inc))

artifact_adj_sdsim <- soren_dice_sim(g_assemblages_inc)
