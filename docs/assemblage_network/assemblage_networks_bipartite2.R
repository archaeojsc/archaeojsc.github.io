require(tidyverse)
require(igraph)
require(ggraph)

# Import data from file ---------------------------------------------------

dat <- read_csv("Catalog_SiteA.csv",
                col_select = c(LEVEL_ID, CODE))

# Create bipartite graph from unique provenience & artifact pairs ---------

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

# Project one-mode graphs -------------------------------------------------

assemblage_projections <-
  bipartite_projection(g_assemblages_bpg, multiplicity = TRUE)

g_assemblage_prov <- assemblage_projections$proj1

g_assemblage_artifact <- assemblage_projections$proj2

plot(
  g_assemblage_prov,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Provenience Network"
)

plot(
  g_assemblage_artifact,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Artifact Network"
)


# Analyzing one-mode graphs -----------------------------------------------

plot(density(degree(g_assemblage_prov)),
     main = "Distribution of Provenience Degree")

plot(density(E(g_assemblage_prov)$weight, bw = 1),
     main = "Distribution of Provenience Edge Weights")

plot(density(degree(g_assemblage_artifact)),
     main = "Distribution of Artifact Degree")

plot(density(E(g_assemblage_artifact)$weight, bw = 1),
     main = "Distribution of Artifact Edge Weights")


# Thinning out weak edges -------------------------------------------------

prov_edgeweight_thresh <- quantile(E(g_assemblage_prov)$weight, 0.95)

artifact_edgeweight_thresh <- quantile(E(g_assemblage_artifact)$weight, 0.95)

g_assemblage_prov_strong <-
  delete.edges(g_assemblage_prov, which(E(g_assemblage_prov)$weight < prov_edgeweight_thresh))

g_assemblage_prov_strong <-
  delete.vertices(g_assemblage_prov_strong, which(degree(g_assemblage_prov_strong) ==
                                                    0))

g_assemblage_artifact_strong <-
  delete.edges(g_assemblage_artifact, which(E(g_assemblage_artifact)$weight < artifact_edgeweight_thresh))

g_assemblage_artifact_strong <-
  delete.vertices(g_assemblage_artifact_strong, which(degree(g_assemblage_artifact_strong) ==
                                                        0))

plot(
  g_assemblage_prov_strong,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Strong Provenience Network"
)

plot(
  g_assemblage_artifact_strong,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Strong Artifact Network"
)


plot(density(degree(g_assemblage_prov_strong)), col = "green",
     main = "Distribution of Strong Provenience Degree")

plot(density(E(g_assemblage_prov_strong)$weight, bw = 1), col = "green",
     main = "Distribution of Strong Provenience Edge Weights")

plot(density(degree(g_assemblage_artifact_strong)), col = "blue",
     main = "Distribution of Strong Artifact Degree")

plot(density(E(g_assemblage_artifact_strong)$weight, bw = 1), col = "blue",
     main = "Distribution of Strong Artifact Edge Weights")

g_assemblages_bpg %>%
  ggraph(layout = "bipartite") +
  geom_edge_link(color = "gray", alpha = 0.25) +
  geom_node_point(aes(color = type), size = 2) +
  scale_color_manual(
    name = "Node Type",
    breaks = c(FALSE, TRUE),
    labels = c("Provenience", "Artifact"),
    values = c("green", "blue")
  ) +
  ggtitle("Bipartite Graph of Provenience and Artifact Type")

g_assemblage_prov %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(color = "green", size = 2) +
  ggtitle("Projected Provenience Associations")

g_assemblage_artifact %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(color = "blue", size = 2) +
  ggtitle("Projected Artifact Type Associations")



