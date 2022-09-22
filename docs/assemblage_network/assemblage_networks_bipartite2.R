require(igraph)
require(ggraph)

# Import data from file ---------------------------------------------------

dat <- read_csv("Catalog_A2.csv")

# Create bipartite graph from unique provenience & artifact pairs ---------

G_assemblages_bpg <-
  graph_from_data_frame(unique.data.frame(select(dat, LEVEL_ID, CODE)),
                        directed = FALSE)

V(G_assemblages_bpg)$type <-
  bipartite_mapping(G_assemblages_bpg)$type

plot(
  G_assemblages_bpg,
  layout = layout.bipartite,
  vertex.size = 5,
  vertex.label.cex = 0.3
)

# Project one-mode graphs -------------------------------------------------

assemblage_projections <-
  bipartite_projection(G_assemblages_bpg, multiplicity = TRUE)

G_assemblage_prov <- assemblage_projections$proj1

G_assemblage_artifact <- assemblage_projections$proj2

plot(
  G_assemblage_prov,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Provenience Network"
)

plot(
  G_assemblage_artifact,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Artifact Network"
)


# Analyzing one-mode graphs -----------------------------------------------

plot(density(degree(G_assemblage_prov)),
     main = "Distribution of Provenience Degree")

plot(density(E(G_assemblage_prov)$weight, bw = 1),
     main = "Distribution of Provenience Edge Weights")

plot(density(degree(G_assemblage_artifact)),
     main = "Distribution of Artifact Degree")

plot(density(E(G_assemblage_artifact)$weight, bw = 1),
     main = "Distribution of Artifact Edge Weights")


# Thinning out weak edges -------------------------------------------------

prov_edgeweight_thresh <- quantile(E(G_assemblage_prov)$weight, 0.95)

artifact_edgeweight_thresh <- quantile(E(G_assemblage_artifact)$weight, 0.95)

G_assemblage_prov_strong <-
  delete.edges(G_assemblage_prov, which(E(G_assemblage_prov)$weight < prov_edgeweight_thresh))

G_assemblage_prov_strong <-
  delete.vertices(G_assemblage_prov_strong, which(degree(G_assemblage_prov_strong) ==
                                                    0))

G_assemblage_artifact_strong <-
  delete.edges(G_assemblage_artifact, which(E(G_assemblage_artifact)$weight < artifact_edgeweight_thresh))

G_assemblage_artifact_strong <-
  delete.vertices(G_assemblage_artifact_strong, which(degree(G_assemblage_artifact_strong) ==
                                                        0))

plot(
  G_assemblage_prov_strong,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Strong Provenience Network"
)

plot(
  G_assemblage_artifact_strong,
  vertex.size = 5,
  vertex.label = NA,
  layout = layout.fruchterman.reingold,
  main = "Assmeblage Strong Artifact Network"
)


plot(density(degree(G_assemblage_prov_strong)), col = "green",
     main = "Distribution of Strong Provenience Degree")

plot(density(E(G_assemblage_prov_strong)$weight, bw = 1), col = "green",
     main = "Distribution of Strong Provenience Edge Weights")

plot(density(degree(G_assemblage_artifact_strong)), col = "blue",
     main = "Distribution of Strong Artifact Degree")

plot(density(E(G_assemblage_artifact_strong)$weight, bw = 1), col = "blue",
     main = "Distribution of Strong Artifact Edge Weights")

G_assemblages_bpg %>%
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

G_assemblage_prov %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(color = "green", size = 2) +
  ggtitle("Projected Provenience Associations")

G_assemblage_artifact %>%
  ggraph(layout = "fr") +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(color = "blue", size = 2) +
  ggtitle("Projected Artifact Type Associations")



