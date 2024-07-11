rm(list = ls())

# Libraries ---------------------------------------------------------------
require(visNetwork, quietly = TRUE)
library(igraph)


# Data Preparation --------------------------------------------------------

# Create datasetnorganisedinto dataframe 
nnodes <- 10
nnedges <- 20


#Nodes
nodes <- data.frame(id = 1:nnodes, 
                    # add labels on nodes
                    label = paste("Node", 1:10),
                    
                  
                    # add groups on nodes 
                    # should be 2,5 or 10 (divider for 10)
                    group = c("Hazard", "Top Event", "Consequences", "Threats", 
                               "Barriers"),
                    # size adding value
                    value = 1:10, 
                    # Control shape of nodes
                    # should be 2,5 or 10 (divider for 10)
                    # String. Default to 'ellipse'. The shape defines what the node looks like. 
                    # There are two types of nodes. One type has the label inside of it 
                    # and the other type has the label underneath it. 
                    # The types with the label inside of it are: 
                    # ellipse, circle, database, box, text. 
                    # The ones with the label outside of it are: 
                    # image, circularImage, diamond, dot, star, triangle, triangleDown, hexagon, square and icon.
                    shape = c("square", "triangle", "box", "circle", "dot", "star",
                              "ellipse", "database", "text", "diamond"),
                    value = 1:10,
                    # tooltip (html or character), when the mouse is above
                    title = paste0("<p><b>", 1:10,"</b><br>Node !</p>"),
                    # color
                    # should be 2,5 or 10 (divider for 10)
                    color = c("darkred", "grey", "orange", "darkblue", "purple"),
                    # shadow
                    shadow = c(FALSE, TRUE, FALSE, TRUE, TRUE))   
#id has to be the same like from and to columns in edges
#nodes$id <- nodes$label

#Edges
edges <- data.frame(from = sample(1:nnodes, nnedges, replace = T),
                    to = sample(1:nnodes, nnedges, replace = T),
                    width = sample(1:nnodes, nnedges, replace = T))
visNetwork(nodes, edges, width = "100%")
# with defaut layout
visNetwork(nodes, edges, height = "500px") %>%
  visIgraphLayout() %>%
  visNodes(size = 10)