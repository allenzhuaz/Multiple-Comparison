library(gMCP)
graphGUI()


# Chapter 1
graph <- BonferroniHolm(3)
gMCP(graph, pvalues = c(0.01,0.07,0.02))

# Chapter 2
m <- rbind(H11=c(0, 0.5, 0, 0.5, 0, 0 ),
           H21=c(1/3, 0, 1/3, 0, 1/3, 0 ),
           H31=c(0, 0.5, 0, 0, 0, 0.5),
           H12=c(0, 1, 0, 0, 0, 0 ),
           H22=c(0.5, 0, 0.5, 0, 0, 0 ),
           H32=c(0, 1, 0, 0, 0, 0 ))
graph <- matrix2graph(m)
graph <- setWeights(graph, c(1/3, 1/3, 1/3, 0, 0, 0))
graph

graph@nodeAttr$X <- c(H11=100, H21=300, H31=500, H12=100, H22=300, H32=500)
graph@nodeAttr$Y <- c(H11=100, H21=100, H31=100, H12=300, H22=300, H32=300)

graph <- placeNodes(graph, nrow=2)

cat(graph2latex(graph))
edgeAttr(graph, "H11", "H21", "labelX") <- 200
edgeAttr(graph, "H11", "H21", "labelY") <- 80

graphGUI("graph")

# Chapter 3
