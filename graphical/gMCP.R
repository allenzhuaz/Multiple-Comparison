library(gMCP)
graphGUI()

###################
#### Chapter 1 ####
###################
graph <- BonferroniHolm(3)
gMCP(graph, pvalues = c(0.01,0.07,0.02))

###################
#### Chapter 2 ####
###################
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

###################
#### Chapter 3 ####
###################

graph <- BretzEtAl2011()
# We can reject a single node:
print(rejectNode(graph, "H11"))

# Or given a vector of pvalues let the function gMCP do all the work:
pvalues <- c(0.1, 0.008, 0.005, 0.15, 0.04, 0.006)
result <- gMCP(graph, pvalues)
print(result)

gMCPReport(result, "Report.tex")

# Estimates:
est <- c("H1"=0.860382, "H2"=0.9161474, "H3"=0.9732953)
# Sample standard deviations:
ssd <- c("H1"=0.8759528, "H2"=1.291310, "H3"=0.8570892)
pval <- c(0.01260, 0.05154, 0.02124)/2
simConfint(BonferroniHolm(3), pvalues=pval,
           confint=function(node, alpha) {
             c(est[node]-qt(1-alpha,df=9)*ssd[node]/sqrt(10), Inf)
           }, estimates=est, alpha=0.025, mu=0, alternative="greater")
## lower bound estimate upper bound
## H1 0.0000 0.8604 Inf
## H2 -0.0076 0.9161 Inf
## H3 0.0000 0.9733 Inf
# Note that the sample standard deviations in the following call
# will be calculated from the pvalues and estimates.
simConfint(BonferroniHolm(3), pvalues=pval,
           confint="t", df=9, estimates=est, alpha=0.025, alternative="greater")
## lower bound estimate upper bound
## [1,] 0.000000 0.8604 Inf
## [2,] -0.007581 0.9161 Inf
## [3,] 0.000000 0.9733 Inf

###################
#### Chapter 5 ####
###################
m <- rbind(H1=c(0, 0, 0.5, 0.5 ),
           H2=c(0, 0, 0.5, 0.5 ),
           H3=c("\\epsilon", 0, 0, "1-\\epsilon"),
           H4=c(0, "\\epsilon", "1-\\epsilon", 0 ))
graph <- matrix2graph(m)
#graph <- improvedParallelGatekeeping()
graph

substituteEps(graph, eps=0.001)

gMCP(graph, pvalues=c(0.02, 0.04, 0.01, 0.02), eps=0.001)

###################
#### Chapter 6 ####
###################
m <- rbind(H1=c(0, 0, 1, 0, 0),
           H2=c(0, 0, 1, 0, 0),
           H3=c(0, 0, 0, 0.9999, 1e-04),
           H4=c(0, 1, 0, 0, 0),
           H5=c(0, 0, 0, 0, 0))
weights <- c(1, 0, 0, 0, 0)
subgraph1 <- new("graphMCP", m=m, weights=weights)
m <- rbind(H1=c(0, 0, 1, 0, 0),
           H2=c(0, 0, 1, 0, 0),
           H3=c(0, 0, 0, 1e-04, 0.9999),
           H4=c(0, 0, 0, 0, 0),
           H5=c(1, 0, 0, 0, 0))
weights <- c(0, 1, 0, 0, 0)
subgraph2 <- new("graphMCP", m=m, weights=weights)
weights <- c(0.5, 0.5)

graph <- new("entangledMCP", subgraphs=list(subgraph1, subgraph2), weights=weights)

###################
#### Chapter 7 ####
###################
cr <- rbind(H1=c(1 , 0.5 , 0.3 , 0.15),
            H2=c(0.5 , 1 , 0.15, 0.3 ),
            H3=c(0.3 , 0.15, 1 , 0.5 ),
            H4=c(0.15, 0.3 , 0.5 , 1 ))

library(bindata)
n1 <- 20
n2 <- 1e4
pvals <- t(replicate(n2,
                     sapply(colSums(rmvbin(n1, margprob = c(0.35, 0.4, 0.25, 0.3), bincorr = cr)),
                            function(x, ...) {binom.test(x, ...)$p.value}, n=n1, alternative="less")
))

graph <- generalSuccessive(gamma=0, delta=0)
out <- graphTest(pvalues=pvals, graph = graph)
extractPower(out)

# Bonferroni adjustment
G <- diag(2)
weights <- c(0.5,0.5)
corMat <- diag(2)+matrix(1,2,2)
theta <- c(1,2)
calcPower(weights, alpha=0.025, G, theta, corMat)

calcPower(weights, alpha=0.025, G, 2*theta, 2*corMat)

graph <- generalSuccessive()
graph



##############################
#### Chapter 9 Case Study ####
##############################
data(hydroquinone)
pvalues <- c()
x <- hydroquinone$micronuclei[hydroquinone$group=="C-"]
for (dose in c("30 mg/kg", "50 mg/kg", "75 mg/kg", "100 mg/kg", "C+")) {
  y <- hydroquinone$micronuclei[hydroquinone$group==dose]
  result <- wilcox.test(x, y, alternative="less", correct=TRUE)
  pvalues <- c(result$p.value, pvalues)
}
pvalues
## [1] 0.004929 0.002634 0.002634 0.004319 0.066255
library(coin, quietly=TRUE)
pvalues <- c()
for (dose in c("30 mg/kg", "50 mg/kg", "75 mg/kg", "100 mg/kg", "C+")) {
  subdata <- droplevels(hydroquinone[hydroquinone$group %in% c("C-", dose),])
  result <- wilcox_test(micronuclei ~ group, data=subdata, distribution="exact")
  pvalues <- c(pvalue(result), pvalues)
}
pvalues
## [1] 0.006061 0.001263 0.001263 0.005051 0.135101
