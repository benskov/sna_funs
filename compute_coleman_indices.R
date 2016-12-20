# This function implements the individual-level Coleman's segregation index formula proposed 
# by Bojanowski and Corten (2014). It uses each ego's local network as the sub-graph to 
# compute ego's segregation index. It extends the index to k groups instead of just two.

# A weakness of this function is that it's black-and-white, i.e., either you're in group i
# or you're not. But groups i and j may be more alike than i and k, and there should
# be a way to reflect this. It's especially true if the grouping variable is actually
# continuous, so not grouping per se.

segregation <- function (G, g, inf.replace = Inf, index = "coleman") { 
    # G is the network object or a square matrix
    # g is a vector with group assignments, preferably named to have unambigious group assignment to each vertix
    # inf.replace is used for isolates; trying to follow the same line of thought as sna::geodist()
    # only Coleman index supported so far
    
    # Savety moves and housekeeping
    stopifnot(!is.null(g), index == "coleman")
    if (is.network(G)) {
        A <- as.sociomatrix.sna(G, force.bipartite = TRUE)
        vertices <- network.vertex.names(G)
        N <- network.size(G)
    } else if (is.matrix(G) & nrow(G) == ncol(G) & !is.null(rownames(G)) & 
               all(rownames(G) == colnames(G))) {
        A <- G
        vertices <- colnames(G)
        N <- nrow(G)
    } else {
        stop("Please, provide valid arguments.")
    }
    diag(A) <- 1 # Now we don't have to include i in its own sub-graph in each loop below
    out_degs <- rowSums(A) - 1 # out_degrees; 1 is subtracted to counter 1-diagonal
    g <- if (!is.null(names(g))) as.factor(g) else setNames(as.factor(g), as.character(vertices))
    n <- table(g) # Look-up table with group frequencies; used in for loop
    
    # Matrix of same-group nominations
    W_group <- matrix(0, nrow = N, ncol = N) # Populated below
    for (i in levels(g[vertices])) {
        same_group <- which(i == g) # Find all nodes in i'th group
        W_group[same_group, same_group] <- 1 # w_ij set to 1 when i and j are in the same group
    }
    
    # Same-group adjacency matrix by element-wise multiplication
    A_group <- A * W_group 
    
    # Calculating indices
    indices <- list()
    for (v in vertices) {
        subgraph <- which(A[v, ] == 1)
        if (sum(out_degs[subgraph]) == 0) {
            indices[[v]] <- inf.replace
            next
        }
        m_exp <- sum(out_degs[subgraph]) * (n[g[v]] - 1)/(N - 1) # i's expected number of same-group nominations 
        m <- sum(A_group[subgraph, subgraph]) # i's actual number of same-group nominations
        indices[[v]] <- (m - m_exp)/ifelse(m >= m_exp, sum(out_degs[subgraph]) - m_exp, m_exp)
    }
    
    # Output 
    data.frame(vertex = network.vertex.names(G), 
               index = unlist(indices[vertices]), 
               group = g[vertices])
}
