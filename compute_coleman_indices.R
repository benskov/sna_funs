# This function implements the individual-level Coleman's segregation index formula proposed 
# by Bojanowski and Corten (2014). It uses each ego's local network as the sub-graph to 
# compute ego's segregation index.

segregation <- function (G, g, index = "coleman") { 
    # G is the network object, g is a vector with group assignments (coerced to factor)
    
    # Safety moves
        stopifnot(!is.null(g), index == "coleman") # only Coleman index supported so far
        
    # Housekeeping
        if (is.network(G)) {
            A <- as.sociomatrix.sna(G, force.bipartite = TRUE)
            vertices <- network.vertex.names(G)
            N <- network.size(G)
        } else if (is.matrix(G) & !is.null(rownames(G)) & all(rownames(G) == colnames(G))) {
            A <- G
            vertices <- colnames(G)
            N <- nrow(G)
        } else {
            stop("Please, provide valid arguments.")
        }
        diag(A) <- 1 # Set diagonals to 1 (this is taken into account when computing out degrees below)
        out_degs <- rowSums(A) - 1 # 1-by-N vector with out-degrees; 1 is subtracted to counter 1-diagonals
        g <- if (is.null(names(g))) setNames(as.factor(g), vertices) else as.factor(g) # Coerce g into factor
        n <- table(g) # Look-up table with frequencies of each group (used in for loop)
    
    # Create "group-weight" matrix with information on same-group nominations
    W_group <- matrix(0, nrow = N, ncol = N) # An empty matrix we populate below
    for (i in levels(g)) {
        same_group <- which(i == g) # Find all nodes in i'th group
        W_group[same_group, same_group] <- 1 # w_ij set to 1 when i and j are in the same group
    }
    # A weakness of this method is that it's black-and-white, i.e., either you're in group i
    # or you're not. But groups i and j may be more alike than i and k, and there should
    # be a way to reflect this. It's especially true if the grouping variable is actually
    # continuous, so not grouping per se.
    
    # Same-group adjacency matrix by element-wise multiplication
    A_group <- A * W_group 
    
    # Generating the data for output
    o <- data.frame(vertex = vertices, m_exp = NA_integer_, m = NA_integer_)
    subgraphs <- list()
    for (v in vertices) {
        subgraphs[[v]] <- unname(A[rownames(A) == v, ] == 1) # Sub-graph of i and i's alters; unnamed so subsetting possible (below)
        subgraph <- subgraphs[[v]]
        vert <- which(o$vertex == v)
        o[vert, "m_exp"] <- sum(out_degs[subgraph]) * (n[g[v]] - 1)/(N - 1) # i's expected number of same-group nominations 
        o[vert, "m"] <- sum(A_group[subgraph, subgraph]) # i's actual number of same-group nominations
    }
    o$score <- ifelse(o$m >= o$m_exp, 
                      (o$m - o$m_exp)/(sum(out_degs[subgraph]) - o$m_exp),
                      (o$m - o$m_exp)/o$m_exp)
    o$group <- g[vertices]
    o[, -c(2, 3)] # Removes irrelevant columns from output
}

