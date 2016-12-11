# This function implements the individual-level Coleman's segregation index formula proposed 
# by Bojanowski and Corten (2014). It uses each ego's local network as the sub-graph to 
# compute ego's segregation index.

coleman_index <- function (G, g) { 
  # G is the network object, g is a factor with group assignments
  
  # Safety move
    stopifnot(is.network(G), is.factor(g)) # Making sure G has correct class, and g is a factor
    
  # Housekeeping
    diag.remove(G, remove.val = 1) # Set diagonals to 1 (this is taken into account when computing out degrees below)
    vertices <- network.vertex.names(G)
    g <- as.factor(g) # Coerce g into factor
    N <- G$gal$n # Network size
    A <- as.sociomatrix.sna(G, force.bipartite = TRUE) # Adjacency matrix, which we use repeatedly in loop
    out_degs <- rowSums(A) - 1 # 1-by-N vector with out-degrees; 1 is subtracted because all diagonals are set to one above
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
  # FIX THIS. HAVE TO HIT THE SACK NOW. CHECK IF IT WORKS
    o <- data.frame(vertex = network.vertex.names(G), m_exp = NA_integer_, m = NA_integer_)
    subgraphs <- list()
    for (v in vertices) {
      subgraphs[[v]] <- unname(A[v, ] == 1) # Sub-graph consists of i and i's alters; unnamed so it can be used for subsetting later on
      # subgraph[i] <- TRUE # Makes i part of its own sub-graph
        # Consider trying out setting diag(A) = 1, so we don't have to do this N times
      # o[i, "m_exp"] <- sum(out_degs[subgraph]) * (n[g[i]] - 1)/(N - 1) # i's expected number of same-group nominations 
      # o[i, "m"] <- sum(A_group[subgraph, subgraph]) # i's actual number of same-group nominations
    }
    o$score <- ifelse(o$m >= o$m_exp, 
                      (o$m - o$m_exp)/(sum(out_degs[subgraphs[o$vertex]]) - o$m_exp,
                      (o$m - o$m_exp)/o$m_exp)
    o$group <- g
    o[, -c("m_exp", "m") # Returns relevant columns of data framE
}

timetest <- data.frame(N = 0, user.self = 0, sys.self = 0, elapsed = 0)
for(N in seq(10, 3000, length.out = 15)) {
  N <- round(N)
  G <- network(rgraph(N)$adj_matrix, directed = TRUE, matrix.type = "adjacency") 
  g <- rep(LETTERS[1:5], N/5+10)[1:N]
  timetest <- rbind(timetest, c(N, system.time(coleman_index(G, g))[1:3]))
}; timetest[-1, ]; plot(timetest$elapsed[-1])

coleman_index() %>% ggplot(aes(x = coleman_segregation_index, fill = group)) + 
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~ group, ncol = 1) + 
  scale_x_continuous(limits = c(-1.1, 1.1))
