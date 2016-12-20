#' Compute Coleman's homophily index
#' 
#' \code{colemanindex} returns the individual-level Coleman's homophily indices 
#' proposed by Bojanowski and Corten (2014).
#' 
#' The number of groups has no upper limits. Increasing the number of groups
#' comes without computational costs.
#' 
#' The current version is binary: either you're in group i or you're not. Groups
#' i and j, however, may be more alike than i and k, which should probably be
#' reflected in future versions. This is especially true with continuous 
#' "grouping" variables.
#' 
#' @param G an \code{network} object of square adjacency matrix with named 
#'   dimensions. Dimnames are considered vertex names.
#' @param g a vector with group assignments. Preferably, this vector is named 
#'   with vertex names. If unnamed, elements must be in the same order as 
#'   vertices in the adjacency matrix, to produce unambigious group assignments.
#'   The length of \code{g} must be equal to the network size.
#' @param inf.replace the homophily index of disconnected nodes; by default, 
#'   this is equal to \code{Inf}.
#'   
#' @return 
#' A tidy data frame with three columns: node name, homophily index, and
#' group assignment.
#' 
#' @references 
#' Bojanowski, M. and Corten, R. (2014) **Measuring segregation in
#' social networks**. Social Networks. 2014;39:14-32.
#' 
#' @examples 
#' nodes <- LETTERS[1:20] # Node names
#' A <- matrix(sample(0:1, 400, TRUE, c(0.75, 0.25)), 20, 20, dimnames = list(nodes, nodes))
#' diag(A) <- 0 # No auto-nominations
#' groups <- setNames(c(rep("smoker", 10), rep("non-smoker", 10)), nodes)
#' colemanindex(A, groups, inf.replace = 0)

colemanindex <- function (G, g, inf.replace = Inf) { 
    # Savety moves and housekeeping
    if (is.network(G)) {
        A <- sna::as.sociomatrix.sna(G)
        vertices <- sna::network.vertex.names(G)
        N <- sna::network.size(G)
    } else if (is.matrix(G) & nrow(G) == ncol(G) & !is.null(rownames(G)) & 
               all(rownames(G) == colnames(G))) {
        A <- G
        vertices <- colnames(G)
        N <- nrow(G)
    } else {
        stop("Please, provide valid arguments.")
    }
    stopifnot(!is.null(g), length(g) == N)
    diag(A) <- 1 # This way, we don't have to include i in its own sub-graph in each loop below
    out_degs <- rowSums(A) - 1 # Out-degrees; 1 is subtracted to counter 1-diagonal
    g <- if (!is.null(names(g))) as.factor(g) else setNames(as.factor(g), as.character(vertices))
    n <- table(g) # Look-up table with group frequencies; used in for loop
    
    # Matrix of same-group nominations
    W_group <- matrix(0, nrow = N, ncol = N) # Populated below
    for (i in levels(g[vertices])) {
        same_group <- which(i == g) # Find all nodes in i'th group
        W_group[same_group, same_group] <- 1 # w_ij set to 1 when i and j are in the same group
    }
    
    # Same-group adjacency matrix
    A_group <- A * W_group 
    
    # Calculating indices
    indices <- list()
    for (v in vertices) {
        subgraph <- which(A[v, ] == 1)
        if (sum(out_degs[subgraph]) == 0) { 
            # We need to test this somewhere, so better here so we can skip the 
            # rest of the loop for isolates
            indices[[v]] <- inf.replace
            next
        }
        m_exp <- sum(out_degs[subgraph]) * (n[g[v]] - 1)/(N - 1) # expected number of same-group nominations in i's local network
        m <- sum(A_group[subgraph, subgraph]) # actual number of same-group nominations in i's local network
        indices[[v]] <- (m - m_exp)/ifelse(m >= m_exp, sum(out_degs[subgraph]) - m_exp, m_exp)
    }
    
    # Output 
    data.frame(node = vertices, 
               index = unlist(indices[vertices]), 
               group = g[vertices],
               row.names = NULL)
}
