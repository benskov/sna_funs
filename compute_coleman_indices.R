#' Compute Coleman's homophily index
#' 
#' \code{colemanindex} returns the individual-level Coleman's homophily indices 
#' proposed by Bojanowski and Corten (2014).
#' 
#' The function supports lists of graph objects. For now, the vector with group assignments is used for all networks, 
#' so it should (a) hold group assignments for \italic{all} nodes, and (b) do so with unique node names. This
#' is obviously suboptimal, and will be fixed later.
#' 
#' The number of groups has no upper limits. Increasing the number of groups
#' comes without real computational costs.
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
#' @param method string; should the function use an edge list (default) or matrix format. The prior should be efficient, the latter faster.
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
#' A <- matrix(sample(0:1, 400, TRUE, c(0.75, 0.25)), 20, 20, dimnames = list(nodes, nodes)) # Toy network
#' diag(A) <- 0 # Remove accidental auto-nominations
#' groups <- setNames(c(rep("smoker", 10), rep("non-smoker", 10)), nodes)
#' colemanindex(A, groups, inf.replace = 0)

colemanindex <- function (G, g, method = "edgelist", inf.replace = Inf) { 
    if ((is.list(G) & !any(class(G) == "network")) & (!(any(class(G) %in% 
        c("network", "matrix.csr", "matrix.csc", "matrix.ssr", "matrix.ssc", 
          "matrix.hb", "data.frame"))))) {
        return(lapply(G, colemanindex, g = g, method = method, inf.replace = inf.replace))
    }
    
    # Utility function
    compute_index <- function(m, sum_subgraph_out_degs, group_freq) {
        if (sum_subgraph_out_degs == 0) return(inf.replace)
        m_exp <- sum_subgraph_out_degs * (group_freq - 1)/(N - 1) # expected number of same-group nominations in i's local network
        (m - m_exp)/ifelse(m >= m_exp, sum_subgraph_out_degs - m_exp, m_exp)
    }
    
    if (method == "edgelist") {
        if (!sna::is.edgelist.sna(G)) G <- sna::as.edgelist.sna(G)
        el <- as.data.frame(G)
        vertices <- attr(G, "vnames")
        stopifnot(!is.null(vertices))
        el[, 1] <- vertices[el[, 1]]
        el[, 2] <- vertices[el[, 2]]
        N <- attr(G, "n")
        stopifnot(!is.null(g), length(g) == N)
        g <- if (!is.null(names(g))) as.factor(g) else setNames(as.factor(g), as.character(vertices))
        n <- table(g) # Look-up table with group frequencies; used in for loop
        
        out_degs <- table(el[, 1])
        el$sg <- g[el[, 1]] == g[el[, 2]] # Are sender and receiver in the same group?
        
        # Calculating indices
        indices <- list()
        for (v in vertices) {
            subgraph <- c(v, el[el[, 1] == v, 2])
            m <- sum(el[el[, 1] %in% subgraph & el[, 2] %in% subgraph, "sg"])
            indices[[v]] <- compute_index(m, sum(out_degs[subgraph], na.rm = TRUE), n[g[v]])
        }
    } else if (method == "matrix") {
        if (!is.matrix(G)) G <- network::as.sociomatrix(G)
        if (nrow(G) == ncol(G) & !is.null(rownames(G)) & all(rownames(G) == colnames(G))) {
            A <- G
            vertices <- colnames(G)
            N <- nrow(G)
        } else {
            stop("Please, provide valid arguments.")
        }
        stopifnot(!is.null(g), length(g) == N)
        out_degs <- rowSums(A) # Out-degrees
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
            subgraph <- which(A[v, ] == 1 | colnames(A) == v)
            m <- sum(A_group[subgraph, subgraph]) # actual number of same-group nominations in i's local network
            indices[[v]] <- compute_index(m, sum(out_degs[subgraph], na.rm = TRUE), n[g[v]])
        }
    } else {
        stop("Choose either edgelist or matrix method.")
    }
    
    # Output 
    data.frame(node = vertices, 
               index = unlist(indices[vertices]), 
               group = g[vertices],
               row.names = NULL)
}
