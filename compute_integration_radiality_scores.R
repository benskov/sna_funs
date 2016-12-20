# DOCUMENTATION NOT PROPERLY FORMATTED TO ACCOUNT FOR INTEGRATION() AND RADIALITY() INVOKING THE ACTUAL FUNCTION. NEEDS TO BE FIXED.

#' Compute integration and radiality scores
#' 
#' Implements the integration and radiality scores proposed by Valente and 
#' Foreman (1998). Integration is a global measure of inward centrality, 
#' radiality of outward centrality.
#' 
#' These functions draw heavily on \code{sna::geodist()} and may be slow for 
#' large, dense networks.
#' 
#' @param score length-one character indicating which score to compute.
#' @param G a \code{network} object of matrix coericble into one.
#' @param rm.autonom logical; should autonomiations be ignored when computing
#'   the scores?
#'   
#' @return 
#' A tidy data frame with four columns: node names, absolute scores,
#' relative scores, and reversed scores (as percentages).
#' 
#' @examples 
#' nodes <- LETTERS[1:20] # Node names
#' A <- matrix(sample(0:1, 400, TRUE, c(0.75, 0.25)), 20, 20, dimnames = list(nodes, nodes))
#' integration(A)
#' radiality(A)

integration <- function(...) {
    compute_int_rad_scores("i", ...)
}

radiality <- function(...) {
    compute_int_rad_scores("r", ...)
}

compute_int_rad_scores <- function (score = "i", G, rm.autonom = TRUE) { 
    # Internal function; invoked by integration() and radiality()
    
    # Savety moves and housekeeping
    if (!is.matrix(G) & !network::is.network(G)) stop("The function requires a network object or a matrix coercible into one.")
    G <- network::as.network(G)
    if (rm.autonom) sna::diag.remove(G, 0)
    N <- network::network.size(G)
    
    geodesics <- sna::geodist(G, inf.replace = 0)$gdist # max(geodesics) meaningless if inf.replace != 0
    max_geodesic <- max(geodesics)
    reversed_distances <- ifelse(geodesics == 0, 0, (1 + max_geodesic) - geodesics)
        # Generates matrix of reversed distances, in which:
        # - Vertices have RD = 0 to themselves, and
        # - Isolates have RD = 0 to all other vertices
    max_reversed_distance <- max(reversed_distances)
    scores <- if (score == "i") { # Integration scores (default)
        colSums(reversed_distances)/(N - 1)
    } else if (score == "r") { # Radiality score
        rowSums(reversed_distances)/(N - 1) 
    } else {
        stop("Incorrect type of score; must be either integration (i) or radiality (r). Please, use the dedicated functions integration() or radiality() to avoid this problem.")
    }
    
    # Build and return output as tidy data frame:
    o <- data.frame(node = network.vertex.names(G), absolute = scores)
    o$relative <- o$absolute/max_reversed_distance
    o$reversed <- (1 - o$relative) * 100 # As percentages
    o
}
