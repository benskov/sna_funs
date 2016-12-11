integration <- function(...) {
  compute_int_rad_scores("i", ...)
}

radiality <- function(...) {
  compute_int_rad_scores("r", ...)
}

compute_int_rad_scores <- function (score = "i", G, rm.autonom = TRUE) { 
  # Internal function; used by integration() and radiality()
  
  # Implements integration and radiality scores, as proposed by Valente and Foreman (1998) in 
  # Social Networks 20 (1998); 89-105
  
  stopifnot(is.network(G)) 
  if (rm.autonom) sna::diag.remove(G, 0) # Removing autonominations
  
  geodesics <- sna::geodist(G, inf.replace = 0)$gdist # max(geodesics) fails if inf.replace != 0
  max_geodesic <- max(geodesics)
  reversed_distances <- ifelse(geodesics == 0, 0, (1 + max_geodesic) - geodesics)
    # Generates matrix of reversed distances, in which:
    # - Vertices have RD = 0 to themselves, and
    # - Isolates have RD = 0 to all other vertices
  max_reversed_distance <- max(reversed_distances)
  scores <- if (score == "r") { # Radiality score
    rowSums(reversed_distances)/(G$gal$n - 1) 
  } else if (score == "i") { # Integration scores (default)
    colSums(reversed_distances)/(G$gal$n - 1)
  } else {
    stop("Incorrect score type; must be either integration (i) or radiality (r).")
  } # G$gal$n = number of vertices
  
  # Build and return output as tidy data frame:
  o <- data.frame(ego = network.vertex.names(G), absolute = scores)
  o$relative <- o$absolute/max_reversed_distance
  o$reversed <- (1 - o$relative) * 100 # As percentages
  o
}

timetest <- data.frame(N = 0, user.self = 0, sys.self = 0, elapsed = 0)
for(N in seq(10, 3000, length.out = 15)) {
  N <- round(N)
  G <- network(rgraph(N)$adj_matrix, directed = TRUE, matrix.type = "adjacency") 
  timetest <- rbind(timetest, c(N, system.time(integration(G))[1:3]))
}; timetest[-1, ]; plot(timetest$elapsed[-1], log = "xy")