d2dgeo <- function(D, knn, ciso = FALSE, is_dist = TRUE) {
  
  ## dist knows how to handle missing values and rows are scaled accordingly!
  message(Sys.time(), ": ... to full distance matrix")
  dx <- as.matrix(D)
  
  if (anyNA(dx)) warning("NAs in the distance matrix")
  
  ## Keep only the knn nearest neighbors in distance matrix, ignore points with
  ## distance zero, later we set these distances to a small value so that
  ## igraph::graph_from_adjacency_matrix will not ignore them
  diag(dx) <- NA
  off_diag_zeros <- dx == 0
  dx[dx == 0] <- NA
  message("Sum off diagonal zeros dx: ", sum(off_diag_zeros, na.rm = TRUE))
  message(Sys.time(), ": Computing knn adjacency matrix")
  
  for (i in 1:ncol(dx)) {
    ri <- rank(dx[, i], na.last = TRUE, ties.method = "first")
    dx[ri > knn, i] <- NA
  }
  
  if (ciso) {
    message(Sys.time(), ": CIso, calculate scaling")
    sqrt_mean_knn <- sqrt(colMeans(dx, na.rm = TRUE))
    
    message("min(sqrt_mean_knn) = ", min(sqrt_mean_knn, na.rm = TRUE))
    
    if (anyNA(sqrt_mean_knn))
      stop("There are NAs in the kNN means.")
    
    dx <- dx / (sqrt_mean_knn %o% sqrt_mean_knn)
    message("min(dx) = ", min(dx, na.rm = TRUE))
    if (anyNA(dx))
      warning("There are NAs in the scaled distance matrix.")
  }
  
  ## igraph::graph_from_adjacency_matrix ignores zero weights
  dx[is.na(dx)] <- 0
  dx[off_diag_zeros] <- 1e-5
  
  message(Sys.time(), ": Creating graph")
  gx <- graph_from_adjacency_matrix(
    dx, weighted = TRUE,
    mode = "undirected",
    diag = FALSE)
  if (!igraph::is_connected(gx))
    stop("kNN Graph is not connected, increase knn")
  
  message(Sys.time(), ": Dijkstra")
  dgx <- igraph::distances(gx, algorithm = "dijkstra")
  if (anyNA(dgx))
    stop("There are NAs in the final distance matrix, ",
         "this should not have happened, something went wrong.")
  
  return(dgx)
} 