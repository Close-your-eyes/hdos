#' Sample points from a n-dimensional rectangle prism (cuboid)
#'
#' @param n_samples
#' @param n_dim
#' @param dimnames
#' @param lengths
#' @param n_bins
#' @param origin
#' @param return_coords_only
#' @param return_3D_radii
#' @param return_3D_angles
#' @param inside
#' @param rot_mat
#'
#' @return
#' @export
#'
#' @examples
cuboid_nd <- function(n_samples = 1000,
                      n_dim = 3,
                      dimnames = paste0("x", 1:n_dim),
                      lengths = rep(1, n_dim),
                      n_bins = 9,
                      origin = rep(0, n_dim),
                      return_coords_only = F,
                      return_3D_radii = T,
                      return_3D_angles = T,
                      inside = T,
                      rot_mat = diag(n_dim)) {

  # Check if lengths match the dimensions
  if (length(lengths) != n_dim) {
    stop("The 'lengths' vector must have the same length as the 'n_dim' parameter.")
  }

  name <- paste0(n_dim, "Dcuboid_", paste(lengths, collapse = ""), "_", paste(origin, collapse = ""), "_", n_samples)

  # hashname for coloring of seperate objects later
  hashname <- rlang::hash(list(n_samples,
                               n_dim,
                               sort(dimnames),
                               lengths,
                               origin,
                               inside))

  # Generate lower and upper bounds for the rectangle in each dimension
  bounds <- matrix(c(rep(0, n_dim), lengths), ncol = 2)

  # Generate random points
  points <- matrix(0, nrow = n_samples, ncol = n_dim)
  colnames(points) <- dimnames

  if (inside) {
    # Sample from the inside
    for (i in 1:n_dim) {
      points[, i] <- runif(n_samples, min = bounds[i, 1], max = bounds[i, 2])
    }
  } else {
    # Sample from the inside
    for (j in 1:n_samples) {
      # Choose a random dimension to "fix" for inside sampling
      fixed_dim <- sample(1:n_dim, 1)
      points[j, ] <- sapply(1:n_dim, function(i) {
        if (i == fixed_dim) {
          sample(c(bounds[i, 1], bounds[i, 2]), 1) # Pick one of the boundary values
        } else {
          runif(1, min = bounds[i, 1], max = bounds[i, 2]) # Pick a random value within bounds
        }
      })
    }
  }

  points <- scale(points, scale = F) # center around zero

  if (return_coords_only) {
    if (any(origin != 0)) {
      points <- sweep(points, 2, origin, FUN = "+")
    }
    if (!identical(rot_mat, diag(3)) && n_dim == 3) {
      points <- rotate_space_3D(points, rot_mat)
    }
    return(list(coord = as.data.frame(points), meta = data.frame(name = rep(name, nrow(points)), hash = rep(hashname, nrow(points)))))
  }

  bins <- apply(points, 2, cut, breaks = n_bins, labels = F, simplify = F)
  names(bins) <- paste0(names(bins), "_bin")

  radii_3d <- NULL
  if (return_3D_radii) {
    radii_3d <- dplyr::bind_cols(get_radii_3d_all(points))
    if (ncol(radii_3d) == 1) {
      names(radii_3d) <- strsplit(names(radii_3d), "_")[[1]][1]
    }
  }

  angle_3d <- NULL
  if (return_3D_angles) {
    angle_3d <- dplyr::bind_cols(get_angle_3d_all(points))
    if (ncol(angle_3d) == 2) {
      names(angle_3d) <- sapply(strsplit(names(angle_3d), "_"), "[", 1)
    }
  }

  # do shifting last, it distorts angle and radii calc
  if (any(origin != 0)) {
    points <- sweep(points, 2, origin, FUN = "+")
  }
  if (!identical(rot_mat, diag(3)) && n_dim == 3) {
    points <- rotate_space_3D(points, rot_mat)
  }

  points <- list(coord = as.data.frame(points),
                 meta = do.call(dplyr::bind_cols, list(
                   as.data.frame(bins),
                   radii_3d,
                   angle_3d,
                   data.frame(name = name, hash = hashname)
                 )))

  return(points)
}
