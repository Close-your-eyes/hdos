#' Sample points from a 3-dimensional torus
#'
#' @param n_samples
#' @param R
#' @param r
#' @param amplitude
#' @param frequency
#' @param n_dim
#' @param dimnames
#' @param scaling
#' @param origin
#' @param inside
#' @param n_bins
#' @param return_3D_angles
#' @param return_coords_only
#' @param rot_mat
#' @param return_3D_radii
#'
#' @return
#' @export
#'
#' @examples
torus_3d <- function(n_samples = 1000,
                     R = 3,
                     r = 1,
                     amplitude = 0,
                     frequency = 2,
                     n_dim = 3,
                     dimnames = paste0("x", 1:n_dim),
                     scaling = rep(1, n_dim),
                     origin = rep(0, n_dim),
                     inside = T,
                     n_bins = 9,
                     return_3D_angles = T,
                     return_3D_radii = T,
                     return_coords_only = F,
                     rot_mat = diag(n_dim)) {

  if (n_dim != 3) {
    stop("n_dim must be 3.")
  }

  if (length(scaling) != n_dim) {
    message("adjusting len of scaling.")
    if (length(scaling) < n_dim) {
      scaling <- c(scaling, rep(1, n_dim-length(scaling)))
    } else {
      scaling <- scaling[1:n_dim]
    }
  }
  if (length(origin) != n_dim) {
    message("adjusting len of origin.")
    if (length(origin) < n_dim) {
      origin <- c(origin, rep(1, n_dim-length(origin)))
    } else {
      origin <- origin[1:n_dim]
    }
  }

  name <- paste0(n_dim, "Dtorus_r", R, r, "_", paste(scaling, collapse = ""), "_", round(amplitude,1), "_", round(frequency, 1), "_", paste(origin, collapse = ""), "_", n_samples)

  # hashname for coloring of separate objects later
  hashname <- rlang::hash(list(n_samples,
                               sort(dimnames),
                               amplitude,
                               frequency,
                               r,
                               R,
                               scaling,
                               origin,
                               inside))


  # Generate random angles
  theta <- runif(n_samples, 0, 2 * pi)  # Around the central circle
  phi <- runif(n_samples, 0, 2 * pi)    # Within the tube cross-section

  # Add sine variation to the z-coordinate
  z_sine <- amplitude * sin(frequency * theta)

  if (!inside) {
    x <- (R + r * cos(phi)) * cos(theta)
    y <- (R + r * cos(phi)) * sin(theta)
    z <- r * sin(phi) + z_sine
  } else {
    # Convert to Cartesian coordinates
    # Random radii within the tube cross-section
    rho <- sqrt(runif(n_samples, 0, 1)) * r  # sqrt ensures uniform distribution in the circle
    x <- (R + rho * cos(phi)) * cos(theta)
    y <- (R + rho * cos(phi)) * sin(theta)
    z <- rho * sin(phi) + z_sine  # Add sine wave to z
  }

  # Combine into a matrix
  points <- matrix(c(x, y, z), ncol = 3)
  colnames(points) <- dimnames
  if (any(scaling != 1)) {
    points <- points %*% diag(scaling) # this (matmult) removes colnames of points
  }
  colnames(points) <- dimnames

  ## correct for minor deviation from origin, center only
  #points <- scale(points, scale = F) # d

  if (return_coords_only) {
    #radii <- sqrt(rowSums(points^2))
    if (any(origin != 0)) {
      points <- sweep(points, 2, origin, FUN = "+")
    }
    if (!identical(rot_mat, diag(3))) {
      points <- rotate_space_3D(points, rot_mat)
    }
    points <- list(coord = as.data.frame(points), meta = data.frame(name = rep(name, nrow(points)), hash = rep(hashname, nrow(points))))
    attributes(points) <- list(type = "torus", center = origin, radius_major = R, radius_minor = r, names = names(points))
    return(points)
  }

  bins <- apply(points, 2, cut, breaks = n_bins, labels = F, simplify = F)
  names(bins) <- paste0(names(bins), "_bin")

  # radius
  #radii <- sqrt(rowSums(points^2))

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
    names(angle_3d) <- sapply(strsplit(names(angle_3d), "_"), "[", 1)
  }

  # do shifting last, it distorts angle and radii calc
  if (any(origin != 0)) {
    points <- sweep(points, 2, origin, FUN = "+")
  }

  if (!identical(rot_mat, diag(3))) {
    points <- rotate_space_3D(points, rot_mat)
  }

  normal_vec <- get_torus_normal_vec_by_pca(points)

  points <- list(coord = as.data.frame(points),
                 meta = do.call(dplyr::bind_cols, list(as.data.frame(bins),
                                                       angle_3d,
                                                       radii_3d,
                                                       data.frame(name = name, hash = hashname)))) # r = radii
  attributes(points) <- list(
    type = "torus",
    center = origin,
    radius_major = R,
    radius_minor = r,
    normal_vector = normal_vec,
    names = names(points)
  )

  return(points)
}
