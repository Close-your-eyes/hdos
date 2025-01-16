#' Sample points from a n-dimensional sphere
#'
#' @param n_samples
#' @param n_dim
#' @param dimnames
#' @param density_fun
#' @param radius
#' @param ellipsoid_scaling
#' @param origin
#' @param inside
#' @param n_bins
#' @param return_3D_angles
#' @param return_3D_radii
#' @param return_coords_only
#'
#' @return
#' @export
#'
#' @examples
sphere_nd <- function(n_samples = 1000,
                      n_dim = 3,
                      dimnames = paste0("x", 1:n_dim),
                      density_fun = function(r) 3*r^2,
                      radius = 1,
                      ellipsoid_scaling = rep(1, n_dim),
                      origin = rep(0, n_dim),
                      inside = T,
                      n_bins = 9,
                      return_3D_angles = T,
                      return_3D_radii = T,
                      return_coords_only = F) {

  if (length(ellipsoid_scaling) != n_dim) {
    message("adjusting len of ellipsoid_scaling.")
    if (length(ellipsoid_scaling) < n_dim) {
      ellipsoid_scaling <- c(ellipsoid_scaling, rep(1, n_dim-length(ellipsoid_scaling)))
    } else {
      ellipsoid_scaling <- ellipsoid_scaling[1:n_dim]
    }
  }
  if (any(ellipsoid_scaling <= 0)) {
    message("ellipsoid_scaling must be greater than zero. setting values to one.")
    ellipsoid_scaling[which(ellipsoid_scaling <= 0)] <- 1
  }
  if (length(origin) != n_dim) {
    message("adjusting len of origin.")
    if (length(origin) < n_dim) {
      origin <- c(origin, rep(1, n_dim-length(origin)))
    } else {
      origin <- origin[1:n_dim]
    }
  }

  # hashname for coloring of separate objects later
  hashname <- rlang::hash(list(n_samples,
                               n_dim,
                               sort(dimnames),
                               density_fun,
                               radius,
                               ellipsoid_scaling,
                               origin,
                               inside))
  hashname <- substr(hashname, 1, 6)

  if (any(ellipsoid_scaling != 1)) {
    name <- paste0(n_dim, "Dellips_", paste(round(ellipsoid_scaling,1), collapse = ""), "_", paste(round(origin,1), collapse = ""), "_", n_samples)
  } else {
    name <- paste0(n_dim, "Dsphere_r", round(radius,1), "_", paste(round(origin,1), collapse = ""), "_", n_samples)
  }

  ## this function calculates everything based on random points generated initially

  # Generate standard normal random points
  points <- matrix(rnorm(n_samples * n_dim), ncol = n_dim)
  colnames(points) <- dimnames

  # Normalize to lie on the unit sphere
  magnitudes <- sqrt(rowSums(points^2))
  points <- points / magnitudes  # Points on the unit sphere


  if (inside) {
    # Scale points to be inside the sphere uniformly
    #radii <- runif(n_samples)^(1 / n_dim) * radius  # Radial distribution for uniformity

    radii <- generate_radii(n_samples, density_fun) # same as sqrt(rowSums(points^2)), also for higher dimensions
    radii <- radii * radius
    points <- points * radii
  } else {
    # Scale to the surface of the sphere
    points <- points * radius
    #radii <- rep(radius, n_samples)
  }

  # scale to ellipsoid
  if (any(ellipsoid_scaling != 1)) {
    points <- points %*% diag(ellipsoid_scaling)
    colnames(points) <- dimnames
    # recalc radii, needed after elipsoid scaling
    # radii <- sqrt(rowSums(points^2)) # without elipsoid scaling, this is the same as radii <- generate_radii(n_samples, density_fun)
  }

  ## correct for minor deviation from origin, center only
  points <- scale(points, scale = F)

  if (return_coords_only) {
    if (any(origin != 0)) {
      points <- sweep(points, 2, origin, FUN = "+")
    }
    points <- list(coord = tibble::as_tibble(points), meta = tibble::tibble(name = rep(name, nrow(points)), hash = rep(hashname, nrow(points))))
    attributes(points) <- list(type = "sphere", center = origin, radius = radius, names = names(points))
    return(points)
  }

  bins <- apply(points, 2, cut, breaks = n_bins, labels = F, simplify = F)
  bins <- lapply(bins, function(x) factor(x, levels = sort(unique(x))))
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

  #dens_3d <- get_density_3d_all(points)

  points <- list(coord = tibble::as_tibble(points),
                 meta = do.call(dplyr::bind_cols, list(tibble::as_tibble(bins),
                                                       #r = radii,
                                                       radii_3d,
                                                       angle_3d,
                                                       tibble::tibble(hash = hashname, name = name))))

  attributes(points) <- list(type = "sphere", center = origin, radius = radius, names = names(points))
  return(points)
}
