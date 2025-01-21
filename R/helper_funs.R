get_density_2d <- function(x, y, ...) {
  if (missing(y) && (is.matrix(x) || is.data.frame(x))) {
    y <- x[,2]
    x <- x[,1]
  }
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(list(raw = dens$z[ii], scaled = dens$z[ii]/max(dens$z[ii])))
}
get_density_3d <- function(x, y, z) {
  if (missing(y) && missing(z) && (is.matrix(x) || is.data.frame(x))) {
    y <- x[,2]
    z <- x[,3]
    x <- x[,1]
  }
  dens <- kde3d(x,y,z)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  iz <- findInterval(z, dens$z)
  ii <- cbind(cbind(ix, iy), iz)
  return(list(raw = dens$d[ii], scaled = dens$d[ii]/max(dens$d[ii])))
}
get_density_4d <- function(x, y, z, w) {
  if (missing(y) && missing(z) && missing(w) &&  (is.matrix(x) || is.data.frame(x))) {
    y <- x[,2]
    z <- x[,3]
    w <- x[,4]
    x <- x[,1]
  }
  dens <- kde4d(x,y,z,w)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  iz <- findInterval(z, dens$z)
  iw <- findInterval(w, dens$w)
  ii <- cbind(cbind(cbind(ix, iy), iz), iw)
  return(list(raw = dens$d[ii], scaled = dens$d[ii]/max(dens$d[ii])))
}

# Function to calculate polar and azimuthal angles
# radial_density_function (rdf)
rdf_uniform <- function(r) 3*r^2
rdf_low_center <- function(r) 5*r^4
rdf_linear_incr_to_center <- function(r) 6*r
rdf_high_center <- function(r) 3*(r-1)^2

general_rdf <- function(k) {
  # higher k --> higher density towards outside
  function(r) {
    (k+1)*r^k
  }
}

# Function to calculate CDF from the density function
calculate_cdf <- function(r, density_function) {
  integrate(density_function, 0, r)$value
}

# Generate inverse CDF using the density function
generate_radii <- function(n, density_function) {
  # Precompute CDF values
  cdf_values <- sapply(seq(0, 1, length.out = 1000), function(r) calculate_cdf(r, density_function))
  r_values <- seq(0, 1, length.out = 1000)

  # Normalize CDF to make it go from 0 to 1
  cdf_values <- cdf_values / max(cdf_values)

  # Map uniform random numbers through inverse CDF
  u <- runif(n, 0, 1)
  approx(cdf_values, r_values, xout = u)$y  # Interpolation
}

kde3d <- function (x, y, z, h, n = 20, lims = c(range(x), range(y), range(z))) {
  # from misc3d::kde3d
  nx <- length(x)
  if (length(y) != nx || length(z) != nx)
    stop("data vectors must be the same length")
  if (missing(h))
    h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y),
           MASS::bandwidth.nrd(z))/6
  else if (length(h) != 3)
    h <- rep(h, length = 3)
  if (length(n) != 3)
    n <- rep(n, length = 3)
  if (length(lims) == 2)
    lims <- rep(lims, length = 6)
  gx <- seq(lims[1], lims[2], length = n[1])
  gy <- seq(lims[3], lims[4], length = n[2])
  gz <- seq(lims[5], lims[6], length = n[3])
  mx <- matrix(outer(gx, x, dnorm, h[1]), n[1], nx)
  my <- matrix(outer(gy, y, dnorm, h[2]), n[2], nx)
  mz <- matrix(outer(gz, z, dnorm, h[3]), n[3], nx)
  v <- array(0, n)
  tmy.nx <- t(my)/nx
  for (k in 1:n[3]) {
    tmy.nz.zk <- tmy.nx * mz[k, ]
    v[, , k] <- mx %*% tmy.nz.zk
  }
  return(list(x = gx, y = gy, z = gz, d = v))
}

kde4d <- function(x, y, z, w, h, n = 20, lims = c(range(x), range(y), range(z), range(w))) {
  # from Chatgpt
  # Ensure all data vectors have the same length
  nx <- length(x)
  if (length(y) != nx || length(z) != nx || length(w) != nx)
    stop("All data vectors must be the same length")

  # Set bandwidths, default to scaled bandwidth estimators
  if (missing(h))
    h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y),
           MASS::bandwidth.nrd(z), MASS::bandwidth.nrd(w)) / 6
  else if (length(h) != 4)
    h <- rep(h, length = 4)

  # Ensure grid size and limits are properly defined
  if (length(n) != 4)
    n <- rep(n, length = 4)
  if (length(lims) == 2)
    lims <- rep(lims, length = 8)

  # Generate the grid for each dimension
  gx <- seq(lims[1], lims[2], length = n[1])
  gy <- seq(lims[3], lims[4], length = n[2])
  gz <- seq(lims[5], lims[6], length = n[3])
  gw <- seq(lims[7], lims[8], length = n[4])

  # Compute the KDE values
  mx <- matrix(outer(gx, x, dnorm, h[1]), n[1], nx)
  my <- matrix(outer(gy, y, dnorm, h[2]), n[2], nx)
  mz <- matrix(outer(gz, z, dnorm, h[3]), n[3], nx)
  mw <- matrix(outer(gw, w, dnorm, h[4]), n[4], nx)

  # Initialize the 4D density array
  v <- array(0, n)

  # Compute the density in 4D space
  tmy.nx <- t(my) / nx
  for (l in 1:n[4]) {
    tmw.nx <- tmy.nx * mw[l, ]
    for (k in 1:n[3]) {
      tmz.nx <- tmw.nx * mz[k, ]
      v[, , k, l] <- mx %*% tmz.nx
    }
  }

  return(list(x = gx, y = gy, z = gz, w = gw, d = v))
}

calculate_3d_angles <- function(x, y, z, r) {

  if (missing(y) && missing(z) && (is.matrix(x) || is.data.frame(x))) {
    y <- x[,2]
    z <- x[,3]
    x <- x[,1]
  }

  if (missing(r)) {
    r <- sqrt(x^2 + y^2 + z^2)
  }

  # Check if the point is at the origin
  if (any(r == 0)) {
    stop("Point is at the origin; angles are undefined.")
  }

  # Calculate polar angle (theta)
  theta <- acos(z / r)

  # Calculate azimuthal angle (phi)
  phi <- atan2(y, x)

  # Return the results as a list
  return(data.frame(phi, theta))
}

#' Title
#'
#' @param data
#' @param x
#' @param y
#' @param z
#' @param color
#' @param showlegend
#'
#' @return
#' @export
#'
#' @examples
plotly3dplot <- function(data,
                         x,
                         y,
                         z,
                         color,
                         showlegend = T) {

  if (is.numeric(data[["meta"]][[color]])) {
    color_mode <- "c"
  } else {
    color_mode <- "d"
  }
  color_mode <- match.arg(color_mode, c("continuous", "discrete"))

  if (is.list(data) && identical(names(data), c("coord", "meta"))) {
    data <- dplyr::bind_cols(data)
  }

  if (missing(x)) {
    x <- names(data)[1]
  }
  if (missing(y)) {
    y <- names(data)[2]
  }
  if (missing(z)) {
    z <- names(data)[3]
  }

  if (missing(color)) {
    plotly::plot_ly(data,
                    x = data[[x]], y = data[[y]], z = data[[z]],
                    type = "scatter3d", mode = "markers",
                    size = 2)
  } else if (color_mode == "discrete") {
    plotly::plot_ly(data,
                    x = data[[x]], y = data[[y]], z = data[[z]],
                    type = "scatter3d", mode = "markers",
                    size = 2,
                    # reverse to match ggplot
                    colors = rev(as.character(fcexpr::col_pal("custom", n = length(unique(data[[color]])), direction = -1))),
                    color = data[[color]],
                    showlegend = showlegend)
  } else if (color_mode == "continuous") {
    plotly::plot_ly(data,
                    x = data[[x]], y = data[[y]], z = data[[z]],
                    type = "scatter3d", mode = "markers",
                    marker = plotly_marker(data = data,
                                           color = color,
                                           size = 2,
                                           colorscale = plotly_colorscale(fcexpr::col_pal("RColorBrewer::Spectral", direction = -1)),
                                           showscale = showlegend))
  }
}

plotly_colorscale <- function(colors) {
  n <- length(colors)
  stops <- seq(0, 1, length.out = n) # Evenly spaced values from 0 to 1
  lapply(seq_along(colors), function(i) list(stops[i], colors[i]))
}

plotly_marker <- function(data, color, ...) {
  list(color = data[[color]],
       cmin = min(data[[color]]),
       cmax = max(data[[color]]),
       ...)

}

get_angle_3d_all <- function(points) {
  dims <- colnames(points)
  combs <- combn(dims, 3, simplify = F)
  names(combs) <- sapply(combs, paste, collapse = "")
  # check if center is 0
  # Using base R
  # abs(x - y) < sqrt(.Machine$double.eps)
  if (!any(dplyr::near(colMeans(points), 0))) {
    message("centering before 3d angle calc.")
    points <- scale(points, scale = F)
  }
  lapply(combs, function(cmb) {
    ret <- calculate_3d_angles(points[, cmb])
    names(ret) <- paste0(names(ret), "_", paste(cmb, collapse = ""))
    return(ret)
  })
}

get_radii_3d_all <- function(points) {
  colcombs <- combn(colnames(points), 3, simplify = F)
  if (!any(dplyr::near(colMeans(points), 0))) {
    message("centering before 3d radii calc.")
    points <- scale(points, scale = F)
  }
  radii_3d <- lapply(colcombs, function(cmb) {
    sqrt(rowSums(points[,cmb]^2)) # scale if sweeping happens before
  })
  names(radii_3d) <- paste0("r_", sapply(colcombs, paste, collapse = ""))
  return(radii_3d)
}

get_density_3d_all <- function(points) {
  dims <- colnames(points)
  combs <- combn(dims, 3, simplify = F)
  names(combs) <- sapply(combs, paste, collapse = "")
  lapply(combs, function(cmb) {
    get_density_3d(points[, cmb])
  })
}

get_density_2d_all <- function(points) {
  dims <- colnames(points)
  combs <- combn(dims, 2, simplify = F) # All 2D combinations
  names(combs) <- sapply(combs, paste, collapse = "")

  lapply(combs, function(cmb) {
    get_density_2d(points[, cmb])
  })
}

#' Title
#'
#' @param data
#' @param extend_axes
#' @param density
#'
#' @return
#' @export
#'
#' @examples
add_noise_points <- function(data,
                             extend_axes = 1.2,
                             density = 0.1) {
  if (is.list(data) && identical(names(data), c("coord", "meta"))) {
    points <- data[["coord"]]
  } else {
    points <- data
  }

  library(mclust)
  clustdensl <- lapply(points, function(x) {
    # estimate densities of existing communities
    res <- mclust::Mclust(x[which(x != 0)]) # exclusion of points at exactly 0 wont change much but removes points who miss this dimension
    x <- split(x, res$classification)
    sapply(x, function(y) {
      diff <- abs(min(y) - max(y))
      dens <- length(y)/diff
      return(dens)
    })
  })
  ranges <- lapply(points, range)
  n_noise <- median(mapply(x = clustdensl, y = ranges, function(x,y) {
    y <- y*extend_axes
    diff <- abs(max(y) - min(y))
    median(x)*diff*density # scaling here
  }))
  noise_points <- do.call(cbind, lapply(ranges, function(y) {
    y <- y*extend_axes
    runif(n_noise, min(y), max(y))
  }))


  if (is.list(data) && identical(names(data), c("coord", "meta"))) {
    data[["coord"]] <- dplyr::bind_rows(data[["coord"]], as.data.frame(noise_points))
    data[["meta"]] <- dplyr::bind_rows(data[["meta"]], data.frame(name = rep("noise", nrow(noise_points))))
    return(data)
  } else {
    noise_mat <- data.frame(as.data.frame(noise_points), name = "noise")
    return(noise_mat)
  }

}


#' Title
#'
#' @param object_list
#' @param missing_dim_fill
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
bind_objects <- function(object_list,
                         missing_dim_fill = c("0", "runif"),
                         ...) {

  missing_dim_fill <- match.arg(missing_dim_fill, c("0", "runif"))

  out <- list(coord = do.call(dplyr::bind_rows, sapply(object_list, "[[", "coord", simplify = F)),
              meta = do.call(dplyr::bind_rows, sapply(object_list, "[[", "meta", simplify = F)))

  if (missing_dim_fill == "0") {
    out[["coord"]][is.na(out[["coord"]])] <- 0
  } else {
    for (i in 1:ncol(out[["coord"]])) {
      out[["coord"]][[i]][is.na(out[["coord"]][[i]])] <- runif(n = sum(is.na(out[["coord"]][[i]])))
    }
  }

  return(out)

}

#' Title
#'
#' @param object
#' @param dimnames
#' @param overwrite
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
add_umap <- function(object,
                     dimnames = c("UMAP1", "UMAP2"),
                     overwrite = T,
                     ...) {
  um <- as.data.frame(uwot::umap(as.matrix(object[["coord"]]), scale = T, verbose = T, ...))

  if (overwrite) {
    col_inds <- which(names(object[["meta"]]) %in% dimnames)
    if (length(col_inds) > 0) {
      object[["meta"]] <- object[["meta"]][, -col_inds]
    }
  } else {
    dimnames <- make.unique(c(names(object[["meta"]]), dimnames))
    dimnames <- dimnames[c(length(dimnames)-1, length(dimnames))]
  }
  names(um) <- dimnames

  object[["meta"]] <- dplyr::bind_cols(object[["meta"]], um)
  return(object)
}

#' Title
#'
#' @param object
#' @param x
#' @param y
#' @param color
#' @param ncol_legend
#'
#' @return
#' @export
#'
#' @examples
ggumap <- function(object,
                   x = "UMAP1",
                   y = "UMAP2",
                   color = "name",
                   ncol_legend = 2) {

  # shuffle
  object[["meta"]] <- object[["meta"]][sample(1:nrow(object[["meta"]])),]

  p <- ggplot(object[["meta"]], aes(!!rlang::sym(x), !!rlang::sym(y), color = !!rlang::sym(color))) +
    geom_point(size = 0.2)

  if (is.numeric(object[["meta"]])) {
    p <- p + scale_color_spectral()
  } else {
    p <-
      p +
      theme(legend.position = "bottom", legend.title = element_blank()) +
      scale_color_custom() +
      guides(color = guide_legend(ncol = ncol_legend, override.aes = c(size = 3)))
  }

  return(p)
}

#' Title
#'
#' @param object
#' @param x
#' @param y
#' @param color
#' @param facet_var
#' @param ncol_legend
#' @param legend.position
#' @param theme
#'
#' @return
#' @export
#'
#' @examples
ggbin <- function(object,
                  x = "x1",
                  y = "x2",
                  color = "name",
                  facet_var = "x3_bin",
                  ncol_legend = 2,
                  legend.position = "bottom",
                  theme = theme_void) {

  if (is.list(object) && identical(names(object), c("coord", "meta"))) {
    object <- dplyr::bind_cols(object)
  }

  p <- ggplot(object, aes(!!rlang::sym(x), !!rlang::sym(y), color = !!rlang::sym(color))) +
    geom_point(size = 0.2) +
    do.call(theme, args = list()) +
    coord_fixed() +
    facet_wrap(vars(!!rlang::sym(facet_var)))

  if (!is.numeric(object[[color]])) {
    # discrete
    p <- p +
      ggplot2::theme(legend.position = legend.position, legend.title = element_blank(), plot.background = element_rect(colour = "black")) +
      scale_color_custom() +
      guides(color = guide_legend(ncol = ncol_legend, override.aes = c(size = 3)))
  } else {
    p <- p +
      scale_color_spectral() +
      ggplot2::theme(legend.position = legend.position, plot.background = element_rect(colour = "black"))
  }

  return(p)
}


scale_color_spectral <- function(colors = fcexpr::col_pal("RColorBrewer::Spectral", direction = -1),
                                 values = NULL,
                                 name = waiver(),
                                 na.value = "grey50",
                                 guide = "colorbar",
                                 ...) {
  if (is.null(values)) {
    values <- seq(0, 1, length.out = length(colors))
  }
  ggplot2::scale_color_gradientn(
    colours = colors,
    values = scales::rescale(values),
    name = name,
    na.value = na.value,
    guide = guide,
    ...
  )
}

scale_color_custom <- function(colors = fcexpr::col_pal("custom"), name = waiver(), na.value = "grey50", ...) {
  ggplot2::scale_color_manual(
    values = colors,
    name = name,
    na.value = na.value,
    ...
  )
}


is_rotation_matrix <- function(mat) {

  # for a rotation matrix 3 properties must be fulfilled:

  if (nrow(mat) != ncol(mat)) {
    stop("mat must be nxn.")
  }

  # Check if the matrix is orthogonal: mat multiplied by its transpose must give identity matrix
  #orthogonal_check <- identical(t(mat) %*% mat, diag(nrow(mat)))
  orthogonal_check <- all(dplyr::near(t(mat) %*% mat, diag(nrow(mat))))

  # Check if the determinant is 1 (no squishing of space)
  determinant_check <- dplyr::near(det(mat), 1)

  return(orthogonal_check && determinant_check)
}

rotate_space_3D <- function(coord, rot_mat, long_format = T) {
  if (!is_rotation_matrix(rot_mat)) {
    message("rot_mat is not a rotation matrix.")
  }
  cnames <- colnames(coord)
  ismat <- is.matrix(coord)
  if (long_format) {
    # transpose forth and back
    coord <- t(as.matrix(coord))
    coord <- rot_mat %*% coord
    coord <- t(coord)
    if (!ismat) {
      coord <- as.data.frame(coord)
    }
    colnames(coord) <- cnames
    return(coord)
  } else {
    coord <- rot_mat %*% as.matrix(coord)
    if (!ismat) {
      coord <- as.data.frame(coord)
    }
    colnames(coord) <- cnames
    return(coord)
  }
}

#' Title
#'
#' @param rx
#' @param ry
#' @param rz
#'
#' @return
#' @export
#'
#' @examples
get_rotation_matrix <- function(rx, ry, rz) {
  #rz, ry, rz are rotations around respective axis in radians

  # Rotation around x-axis
  Rx <- matrix(c(1, 0, 0,
                 0, cos(rx), -sin(rx),
                 0, sin(rx), cos(rx)),
               nrow = 3, byrow = TRUE)

  # Rotation around y-axis
  Ry <- matrix(c(cos(ry), 0, sin(ry),
                 0, 1, 0,
                 -sin(ry), 0, cos(ry)),
               nrow = 3, byrow = TRUE)

  # Rotation around z-axis
  Rz <- matrix(c(cos(rz), -sin(rz), 0,
                 sin(rz), cos(rz), 0,
                 0, 0, 1),
               nrow = 3, byrow = TRUE)

  ## matmult of rotation matrices in 3D is non-commutative
  ## so the order of application matters, as the first rotation changes the axis orientations
  R <- Rz %*% Ry %*% Rx
  return(R)
}

get_torus_normal_vec_by_pca <- function(coord, rounding = 1) {
  # rationale: lowest variance in normal direction from torus' central plane
  torus_pca <- prcomp(coord, center = T, scale. = F)
  normal_vector <- torus_pca$rotation[, 3] # third column: last PC --> lowest variance
  unit_normal_vector <- normal_vector / sqrt(sum(normal_vector^2))
  return(abs(round(unit_normal_vector,rounding))) # round since PC is approximation; abs just to make it positive always
}

#' Title
#'
#' @param object
#' @param color
#' @param showlegend
#'
#' @return
#' @export
#'
#' @examples
rgl3dplot <- function(object,
                      x = "x1",
                      y = "x2",
                      z = "x3",
                      color = "name",
                      showlegend = T) {

  type <- "discrete"
  if (is.factor(object[["meta"]][[color]])) {
    inds <- as.numeric(factor(object[["meta"]][[color]]))
    col_vec <- fcexpr::col_pal("custom", n = length(unique(inds)), direction = 1)[inds]
  } else if (is.character(object[["meta"]][[color]])) {
    object[["meta"]][[color]] <- factor(object[["meta"]][[color]], levels = unique(object[["meta"]][[color]])) #sort?
    inds <- as.numeric(object[["meta"]][[color]])
    col_vec <- fcexpr::col_pal("custom", n = length(unique(inds)), direction = 1)[inds]
  } else {
    type <- "continuous"
    # Map values to the color palette
    col_vec <- suppressMessages(fcexpr::col_pal("RColorBrewer::Spectral", n = 100, direction = -1))
    #normvals <- (object[["meta"]][[color]] - min(object[["meta"]][[color]])) / (max(object[["meta"]][[color]]) - min(object[["meta"]][[color]]))
    #col_vec <- col_vec[round(normvals * (length(col_vec) - 1)) + 1]
    normvals <- scales::rescale(object[["meta"]][[color]])
    col_vec <- col_vec[normvals*99+1]
  }

  rgl::open3d()
  rgl::plot3d(x = object[["coord"]][[x]],
              y = object[["coord"]][[y]],
              z = object[["coord"]][[z]],
              xlab = x,
              ylab = y,
              zlab = z,
              col = col_vec,
              aspect = F)
  if (type == "discrete" && showlegend) {
    # add legend
    legend_vec <- stats::setNames(unique(col_vec), unique(object[["meta"]][[color]]))[levels(object[["meta"]][[color]])]
    rgl::legend3d(
      "topright",
      legend = names(legend_vec),
      pch = 16,
      col = legend_vec,
      cex = 1,
      inset = c(0.02)
    )
  }
  rgl::rglwidget()
}
