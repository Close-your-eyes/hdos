fun_nd <- function(n_samples = 1000,
                   n_dim = 3,
                   dimnames = paste0("x", 1:n_dim),
                   fun = function(x1, x2) sin(x1*x2),
                   range1 = c(-5,5),
                   range2 = c(-5,5),
                   R2 = 0.9) {


  x <- runif(n_samples, range1[1], range1[2])  # Random values for x
  y <- runif(n_samples, range2[1], range2[2])  # Random values for y

  # compute the true output
  z_true <- fun(x, y)

  # add noise with controlled variance
  total_variance <- var(z_true)
  explained_variance <- total_variance * R2
  noise_variance <- total_variance - explained_variance
  noise_sd <- sqrt(noise_variance)

  # Generate noisy z values
  z <- z_true + rnorm(length(x), mean = 0, sd = noise_sd)

  points <- data.frame(x,y,z)
  names(points) <- dimnames

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
    names(angle_3d) <- sapply(strsplit(names(angle_3d), "_"), "[", 1)
  }


  points <- list(coord = as.data.frame(points),
                 meta = do.call(dplyr::bind_cols, list(as.data.frame(bins),
                                                       radii_3d,
                                                       angle_3d
                                                       #data.frame(hash = hashname, name = name)
                                                       )))
  return(points)
}


fun_nd2 <- function() {

  # Step 1: Define the helical path
  t <- seq(0, 4 * pi, length.out = 200)  # Parameter for the central curve
  x_center <- sin(t)  # Helical x-component
  y_center <- cos(t)  # Helical y-component
  z_center <- sin(t)       # Helical z-component (vertical rise)

  # Step 2: Define the tube's radius and circular cross-sections
  radius <- 2  # Radius of the tube
  num_points <- 50  # Points per circle (resolution of tube)

  # Generate points around the circular cross-sections
  theta <- seq(0, 2 * pi, length.out = num_points)  # Circle parameter
  circle_x <- radius * cos(theta)  # Circle x-coordinates
  circle_y <- radius * sin(theta)  # Circle y-coordinates

  # Step 3: Generate the tube's surface
  x <- c()
  y <- c()
  z <- c()

  for (i in seq_along(t)) {
    # Rotate the circle to align with the helix direction
    x <- c(x, x_center[i] + circle_x)
    y <- c(y, y_center[i] + circle_y)
    z <- c(z, rep(z_center[i], num_points))
  }

  # Step 4: Visualize the tube
  plot_ly(x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers",
          marker = list(size = 2, color = z, colorscale = "Viridis")) %>%
    layout(scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    ))

}




