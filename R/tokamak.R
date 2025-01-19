#' Title
#'
#' @param n_samples
#' @param R_plasma
#' @param r_plasma
#' @param R_magnet
#' @param r_magnet
#' @param RR_magnet
#' @param magnet_freq
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
tokamak <- function(n_samples = 1000,
                    R_plasma = 3,
                    r_plasma = 0.05,
                    R_magnet = 2,
                    r_magnet = 0.1,
                    RR_magnet = 3,
                    magnet_freq = pi/8,
                    ...) {
  tm <- lapply(seq(0, 2*pi, magnet_freq), function(x) {
    torus_3d(
      n_samples,
      R_magnet,
      r_magnet,
      origin = c(RR_magnet*cos(x),RR_magnet*sin(x),0),
      rot_mat = get_rotation_matrix(pi/2,0,x),
      ...
    )
  }, ...)
  tp <- torus_3d(n_samples, R_plasma, r_plasma, ...)
  return(bind_objects(c(list(tp), tm)))
}
