g_acc <- 9.81 # gravitational acceleration (m/s^2)

# Wave energy constant (using density of water at 20 C i.e. 998 kg/m^3)
we_const <- 0.998 * g_acc^2 / (64 * pi)


#' Calculate the wave energy flux
#'
#' Calculates the wave energy flux (power per meter of wave crest) given
#' either (1) the significant wave height and peak period or (2) the wind speed
#' at 10m, fetch length and (optionally) water depth.
#'
#' Given the significant height (\emph{H}) and peak period (\emph{T}),
#' the wave energy flux is calculated as: \deqn{\frac{\rho g^2}{64 \pi} H^2 T},
#' where \eqn{\rho} is the density of water (998 kg/m^3) and \emph{g} is the
#' acceleration of gravity (9.81 m/s^2).
#'
#' If both \code{height} and \code{period} are missing, they are estimated from
#' on the wind speed at 10m (\eqn{U_{10}}) and the fetch length (\emph{F}) as
#' described in Resio et al. (2003):
#' \deqn{{U_f}^2 = 0.001 (1.1 + 0.035 U_{10}) {U_{10}}^2} (friction velocity)
#' \deqn{\frac{g H}{{U_f}^2} = \min (0.0413 \sqrt{\frac{g F}{{U_f}^2}}, 211.5)}
#' \deqn{\frac{g T}{U_f} = \min (0.651 (\frac{g F}{{U_f}^2})^{1/3}, 239.8)}
#' If the depth (\emph{d}) is specified, it imposes a limit on the peak period:
#' \deqn{T_{max} = 9.78 \sqrt{\frac{d}{g}}} (in seconds)
#'
#' @param height Significant wave height, in meters.
#' @param period Peak wave period, in seconds.
#' @param wind Wind speed at 10m, in m/s.
#' @param fetch Fetch length, in meters.
#' @param depth Water depth, in meters.
#' @return Wave energy flux, in kW/m.
#' @references Resio, D.T., Bratos, S.M., and Thompson, E.F. (2003). Meteorology
#'  and Wave Climate, Chapter II-2. Coastal Engineering Manual.
#'  US Army Corps of Engineers, Washington DC, 72pp.
#' @examples
#'  # With height and period arguments
#'  wave_energy(8, 1)
#'
#'  # With wind, fetch and depth arguments (must be named)
#'  wave_energy(wind = 12, fetch = 15000, depth = 10)
#' @export
wave_energy <- function(height = NA, period = NA,
                        wind = NA, fetch = NA, depth = NA) {
    if (all(is.na(height)) && all(is.na(period))) {
        # Friction velocity squared
        uf2 <- 0.001 * (1.1 + 0.035 * wind) * wind^2
        # Non-dimensional fetch, height and period
        fetch_nd <- fetch * g_acc / uf2
        height_nd <- pmin(0.0413 * sqrt(fetch_nd), 211.5)
        period_nd <- pmin(0.651 * fetch_nd^(1/3), 239.8)
        # Calculate height and period (latter with depth limitation, if available)
        #  then wave energy
        height <- height_nd * uf2 / g_acc
        period <- period_nd * sqrt(uf2) / g_acc
        if (any(!is.na(depth)))
            period <- pmin(period, 9.78 * sqrt(depth / g_acc))
    } else {
        if (any(!is.na(wind)) || any(!is.na(fetch)) || any(!is.na(depth)))
            warning("height or period provided; ignoring other arguments to wave_energy.")
    }
    we_const * height^2 * period
}
