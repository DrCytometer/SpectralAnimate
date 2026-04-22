# unmixing_matrix_animate.R

utils::globalVariables( c(
  "Detector", "Fluorophore", "value", "state"
) )

#' @title Unmixing Matrix Animation
#'
#' @description
#' Animates the transition from the spectral spillover (mixing) matrix to the
#' unmixing matrix (Moore-Penrose pseudoinverse via OLS). Both matrices are
#' displayed in a consistent fluorophore x detector orientation — the unmixing
#' matrix is transposed before display so axes remain stable across the
#' transition.
#'
#' Each matrix uses its own independent colour scale so that the dynamic range
#' of peaks and troughs is fully visible in both the mixing and unmixing
#' heatmaps. This is achieved by rescaling each matrix to \[0, 1\] before
#' combining into a single animated `ggplot`, then using a shared
#' `scale_fill_viridis_c` with `limits = c(0, 1)`.
#'
#' The mixing matrix (normalised spectral signatures) has fluorophores as rows
#' and detectors as columns. The pseudoinverse returned by `unmix.ols()` has the
#' opposite orientation (detectors x fluorophores); it is transposed internally
#' so that both heatmaps share the same axis layout.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_classic
#' @importFrom ggplot2 coord_fixed element_text labs theme
#' @importFrom gganimate transition_states animate gifski_renderer
#'
#' @param spectra Normalised spectral signatures (mixing matrix). Fluorophores
#'   in rows, detectors in columns. Both row and column names are required.
#' @param color.palette Viridis palette name for the heatmap fill. Default is
#'   `"viridis"`. Options: `"magma"`, `"inferno"`, `"plasma"`, `"viridis"`,
#'   `"cividis"`, `"rocket"`, `"mako"`, `"turbo"`.
#' @param title Optional character string prepended to the animation filename.
#'   Default is `NULL`.
#' @param figure.dir Output directory for the saved GIF. Default is `"."`.
#' @param animation.filename Base filename for the GIF output. Default is
#'   `"unmixing_matrix_transition.gif"`.
#' @param show.legend Logical. Whether to show the fill legend. Default is
#'   `FALSE`.
#' @param plot.width Width of each frame in inches. Default is `NULL`
#'   (auto-scaled by number of detectors).
#' @param plot.height Height of each frame in inches. Default is `NULL`
#'   (auto-scaled by number of fluorophores).
#' @param nframes Total number of animation frames. Default is `80`.
#' @param fps Frames per second for the output GIF. Default is `15`.
#' @param transition.length Relative length of the morphing transition between
#'   states (passed to `transition_states()`). Default is `3`.
#' @param state.length Relative length of the pause on each state (passed to
#'   `transition_states()`). Default is `2`.
#' @param res Resolution in ppi for the rendered frames. Default is `150`.
#' @param save Logical. If `TRUE`, saves the GIF to `figure.dir`. If `FALSE`,
#'   returns the rendered animation object invisibly. Default is `TRUE`.
#'
#' @return The rendered `gganimate` animation object, invisibly.
#'
#' @export

unmixing.matrix.animate <- function(
    spectra,
    color.palette      = "viridis",
    title              = NULL,
    figure.dir         = ".",
    animation.filename = "unmixing_matrix_transition.gif",
    show.legend        = FALSE,
    plot.width         = NULL,
    plot.height        = NULL,
    nframes            = 80,
    fps                = 15,
    transition.length  = 3,
    state.length       = 2,
    res                = 150,
    save               = TRUE
) {

  # ── Input checks ──────────────────────────────────────────────────────────
  if ( is.null( rownames( spectra ) ) || is.null( colnames( spectra ) ) )
    stop( "`spectra` must have both row names (fluorophores) and column names (detectors)." )

  # ── Compute unmixing matrix (pseudoinverse via SVD) ───────────────────────
  sv              <- svd( t( spectra ) )
  unmixing.matrix <- sv$v %*% ( t( sv$u ) / sv$d )
  rownames( unmixing.matrix ) <- rownames( spectra )
  colnames( unmixing.matrix ) <- colnames( spectra )

  # ── Auto-scale dimensions ─────────────────────────────────────────────────
  n.fluor    <- nrow( spectra )
  n.detector <- ncol( spectra )

  if ( is.null( plot.width ) )
    plot.width  <- max( ( n.detector / 64 ) * 12, 4 )
  if ( is.null( plot.height ) )
    plot.height <- max( 5 + round( n.fluor / 8, 0 ), 3 )

  # ── Helper: matrix → long data.frame with per-matrix [0,1] rescaling ─────
  rescale01 <- function( x ) {
    rng <- range( x, na.rm = TRUE )
    if ( diff( rng ) == 0 ) return( rep( 0, length( x ) ) )
    ( x - rng[1] ) / diff( rng )
  }

  mat_to_long <- function( mat, state.label ) {
    raw.vals <- as.vector( mat )
    data.frame(
      Fluorophore = factor(
        rep( rownames( mat ), times = n.detector ),
        levels = rev( rownames( mat ) )
      ),
      Detector    = factor(
        rep( colnames( mat ), each = n.fluor ),
        levels = colnames( mat )
      ),
      value       = rescale01( raw.vals ),
      state       = state.label,
      stringsAsFactors = FALSE
    )
  }

  mixing.long   <- mat_to_long( spectra,          "Spectral signatures" )
  unmixing.long <- mat_to_long( unmixing.matrix,  "Unmixing matrix"     )

  anim.data       <- rbind( mixing.long, unmixing.long )
  anim.data$state <- factor( anim.data$state, levels = unique( anim.data$state ) )

  # ── Build animated plot ───────────────────────────────────────────────────
  heatmap.anim <- ggplot2::ggplot(
    anim.data,
    ggplot2::aes( Detector, Fluorophore, fill = value )
  ) +
    ggplot2::geom_tile( colour = "white", linewidth = 0.15 ) +
    ggplot2::scale_fill_viridis_c(
      option = color.palette,
      name   = "Scaled\nvalue",
      limits = c( 0, 1 )
    ) +
    ggplot2::coord_fixed( ratio = 1 ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 ),
      plot.title  = ggplot2::element_text( size = 14, face = "bold", hjust = 0.5 )
    ) +
    ggplot2::labs(
      x     = "Detector",
      y     = NULL,
      title = "{closest_state}"
    ) +
    gganimate::transition_states(
      state,
      transition_length = transition.length,
      state_length      = state.length,
      wrap              = FALSE
    )

  if ( !show.legend )
    heatmap.anim <- heatmap.anim +
    ggplot2::theme( legend.position = "none" )

  # ── Render ─────────────────────────────────────────────────────────────────
  if ( !dir.exists( figure.dir ) ) dir.create( figure.dir, recursive = TRUE )

  anim.file <- if ( !is.null( title ) )
    paste0( title, "_", animation.filename )
  else
    animation.filename

  rendered <- gganimate::animate(
    heatmap.anim,
    nframes  = nframes,
    fps      = fps,
    width    = plot.width,
    height   = plot.height,
    units    = "in",
    res      = res,
    device   = "ragg_png",
    renderer = gganimate::gifski_renderer( file.path( figure.dir, anim.file ) )
  )

  if ( save )
    message( "Animation saved to: ", file.path( figure.dir, anim.file ) )

  return( invisible( rendered ) )
}
