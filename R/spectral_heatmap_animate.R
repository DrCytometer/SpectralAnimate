# spectral_heatmap_animate.R

utils::globalVariables( c(
  "Detector", "Fluorophore", "value", "reveal.frame"
) )

#' @title Animated Spectral Heatmap
#'
#' @description
#' Animates the spectral matrix heatmap using `gganimate`, revealing one row of
#' spectra (fluorophore) at a time from bottom to top, then pausing on the
#' completed heatmap. Saves the result as a GIF file.
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_classic
#' @importFrom ggplot2 coord_fixed element_text labs theme
#' @importFrom gganimate animate anim_save transition_manual gifski_renderer
#'
#' @param spectra Matrix or dataframe containing spectral data;
#'   format: fluorophores x detectors.
#' @param title Optional prefix for the output animation filename.
#' @param legend.label Character string that will appear on the heatmap legend.
#'   Default is \code{"Intensity"}.
#' @param plot.dir Optional output directory. Default is \code{NULL}, in which
#'   case \code{getwd()} is used.
#' @param color.palette Optional character string defining the viridis color
#'   palette. Default is \code{"viridis"}. Options: \code{"magma"},
#'   \code{"inferno"}, \code{"plasma"}, \code{"viridis"}, \code{"cividis"},
#'   \code{"rocket"}, \code{"mako"}, \code{"turbo"}.
#' @param show.legend Logical. If \code{TRUE}, the fill legend is included.
#'   Default is \code{TRUE}.
#' @param plot.width Width of the output animation in pixels. Default is
#'   \code{NULL}, which auto-scales based on the number of detectors.
#' @param plot.height Height of the output animation in pixels. Default is
#'   \code{NULL}, which auto-scales based on the number of fluorophores.
#' @param fps Frames per second for the output GIF. Default is \code{10}.
#' @param pause.seconds How long (in seconds) to hold on the completed heatmap
#'   at the end of the animation. Default is \code{2}.
#' @param save Logical. If \code{TRUE} (default), saves a GIF to
#'   \code{plot.dir}. If \code{FALSE}, returns the \code{gganimate} object
#'   invisibly for further use.
#'
#' @return If \code{save = TRUE}, writes a GIF and returns the file path
#'   invisibly. If \code{save = FALSE}, returns the animated plot object.
#'
#' @export

spectral.heatmap.animate <- function(
    spectra,
    title         = NULL,
    plot.dir      = NULL,
    legend.label  = "Intensity",
    color.palette = "viridis",
    show.legend   = FALSE,
    plot.width    = NULL,
    plot.height   = NULL,
    fps           = 2,
    pause.seconds = 2,
    save          = TRUE
) {

  # ── dependencies ────────────────────────────────────────────────────────────
  if ( !requireNamespace( "ggplot2",   quietly = TRUE ) ) stop( "Package 'ggplot2' is required."   )
  if ( !requireNamespace( "gganimate", quietly = TRUE ) ) stop( "Package 'gganimate' is required." )
  if ( !requireNamespace( "gifski",    quietly = TRUE ) ) stop( "Package 'gifski' is required (GIF renderer)." )

  # ── file paths ───────────────────────────────────────────────────────────────
  if ( is.null( plot.dir ) )
    plot.dir <- getwd()

  anim.filename <- if ( !is.null( title ) )
    paste0( title, "_spectral_heatmap_animated.gif" )
  else
    "spectral_heatmap_animated.gif"

  # ── prepare data ─────────────────────────────────────────────────────────────
  heatmap.df             <- data.frame( spectra, check.names = FALSE )
  row.levels             <- rownames( heatmap.df )   # fluorophore order (top → bottom)
  col.levels             <- colnames( heatmap.df )   # detector order

  # auto-scale pixel dimensions
  n.rows <- length( row.levels )
  n.cols <- length( col.levels )

  if ( is.null( plot.width ) )
    plot.width  <- as.integer( max( ( ( n.cols - 1 ) / 64 * 12 ) * 96, 400 ) )
  if ( is.null( plot.height ) )
    plot.height <- as.integer( ( 5 + round( n.rows / 8, 0 ) ) * 96 )

  # pivot to long format (identical to spectral.heatmap)
  heatmap.df$Fluorophore <- row.levels
  heatmap.long <- data.frame(
    Fluorophore = rep( row.levels, times = n.cols ),
    Detector    = rep( col.levels, each  = n.rows ),
    value       = as.vector(
      as.matrix( heatmap.df[ , colnames( heatmap.df ) != "Fluorophore" ] )
    ),
    stringsAsFactors = FALSE
  )

  # factor levels: Fluorophore reversed so row 1 appears at the top of the plot
  heatmap.long$Fluorophore <- factor( heatmap.long$Fluorophore, levels = rev( row.levels ) )
  heatmap.long$Detector    <- factor( heatmap.long$Detector,    levels = col.levels )

  # ── animation state column ───────────────────────────────────────────────────
  # Each fluorophore (row) gets its own reveal frame.
  # row.levels[1] is the top row visually; rev() makes it appear last so the
  # heatmap builds upward (bottom-to-top reveal).

  heatmap.long$reveal.frame <- match(
    as.character( heatmap.long$Fluorophore ),
    row.levels
  )

  # total animation frames = one per row + pause frames at end
  n.rows       <- length( row.levels )
  pause.frames <- round( pause.seconds * fps )
  total.frames <- n.rows + pause.frames

  # ── build plot ───────────────────────────────────────────────────────────────
  heatmap.plot <- ggplot2::ggplot(
    heatmap.long,
    ggplot2::aes( Detector, Fluorophore, fill = value )
  ) +
    ggplot2::geom_tile() +
    ggplot2::theme_classic() +
    ggplot2::coord_fixed( ratio = 1 ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 )
    ) +
    ggplot2::labs( x = NULL, y = NULL, fill = legend.label ) +
    ggplot2::scale_fill_viridis_c( option = color.palette )

  if ( !show.legend )
    heatmap.plot <- heatmap.plot + ggplot2::theme( legend.position = "none" )

  # ── gganimate layer ──────────────────────────────────────────────────────────
  # transition_manual steps through reveal.frame one at a time.
  # shadow_mark(past = TRUE) keeps previously revealed rows visible.
  heatmap.anim <- heatmap.plot +
    gganimate::transition_manual( frames = reveal.frame, cumulative = TRUE )

  # ── render ───────────────────────────────────────────────────────────────────
  rendered <- gganimate::animate(
    heatmap.anim,
    nframes   = total.frames,
    fps       = fps,
    width     = plot.width,
    height    = plot.height,
    renderer  = gganimate::gifski_renderer(),
    # end_pause holds the final frame for `pause.frames` ticks
    end_pause = pause.frames
  )

  if ( save ) {
    out.path <- file.path( plot.dir, anim.filename )
    gganimate::anim_save( out.path, animation = rendered )
    message( "Animation saved to: ", out.path )
    return( invisible( out.path ) )
  } else {
    return( rendered )
  }
}
