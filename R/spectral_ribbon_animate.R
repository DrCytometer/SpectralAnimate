# spectral_ribbon_animate.R

utils::globalVariables( c(
  "channel", "value", "channel.num"
) )

#' @title Spectral Ribbon Animation
#'
#' @description
#' Generates an animated spectral ribbon plot that builds the density heatmap
#' detector by detector using `gganimate`. Each frame adds the next detector
#' column cumulatively, giving a left-to-right reveal of the full ribbon. An
#' optional static JPEG is also saved alongside the GIF.
#'
#' The function is called internally by `clean.controls()`. To use it directly,
#' supply data via `ribbon.plot.data` (a matrix or data.frame, events x
#' detectors) and a valid `asp` parameter list.
#'
#' @importFrom ggplot2 ggplot aes scale_y_continuous geom_bin2d labs
#' @importFrom ggplot2 theme_minimal theme element_text element_blank
#' @importFrom ggplot2 scale_fill_gradientn scale_fill_viridis_c ggsave
#' @importFrom utils stack
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom ragg agg_jpeg
#' @importFrom gganimate transition_manual animate anim_save gifski_renderer
#'
#' @param ribbon.plot.data A matrix or dataframe containing the raw spectral
#'   data. Rows are events, columns are detectors.
#' @param asp The AutoSpectral parameter list. Used for biexponential
#'   transformation parameters, axis breaks, limits, and plot styling.
#' @param title An optional character string prepended to both the plot and
#'   animation filenames. Default is `NULL`.
#' @param figure.dir Output directory for the saved files. If `NULL` (default),
#'   `asp$figure.spectral.ribbon.dir` is used.
#' @param save Logical. If `TRUE` (default), the plot and animation are saved to
#'   `figure.dir`. If `FALSE`, the rendered animation is returned invisibly
#'   without saving.
#' @param color.palette Character string defining the color palette. Default is
#'   `"rainbow"`, mimicking a FlowJo colour scheme. Viridis options are also
#'   accepted: `"magma"`, `"inferno"`, `"plasma"`, `"viridis"`, `"cividis"`,
#'   `"rocket"`, `"mako"`, `"turbo"`.
#' @param plot.width Width of the saved plot in inches. Default is `15`.
#' @param plot.height Height of the saved plot in inches. Default is `10`.
#' @param plot.filename Base filename for the static JPEG. Default is
#'   `"spectral_ribbon_plot.jpg"`.
#' @param animation.filename Base filename for the GIF. Default is
#'   `"spectral_animation.gif"`.
#'
#' @return The rendered `gganimate` animation object, invisibly. When
#'   `save = TRUE`, the GIF and a static JPEG are also written to `figure.dir`.
#'
#' @export

spectral.ribbon.animate <- function(
    ribbon.plot.data,
    asp,
    title              = NULL,
    figure.dir         = NULL,
    save               = TRUE,
    color.palette      = "rainbow",
    plot.width         = 15,
    plot.height        = 10,
    plot.filename      = "spectral_ribbon_plot.jpg",
    animation.filename = "spectral_animation.gif"
) {

  if ( is.null( figure.dir ) ) figure.dir <- asp$figure.spectral.ribbon.dir
  if ( !dir.exists( figure.dir ) ) dir.create( figure.dir )

  # shift to long format for plotting
  channels         <- colnames( ribbon.plot.data )
  ribbon.plot.data <- as.data.frame( ribbon.plot.data[ , channels ] )
  ribbon.plot.long <- utils::stack( ribbon.plot.data )
  names( ribbon.plot.long ) <- c( "value", "channel" )
  ribbon.plot.long$channel  <- factor( ribbon.plot.long$channel, levels = channels )

  # setting scales and transformation
  ribbon.breaks <- asp$ribbon.breaks
  ribbon.labels <- sapply( ribbon.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  ribbon.limits <- c( asp$ribbon.plot.min, asp$expr.data.max )

  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )

  # animation transition column
  ribbon.plot.long$channel.num <- as.numeric( ribbon.plot.long$channel )

  # create plot
  ribbon.plot <- ggplot2::ggplot(
    ribbon.plot.long,
    ggplot2::aes( channel, biexp.transform( value ) )
  ) +
    ggplot2::scale_y_continuous(
      limits = biexp.transform( ribbon.limits ),
      breaks = biexp.transform( ribbon.breaks ),
      labels = ribbon.labels
    ) +
    ggplot2::geom_bin2d(
      bins     = c( length( channels ), asp$ribbon.bins ),
      boundary = 0.5,
      na.rm    = TRUE
    ) +
    ggplot2::labs( title = title, x = "Detector", y = "Intensity" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(
        angle = asp$ribbon.plot.axis.text.angle,
        vjust = 1,
        hjust = 1
      ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none",
      strip.text       = ggplot2::element_text(
        size = asp$ribbon.plot.strip.text.size,
        face = asp$ribbon.plot.strip.text.face
      )
    )

  # color options
  viridis.colors <- c(
    "magma", "inferno", "plasma", "viridis",
    "cividis", "rocket", "mako", "turbo"
  )
  if ( color.palette %in% viridis.colors ) {
    ribbon.plot <- ribbon.plot +
      ggplot2::scale_fill_viridis_c( option = color.palette )
  } else {
    ribbon.plot <- ribbon.plot +
      ggplot2::scale_fill_gradientn(
        colours = asp$density.palette.base.color,
        values  = asp$ribbon.scale.values
      )
  }

  animation.file.name <- if ( !is.null( title ) )
    paste( title, animation.filename, sep = "_" )
  else
    animation.filename

  anim <- ribbon.plot +
    gganimate::transition_manual( channel.num, cumulative = TRUE )

  rendered.anim <- gganimate::animate(
    anim,
    nframes  = length( channels ),
    width    = plot.width  * 72,
    height   = plot.height * 72,
    renderer = gganimate::gifski_renderer(
      file.path( figure.dir, animation.file.name )
    ),
    device   = "ragg_png"
  )

  # save or return
  if ( save ) {
    ribbon.plot.filename <- if ( !is.null( title ) )
      paste( title, plot.filename, sep = "_" )
    else
      plot.filename

    ggplot2::ggsave(
      ribbon.plot.filename,
      plot      = ribbon.plot,
      device    = ragg::agg_jpeg,
      path      = figure.dir,
      width     = plot.width,
      height    = plot.height,
      limitsize = FALSE,
      create.dir = TRUE
    )

    gganimate::anim_save(
      filename  = animation.file.name,
      animation = rendered.anim,
      path      = figure.dir
    )
  }

  return( invisible( rendered.anim ) )
}
