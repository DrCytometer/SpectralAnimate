# spectral_normalization_animate.R

utils::globalVariables( c(
  "channel", "y_val", "event_id", "point_alpha", "state_label"
) )

#' @title Spectral Normalization Animation
#'
#' @description
#' Animates the transition of raw spectral data to a normalized signature,
#' helping visualize how data is rescaled based on median signal intensities.
#' The animation transitions between two states: raw biexponentially-transformed
#' event intensities, and the per-channel normalized spectral signature derived
#' from channel medians.
#'
#' @importFrom utils stack
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom ggplot2 ggplot aes scale_y_continuous geom_point geom_line
#' @importFrom ggplot2 theme_minimal labs theme element_text element_blank scale_alpha_identity
#' @importFrom gganimate transition_states view_follow animate anim_save gifski_renderer
#' @importFrom scales trans_new
#' @importFrom stats median
#'
#' @param ribbon.plot.data Dataframe or matrix containing the raw spectral event
#'   intensities. Rows are events, columns are spectral detectors.
#' @param asp The AutoSpectral parameter list. Used to extract biexponential
#'   transformation parameters, axis breaks, limits, and plot styling.
#' @param figure.dir Output directory for the saved GIF. Default is
#'   `"./figure_spectral_ribbon"`.
#' @param title Animation filename prefix. Default is
#'   `"Normalization_Transition"`.
#' @param downsample.n Number of events to sample prior to rendering for
#'   performance. Default is `500`.
#' @param plot.width Width of the output animation in inches. Default is `16`.
#' @param plot.height Height of the output animation in inches. Default is `9`.
#'
#' @return The rendered `gganimate` animation object, invisibly. The GIF is also
#'   written to `figure.dir`.
#'
#' @export

spectral.normalization.animate <- function(
    ribbon.plot.data,
    asp,
    figure.dir   = "./figure_spectral_ribbon",
    title        = "Normalization_Transition",
    downsample.n = 500,
    plot.width   = 16,
    plot.height  = 9
) {

  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )

  biexp.inverse <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = TRUE
  )

  biexp_scale_trans <- scales::trans_new(
    name      = "biexp",
    transform = biexp.transform,
    inverse   = biexp.inverse
  )

  if ( !dir.exists( figure.dir ) ) dir.create( figure.dir )

  # downsample
  if ( nrow( ribbon.plot.data ) > downsample.n ) {
    set.seed( 42 )
    row.idx          <- sample( nrow( ribbon.plot.data ), downsample.n )
    ribbon.plot.data <- ribbon.plot.data[ row.idx, , drop = FALSE ]
  }

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

  ribbon.plot.long$event_id <- rep(
    seq_len( nrow( ribbon.plot.data ) ),
    times = length( channels )
  )

  state1 <- ribbon.plot.long
  state1$y_val       <- biexp.transform( state1$value )
  state1$state_label <- "Raw Spectral Data"
  state1$point_alpha <- 0.4

  # channel medians for normalization
  medians       <- tapply( ribbon.plot.long$value, ribbon.plot.long$channel, median, na.rm = TRUE )
  median.vals   <- as.numeric( medians )
  norm.range    <- range( median.vals, na.rm = TRUE )
  norm.m        <- ( median.vals - norm.range[1] ) / diff( norm.range )
  names( norm.m ) <- names( medians )

  state2 <- ribbon.plot.long
  state2$y_val       <- norm.m[ as.character( state2$channel ) ] *
    max( biexp.transform( ribbon.limits ) )
  state2$state_label <- "Normalized Spectral Signature"
  state2$point_alpha <- 1.0

  animation_data             <- rbind( state1, state2 )
  animation_data$state_label <- factor(
    animation_data$state_label,
    levels = c( "Raw Spectral Data", "Normalized Spectral Signature" )
  )

  anim.plot <- ggplot2::ggplot(
    animation_data,
    ggplot2::aes( x = channel, y = y_val, group = event_id )
  ) +
    ggplot2::scale_y_continuous(
      trans  = biexp_scale_trans,
      limits = ribbon.limits,
      breaks = ribbon.breaks,
      labels = ribbon.labels
    ) +
    ggplot2::geom_point(
      ggplot2::aes( alpha = point_alpha ),
      color = "royalblue",
      size  = 0.8
    ) +
    ggplot2::geom_line(
      data = state2,
      ggplot2::aes( group = 1 ),
      color     = "black",
      linewidth = 1
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme_minimal() +
    ggplot2::labs( title = "{closest_state}", x = "Detector", y = "Intensity" ) +
    ggplot2::theme(
      axis.text.x       = ggplot2::element_text(
        angle = asp$ribbon.plot.axis.text.angle,
        vjust = 1,
        hjust = 1
      ),
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "none",
      strip.text        = ggplot2::element_text(
        size = asp$ribbon.plot.strip.text.size,
        face = asp$ribbon.plot.strip.text.face
      )
    ) +
    gganimate::view_follow( fixed_x = TRUE ) +
    gganimate::transition_states(
      state_label,
      transition_length = 3,
      state_length      = 2,
      wrap              = FALSE
    )

  anim_rendered <- gganimate::animate(
    anim.plot,
    nframes   = 100,
    device    = "ragg_png",
    renderer  = gganimate::gifski_renderer(),
    end_pause = 30,
    width     = plot.width,
    height    = plot.height,
    units     = "in",
    res       = 150
  )

  animation.file.name <- paste0( title, "_normalization.gif" )

  gganimate::anim_save(
    filename  = animation.file.name,
    animation = anim_rendered,
    path      = figure.dir
  )

  return( invisible( anim_rendered ) )
}
