# spectral_variant_animate.R

utils::globalVariables( c(
  "detector_idx", "intensity", "id", "color_val", "min_val", "max_val", "frame"
) )

#' @title Spectral Variant Animation
#'
#' @description
#' Visualizes the variability of spectral signatures by animating individual
#' variant lines sequentially, followed by a dynamic ribbon representing the
#' full range of variability across all variants. The first few variants appear
#' slowly (one every ten frames) to allow inspection; subsequent variants are
#' added rapidly (one per frame). A shaded ribbon summarising the min/max range
#' across all variants is then revealed at the end.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon scale_x_continuous
#' @importFrom ggplot2 scale_color_identity labs theme_minimal theme element_text
#' @importFrom stats reshape
#' @importFrom gganimate transition_manual animate anim_save gifski_renderer
#'
#' @param spectra.variants Matrix of spectral variants; rows are variants,
#'   columns are detectors. Row names are optional; column names are used as
#'   x-axis labels.
#' @param title Filename prefix for the output GIF. Default is
#'   `"Spectral_Variant_Animation"`.
#' @param save Logical. If `TRUE` (default), saves the GIF to `plot.dir`.
#' @param plot.dir Output directory for the saved GIF. Default is
#'   `"./figure_spectral_variants"`.
#' @param variant.fill.color Fill color for the variability ribbon shown at the
#'   end of the animation. Default is `"red"`.
#' @param variant.fill.alpha Alpha transparency for the variability ribbon.
#'   Default is `0.3`.
#' @param variant.linewidth Line thickness for the individual variant traces.
#'   Default is `1`.
#' @param plot.width Width of the output animation in inches. Default is `16`.
#' @param plot.height Height of the output animation in inches. Default is `9`.
#'
#' @return The rendered `gganimate` animation object, invisibly. When
#'   `save = TRUE`, the GIF is also written to `plot.dir`.
#'
#' @export

spectral.variant.animate <- function(
    spectra.variants,
    title              = "Spectral_Variant_Animation",
    save               = TRUE,
    plot.dir           = "./figure_spectral_variants",
    variant.fill.color = "red",
    variant.fill.alpha = 0.3,
    variant.linewidth  = 1,
    plot.width         = 16,
    plot.height        = 9
) {

  detector.names <- colnames( spectra.variants )
  num_variants   <- nrow( spectra.variants )
  detector.n     <- seq_len( ncol( spectra.variants ) )

  # 1. Calculate Frame Timing
  # IDs 1-4: increments of 10 (Slow) | IDs 5+: increments of 1 (Fast)
  line_starts  <- c( seq( 1, 40, by = 10 ), 41:( 41 + ( num_variants - 5 ) ) )
  ribbon_start <- max( line_starts ) + 5
  end_frame    <- ribbon_start + 20

  # 2. Build explicit frame data for lines
  long_frames <- do.call( rbind, lapply( seq_len( end_frame ), function( f ) {
    active_ids <- which( line_starts <= f )
    if ( length( active_ids ) == 0 ) return( NULL )

    df       <- as.data.frame( spectra.variants[ active_ids, , drop = FALSE ] )
    df$id    <- active_ids
    df$frame <- f

    frame_long <- stats::reshape(
      df,
      direction = "long",
      varying   = detector.names,
      v.names   = "intensity",
      timevar   = "detector_idx",
      times     = detector.n
    )

    # Color Logic: Black for the first 5 frames of a line's life, then variant colour
    frame_long$color_val <- ifelse(
      f < ( line_starts[ frame_long$id ] + 5 ),
      "black",
      variant.fill.color
    )
    frame_long
  } ) )

  # 3. Build explicit frame data for the Ribbon
  ribbon_data <- do.call( rbind, lapply( seq( ribbon_start, end_frame ), function( f ) {
    data.frame(
      detector_idx = detector.n,
      min_val      = apply( spectra.variants, 2, min, na.rm = TRUE ),
      max_val      = apply( spectra.variants, 2, max, na.rm = TRUE ),
      frame        = f
    )
  } ) )

  # 4. Plotting
  anim <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = long_frames,
      ggplot2::aes(
        x     = detector_idx,
        y     = intensity,
        group = id,
        color = color_val
      ),
      alpha     = 0.4,
      linewidth = variant.linewidth
    ) +
    ggplot2::geom_ribbon(
      data = ribbon_data,
      ggplot2::aes( x = detector_idx, ymin = min_val, ymax = max_val ),
      fill  = variant.fill.color,
      alpha = variant.fill.alpha
    ) +
    ggplot2::scale_x_continuous(
      breaks = detector.n,
      labels = detector.names
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::labs( x = "Detector", y = "Intensity", title = title ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text( angle = 45, hjust = 1 )
    ) +
    gganimate::transition_manual( frame )

  # 5. Rendering
  anim_rendered <- gganimate::animate(
    anim,
    nframes   = end_frame,
    fps       = 20,
    end_pause = 40,
    width     = plot.width,
    height    = plot.height,
    units     = "in",
    res       = 150,
    renderer  = gganimate::gifski_renderer()
  )

  if ( save ) {
    if ( !dir.exists( plot.dir ) ) dir.create( plot.dir )
    gganimate::anim_save(
      filename  = paste0( title, ".gif" ),
      animation = anim_rendered,
      path      = plot.dir
    )
  }

  return( invisible( anim_rendered ) )
}
