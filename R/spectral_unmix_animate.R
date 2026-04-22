# spectral_unmix_animate.R

utils::globalVariables( c(
  "channel", "value", "fluor", "Fluorophore", "Detector",
  "value_raw", "event_id", "channel_raw", "y", "fluor_idx",
  "value_resid", "value_unmixed", "state", "point_alpha"
) )

#' @title Spectral Unmixing Animation
#'
#' @description
#' Produces a three-panel animation illustrating spectral unmixing of flow
#' cytometry data. The layout is:
#'
#'   Raw spectral scatter | Unmixing matrix heatmap | Unmixed scatter
#'
#' For each single-colour control (SCC) provided in `data.list`, the animation
#' shows all fluorophore layers as coloured scatter points on the raw (left)
#' panel, transitions them to their residual positions (`residuals = raw -
#' unmixed %*% spectra`), and simultaneously reveals the unmixed data on the
#' right panel. The central panel shows a static heatmap of the unmixing matrix.
#' Panels are stitched per-frame using `magick`.
#'
#' @section Animation stages per layer:
#' Each SCC cycles through five named states that `gganimate::transition_states()`
#' interpolates between:
#' \describe{
#'   \item{`"raw"`}{All layers visible as a density ribbon; active layer
#'     highlighted as coloured points on top.}
#'   \item{`"peel"`}{Active layer lifts off; remaining points shift to residuals.}
#'   \item{`"transit"`}{Active-layer cloud moves across the centre panel.}
#'   \item{`"land"`}{Cloud arrives at the unmixed channel on the right panel.}
#'   \item{`"settle"`}{All previously unmixed layers visible; new layer settled.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_tile scale_y_continuous
#' @importFrom ggplot2 scale_x_continuous scale_x_discrete scale_fill_gradientn
#' @importFrom ggplot2 scale_fill_viridis_c scale_colour_manual scale_alpha_identity
#' @importFrom ggplot2 theme_minimal theme_classic theme element_text element_blank
#' @importFrom ggplot2 labs coord_fixed ggsave
#' @importFrom gganimate transition_states animate anim_save gifski_renderer
#' @importFrom gganimate shadow_mark ease_aes
#' @importFrom cowplot plot_grid
#' @importFrom flowWorkspace flowjo_biexp
#' @importFrom dplyr bind_rows mutate select all_of
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom scales hue_pal
#' @importFrom magick image_read image_append image_animate image_join image_write
#' @importFrom magick image_scale
#'
#' @param data.list A **named** list of matrices or data.frames. Each element
#'   is a single-colour control: rows are events, columns are spectral detectors.
#'   Names are used as fluorophore labels. Column names must match
#'   `spectral.channel`.
#' @param spectra Normalised spectral signatures (mixing matrix), fluorophores x
#'   detectors. Row names must match `names(data.list)`, column names must match
#'   `spectral.channel`.
#' @param spectral.channel Character vector of spectral detector names (the raw
#'   x-axis). Must match column names in each element of `data.list` and in
#'   `spectra`.
#' @param asp The AutoSpectral parameter list. Used for biexponential transform,
#'   axis breaks, limits, and palette parameters.
#' @param fluor.colors Optional named character vector mapping fluorophore names
#'   to colours. If `NULL` (default), colours are drawn automatically from
#'   `scales::hue_pal()`.
#' @param max.points Maximum number of events to sample **per element of
#'   `data.list`** for the scatter plot panels. Default is `500`. Sampling is
#'   applied independently to each fluorophore's data to keep rendering fast.
#' @param color.palette Viridis palette name for the central heatmap fill.
#'   Default is `"viridis"`.
#' @param ribbon.color.palette Palette for the raw density ribbon fill. Default
#'   is `"rainbow"`, matching `spectral.ribbon.plot()`. Set to a viridis name to
#'   override.
#' @param title Optional character string prepended to the animation filename.
#' @param figure.dir Output directory for the saved GIF. Default is `"."`.
#' @param animation.filename Base filename for the GIF. Default is
#'   `"spectral_unmix.gif"`.
#' @param plot.width Width of the combined figure in inches. Default is `18`.
#' @param plot.height Height of the combined figure in inches. Default is `6`.
#' @param nframes Total animation frames. Default is `120`.
#' @param fps Frames per second. Default is `15`.
#' @param res Resolution in ppi. Default is `150`.
#' @param transition.length Relative transition length per state pair. Default
#'   is `2`.
#' @param state.length Relative pause length per state. Default is `1`.
#' @param save Logical. Save GIF if `TRUE`. Default is `TRUE`.
#'
#' @return Invisibly returns a list with elements `raw.plot`, `heatmap.plot`,
#'   `right.plot`, and `animation` (the rendered `magick` GIF object).
#'
#' @export
spectral.unmix.animate <- function(
    data.list,
    spectra,
    spectral.channel,
    asp,
    fluor.colors        = NULL,
    max.points          = 500,
    color.palette       = "viridis",
    ribbon.color.palette = "rainbow",
    title               = NULL,
    figure.dir          = ".",
    animation.filename  = "spectral_unmix.gif",
    plot.width          = 18,
    plot.height         = 6,
    nframes             = 120,
    fps                 = 100,
    res                 = 150,
    transition.length   = 3,
    state.length        = 2,
    save                = TRUE
) {

  # -- Input checks ------------------------------------------------------------
  if ( is.null( names( data.list ) ) )
    stop( "`data.list` must be a *named* list; names are used as fluorophore labels." )
  if ( !all( names( data.list ) %in% rownames( spectra ) ) )
    stop( "All names in `data.list` must appear as row names in `spectra`." )

  fluor.names <- names( data.list )
  n.fluor     <- length( fluor.names )
  n.det <- length( spectral.channel )

  # -- Colours ------------------------------------------------------------------
  if ( is.null( fluor.colors ) ) {
    palette.fn  <- scales::hue_pal()
    auto.colors <- palette.fn( n.fluor )
    fluor.colors <- stats::setNames( auto.colors, fluor.names )
  }

  # -- Biexponential transform -------------------------------------------------
  biexp.transform <- flowWorkspace::flowjo_biexp(
    channelRange = asp$default.transformation.param$length,
    maxValue     = asp$default.transformation.param$max.range,
    pos          = asp$default.transformation.param$pos,
    neg          = asp$default.transformation.param$neg,
    widthBasis   = asp$default.transformation.param$width,
    inverse      = FALSE
  )

  ribbon.breaks <- asp$ribbon.breaks
  ribbon.labels <- sapply( ribbon.breaks, function( x ) {
    if ( x == 0 ) "0" else parse( text = paste0( "10^", log10( abs( x ) ) ) )
  } )
  ribbon.limits <- c( asp$ribbon.plot.min, asp$expr.data.max )

  # -- Compute unmixing matrix (pseudoinverse) ----------------------------------
  sv               <- svd( t( spectra ) )
  unmixing.matrix  <- sv$v %*% ( t( sv$u ) / sv$d )

  # For the static central heatmap we keep the unmixing matrix as
  # fluorophores ?? detectors (transposed relative to its mathematical form)
  # so the heatmap axes match the mixing matrix layout visually.
  rownames( unmixing.matrix ) <- rownames( spectra )
  colnames( unmixing.matrix ) <- colnames( spectra )

  fluor.channel <- rownames( spectra )   # unmixed axis labels

  # -- Build sampled raw data for the scatter left panel -----------------------
  # Downsample each element of data.list to max.points for scatter display.
  # Each layer gets its own color matching the unmixed-data panel.
  all.raw.sampled <- dplyr::bind_rows( lapply( fluor.names, function( nm ) {
    df <- as.data.frame( data.list[[ nm ]] )[ , spectral.channel, drop = FALSE ]
    if ( nrow( df ) > max.points ) {
      set.seed( 42 )
      df <- df[ sample( nrow( df ), max.points ), , drop = FALSE ]
    }
    df |>
      dplyr::mutate( fluor = nm ) |>
      tidyr::pivot_longer(
        cols      = dplyr::all_of( spectral.channel ),
        names_to  = "channel",
        values_to = "value"
      )
  } ) )
  all.raw.sampled$channel <- factor( all.raw.sampled$channel, levels = spectral.channel )

  raw.density.plot <- ggplot2::ggplot(
    all.raw.sampled,
    ggplot2::aes( channel, biexp.transform( value ), colour = fluor )
  ) +
    ggplot2::geom_point( size = 0.8, alpha = 0.4, na.rm = TRUE ) +
    ggplot2::scale_colour_manual( values = fluor.colors, guide = "none" ) +
    ggplot2::scale_y_continuous(
      limits = biexp.transform( ribbon.limits ),
      breaks = biexp.transform( ribbon.breaks ),
      labels = ribbon.labels
    ) +
    ggplot2::labs( title = "Raw spectral data", x = "Detector", y = "Intensity" ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(
        angle = asp$ribbon.plot.axis.text.angle, vjust = 1, hjust = 1 ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none",
      plot.title       = ggplot2::element_text( size = 14, face = "bold", hjust = 0.5 )
    )

  # -- Static central heatmap (unmixing matrix) ---------------------------------
  hm.long <- as.data.frame( unmixing.matrix ) |>
    tibble::rownames_to_column( "Fluorophore" ) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of( spectral.channel ),
      names_to  = "Detector",
      values_to = "value"
    ) |>
    dplyr::mutate(
      Fluorophore = factor( Fluorophore, levels = fluor.names ),
      Detector    = factor( Detector,    levels = rev( spectral.channel ) )
    )

  # Width scaled tightly so it looks like a bridge between the two ribbons
  hm.width.ratio <- n.det / ( n.det + n.fluor )

  # Rotated 90 degrees: Fluorophore on X-axis, Detector on Y-axis
  heatmap.plot <- ggplot2::ggplot(
    hm.long,
    ggplot2::aes( Fluorophore, Detector, fill = value )
  ) +
    ggplot2::geom_tile( colour = "white", linewidth = 0.1 ) +
    ggplot2::scale_fill_viridis_c( option = color.palette, name = "Weight" ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text( angle = 45, hjust = 1, size = 7 ),
      axis.text.y  = ggplot2::element_text( size = 7 ),
      plot.title   = ggplot2::element_text( size = 14, face = "bold", hjust = 0.5 ),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = "Unmixing matrix",
      x = "Fluorophore  ->",
      y = "??? Detector"
    ) +
    ggplot2::coord_fixed()

  # -- Build animated point data: three global states, all fluorophores at once --
  # State 1 "raw"      -- all points at raw spectral positions  (left panel only)
  # State 2 "residual" -- all points transition to residual positions (left);
  #                      unmixed data appears on right panel simultaneously
  # State 3 "pause"    -- hold the final state of state 2 (both panels)
  #
  # Both panels share the same three state labels, guaranteeing identical frame
  # counts from gganimate and eliminating stitching mismatches.

  per.fluor <- lapply( seq_along( fluor.names ), function( fi ) {
    nm <- fluor.names[ fi ]

    raw.df <- as.data.frame( data.list[[ nm ]] )[ , spectral.channel, drop = FALSE ]

    # downsample
    if ( nrow( raw.df ) > max.points ) {
      set.seed( 42 + fi )
      raw.df <- raw.df[ sample( nrow( raw.df ), max.points ), , drop = FALSE ]
    }
    n.ev <- nrow( raw.df )

    # unmix
    unmixed.df <- as.data.frame( as.matrix( raw.df ) %*% t( unmixing.matrix ) )
    colnames( unmixed.df ) <- fluor.channel

    # residuals: raw minus full spectral reconstruction
    residual.mat <- as.matrix( raw.df ) - as.matrix( unmixed.df ) %*% spectra

    # globally unique event ids across fluorophores
    global.event.id <- paste0( fi, "_", seq_len( n.ev ) )

    # -- raw long --------------------------------------------------------------
    raw.long <- raw.df |>
      dplyr::mutate( event_id = global.event.id, fluor = nm, fluor_idx = fi ) |>
      tidyr::pivot_longer(
        cols      = dplyr::all_of( spectral.channel ),
        names_to  = "channel_raw",
        values_to = "value_raw"
      ) |>
      dplyr::mutate( y = biexp.transform( value_raw ) ) |>
      dplyr::select( event_id, channel_raw, y, fluor, fluor_idx )

    # -- residual long ---------------------------------------------------------
    residual.long <- as.data.frame( residual.mat ) |>
      dplyr::mutate( event_id = global.event.id, fluor = nm, fluor_idx = fi ) |>
      tidyr::pivot_longer(
        cols      = dplyr::all_of( spectral.channel ),
        names_to  = "channel_raw",
        values_to = "value_resid"
      ) |>
      dplyr::mutate( y = biexp.transform( value_resid ) ) |>
      dplyr::select( event_id, channel_raw, y, fluor, fluor_idx )

    # -- unmixed long (right panel) ??? ALL fluorophore columns -----------------
    # Each event appears in every unmixed column, not just its own fluorophore,
    # so the full unmixing result is shown. Points retain their source colour.
    unmixed.long <- unmixed.df |>
      dplyr::mutate( event_id = global.event.id, fluor = nm, fluor_idx = fi ) |>
      tidyr::pivot_longer(
        cols      = dplyr::all_of( fluor.channel ),
        names_to  = "channel_raw",
        values_to = "value_unmixed"
      ) |>
      dplyr::mutate( y = biexp.transform( value_unmixed ) ) |>
      dplyr::select( event_id, channel_raw, y, fluor, fluor_idx )

    list(
      raw.long      = raw.long,
      residual.long = residual.long,
      unmixed.long  = unmixed.long
    )
  } )

  # -- Combine across all fluorophores ------------------------------------------
  all.raw.long      <- dplyr::bind_rows( lapply( per.fluor, `[[`, "raw.long"      ) )
  all.residual.long <- dplyr::bind_rows( lapply( per.fluor, `[[`, "residual.long" ) )
  all.unmixed.long  <- dplyr::bind_rows( lapply( per.fluor, `[[`, "unmixed.long"  ) )

  # Fix channel_raw factor levels so scale_x_discrete works correctly
  all.raw.long$channel_raw      <- factor( all.raw.long$channel_raw,      levels = spectral.channel )
  all.residual.long$channel_raw <- factor( all.residual.long$channel_raw, levels = spectral.channel )
  all.unmixed.long$channel_raw  <- factor( all.unmixed.long$channel_raw,  levels = fluor.channel )

  # -- Build left-panel animation data (3 states) -------------------------------
  left.s1       <- all.raw.long;      left.s1$state <- "1_raw"
  left.s2       <- all.residual.long; left.s2$state <- "2_residual"
  left.s3       <- left.s2;           left.s3$state <- "3_pause"

  left.data       <- dplyr::bind_rows( left.s1, left.s2, left.s3 )
  left.data$state <- factor( left.data$state, levels = c( "1_raw", "2_residual", "3_pause" ) )

  # -- Build right-panel animation data (3 states) ------------------------------
  right.s1             <- all.unmixed.long; right.s1$state <- "1_raw";       right.s1$point_alpha <- 0.0
  right.s2             <- all.unmixed.long; right.s2$state <- "2_residual";  right.s2$point_alpha <- 0.6
  right.s3             <- right.s2;         right.s3$state <- "3_pause"

  right.data       <- dplyr::bind_rows( right.s1, right.s2, right.s3 )
  right.data$state <- factor( right.data$state, levels = c( "1_raw", "2_residual", "3_pause" ) )

  # -- Left panel animated plot -------------------------------------------------
  left.anim.plot <- ggplot2::ggplot(
    left.data,
    ggplot2::aes(
      x      = channel_raw,
      y      = y,
      colour = fluor,
      group  = interaction( fluor_idx, event_id )
    )
  ) +
    ggplot2::geom_point( size = 0.8, alpha = 0.5, na.rm = TRUE ) +
    ggplot2::scale_colour_manual( values = fluor.colors, guide = "none" ) +
    ggplot2::scale_x_discrete( limits = spectral.channel ) +
    ggplot2::scale_y_continuous(
      limits = biexp.transform( ribbon.limits ),
      breaks = biexp.transform( ribbon.breaks ),
      labels = ribbon.labels
    ) +
    ggplot2::labs(
      title = "Raw spectral data",
      x     = "Detector",
      y     = "Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(
        angle = asp$ribbon.plot.axis.text.angle, vjust = 1, hjust = 1 ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none",
      plot.title       = ggplot2::element_text( size = 14, face = "bold", hjust = 0.5 )
    ) +
    gganimate::transition_states(
      state,
      transition_length = c( state.length, transition.length, state.length ),
      state_length      = c( state.length, state.length,      state.length * 3 ),
      wrap              = FALSE
    ) +
    gganimate::ease_aes( "cubic-in-out" )

  # -- Right panel animated plot ------------------------------------------------
  right.anim.plot <- ggplot2::ggplot(
    right.data,
    ggplot2::aes(
      x      = channel_raw,
      y      = y,
      colour = fluor,
      alpha  = point_alpha,
      group  = interaction( fluor_idx, event_id )
    )
  ) +
    ggplot2::geom_point( size = 0.8, na.rm = TRUE ) +
    ggplot2::scale_colour_manual( values = fluor.colors, guide = "none" ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_x_discrete( limits = fluor.channel ) +
    ggplot2::scale_y_continuous(
      limits = biexp.transform( ribbon.limits ),
      breaks = biexp.transform( ribbon.breaks ),
      labels = ribbon.labels
    ) +
    ggplot2::labs(
      title = "Unmixed data",
      x     = "Fluorophore",
      y     = "Intensity"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(
        angle = asp$ribbon.plot.axis.text.angle, vjust = 1, hjust = 1 ),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "none",
      plot.title       = ggplot2::element_text( size = 14, face = "bold", hjust = 0.5 )
    ) +
    gganimate::transition_states(
      state,
      transition_length = c( state.length, transition.length, state.length ),
      state_length      = c( state.length, state.length,      state.length * 3 ),
      wrap              = FALSE
    ) +
    gganimate::ease_aes( "cubic-in-out" )

  # -- Render both panels --------------------------------------------------------
  if ( !dir.exists( figure.dir ) ) dir.create( figure.dir, recursive = TRUE )

  total.frames <- nframes
  px.w         <- round( plot.width  * res / 3 )   # per-panel width
  px.h         <- round( plot.height * res )

  message( "Rendering left panel animation..." )
  left.rendered <- gganimate::animate(
    left.anim.plot,
    nframes  = total.frames,
    fps      = fps,
    width    = px.w,
    height   = px.h,
    device   = "ragg_png",
    renderer = gganimate::gifski_renderer()
  )

  message( "Rendering right panel animation..." )
  right.rendered <- gganimate::animate(
    right.anim.plot,
    nframes  = total.frames,
    fps      = fps,
    width    = px.w,
    height   = px.h,
    device   = "ragg_png",
    renderer = gganimate::gifski_renderer()
  )

  # -- Stitch frames together with static heatmap in the centre -----------------
  # Save static heatmap as PNG, then read each frame pair and combine with magick
  message( "Stitching panels with static heatmap..." )

  heatmap.tmp <- tempfile( fileext = ".png" )
  ggplot2::ggsave(
    heatmap.tmp,
    plot   = heatmap.plot,
    width  = px.w / res * 0.4,
    height = px.h / res,
    dpi    = res,
    device = "png"
  )
  hm.img <- magick::image_read( heatmap.tmp )
  hm.img <- magick::image_scale( hm.img, paste0( px.w, "x", px.h ) )

  left.frames  <- magick::image_read( left.rendered )
  right.frames <- magick::image_read( right.rendered )

  # gganimate may produce fewer frames than requested (e.g. when nframes is
  # not evenly divisible across states/transitions). Derive the frame count
  # from the actual rendered objects to avoid subscript-out-of-bounds.
  n.left  <- length(left.frames)
  n.right <- length(right.frames)
  if ( n.left != n.right ) {
    warning(
      "Left panel has ", n.left, " frames but right panel has ", n.right,
      " frames. Trimming to the shorter of the two."
    )
  }
  actual.frames <- min( n.left, n.right )

  combined.frames <- lapply( seq_len( actual.frames ), function( i ) {
    magick::image_append(
      c( left.frames[i], hm.img[1], right.frames[i] ),
      stack = FALSE
    )
  } )

  combined.gif <- magick::image_animate(
    magick::image_join( combined.frames ),
    fps = fps
  )

  anim.file <- if ( !is.null( title ) )
    paste0( title, "_", animation.filename )
  else
    animation.filename

  out.path <- file.path( figure.dir, anim.file )

  if ( save ) {
    magick::image_write( combined.gif, path = out.path )
    message( "Animation saved to: ", out.path )
  }

  return( invisible( list(
    raw.plot     = raw.density.plot,
    heatmap.plot = heatmap.plot,
    right.plot   = right.anim.plot,
    animation    = combined.gif
  ) ) )
}
