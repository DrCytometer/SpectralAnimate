# spectral_match_animate.R

utils::globalVariables( c(
  "Detector", "Intensity", "AF_Type", "status", "size", "alpha", "state"
) )

#' @title Spectral Signature Matching Animation
#'
#' @description
#' Produces an animation showing the automated matching of a cell's spectral
#' trace against a library of autofluorescence (AF) signatures using cosine
#' similarity.
#'
#' @section Animation states:
#' \describe{
#'   \item{Initial}{Displaying the library of AF candidates.}
#'   \item{Cell Reveal}{The unknown cell trace is overlaid as a dashed black line.}
#'   \item{Testing}{Each AF candidate is highlighted in orange sequentially.}
#'   \item{Winner}{The best match is identified and highlighted in green.}
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual scale_linewidth_identity
#' @importFrom ggplot2 scale_alpha_identity theme_minimal theme element_text labs
#' @importFrom gganimate transition_states shadow_mark animate anim_save
#'
#' @param af.matrix Matrix of AF signatures; rows are AF types, columns are
#'   detectors.
#' @param cell.trace Numeric vector representing the cell spectral signature to
#'   match. Must have the same length as the number of columns in `af.matrix`.
#' @param figure.dir Output directory for the saved GIF. Default is
#'   `"./af_matching"`.
#' @param title Optional character string prepended to the animation filename.
#'   Default is `"AF_Signature_Matching"`.
#' @param plot.width Width of the figure in inches. Default is `16`.
#' @param plot.height Height of the figure in inches. Default is `9`.
#'
#' @return The rendered `gganimate` animation object, invisibly. The GIF is also
#'   written to `figure.dir`.
#'
#' @export

spectral.match.animate <- function(
    af.matrix,
    cell.trace,
    figure.dir = "./af_matching",
    title = "AF_Signature_Matching",
    plot.width = 16,
    plot.height = 9
) {

  if ( !dir.exists( figure.dir ) ) dir.create( figure.dir )

  # Calculate Cosine Similarity to find the winner
  detectors <- colnames( af.matrix )
  calc_cosine <- function( a, b ) {
    sum( a * b ) / ( sqrt( sum( a^2 ) ) * sqrt( sum( b^2 ) ) )
  }

  similarities <- apply( af.matrix, 1, calc_cosine, b = cell.trace )
  winner_name  <- names( which.max( similarities ) )
  af_names     <- rownames( af.matrix )

  # Build long-format AF data using base R
  af_long <- data.frame(
    AF_Type   = rep( af_names, each = length( detectors ) ),
    Detector  = rep( detectors, times = length( af_names ) ),
    Intensity = as.vector( t( af.matrix ) ),
    stringsAsFactors = FALSE
  )
  af_long$Detector <- factor( af_long$Detector, levels = detectors )

  animation_list <- list()
  curr_state     <- 1

  # Step A: Initial Intro (AF only)
  for ( i in 1:2 ) {
    tmp <- af_long
    tmp$state         <- curr_state
    tmp$status        <- "Candidate"
    tmp$size          <- 0.7
    tmp$alpha         <- 0.6
    tmp$show_cell     <- FALSE
    tmp$current_label <- "Library Overview"
    animation_list[[ curr_state ]] <- tmp
    curr_state <- curr_state + 1
  }

  # Step B: Cell Trace Reveal
  for ( i in 3:4 ) {
    tmp <- af_long
    tmp$state         <- curr_state
    tmp$status        <- "Candidate"
    tmp$size          <- 0.7
    tmp$alpha         <- 0.6
    tmp$show_cell     <- TRUE
    tmp$current_label <- "Cell Background"
    animation_list[[ curr_state ]] <- tmp
    curr_state <- curr_state + 1
  }

  # Step C: Testing loop with FLASH
  for ( i in seq_along( af_names ) ) {
    current_af <- af_names[ i ]

    # Sub-step 1: The Highlight/Flash
    tmp <- af_long
    tmp$state  <- curr_state
    tmp$status <- ifelse(
      tmp$AF_Type == current_af, "Testing",
      ifelse( match( tmp$AF_Type, af_names ) < i, "Rejected", "Candidate" )
    )
    tmp$size          <- ifelse( tmp$AF_Type == current_af, 4, 0.7 )
    tmp$alpha         <- ifelse( tmp$AF_Type == current_af, 1, 0.2 )
    tmp$show_cell     <- TRUE
    tmp$current_label <- paste( "Testing:", current_af )
    animation_list[[ curr_state ]] <- tmp
    curr_state <- curr_state + 1

    # Sub-step 2: Settle/Observe
    tmp$state  <- curr_state
    tmp$size   <- ifelse( tmp$AF_Type == current_af, 2, 0.7 )
    animation_list[[ curr_state ]] <- tmp
    curr_state <- curr_state + 1
  }

  # Step D: The Winner
  tmp <- af_long
  tmp$state         <- curr_state
  tmp$status        <- ifelse( tmp$AF_Type == winner_name, "Winner", "Rejected" )
  tmp$size          <- ifelse( tmp$AF_Type == winner_name, 3, 0.7 )
  tmp$alpha         <- ifelse( tmp$AF_Type == winner_name, 1, 0.2 )
  tmp$show_cell     <- TRUE
  tmp$current_label <- paste( "Best Match Found:", winner_name )
  animation_list[[ curr_state ]] <- tmp

  plot_data <- do.call( rbind, animation_list )

  state_labels <- unique( plot_data[ order( plot_data$state ), c( "state", "current_label" ) ] )
  state_labels <- state_labels$current_label

  # Build the Plot
  cell.df <- plot_data[ plot_data$show_cell == TRUE, ]
  # one cell-trace row per (state x detector)
  cell.df <- cell.df[ !duplicated( cell.df[ , c( "state", "Detector" ) ] ), ]
  cell.df$Intensity <- rep( as.numeric( cell.trace ), length.out = nrow( cell.df ) )

  anim.plot <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = plot_data,
      ggplot2::aes(
        x         = Detector,
        y         = Intensity,
        group     = AF_Type,
        color     = status,
        linewidth = size,
        alpha     = alpha
      )
    ) +
    ggplot2::geom_line(
      data = cell.df,
      ggplot2::aes( x = Detector, y = Intensity, group = 1 ),
      color     = "black",
      linewidth = 1.5,
      linetype  = "dashed"
    ) +
    ggplot2::scale_color_manual( values = c(
      "Candidate" = "steelblue",
      "Testing"   = "orange",
      "Rejected"  = "gray80",
      "Winner"    = "forestgreen"
    ) ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x    = ggplot2::element_text( angle = 45, hjust = 1 ),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title    = "{state_labels[as.integer(closest_state)]}",
      subtitle = "Dashed line = Cell Background"
    ) +
    gganimate::transition_states(
      state,
      transition_length = 1,
      state_length      = 2,
      wrap              = FALSE
    ) +
    gganimate::shadow_mark( past = FALSE, future = FALSE )

  # Render
  rendered <- gganimate::animate(
    anim.plot,
    nframes   = curr_state * 8,
    end_pause = 40,
    width     = plot.width,
    height    = plot.height,
    units     = "in",
    res       = 150
  )

  animation.file.name <- paste0( title, ".gif" )

  gganimate::anim_save(
    filename  = animation.file.name,
    animation = rendered,
    path      = figure.dir
  )

  return( invisible( rendered ) )
}
