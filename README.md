
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SpectralAnimate

<!-- badges: start -->

<!-- badges: end -->

SpectralAnimate is just a bit of fun to animate transitions between some
of the graphics generated in the
[AutoSpectral](https://github.com/DrCytometer/AutoSpectral) package. The
rendering is slow, although I’ve done a little work to speed it up.
Anyway, can be fun to see transitions between matrices, or layering of
curves. These functions were used to build the CytoBytes YouTube video.
Some functions require parameters from AutoSpectral on a
cytometer-specific basis.

## Installation

You can install the development version of SpectralAnimate from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("DrCytometer/SpectralAnimate")
```

## Example

To animate the transition between the spectral profiles heatmap (mixing
matrix) and the unmixing matrix heatmap:

``` r
library(SpectralAnimate)
unmixing.matrix.animate(spectra,
                        color.palette = "viridis")
```

<figure>
<img src="man/figures/unmixing_matrix_transition.gif"
alt="Unmixing matrix transition" />
<figcaption aria-hidden="true">Unmixing matrix transition</figcaption>
</figure>
