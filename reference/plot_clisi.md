# Making the plot for local enrichment

Making the plot for local enrichment

## Usage

``` r
plot_clisi(
  local_1,
  local_2,
  col_vec = (scales::hue_pal())(nrow(local_1$common_clisi$membership_info)),
  l_bg = 75,
  c_bg = 50,
  alpha_bg = 0.5,
  xlab1 = "Distinct enrichment",
  xlab2 = "Distinct enrichment",
  ylab = "Common enrichment",
  main1 = "Modality 1",
  main2 = "Modality 2",
  ...
)
```

## Arguments

- local_1:

  output of `clisi_information` on one Modality

- local_2:

  output of `clisi_information` on the other Modality

- col_vec:

  vector of colors

- l_bg:

  `l` parameter (luminosity, i.e., brightness) for individual cells

- c_bg:

  `c` parameter (chroma, i.e., color intensity) for individual cells

- alpha_bg:

  `alpha` parameter (color) for individual cells

- xlab1:

  `xlab` for Modality 1

- xlab2:

  `xlab` for Modality 2

- ylab:

  `ylab` for the shared modality

- main1:

  Title for plot corresponding to Modality 1

- main2:

  Title for plot corresponding to Modality 2

- ...:

  extra parameters for
  [`ggrepel::geom_text_repel`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)

## Value

List of two `gg` objects
