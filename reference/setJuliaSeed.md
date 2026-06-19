# Set Seed for julia's Random Number Generator

Sets the seed for julia's random number generator to ensure
reproducibility.

## Usage

``` r
setJuliaSeed(julia_seed)
```

## Arguments

- julia_seed:

  An integer seed value to be passed to julia's random number generator.

## Details

This function initializes the necessary julia functions and sets the
random seed for julia. If the provided seed is NULL, the function does
nothing.

## Author

Manuel Huth
