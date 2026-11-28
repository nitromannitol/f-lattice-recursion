# f-lattice-recursion

This repository contains Julia code to produce tile odometers and tiles for the Abelian sandpile on the F-lattice, following the algorithm in

> Ahmed Bou-Rabee, *Integer superharmonic matrices on the F-lattice*, [arXiv:2110.07556](https://arxiv.org/abs/2110.07556). [PDF](https://nitromannitol.github.io/manuscripts/farey_recursion.pdf)

All heights are computed with exact integer arithmetic. The code was updated in November 2026, when Claude was used to repair a minor error in the original version.

## Instructions

1. Install Julia 1.7 or later: https://julialang.org/downloads/. No other packages are needed.

2. From this folder, draw the standard and alternate tiles for a word, for example `3312`:

   ```sh
   julia --project=. run_code.jl 3312 patterns.svg
   ```

   Words use the digits 1, 2 and 3. Each tile is colored by the Laplacian of its odometer: black is −1, gray is 0, and white is outside the tile. This is what the output looks like: [examples/3312.svg](examples/3312.svg)

3. To work with the tiles directly, start Julia with `julia --project=.` and run, for example

   ```julia
   using FLatticeExact
   str = "31"
   mats, alt_mats = get_odometers(str)
   p1, q1, p0, q0 = getFracs(str)
   ```

   `mats` holds the four standard tiles and `alt_mats` the four alternate tiles, in the order odd, even, odd parent, even parent. `getFracs` returns their four slopes in the same order. Each tile is a dictionary from lattice points, written as complex integers `x + y*im`, to exact integer heights. Points outside the tile are not in the dictionary.

4. To get the Laplacian of a tile at a point:

   ```julia
   field = periodic_pattern(mats[1], p1)
   pattern_value(field, 4im)     # Laplacian at (0, 4)
   periodic_height(field, 4im)   # height at (0, 4), extended periodically
   ```

   `periodic_pattern` throws an error if the tile is not consistent with its period.

---

- `include("init.jl")` instead of `using FLatticeExact` gives the same tiles as matrices. The first index is the x-coordinate, and points outside the tile are `NaN`, as is a one-cell border.
- To save a tile as a table of `x y height`, run for example `julia --project=. examples/export.jl 23 alt odd > word23.tsv`. Use `std` or `alt`, then `odd` or `even`.
- To run the tests: `julia --project=. test/runtests.jl`
