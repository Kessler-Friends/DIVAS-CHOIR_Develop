# Software source ownership

`pkg/` is the R package built by CI and loaded by the hub's analysis scripts.
Edit implementation files in `pkg/R/`.

There is no top-level `R/` mirror. Keep a single implementation in `pkg/R/`.

Categorical signal extraction uses the multinomial parametric bootstrap in
`pkg/R/DIVASCHOIRSignalExtract.R`, with category probabilities derived from
the empirical proportions and fitted signal.
