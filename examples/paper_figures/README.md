# Figures 3 and 6 of the paper

Run from the root of the repository with soapff importable (`PYTHONPATH=src`).

- `pca.py`: Fig. 3, the SOAP vectors of the O rows of the `clayff` table (Na montmorillonite) and of the `oplsaa` table (the
  O of the 88 references) on their first two principal components, coloured by type (the type is not used to find the
  axes). Writes `figs/pca.png` and `figdata/pca/numbers.json`.
- `soap_params.py`: Fig. 6, the typing with `r_cut`, `n_max` and `l_max` changed one at a time from those of the bundled
  libraries (`sigma` and the depth of 4 bonds kept). For each setting, a library that extends the bundled one with this
  `[soap]` (its table computed again), and the test structures typed with it:
  - `clayff`: the 10 (001) interfaces with water (`../interfaces_001_water`) and the 6 water shapes (`../interface_shapes`),
    structures made of a mineral and water although the library has the minerals and the water only apart;
  - `clayff_edge`: the particle with edges, the particle twice and the particle with water and ions (`../montmorillonite_edge`);
  - `oplsaa`: the 18 test structures of `../oplsaa/test` that are not MD frames.

```
SOAPFF_CACHE=WORK/cache PYTHONPATH=src python examples/paper_figures/soap_params.py run clayff WORK   # also clayff_edge, oplsaa
python examples/paper_figures/soap_params.py summary WORK     # -> figdata/soap_params/results.json
python examples/paper_figures/soap_params.py figure           # -> figs/soap_params.png
```

`figdata/soap_params/results.json` is the summary used for the paper: for every setting and structure, the atoms with
another type or charge than the default setting (clay) or Foyer (polymers), the atoms below the threshold 0.9, the lowest
similarity and the smallest margin.
