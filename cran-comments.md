## R CMD check results

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Kevin Berg <Kevin.Berg@informatik.uni-regensburg.de>'

Resubmission
* checking for future file timestamps ... NOTE
unable to verify current time

0 errors | 0 warnings | 2 notes

## Resubmission

This is a resubmission of version 0.1.1.

In the previous CRAN check, examples were executed with `--run-donttest` and failed because they rely on external vignette objects (for example `data`, `trajectories`, `D.list`) that are not available during checks.

I fixed this by changing affected examples from `\donttest{}` to `\dontrun{}` in:

* `Hetseq`, `HetseqClassify`, `HetseqDoubleML`
* `PlotClassify`, `PlotDoubleML`
* `distmat`, `prune`, `mincostflow`

and synchronized the corresponding `.Rd` files.
