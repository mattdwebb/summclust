# School-program teaching example

These **fictional data** let you run the same analysis in Stata, R, or Python.
They do not contain real students or schools and are not evidence about a
real program.

The file [school_program.csv](school_program.csv) has 1,710 student observations
in 24 schools over 2019–2023. School sizes differ: most contribute 40–100
observations across all five years, while School 24 contributes 300. Schools
introduce the program at different dates; some never introduce it during
the observed period. The rows represent different student observations,
not a panel following the same student.

| Variable | Meaning |
|---|---|
| `school` | String school identifier, from `School 01` to `School 24`. |
| `year` | Year of the observation. |
| `program` | 1 if the program was available in that school and year, 0 otherwise. |
| `baseline_score` | Student's score before the outcome measurement. |
| `test_score` | Student's outcome score, in points. |

The main example regresses `test_score` on `program` and `baseline_score`,
includes year indicators, and clusters by school. It reports an adjusted
program difference and examines how the result depends on individual
schools. The regression specification alone does not establish a causal
effect.

## Run an example

Start in the top-level repository folder, not in this `examples` folder.
The language guides explain setup and interpret the displayed output.

| Language | Command | Guide |
|---|---|---|
| Stata | `do examples/school_program.do` | [Stata](../Stata/README.md) |
| R | `source("examples/school_program.R")` | [R](../R/README.md) |
| Python | `python examples/school_program.py` | [Python](../python/README.md) |

The Stata script also writes a log and plot to `output/`. The checked-in
[plot](school_program.png) was exported by Stata from this example. The R
and Python scripts print the coefficient table and the five schools with
the largest partial leverage. Output shown in the guides comes from these
scripts; the last displayed digit can differ across systems.
