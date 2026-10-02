# R_work – current code

Czech survey on reproducibility (AV ČR, UK, MU, JČU; 2019) compared with Baker (2016, *Nature*).
Task: P. Jedlička, memo of 21 Sept 2026 ("Data a komparace – Martin").

## Pipeline

1. **`preprocess.Rmd`** – reads the four institutional exports from `../data_and_R_src/objektivita{AV,UK,MU,JU}.xls`,
   gives the variables short names and merges them. Removes one test response (before 11 Apr 2019) and
   saves **`dataAll.RData`** (raw answers, N = 1,021; `old_names` holds the original question texts).
   The rest of the script is F. Kalvas's objectivity analysis; its recoded data are kept in memory only.
2. **`common.R`** – shared preparation, sourced by both reports: Czech filters and respondent flow,
   Baker data (`../data_and_R_src/dataNature.xlsx`, columns checked against question text), harmonised
   answers, field weights and helper functions.
3. **`sample.Rmd`** → `sample.html` – Czech-only description: Q1–Q3, Q38–Q49, and Q9d (no Baker counterpart).
4. **`replication.Rmd`** → `replication.html` – comparison with Baker:
   - **Part 1, main results:** Cramér's V per question, miniature comparisons, Q6 by field, and Q18–Q20
     'yes' share by field with Kruskal–Wallis tests across fields (per survey).
     Only answer categories offered in both surveys ("don't know" kept in Q6/Q16, dropped in Q8/Q9;
     Q18–Q20 without non-experimenters).
   - **Part 2, detailed results (supplementary):** each question in full, field-adjusted Czech results,
     Q8 means, by-field tables.
   - **Part 3, sensitivity:** the same with all answer categories.
   - **Appendix:** Baker reference values (must match the published figures).

Question numbers follow Petr's memo. Both HTML reports are self-contained.

Packages: dplyr, tidyr, readxl, readr, forcats, stringr, ggplot2, knitr, rmarkdown, effectsize.

## Switches (YAML `params` in both Rmds)

| Parameter | Default | Meaning |
|---|---|---|
| `time_cutoff_min` | 0 | Minimum completion time in minutes (0 = none; 10/12/15 for sensitivity) |
| `basic_only` | FALSE | Keep only respondents doing basic research (Q3) |
| `drop_math` | TRUE | Drop mathematicians (no Baker counterpart) |

Set the same values in both Rmds before knitting.

## Other files

- `nature.Rmd`, `nature.html` – F. Kalvas's original comparison (superseded by `replication.Rmd`; kept for
  reference; the Q9 validity item mapping is fixed). It still loads `dataProcessed.RData`, which is no
  longer produced, so it does not knit as is.
- `preprocess.html` – rendered `preprocess.Rmd`.
- `east_and_west_europe.csv` – used by `nature.Rmd` only.

## Open decisions

Time cutoff, basic vs applied research, mathematicians, field grouping, Q6 missing answers (all outside AV),
Q9d, Q16 domain grouping – see the meeting log.
