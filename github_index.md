
# GitHub Repositories

## [Machine Learning from "Scratch"](https://github.com/archaeojsc/ML_from_scratch) &larrhk;

Building machine learning methods from basic principles is a great way to
understand the mathematical and algorithmic intuitions behind the methods. Yes
`scikit-learn` and such are quicker and easier (and generally faster), but it's
always good to have an understanding of what is going on "under the hood".

This repository will be my examples of a variety of machine learning methods
and algorithms from "scratch" (i.e., using minimal or base/common libraries).
Check back, as I'll be expanding the list as I have time and as I'm exploring
more on my own.

Currently available:  

* k-means clustering

Coming soon:  

* Isomap
* Random Forest

---

## [Data Science in Archaeology](https://github.com/archaeojsc/Archaeo_DS) &larrhk;

A "day in the life" repository of using data science as an archaeologist...
samples from some of my work projects and experiments (mainly in R). Much of
this work revolves around simple ETL and EDA, summarization and descriptive
statistics, manipulating and summarizing spatial data, and a few more "advanced"
experiments with probability distribution modeling.

The code in this repository is very rough for now. I'll clean them up into
formal scripts as I have time.

Currently available:

* Summarization utility functions for artifact counts, find depths, and
  temporally diagnostic dates using `tidyverse`
* Experiments with creating mixture distributions (using `mixdist`) for
  temporally diagnostic artifact dates comprised of the independent uniform
  distributions for each artifact type. Working on bootstrap estimates.
* Example automation code for generating standardized summaries for different
  artifact categories and spatial units for use in reports, piped to markdown
  tables using `tidyverse`. Incorporate spatial data pulled form ESRI shapefiles
  with `terra` and `sf`.

