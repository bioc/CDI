% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotFunctions.R
\name{contingency_heatmap}
\alias{contingency_heatmap}
\title{Heatmap of contingency table}
\usage{
contingency_heatmap(
  benchmark_label,
  candidate_label,
  proportion_size = 4,
  show_axis_names = TRUE,
  show_legend = TRUE
)
}
\arguments{
\item{benchmark_label}{A vector of characters indicating the benchmark 
cell type label set of cells.}

\item{candidate_label}{A vector of characters indicating the candidate 
clustering label set of cells.}

\item{proportion_size}{A number indicating the label size of proportion 
values inside each rectangle. 
The label will not be shown if this parameter is set to be FALSE.}

\item{show_axis_names}{A bool value indicating whether the axis names 
should be shown or not in the plot.}

\item{show_legend}{A bool value indicating whether the legend of methods 
should be shown or not in the plot.}
}
\value{
A ggplot object.
}
\description{
If the benchmark cell type label set is available, we can also compare one 
candidate label set 
(e.g. the optimal clustering label set selected by CDI) with the benchmark 
cell type labels. 
Here we provide the heatmap of contingency table for comparison.
Each row represents a cell type in the benchmark label set, 
and each column represents a cluster in the clustering label set. 
Each rectangle is color-scaled by the proportion of the cells in the given 
cluster coming from the benchmark types. 
Each column sums to 1.
}
\examples{
contingency_heatmap(benchmark_label = c(rep("type_a", 100), 
rep("type_b", 100)), 
candidate_label <- c(rep(1, 150), rep(2, 50)))
}