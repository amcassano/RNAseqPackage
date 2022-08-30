#' Plot Gene Counts
#'
#' @description Plot the gene counts of a certain gene of interest
#'
#' @param geneid string, Ensembl Gene ID
#' @param dseq_obj DESeq object
#' @param colors list, colors to use for each condition
#' @param shapes list, shapes to use for each condition
#' @param labels list, labels to use for each condition
#' @param plot_title string, plot title, defaults to "Gene Counts"
#' @param plot_subtitle string, subtitle for plot, defaults to empty string
#' @param groupby string, interest group, defaults to "Condition"
#' @param label_samples boolean, indicates whether samples should be labeled on plot, defaults to TRUE
#'
#' @return image
#' @export
#'
#' @examples
plot_genecounts <- function(geneid, dseq_obj,
                            plot_colors,
                            plot_shapes,
                            plot_labels,
                            plot_title = "Gene Counts",
                            plot_subtitle = "",
                            groupby = "Condition",
                            label_samples = TRUE)
{
  # save the data from the plotCounts function
  plot_data <- DESeq2::plotCounts(dseq_obj,
                                  gene = geneid,
                                  intgroup = groupby,
                                  returnData = TRUE)
  # plot the counts
  countPlot <- ggplot2::ggplot(plot_data,
                               ggplot2::aes(x = Condition,
                                            y = count,
                                            color = Condition,
                                            shape = Condition)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::ggtitle(plot_title,
            subtitle = plot_subtitle) +
    ggplot2::scale_color_manual(labels = plot_labels,
                       values = plot_colors,
                       name = "Condition:") +
    ggplot2::scale_shape_manual(labels = plot_labels,
                       values = plot_shapes,
                       name = "Condition:") +
    ggplot2::xlab(paste("Condition")) +
    ggplot2::ylab(paste0("Normalized Count"))

  #display plot with global theme and labels
  if (label_samples) {
    countPlot + global_theme() + add_labels(plot_data)
  }
  else {
    countPlot + global_theme()
  }
}
