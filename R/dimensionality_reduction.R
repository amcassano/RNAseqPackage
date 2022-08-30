#' PCA Analysis
#'
#' Illustrate PCA plot to show sample relationships
#'
#' @param dseq_transform DESeqTransform object, either VST or RLog transformed
#' @param labels list, conditions being plotted
#' @param colors list, colors for each group
#' @param shapes list, shapes for each group
#' @param plot_title string, defaults to "Principal Component Analysis"
#' @param label_samples boolean, if true, labels for each sample will be added, defaults to TRUE
#'
#' @return PCA plot
#' @export
#'
#' @examples
pca_analysis <- function(dseq_transform,
                         labels,
                         colors,
                         shapes,
                         plot_title = "Principal Component Analysis",
                         label_samples = TRUE)
{
  # save the results of the PCA plot as an object so that we can plot the data using our good friend ggplot2
  pca_data <- DESeq2::plotPCA(dseq_transform, intgroup = "condition", returnData = TRUE)

  # get the percent variation of the PCA data so that we can add that to the plot
  percent_variation <- round(100 * attr(pca_data, "percentVar"))

  # create plot
  pca_plot <-
    ggplot2::ggplot(pca_data, aes(x = PC1,
                         y = PC2,
                         color = Condition,
                         shape = Condition)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(labels = labels,
                       values = colors,
                       name = "Condition:") +
    ggplot2::scale_shape_manual(labels = labels,
                       values = shapes,
                       name = "Condition:") +
    ggplot2::xlab(paste0("PC1: ", percent_variation[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_variation[2], "% variance")) +
    ggplot2::labs(title = plot_title)

  # display plot with global theme +/- sample labels
  if (label_samples) {
    pca_plot + global_theme() + add_labels(pca_data)
  }
  else {
    pca_plot + global_theme()
  }
}


#' UMAP Analysis
#'
#' plot umap for all samples, reducing the number of variables and showing relationships between samples
#'
#' @param dseq_transform DESeqTransform object, either VST or RLog transformed
#' @param neighbors integer, minimum of 2
#' @param m_dist integer
#' @param cond_list list, conditions for each group
#' @param plot_colors list (strings), colors for each group
#' @param plot_shapes list (integers), shapes for each group
#' @param plot_title string, plot title, defaults to "UMAP"
#' @param label_samples boolean, if sample labels should be added to plot, defaults to true
#' @param batch boolean, defaults to false, set to true & add set seed for reproducible results
#'
#' @return UMAP plot
#' @export
#'
#' @examples
umap_analysis <- function(dseq_transform,
                          neighbors,
                          m_dist,
                          cond_list,
                          plot_colors,
                          plot_shapes,
                          plot_title = "UMAP",
                          label_samples = TRUE,
                          batch = FALSE) {

  # data must be transposed before being able to be used as input for UMAP
  transposed_normDEG <- as.data.frame(SummarizedExperiment::assay(dseq_transform))
  transposed_normDEG <- t(transposed_normDEG)

  # condense the list
  labels <- unique(cond_list)

  # store the umap data in its own data frame
  zmap <- as.data.frame(uwot::umap(transposed_normDEG,
                             n_neighbors = neighbors,
                             min_dist = m_dist,
                             batch = batch))
  zmap <- zmap %>%
    dplyr::rename("UMAP1" = "V1") %>%
    dplyr::rename("UMAP2" = "V2")

  sample_names <- row.names(zmap)

  # create a umap table combining the condition, sample name, and zmap data
  umap_data <- cbind(cond_list, sample_names, zmap)
  umap_data <- umap_data %>%
    dplyr::rename("Group" = "cond_list") %>%
    dplyr::rename("SampleName" = "sample_names")

  umap_data$Group <- factor(umap_data$Group, levels = labels)
  head(umap_data)

  # plot the UMAP
  umap_plot <-
    ggplot2::ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = Group, shape = Group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(labels = cond_list,
                       values = plot_colors,
                       name = "Condition") +
    ggplot2::scale_shape_manual(labels = cond_list,
                       values = plot_shapes,
                       name = "Condition") +
    ggplot2::xlab("UMAP1") + ggplot2::ylab("UMAP2") +
    ggplot2::labs(title = plot_title)

  # display plot
  if (label_samples) {
    umap_plot + global_theme() + add_labels(umap_data)
  }
  else {
    umap_plot + global_theme()
  }
}
