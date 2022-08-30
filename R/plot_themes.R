#' Add Labels
#'
#' adds labels to each point on a ggplot
#'
#' @param dat data frame
#'
#' @return geom label repel object
#' @export
#'
#' @examples
add_labels <- function(dat){
  return(ggrepel::geom_label_repel(data = dat,
                          ggplot2::aes(label = rownames(dat)),
                          min.segment.length = 0,
                          size = 3.5,
                          point.padding = 0.8,
                          label.padding = 0.2,
                          box.padding = 1,
                          force_pull = 3,
                          show.legend = FALSE))
}

#' Global Theme
#'
#' Sets the global GGplot theme
#'
#' @return ggplot theme
#' @export
#'
#' @examples
global_theme <- function(){
  return(ggplot2::theme(panel.background = ggplot2::element_rect(fill = "transparent"),
               plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
               #panel.grid.major = element_blank(), # get rid of major grid
               panel.grid.minor = ggplot2::element_blank(), # get rid of minor grid
               axis.line = ggplot2::element_line("black", 0.3),
               aspect.ratio = 1,
               plot.title = ggplot2::element_text(hjust = 0.5),
               plot.subtitle = ggplot2::element_text(hjust = 0.5)))
}
