#'  Draw a manhattan plot from data
#'
#' \code{manhattan_plot()} Create a manhatan plot
#' 
#' @import ggplot2 
#' @export 
#'
#' @param data.table: has to have the columns CHR, POS, p.value 
#' @param fig_name: figure name with extension (.png, .eps, .pdf): Defaut NULL 
#' @param facet_key  
#'
#' @return ggplot  


manhattan_plot <-
  function(DT_in, 
           fig_name = NULL,
           facet_key = NULL){

    DT <- copy(DT_in)

    DT[, position := .GRP, keyby = .(CHR, POS)]

    NCHR = DT[.N, CHR]

    pos_shift <- 
      DT[, round(0.2 * .N / NCHR)]

    DT[, position := position + (CHR - 1) * pos_shift]

    b_info <-
      DT[, .(breaks = 0.5 * (first(position) + last(position))), by = CHR]


    gm <- 
      ggplot(DT, aes(x = position, y = p.value) ) + 
        geom_point(size = 0.5) +
        theme_bw(base_family = "serif", base_size = 11) + 
        theme(panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.line = element_line(color = "black", size = 0.3)) + 
        scale_x_continuous(breaks = unlist(b_info[, breaks]), 
                           label = unlist(b_info[, CHR])) 

    if(!is.null(facet_key))
      gm <- gm +
        facet_wrap(as.formula(paste('~',facet_key)),
                   ncol = 1, 
                   scales = 'free')

    if(!is.null(fig_name))
      ggsave(fig_name)

    return(gm)
  }
