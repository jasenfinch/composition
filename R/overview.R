#' plotStructuralOverview
#' @description plot a matrix structural overview
#' @param structural_classifications should be a tibble containning structural classifications. See details below
#' @details The tibble columns should firstly have MF, Adduct then followed by the structural classification taxonomic levels.
#' @importFrom plotly plot_ly
#' @importFrom dplyr select bind_rows group_by summarise n mutate group_by_all
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom magrittr %>% set_names
#' @export

plotStructuralOverview <- function(structural_classifications){
  con <- structural_classifications %>%
    select(kingdom:names(.)[ncol(.) - 1]) 
  
  tree <- 1:ncol(con) %>%
    map(~{
      i <- .
      if (i == 1) {
        tibble(level = 'kingdom',
               label = con$kingdom,
               parent = '')
      } else {
        con %>%
          select(i,i-1) %>%
          set_names(c('label','parent')) %>%
          mutate(level = names(con)[i])
      }
    }) %>%
    bind_rows() %>%
    group_by_all() %>%
    summarise(N = n())
  
  plot_ly(tree,labels = ~label,parents = ~parent,values = ~N,type = 'sunburst')
}