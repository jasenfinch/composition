#' @importFrom KEGGREST keggLink
#' @importFrom dplyr distinct
#' @importFrom tibble deframe
#' @importFrom stringr str_remove_all

keggCompounds <- function(organism = character()){
  
  if (length(organism) == 0) {
    compounds <- metabolites %>%
      getAccessions() %>%
      .$ACCESSION_ID
  } else {
    enzymes <- keggLink(organism,'enzyme') %>%
      names()
    compounds <- keggLink('compound','enzyme') %>%
      {tibble(Enzyme = names(.),Compound = .)} %>%
      filter(Enzyme %in% enzymes) %>%
      select(Compound) %>%
      distinct() %>%
      deframe() %>%
      str_remove_all('cpd:')
  }
  return(compounds)
}

#' keggPIPs
#' @description Molecular formula putative ionisation products from KEGG.
#' @param MFs Molecular formulas and adducts to search. Should be a tibble containing two character columns named MF and Adduct.
#' @param organism KEGG organism ID.
#' @param adductRules table containing adduct formation rules. Defaults to mzAnnotation::adducts().
#' @importFrom mzAnnotation filterACCESSIONS filterMF filterIP getAccessions adducts
#' @importFrom dplyr everything filter
#' @importFrom tibble is_tibble
#' @export

keggPIPs <- function(MFs, organism = character(), adductRules = adducts()){
  
  if (!is.data.frame(MFs)) {
    stop('MFs should be a tibble or data frame.')
  }
  
  if (!identical(colnames(MFs),c('MF','Adduct'))) {
    stop('MFs should contain the columns MF and Adduct with molecular formulas and their putative adducts respectively.')
  }
  
  compounds <- keggCompounds(organism)
  
  met <- metabolites %>%
    filterACCESSIONS(compounds)
  
  MFs %>%
    split(1:nrow(.)) %>%
    map(~{
      rule <- adductRules %>%
        filter(Name == .x$Adduct)
      met %>%
        filterMF(.x$MF) %>%
        filterIP(rule = rule$Rule) %>%
        getAccessions() %>%
        mutate(MF = .x$MF,Adduct = .x$Adduct) %>%
        select(ACCESSION_ID,MF,Adduct,everything())
    }) %>%
    bind_rows()
}
