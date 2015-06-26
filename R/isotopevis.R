#' Standard Species
#' 
#' Takes a dataframe and identifies the species column (either "species", "Species", 
#' "taxon" or "Taxon") and matches the listed species against a dictionary list of 
#' standardized English names in the "species_lookup_table" The standardized names 
#' are then looked up in a standardized plotting table (plot_lookup_table) to 
#' find the plotting parameters. 
#' 
#' @param df The data.frame in which the isotope ratios 
#'  for all species are stored. 
#' @return A data.frame with the same number of rows as the original data frame,
#'    with the standardized English species name and other standardized grouping 
#'    information for each entry
#'  @examples
#'  out <- standard_species(df)

standard_species <- function(df){
  names <- colnames(df)
  species_names <-  df[,grepl("species", names) |
                         grepl("Species", names) |
                         grepl("taxon", names) |
                         grepl("Taxon", names)]
  species_lookup <- species
  plot_params <- plot_params
  matched_species <- species_lookup[match(
      species_names, species_lookup[, 2]), 1]
  matched_plot_params <- plot_params[match(
      matched_species, plot_params[, 1]), 1:ncol(plot_params)]
  standard <- cbind(df, matched_plot_params)
  data.frame(standard[, 2:ncol(standard)], row.names=1:nrow(standard))
}