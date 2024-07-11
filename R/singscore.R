#' Calculate singscore for specified genes
#'
#' This function calculates the singscore for a given set of genes across specified tissues and cell types.
#'
#' @param genes A character vector of gene names.
#' @param tissue A character vector of tissue types to include in the analysis. Default is all tissues in the metadata.
#' @param cell_type A character vector of cell types to include in the analysis. Default is all cell types in the metadata.
#' @param filter_condition An optional filtering condition to apply to the metadata.
#' @return A data frame with singscore results.
#'
#' @examples
#' singscoreHCA(c("CD3G", "CD8A"), tissue = "lung", cell_type = "t_nk")
#'
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#' @importFrom SummarizedExperiment assay
#' @importFrom CuratedAtlasQueryR get_metadata
#' @importFrom CuratedAtlasQueryR get_pseudobulk
#'
#' @export
singscoreHCA <- function(
    upSet,
    downSet = NULL,
    tissue = get_metadata() |> distinct(tissue_harmonised) |> pull(),
    cell_type = get_metadata() |> distinct(cell_type_harmonised) |> pull(),
    filter_condition = NULL
) {
  # Retrieve metadata and filter by specified tissue and cell type
  my_metadata <- get_metadata() |>
    filter(
      tissue_harmonised %in% !!tissue,
      cell_type_harmonised %in% !!cell_type
    )

  # Apply additional filtering condition if provided
  if (!is.null(filter_condition)) {
    my_metadata <- my_metadata |>
      filter(filter_condition)
  }

  # Generate pseudobulk data and calculate singscore
  my_metadata |>
    get_pseudobulk(assays = "counts") |>
    SE_to_singscore(upSet, downSet)
}

#' Convert SummarizedExperiment to singscore
#'
#' This function converts a SummarizedExperiment object to a singscore object for a given set of genes.
#'
#' @param se A SummarizedExperiment object containing the count data.
#' @param genes A character vector of gene names.
#' @return A data frame with singscore results.
#' @examples
#' SE_to_singscore(se, genes = c("gene1", "gene2"))
#' @importFrom SummarizedExperiment assay
#' @importFrom singscore simpleScore
#' @importFrom singscore rankGenes
#' @importFrom DelayedArray rowSums
#' @import tidySummarizedExperiment
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment colData
#'
#' @noRd
SE_to_singscore <- function(se, upSet, downSet) {

  # Extract count matrix from SummarizedExperiment object
  my_matrix <- assay(se, "counts", withDimnames = TRUE)

  # # Set row and column names (assuming sig_calc_pseudobulk is a predefined object)
  # rownames(my_matrix) <- rownames(se)
  # colnames(my_matrix) <- colnames(se)

  # Filter matrix and calculate singscore

  #my_matrix[DelayedArray::rowSums(my_matrix) != 0 | rownames(my_matrix) %in% genes, , drop = FALSE] |>

  # Ths is possibly a bug
  if(downSet |> is.null())
    my_score =
      my_matrix |>
      as.matrix() |>
      rankGenes() |>
      simpleScore(
        upSet = upSet,
        knownDirection = TRUE,
        centerScore = FALSE
      ) |>
      as_tibble(rownames = ".sample")
  else
    my_score =
      my_matrix |>
      as.matrix() |>
      rankGenes() |>
      simpleScore(
        upSet = upSet,
        downSet = downSet,
        knownDirection = TRUE,
        centerScore = FALSE
      ) |>
      as_tibble(rownames = ".sample")

  se |>
    colData() |>
    as_tibble(rownames = ".sample") |>
    left_join(my_score) |>
    select(
      .sample,
      cell_type_harmonised,
      tissue_harmonised,
      TotalScore,
      TotalDispersion,
      sex, age_days, ethnicity, assay,
      everything()
    )
}
