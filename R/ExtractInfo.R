## -----------------------------------------------------------------------------
##          Internal helpers for extracting counts / batch labels
## -----------------------------------------------------------------------------

#' Extract raw UMI count matrix from input X
#'
#' Internal utility used by exported functions (e.g., \code{size_factor()},
#' \code{calculate_CDI()}). It standardizes how counts are retrieved from
#' different object classes.
#'
#' Supported input types:
#' \itemize{
#'   \item \strong{matrix}: returned as-is
#'   \item \strong{SingleCellExperiment}: retrieved from \code{assays(X)[[count_slot]]}
#'   \item \strong{Seurat}: retrieved via \code{SeuratObject::GetAssayData()}
#' }
#'
#' Notes for Seurat:
#' \itemize{
#'   \item In SeuratObject v5, counts are stored in \emph{layers} (argument \code{layer=})
#'   \item In earlier versions, counts are stored in \emph{slots} (argument \code{slot=})
#'   \item This function accepts \code{count_slot="counts"} (recommended), which is
#'         interpreted as a layer/slot name in the default assay (typically "RNA").
#'   \item For backward compatibility, if \code{count_slot} matches an assay name
#'         (e.g. "RNA"), it is interpreted as the assay name and counts are taken
#'         from the "counts" layer/slot.
#' }
#'
#' @param X A \code{matrix}, \code{Seurat} object, or \code{SingleCellExperiment} object.
#' @param count_slot For \code{SingleCellExperiment}, the assay name in \code{assays(X)}.
#'   For \code{Seurat}, either a layer/slot name (e.g. "counts") or (compatibility) an
#'   assay name (e.g. "RNA").
#'
#' @return A numeric matrix with genes in rows and cells in columns.
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment assays
#' @importFrom SeuratObject GetAssayData Assays DefaultAssay
#' @importFrom Seurat FetchData
#' @noRd
extract_count <- function(X, count_slot = NULL) {
  if (!is(X, "Seurat") && !is(X, "matrix") && !is(X, "SingleCellExperiment")) {
    stop("X must be a matrix, Seurat object, or SingleCellExperiment object.", call. = FALSE)
  }

  ## (1) matrix input: return as-is
  if (is(X, "matrix")) {
    return(as.matrix(X))
  }

  ## (2) SingleCellExperiment: counts are stored in assays(X)[[count_slot]]
  if (is(X, "SingleCellExperiment")) {
    if (is.null(count_slot)) {
      stop("count_slot must be specified for SingleCellExperiment input.", call. = FALSE)
    }
    a <- SummarizedExperiment::assays(X)[[count_slot]]
    if (is.null(a)) {
      stop("count_slot '", count_slot, "' not found in assays(X) for SingleCellExperiment input.", call. = FALSE)
    }
    return(as.matrix(a))
  }

  ## (3) Seurat: support SeuratObject v4 (slot=) and v5 (layer=)
  if (is.null(count_slot)) {
    stop("count_slot must be specified for Seurat input (e.g. 'counts').", call. = FALSE)
  }

  assay_names <- SeuratObject::Assays(X)
  default_assay <- SeuratObject::DefaultAssay(X)

  # Interpret count_slot:
  # - If it matches an assay name (e.g. "RNA"), treat it as assay and use "counts" as layer/slot.
  # - Otherwise treat it as layer/slot name within default assay.
  if (count_slot %in% assay_names) {
    assay_use <- count_slot
    layer_or_slot <- "counts"
  } else {
    assay_use <- default_assay
    layer_or_slot <- count_slot
  }

  ga_formals <- names(formals(SeuratObject::GetAssayData))

  gcmat <- tryCatch(
    {
      if ("layer" %in% ga_formals) {
        # SeuratObject v5
        SeuratObject::GetAssayData(object = X, assay = assay_use, layer = layer_or_slot)
      } else {
        # SeuratObject v4 and earlier
        SeuratObject::GetAssayData(object = X, assay = assay_use, slot = layer_or_slot)
      }
    },
    error = function(e) {
      stop(
        "Failed to extract counts from Seurat object. ",
        "Tried assay='", assay_use, "' and ",
        if ("layer" %in% ga_formals) "layer" else "slot",
        "='", layer_or_slot, "'. ",
        "Original error: ", conditionMessage(e),
        call. = FALSE
      )
    }
  )

  return(as.matrix(gcmat))
}



## -----------------------------------------------------------------------------
##          Internal function -- extract_batch
## -----------------------------------------------------------------------------

#' Extract batch labels from input X (or use provided batch_label)
#'
#' Internal utility to standardize retrieval of batch labels.
#'
#' Behavior:
#' \itemize{
#'   \item If \code{batch_label} is provided (non-NULL), it is returned as a vector.
#'   \item Otherwise, if \code{batch_slot} is provided:
#'     \itemize{
#'       \item For \code{Seurat}, batch is retrieved from meta.data via \code{FetchData(vars=batch_slot)}
#'       \item For \code{SingleCellExperiment}, batch is retrieved from \code{colData(X)[[batch_slot]]}
#'     }
#'   \item If both \code{batch_label} and \code{batch_slot} are NULL, returns NULL.
#' }
#'
#' @param X A \code{matrix}, \code{Seurat} object, or \code{SingleCellExperiment} object.
#' @param batch_label Optional vector of batch labels. If provided, \code{batch_slot} is ignored.
#' @param batch_slot If \code{batch_label} is NULL, the name of the batch field to retrieve
#'   from \code{Seurat} meta.data or \code{SingleCellExperiment} colData.
#'
#' @return A vector of batch labels, or NULL if no batch information is available.
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment colData
#' @importFrom Seurat FetchData
#' @noRd
extract_batch <- function(X, batch_label = NULL, batch_slot = NULL) {
  if (!is(X, "Seurat") && !is(X, "matrix") && !is(X, "SingleCellExperiment")) {
    stop("X must be a matrix, Seurat object, or SingleCellExperiment object.", call. = FALSE)
  }

  # If user supplies batch_label explicitly, use it.
  if (!is.null(batch_label)) {
    return(as.vector(unlist(batch_label)))
  }

  # If no slot specified, there is no batch information.
  if (is.null(batch_slot)) {
    return(NULL)
  }

  # Seurat: batch stored in meta.data; FetchData returns a data.frame
  if (is(X, "Seurat")) {
    df <- tryCatch(
      Seurat::FetchData(object = X, vars = batch_slot),
      error = function(e) {
        stop("batch_slot '", batch_slot, "' not found in meta.data of Seurat object.", call. = FALSE)
      }
    )
    # FetchData returns a data.frame with one column (typically). Use the first column.
    return(as.vector(df[[1]]))
  }

  # SingleCellExperiment: batch stored in colData
  if (is(X, "SingleCellExperiment")) {
    cd <- SummarizedExperiment::colData(X)
    if (!(batch_slot %in% colnames(cd))) {
      stop("batch_slot '", batch_slot, "' not found in colData of SingleCellExperiment object.", call. = FALSE)
    }
    return(as.vector(cd[[batch_slot]]))
  }

  # matrix input: no way to infer batch labels
  return(NULL)
}
