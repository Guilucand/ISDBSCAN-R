#' @importFrom Rcpp sourceCpp
#' @useDynLib ISDBSCAN, .registration = TRUE
#'
NULL


#' @title Density based algorithm (ISDBSCAN) for large single cell sequencing data
#'
#' @description This is an implementation of the ISDBSCAN algorithm of
#'   Cassisi et al (2013) adapted to manage large single cell sequencing data.
#'
#' @details The contribution of this package is to provide support for on-disk
#'   data representations such as HDF5, through the use of \code{DelayedMatrix}
#'   and \code{HDF5Matrix} objects, as well as for sparse data representation
#'   through the classes of the \code{Matrix} package.



#'@rdname ISDBSCAN
#'@importFrom methods is
#'@export
#'@importClassesFrom DelayedArray DelayedMatrix
#'@param x the input data
#'@param k the number of neighbours
#'@param batch_size the size of the readed chunk
#'@param stratif either TRUE or FALSE, TRUE if we want to stratifie the data, FALSE otherwise.
#'@return a list with the following attributes: clusters, layer, border
#'@details This function performs ISDBSCAN clustering.
#'
#'
setMethod(
    "ISDBSCAN",
    signature = signature(x ="ANY"),
    definition = function(x, k=3, batch_size = blocksize(x), stratif = FALSE )
    {

        if(!is(x, "matrix") & !is(x, "Matrix") & !is(x, "HDF5Matrix") &
           !is(x, "DelayedMatrix")) {

            stop("x is not of a supported type")

        } else {

            fit <- cIsdbscan(x, k, batch_size, stratif)

        }

        return(fit)
    })

