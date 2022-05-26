#' @title Breast cancer data sets used in Royston and Altman (2013)
#' @description The \code{gbsg} data set contains patient records from a 1984-1989 trial conducted by the German Breast Cancer Study Group (GBSG) of 720 patients with node positive breast cancer; it retains the 686 patients with complete data for the prognostic variables.
#' @format A data frame with 686 rows and 11 variables:
#' \describe{
#'   \item{\code{pid}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{age}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{meno}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{size}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{grade}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{nodes}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{pgr}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{er}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{hormon}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{rfstime}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{status}}{integer COLUMN_DESCRIPTION}
#'}
#' @details These data sets are used in the paper by Royston and Altman. The Rotterdam data is used to create a fitted model, and the GBSG data for validation of the model. The paper gives references for the data source.
#' @references Patrick Royston and Douglas Altman, External validation of a Cox prognostic model: principles and methods.  BMC Medical Research Methodology 2013, 13:33
#'
#' @source \code{\link[survival]{gbsg}}
"gbsg"