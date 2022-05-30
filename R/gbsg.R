#' @title Breast cancer data sets used in Royston and Altman (2013)
#'
#' @description The \code{gbsg} data set contains patient records from a 1984-1989 trial conducted by the German Breast Cancer Study Group (GBSG) of 720 patients with node positive breast cancer; it retains the 686 patients with complete data for the prognostic variables.
#'
#' @usage gbsg
#'
#' @format A data frame with 686 rows and 11 variables:
#' \describe{
#'   \item{\code{pid}}{Patient identifier}
#'   \item{\code{age}}{Age, years}
#'   \item{\code{meno}}{Menopausal status (0= premenopausal, 1= postmenopausal)}
#'   \item{\code{size}}{Tumor size, mm}
#'   \item{\code{grade}}{Tumor grade}
#'   \item{\code{nodes}}{The number of positive lymph nodes}
#'   \item{\code{pgr}}{Progesterone receptors (fmol/l)}
#'   \item{\code{er}}{Estrogen receptors (fmol/l)}
#'   \item{\code{hormon}}{Hormonal therapy, 0= no, 1= yes}
#'   \item{\code{rfstime}}{Recurrence free survival time; days to first of reccurence, death or last follow-up}
#'   \item{\code{status}}{Censoring indicator, 0= alive without recurrence, 1= recurrence or death}
#'}
#'
#' @details These data sets are used in the paper by Royston and Altman (2013). The Rotterdam data is used to create a fitted model, and the GBSG data for validation of the model. The paper gives references for the data source.
#'
#' @references Patrick Royston and Douglas Altman, External validation of a Cox prognostic model: principles and methods. \emph{BMC Medical Research Methodology} 2013, 13:33
#'
#' @source \code{\link[survival]{gbsg}} in \code{survival} package
"gbsg"
