# R CMD check does not like the variable used by foreach::foreach()
#
# * checking R code for possible problems ... NOTE
# cycle_npreg_mstep: no visible binding for global variable 'iii'
# data_transform_quantile: no visible binding for global variable 'iii'
# fit_cyclical_many: no visible binding for global variable 'iii'
# Undefined global functions or variables:
#   iii
#
# Hack to make the NOTE go away:
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-
# visible-binding-for-global-variable-notes-when
#' @importFrom utils globalVariables
if(getRversion() >= "2.15.1") utils::globalVariables("iii")
