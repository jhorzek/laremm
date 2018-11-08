#' compare regsem and laremm
regsem_v_regctsem <- function(mx_regmodel, lav_regmodel, Jacob_equivalent = F){
  if(Jacob_equivalent == F){
    mx_FML <- getFML(mx_regmodel) # compute FML for mxModel
  }

  lav_FML_pen <- lav_regmodel$optim_fit
    # .5*penalty + FML; not identical with formula reported in Jacobucci 2016
  return_value <- matrix(c(mx_FML[3],lav_FML_pen ), ncol = 2)
  colnames(return_value) <- c("mx FML", "lavaan FML")
  return(return_value)
}

