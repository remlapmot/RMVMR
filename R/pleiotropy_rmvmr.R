#' pleiotropy_rmvmr
#'
#' Generates Q-statistics quantifying the degree of heterogeneity in univariate Radial MR analyses applying a correction using the
#' output from [`ivw_rmvmr`]. The function returns two data frames. The first data frame includes the global Q-statistic for each exposure after applying
#' a correction, as well as a corresponding p-value. The second data frame contains the individual Q-statistic for each SNP in the corrected univariate
#' analyses, relative to the exposure given in column \code{exposure}.
#'
#' @param r_input A formatted data frame using the [`format_rmvmr`] function or an object of class `MRMVInput` from [`MendelianRandomization::mr_mvinput`]
#' @param rmvmr An object containing the output from the [`ivw_rmvmr`] function of class \code{IVW_RMVMR}.
#'
#' @return An object of class \code{"RMVMR_Q"} containing the following components:
#' \describe{
#' \item{\code{gq}}{A data frame containing the global Q-statistic and p-value after applying a correction for each exposure}
#' \item{\code{qdat}}{A data frame containing the individual Q-statistic and p-value for each SNP after applying a correction for each exposure}
#' }
#'
#' @author Wes Spiller; Eleanor Sanderson; Jack Bowden.
#' @references Spiller, W., et al., Estimating and visualising multivariable Mendelian randomization analyses within a radial framework. Forthcoming.
#' @export
#' @examples
#' f.data <- format_rmvmr(
#'     BXGs = rawdat_rmvmr[,c("ldl_beta","hdl_beta","tg_beta")],
#'     BYG = rawdat_rmvmr$sbp_beta,
#'     seBXGs = rawdat_rmvmr[,c("ldl_se","hdl_se","tg_se")],
#'     seBYG = rawdat_rmvmr$sbp_se,
#'     RSID = rawdat_rmvmr$snp)
#' rmvmr_output <- ivw_rmvmr(f.data, FALSE)
#' q_object <- pleiotropy_rmvmr(f.data, rmvmr_output)
#' q_object$gq
#' head(q_object$qdat)
pleiotropy_rmvmr <- function(r_input, rmvmr) {
  # convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_rmvmr_format(r_input)
  }

  # Perform check that r_input has been formatted using format_rmvmr function
  if (
    !("rmvmr_format" %in%
      class(r_input))
  ) {
    stop(
      'The class of the data object must be "rmvmr_format", please resave the object with the output of format_rmvmr().'
    )
  }

  # Extract MVMR estimates
  rmvmr <- rmvmr$coef

  #Define number of exposures included in MVMR model
  exp.number <- length(names(r_input)[-c(1, 2, 3)]) / 2

  #Obtain univariate MR data for each exposure (F>10 SNPs)
  p.dat <- rmvmr_univariate_radial(r_input, exp.number)

  #Apply the MVMR correction and compute corrected Q-statistics
  cor.out <- rmvmr_correction(r_input, rmvmr, exp.number, p.dat)

  multi_return <- function() {
    Out_list <- list("gq" = cor.out$gq, "qdat" = cor.out$qdat)
    class(Out_list) <- "RMVMR_Q"

    return(Out_list)
  }

  OUT <- multi_return()
}
