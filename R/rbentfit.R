#'
#' rank estimation for bent line regression
#'
#' This function use Wilcoxon score functions for
#' fitting the bent line regression model.
#' @rdname rbentfit
#' @param x A numeric variable with change point
#' @param z A vector of covariates
#' @param y A vector of response
#' @param bet.ini  A initial vector of regression coefficients
#' @param tau.ini  A initial value of change point
#' @param tol  tolerance value, 1e-4 for default
#' @param max.iter the maximum iteration steps
#'
#' @return A list with the elements
#' \item{est}{The estimated regression coefficients with intercept.}
#' \item{bp}{The estimated change point.}
#' \item{est.se}{The estimated standard error of the regression coefficients.}
#' \item{bp.est}{The estimated standard error of the change point.}
#' \item{iter}{The iteration steps.}
#'
#' @author Feipeng Zhang
#' @keywords rbentfit
#'
#' @import Rfit
#' @export
#'
#' @examples
#' n <- 150
#' x <- runif(n, 0, 4)
#' z <- rnorm(n, 1, 1)
#' y <- 1 + 0.5*z + 1.5*x  - 3 *pmax(x-2, 0)  + rt(n, 2)
#' rbentfit(y, cbind(1,z), x, bet.ini = c(0, 1, 1, -2), tau.ini = 1)
#'
#' # for the example of  MRS data
#' data(data_mrs)
#' x <- log(data_mrs$mass)
#' y <- log(data_mrs$speed)
#' z <- data_mrs$hopper
#' tau.ini <- 3
#' dat.new <- data.frame(y=y, z1=z, z2 = x, z3=pmax(x-tau.ini,0))
#' library(Rfit)
#' fit.ini <- rfit(y~ z1 + z2 +z3, data= dat.new)   # with intercept
#' bet.ini <- fit.ini$coef
#' fit.rank <- rbentfit(y, cbind(1,z), x, bet.ini, tau.ini)






rbentfit = function(y, z, x, bet.ini, tau.ini, tol=1e-4, max.iter=50)
{
  bet0 = c(bet.ini, 0.005)
  tau0 = tau.ini

  iter = 1
  while(iter <= max.iter) {
    #print(c(iter, bet0, tau0))

    #dat<- data.frame(y=y, z1=z, z2=x, z3=pmax(x-tau0,0), z4= (-1)* ifelse(x>tau0,1,0))
    U <- pmax(x-tau0,0)
    V <- (-1) * ifelse(x>tau0,1,0)
    nz <- ncol(as.matrix(z))
    XX <- cbind(z, x, U, V)
    colnames(XX)[1:nz] <- c(paste("z", 1:nz, sep = ""))
    colnames(XX)[(nz+1):ncol(XX)] <- c("x", "U", "V")
    fit1<- rfit(y~XX-1, data = data.frame(y, XX))
    bet1<- fit1$coef
    bet1.se <- summary(fit1)$coef[,2]

    m <- ncol(XX)
    tau1 =tau0 + 0.6*bet1[m]/(bet1[m-1]+1e-10)
    tau1.se <- bet1.se[m]/(abs(bet1[m-1])+1e-10)

    if( abs(bet1[m])< tol & max(abs(bet1-bet0))<tol) break

    # check if psi1 is admissible
    a<- (max(x)<=tau1)
    b<- (min(x)>=tau1)
    isErr<- sum(a+b)!=0|| is.na(sum(a+b))
    if(isErr){stop("estimated tau out of its range")}

    bet0 = bet1
    tau0 = tau1
    iter = iter+1  # increment index
  } # end while

  return(list(est=bet1[-m], bp=tau1, est.se=bet1.se[-m],
              bp.se=tau1.se, iter=iter))
}  # end function





