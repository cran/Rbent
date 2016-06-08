#'
#' test the existence of change point in the bent line regression
#'
#' This function use Wilcoxon score functions for  calculating
#' the test statistics and p-value by wild bootstrap.
#'
#' @rdname rbenttest
#' @param x A numeric variable with change point
#' @param z A vector of covariates
#' @param y A vector of response
#' @param NB resampling times
#' @param myseed set seed

#'
#' @return A list with the elements
#' \item{Tn}{The statistic based on original data.}
#' \item{Tn.NB}{The statistics by wild bootstrap.}
#' \item{p.value}{The p-value by wild bootstrap.}
#'
#' @author Feipeng Zhang
#' @keywords rbenttest
#'
#' @import Rfit
#' @importFrom stats approxfun density ecdf residuals rnorm runif sd
#' @export
#'
#' @examples
#' # for the example of  MRS data
#' data(data_mrs)
#' x <- log(data_mrs$mass)
#' y <- log(data_mrs$speed)
#' z <- data_mrs$hopper
#' p.value <- rbenttest(y, cbind(1, z), x, NB = 50)$p.value
#'
#' # for the example of bedload transport data
#' data(data_transport)
#' x <- data_transport$x
#' y <- data_transport$y
#' p.value <- rbenttest(y, 1, x, NB = 50)$p.value


rbenttest <- function(y, z, x, NB = 1000, myseed = 1 ){

  ## rank-based dispersion function
  wilcoxon <- function(t){
    sqrt(12)*(t-0.5)
  }

  ## test statistic based on origin data
  testFun <- function(tau){

    ## under H0

    n <- length(y)

    ## by rfit
    xz <- cbind(z, x)
    fit.rank <- rfit(y~ xz-1)
    res.rank <- as.vector(residuals(fit.rank))
    phi.res <- wilcoxon(rank(res.rank)/(n+1))

    Rn.rank <- rep(0, length(tau))
    for (kk in 1:length(tau)){

      Rn.rank[kk] <- 1/sqrt(n)*sum(
        phi.res*(x-tau[kk])*ifelse(x<=tau[kk], 1, 0)
      )
    }

    Tn.rank <- max(abs(Rn.rank))

    return(Tn.rank)
  }

  ## perturbed method to calculate the p-value
  testFun.resample <- function(tau){

    #########################
    ## permutation random errors
    n<- length(y)

    e <- rnorm(n, 0, 1)
    B <- runif(n, 0, 1)
    #v <- (1-sqrt(5))/2*(B<= (sqrt(5)+1)/(2*sqrt(5))) +
    #       (1+sqrt(5))/2*(B> (sqrt(5)+1)/(2*sqrt(5)))

    v <- (-1)*(B<=0.5) +1*(B>0.5)
    ve <- v*e


    #########################
    ## Sn
    xz <-  cbind(z, x)
    Sn <- (t(xz)%*%xz)/n

    xc <- xz - mean(xz)
    Sn.rank <- mean(xc^2)

    ## residual by rfit
    fit.rank <- rfit(y~ xz-1)
    res.rank <- as.vector(residuals(fit.rank))

    #u <- sample(res.rank, n, replace =T)

    ## estimation for density value at residuals
    h <- 1.06* n^(-1/5)* sd(res.rank)
    pdf<- approxfun(density(res.rank, kernel= "epanechnikov", bw=h))
    fe.rank <- pdf(res.rank)

    cdf <- ecdf(res.rank)
    Fe.rank <- cdf(res.rank)

    ## estimation for scales
    taus <- fit.rank$taushat
    tauphi <- fit.rank$tauhat

    ## under H0
    Rn.rank <- rep(0, length(tau))
    for (kk in 1:length(tau)){


      # rfit
      Sn.rank.tau <- apply(
        sqrt(12)*fe.rank*xz*(x-tau[kk])*ifelse(x<= tau[kk], 1, 0), 2, mean
      )
      Rn.rank[kk] <- 1/sqrt(n)*sum(
        ve*(
          wilcoxon(Fe.rank)*(x-tau[kk])*ifelse(x<= tau[kk],1,0)-
            wilcoxon(Fe.rank)* tauphi* xz%*% solve(Sn) %*% Sn.rank.tau
        )
      )
    }
    Tn.rank <- max(abs(Rn.rank))

    return(Tn.rank)
  }

  #######################################################
  ###  calculate the p-value by wild bootstrap
  set.seed(myseed)

  tau <- seq(min(x)+0.01, max(x)-0.01, length=100)
  Tn <-  testFun(tau)
  Tn.NB <- replicate(NB, testFun.resample(tau))

  pv <- mean(Tn.NB >Tn,  na.rm=TRUE)

  return(list(Tn = Tn, Tn.NB = Tn.NB, p.value = pv))


}
