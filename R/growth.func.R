#' Defines pollinator population growth function
#'
#' Generic growth function for bees population according to the available resources, using the inverse logit function
#' @param Nindiv : the current population size as a matrix
#' @param R : the resources as a matrix
#' @param a : the inflexion point in the logistic curve (i.e the median)
#' @param b : the slope of the growth
#' @param max : the asymptote
#' @references Haeussler J, Sahlin U, Baey C, Smith HG, Clough Y (2017) Predicting pollinator population size and pollination ecosystem service responses to enhancing floral and nesting resources. Ecology and Evolution, 7: 1898-1908.\url{http://dx.doi.org/10.1002/ece3.2765}

growth.func <- function(R,a,b,max)
{
  mu <- log(a)
  sigma <- sqrt(log(0.5+0.5*sqrt(1+4*b^2/a^2)))

  newN <- max * plnorm(R,mu,sigma)
  newN[!is.finite(newN)]= 0

  return(newN)
}
