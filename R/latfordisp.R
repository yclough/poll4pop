#' Apply kernels to model foraging or dispersal in space
#' 
#' For foragers occupying a or more cells of a matrix, applies a spatial kernel
#' weighted by the available resources such that foragers appear in cells that are
#' both close-by and have higher resources.
#' 
#' @param N the matrix of foragers (source units)
#' @param alpha the matrix of resources
#' @param decay the dispersal kernel, for example produced by kerncalc
#' @return a list containing :
#' \describe{
#'   \item{Racc_per_forager}{the number of queens at the end of each period and at the beginning of the next season}
#'   \item{visitation_rate}{the visitation rates per cell}
#'   \item{distance_weighted_resources}{matrix of resources weighted by the kernal}
#'  } 
#' 
latfordisp <- function(N,
                       alpha,
                       decay)
{
  alpha[is.na(alpha)] <- 0
  N[is.na(N)] <- 0
  
  # cell-specific sum of the distance-weighted resources
  distWeigthedResources = EBImage::filter2(alpha, decay) 
  distWeigthedResources[distWeigthedResources < 1e-15] <- 0 # replace values that are smaller than the zero-machine by zero
  
  # number of resources per individual, a.k.a repartition of the foragers on the resources
  # e.g. if there are 2 flowers and only 1 bee, it would means that on average there is 1/2 bee per flower
  # on the other hand if there are 6 bees and 2 flowers, there would be approx. 3 bees per flower
  relativeResources <- N*0 # initializing the matrix to have the same size as N
  relativeResources[distWeigthedResources > 0] = N[distWeigthedResources > 0] /  distWeigthedResources[distWeigthedResources > 0]
  relativeResources[distWeigthedResources < 0] = 0 ## the total number of foragers are not preserved here
  
  # dispersion of relative resources
  distWeightedRelativeResources = EBImage::filter2(relativeResources, decay)
  distWeightedRelativeResources[distWeightedRelativeResources < 1e-15] <- 0
  
  visitationRate = distWeightedRelativeResources * alpha
  # the resources accessed by foragers are the distance-weighted resources but only at cells where foragers are actually nesting
  Racc_per_forager = distWeigthedResources*(N>0) 
  
  return(list(Racc_per_forager = Racc_per_forager,
              visitation_rate = visitationRate,
              distance_weighted_resources = distWeigthedResources))
}
