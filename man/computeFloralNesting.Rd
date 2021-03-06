% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/computeFloralNesting.R
\name{computeFloralNesting}
\alias{computeFloralNesting}
\title{Compute floral and nesting resource values for flower-visiting insects}
\usage{
computeFloralNesting(
  landuseMap,
  edgesMap,
  unitEdges = "m",
  widthEdges,
  landuseCodes,
  bees,
  num_floral,
  florNestInfo,
  codeEdges,
  cell.size = 25,
  paramList
)
}
\arguments{
\item{landuseMap}{the map of the landscape containing the landuse categories, given as a raster}

\item{edgesMap}{the map containing the length of field edges in each cell, given as a raster or a rasterstack in the case of multiple edge types}

\item{landuseCodes}{the matching between landuse codes and landuse category names}

\item{bees}{the vector of bee species}

\item{num_floral}{the number of floral periods}

\item{florNestInfo}{a list containing the floral and nesting information, with 4 elements: 1) the floral coverage, 2) bumblee bee info, i.e. floral value and attractiveness, 3) honeybee info and 4) solitary bee info}

\item{codeEdges}{a vector containing the landuse codes for the edge types (e.g. grassy field edge or sown wildflower strip)}

\item{cell.size}{a number in meters corresponding to the width and breadth of the input raster cell}

\item{paramList}{the list of parameters needed to run the model}
}
\value{
a list containing two arrays: one of dimension nb of bees containing the nesting qualities, and one of dimension nb of bees * nb of periods containing the floral values for each species and each period
}
\description{
Function that computes the floral and nesting arrays according to the bees that are accounted for
}
\examples{

nf<-computeFloralNesting(landuseMap=landuse, edgesMap=stack(grassmargins,flowermargins), unitEdges = "sqm", widthEdges=1,
                         landuseCodes, bees=c("GroundNestingBumblebees", "GroundNestingSolitaryBees"), num_floral=3,
                         florNestInfo=parameters$florNestInfo, codeEdges=c(11,21), cell.size = 10,paramList=parameters)
poll<-runpoll_3seasons(M_poll0 = numeric(0), firstyear=TRUE, firstyearfactor = c(1, 1),
            bees = c("GroundNestingBumblebees", "GroundNestingSolitaryBees"), cell.size = 10, paramList=parameters, nest=nf$nest,
            floral=nf$floral, cutoff = 0.99, loc_managed)
}
\references{
Haeussler J, Sahlin U, Baey C, Smith HG, Clough Y (2017) Predicting pollinator population size and pollination ecosystem service responses
to enhancing floral and nesting resources. Ecology and Evolution, 7: 1898-1908.\url{http://dx.doi.org/10.1002/ece3.2765}
Expansion to 3 seasons first used in:
Gardner E, et al. (2020) Reliably predicting pollinator abundance: Challenges of calibrating
process-based ecological models. Methods in Ecology and Evolution. \url{}
}
