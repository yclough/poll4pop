#' Pollinator model function
#'
#' Uses spatial data on floral and nesting quality to compute flower visitation and population dynamics for wild pollinators for a given year.
#' The function can be used for mutliple years by running a loop and using as input the output of the previous model (see example).
#' @param M_poll0 the initial population sizes to be used only if we are not in the first year
#' @param firstyear a boolean identifying whether we run the model for the first time (in which case we need to initialize the model)
#' @param firstyearfactor a factor to control the initial population sizes, defaults to 1.
#' @param bees a vector of bees identifiers corresponding to those in the parameter files
#' @param cell.size the size of a cell in the landscape matrices, in meters (default to 25)
#' @param paramList the list of parameters needed to run the model
#' @param nest the rasterstack containing the nesting qualities (dimensions should be consistent with the number of bees in the model)
#' @param floral the rasterstack containing the floral values (dimensions should be consistent with the number of bees in the model)
#' @param cutoff, the percentile at which the dispersal kernel is cut off
#' @param loc_managed, the location at which managed pollinators are placed (e.g. apiaries)
#' @return a list containing :
#' \describe{
#'   \item{M_poll}{the number of queens at the end of each period and at the beginning of the next season}
#'   \item{N_poll}{the number of foraging bees in each period}
#'   \item{flowvis}{the visitation rates per cell in each period}
#'   \item{floral and nest}{the floral and nesting values (the same as in the inputs, but stored in the results for convenience)}
#'   \item{pollres}{the resources gathered by the bees}
#'  }
#'
#' @references Haeussler J, Sahlin U, Baey C, Smith HG, Clough Y (2017) Predicting pollinator population size and pollination ecosystem service responses
#' to enhancing floral and nesting resources. Ecology and Evolution, 7: 1898-1908.\url{http://dx.doi.org/10.1002/ece3.2765}
#' Expansion to 3 seasons first used in:
#' Gardner E, et al. (2020) Reliably predicting pollinator abundance: Challenges of calibrating
#' process-based ecological models. Methods in Ecology and Evolution. \url{}
#' @examples
#'
#' nf<-computeFloralNesting(landuseMap=landuse, edgesMap=stack(grassmargins,flowermargins), unitEdges = "sqm", widthEdges=1,
#'                          landuseCodes, bees=c("GroundNestingBumblebees", "GroundNestingSolitaryBees"), num_floral=3,
#'                          florNestInfo=parameters$florNestInfo, codeEdges=c(11,21), cell.size = 10,paramList=parameters)
#' poll<-runpoll_3seasons(M_poll0 = numeric(0), firstyear=TRUE, firstyearfactor = c(1, 1),
#'             bees = c("GroundNestingBumblebees", "GroundNestingSolitaryBees"), cell.size = 10, paramList=parameters, nest=nf$nest,
#'             floral=nf$floral, cutoff = 0.99, loc_managed)



runpoll_3seasons<-function(M_poll0=numeric(0),
                  firstyear,
                  firstyearfactor=c(1,1),
		            	bees = c("Bombus spp.","Solitary spp."),
                  cell.size=25,
                  paramList,
                  nest,
                  floral,
            			cutoff=0.99,
                  loc_managed){

  require(raster)
  nr <- nrow(nest[[1]])
  nc <- ncol(nest[[1]])
  ncells <- ncell(nest[[1]])

  hab_names<-as.character(unique(paramList$lfn[,'lu']))
  bee_names<-paramList$poll_names$species_name

  bsel <- match(bees,as.character(paramList$poll_names$species_name))

  num_hab<-length(unique(paramList$lfn[,'code']))
  num_bees=length(bsel)
  num_floral<- nlayers(floral[[1]]) # takes the info from the first species

  # proportion of foraging workers
  pw <- paramList$growth[paramList$growth$Short.notation=='pw',][bsel,3]

  av<-array(paramList$av[ bsel,'best.guess'],c(num_bees))
  solitary_social<-paramList$poll_names$solitary_social[bsel]
  wild_managed<-paramList$poll_names$wild_managed[bsel]
  foraging_period <- paramList$poll_names[bsel,c('foraging_P1','foraging_P2','foraging_P3')] #Edited!!!
  growth <- array(0,c(3,2,num_bees)) #One column for q and one for w

  # Derive growth parameters for the different species #

  for(s in 1:num_bees){
    # identify subset of growth parameters for s
    selec=which(paramList$growth$species==bsel[s])
    w=1 # first period, i.e. queens
    ind<-c(grep("aq",paramList$gr[selec,'Short.notation']),grep("bq",paramList$gr[selec,'Short.notation']),grep("QEmax",paramList$gr[selec,'Short.notation']))
    growth[,w,s]<-paramList$gr[selec,][ind,3] #queens
    w=2
    ind<-c(grep("aw",paramList$gr[selec,'Short.notation']),grep("bw",paramList$gr[selec,'Short.notation']),grep("Wmax",paramList$gr[selec,'Short.notation']))
    growth[,w,s]<-paramList$gr[selec,][ind,3] #workers
  }

  # Define pollinator population size, flower visitation and pollinator resource storage arrays #

  emptyraster<-nest[[1]]
  values(emptyraster)<-NA

  # M_poll is the number of queens and N_poll the number of foragers (e.g. workers)
  M_poll <- stack(mget(rep("emptyraster",num_floral+1)))
 	 names(M_poll)=paste("period",1:(num_floral+1),sep=" ")
 	 M_poll <- mget(rep("M_poll",num_bees))
	names(M_poll)<-bees
  N_poll <- stack(mget(rep("emptyraster",num_floral)))
	  names(N_poll)=paste("period",1:(num_floral),sep=" ")
 	 N_poll <- mget(rep("N_poll",num_bees))
	names(N_poll)<-bees
  floral_periods <- stack(mget(rep("emptyraster",num_floral)))
 	 names(floral_periods)=paste("floral period",1:num_floral,sep=" ")
 	 flowvis <- mget(rep("floral_periods",num_bees))
 	 names(flowvis)=bees
  pollres<-mget(rep("floral_periods",num_bees)) # resource availability in the surrounding landscape for pollinators nesting in the focal cell
  	names(pollres)=bees

  ### ------------------------------------------###
  ###         Kernel computation                ###
  ### ------------------------------------------###
  ##print("Kernel computation ...")
  # The kernels are computed at the beginning and used throughout the simulation
  kernelFor <- list()
  kernelNest <- list()
  for (s in 1:num_bees)
  {
    kernelFor[[s]] <- kerncalc(
        paramList$distance[paramList$distance$species==bsel[s]&paramList$distance$activity=="foraging","best_guess"], cell.size, maxSize = floor(min((nr-1)/2,(nc-1)/2)), decaycut=cutoff)$decay
    kernelNest[[s]] <- kerncalc(
        paramList$distance[paramList$distance$species==bsel[s]&paramList$distance$activity=="nesting","best_guess"], cell.size, maxSize = floor(min((nr-1)/2,(nc-1)/2)), decaycut=cutoff)$decay
  }

  # Initialize population size (1st year) #
    if(firstyear)
    {
      ##print("Initiation ...")
      for(s in 1:num_bees)
      {
        if(wild_managed[s]=="wild"){
          values(M_poll[[s]][[1]]) <- firstyearfactor[s]*values(nest[[s]])*av[s]*cell.size^2/(100^2)
          #print(mean(values(M_poll[[s]][[1]])))
          values(M_poll[[s]][[1]]) <- rpois(length(values(M_poll[[s]][[1]])), values(M_poll[[s]][[1]]))
        }
        if(wild_managed[s]=="managed"){
          values(M_poll[[s]][[1]]) <- values(loc_managed)
        }
      }
    }else{
      for(s in 1:num_bees)
      {
        values(M_poll[[s]][[1]])<-values(M_poll0[[s]][[(num_floral+1)]]) #Edited!!!
      }
    }

  for(s in 1:num_bees)
  {
  	##print(bees[s]);flush.console()
    ##print(paste("nb females start : ",  round(sum(values(M_poll[[s]][[1]])))  ));flush.console()

    ## -----------------  floral period 1  ----------------- ##
    v=1
    ##print("Floral period 1 ...");flush.console()
    ## M_poll (number of queens) aggregate resources (pollres v=1) --> VR1

    # for all species foraging in period 1
    if(foraging_period[s,1]!="none")
    {
      foraging <- latfordisp(N=raster::as.matrix(M_poll[[s]][[1]]),
                           alpha = raster::as.matrix(floral[[s]][[v]]),
                           decay = kernelFor[[s]])
      values(flowvis[[s]][[v]]) <- as.vector(t(foraging$visitation_rate))  			 # flower visitation in the cells surrounding the focal cell
      values(pollres[[s]][[v]]) <- as.vector(t(foraging$Racc_per_forager))	 # aggregated resources experienced by the pollinator nesting in the focal cell
    }

    # for social species with queens foraging in period 1
    if(solitary_social[s]=="social" & foraging_period[s,v]=="q")
    {
      ## ------ # Population growth WORKERS # ------- ##
      ## N_poll, number of workers produced, depending on the number of queens and resources in period 1
      values(N_poll[[s]][[1]]) <- values(M_poll[[s]][[1]]) * growth.func(R = values(pollres[[s]][[v]]),
                                                                         a = growth[1,2,s],
                                                                         b = growth[2,2,s],
                                                                         max = growth[3,2,s])
      ##print("growth func ran ok");flush.console()
  	}

    # for solitary species foraging in period 1
    if(solitary_social[s]=="solitary" & foraging_period[s,v]=="q")
    {
     	values(M_poll[[s]][[2]]) <- values(M_poll[[s]][[1]]) * growth.func(R = values(pollres[[s]][[v]]),
     	                                                                   a = growth[1,1,s],
     	                                                                   b = growth[2,1,s],
     	                                                                   max = growth[3,1,s])
     	M_pollend = M_poll[[s]][[2]] #Added!!!
  	}


    ## -----------------  floral period 2  ----------------- ##
    v=2
    ##print("Floral period 2 ...");flush.console()

    if(solitary_social[s]=="social" & foraging_period[s,v]=="w")
 	  {
   	  ## N_poll (workers) aggregate resources (pollres v=2) --> VR2
  	  ## Only a proportion pw of the workers is actually foraging
  	  foraging <- latfordisp(N = pw[s]*raster::as.matrix(N_poll[[s]][[1]]), alpha = raster::as.matrix(floral[[s]][[v]]), decay = kernelFor[[s]])

  	  values(flowvis[[s]][[v]]) <- as.vector(t(foraging$visitation_rate))		 # flower visitation in the cells surrounding the focal cell
  	  values(pollres[[s]][[v]]) <- as.vector(t(foraging$Racc_per_forager))	 #

   	  ## ------ # Population growth QUEENS # ------- ##
   	  ## M_poll produced at the end of the season depending on M_poll and resources in period 2
   	  # growth[3,2,s] is the maximal number of workers
   	  # N_poll is the number of foraging bees
      #values(M_poll[[s]][[2]]) <- values(M_poll[[s]][[1]]) * growth.func(R = values(pollres[[s]][[v]])*pw[s]*values(N_poll[[s]][[1]])/values(M_poll[[s]][[1]]), a = growth[1,1,s], b = growth[2,1,s], max = growth[3,1,s]) #Replaced!!!
  	  M_poll[[s]][[2]] = M_poll[[s]][[1]] #Added!!!
  	  values(N_poll[[s]][[2]]) <- values(N_poll[[s]][[1]])+values(M_poll[[s]][[1]]) * growth.func(R = values(pollres[[s]][[v]])*pw[s]*values(N_poll[[s]][[1]])/values(M_poll[[s]][[1]]),
                                                                         a = growth[1,2,s],
                                                                         b = growth[2,2,s],
                                                                         max = growth[3,2,s]) #Added!!!
   	 }

    if(solitary_social[s]=="solitary" & foraging_period[s,v]=="q")
    {
      foraging <- latfordisp(N=raster::as.matrix(M_poll[[s]][[1]]),
                             alpha = raster::as.matrix(floral[[s]][[v]]),
                             decay = kernelFor[[s]])

      values(flowvis[[s]][[v]]) <- as.vector(t(foraging$visitation_rate))  			 # flower visitation in the cells surrounding the focal cell
      values(pollres[[s]][[v]]) <- as.vector(t(foraging$Racc_per_forager))	 # aggregated resources experienced by the pollinator nesting in the focal cell

     	values(M_poll[[s]][[2]]) <- values(M_poll[[s]][[1]]) #Added!!!
     	values(M_poll[[s]][[3]]) <- values(M_poll[[s]][[2]]) * growth.func(R = values(pollres[[s]][[v]]),
                                               a = growth[1,1,s],
                                               b = growth[2,1,s],
                                               max = growth[3,1,s]) #Edited!!!
      values(M_poll[[s]][[1]]) <- NA #Added!!!
     	M_pollend = M_poll[[s]][[3]] #Added!!!

    }

  ## -----------------  floral period 3  ----------------- ##Added!!!
    v=3
    ##print("Floral period 3 ...");flush.console()

    if(solitary_social[s]=="social" & foraging_period[s,v]=="w")
 	  {
   	  ## N_poll (workers) aggregate resources (pollres v=2) --> VR2
  	  ## Only a proportion pw of the workers is actually foraging
  	  foraging <- latfordisp(N = pw[s]*raster::as.matrix(N_poll[[s]][[2]]), alpha = raster::as.matrix(floral[[s]][[v]]), decay = kernelFor[[s]])

  	  values(flowvis[[s]][[v]]) <- as.vector(t(foraging$visitation_rate))		 # flower visitation in the cells surrounding the focal cell
  	  values(pollres[[s]][[v]]) <- as.vector(t(foraging$Racc_per_forager))	 #

   	  ## ------ # Population growth QUEENS # ------- ##
   	  ## M_poll produced at the end of the season depending on M_poll and resources in period 2
   	  # growth[3,2,s] is the maximal number of workers
   	  # N_poll is the number of foraging bees
      values(M_poll[[s]][[3]]) <- values(M_poll[[s]][[1]]) * growth.func(R = values(pollres[[s]][[v]])*pw[s]*values(N_poll[[s]][[2]])/values(M_poll[[s]][[1]]), a = growth[1,1,s], b = growth[2,1,s], max = growth[3,1,s]) #Edited!!!
      M_pollend = M_poll[[s]][[3]] #Added!!!
      ##print("growth func ran ok");flush.console()
  	}

    if(solitary_social[s]=="solitary" & foraging_period[s,v]=="q")
    {
      print("Don't allow Season 3 only active solitary bees") #Added!!!

  	}
    ##print(paste("nb new females : ",  round(sum(values(M_pollend)))  ));flush.console() #Edited!!!


    ## M_poll at the beginning of the next season depending on winter survival of M_poll at the end of the season
    # non-honeybees disperse to their new nesting sites #
		if(wild_managed[s]=="wild")
		{
      poll.out.start <- latfordisp(N = raster::as.matrix(M_pollend),
                                   alpha = raster::as.matrix(nest[[s]]),
                                   decay = kernelNest[[s]]) #Edited!!!

      values(M_poll[[s]][[(num_floral+1)]]) <- pmin(as.vector(t(poll.out.start$visitation_rate)),(av[s]*cell.size^2/10000)*values(nest[[s]])) #Edited!!!
	    ##print(paste("nb females surv : ",  round(sum(values(M_poll[[s]][[(num_floral+1)]])))  ));flush.console() #Edited!!!
		}

  } #end species

  return(list(M_poll=M_poll,
              N_poll=N_poll,
              flowvis=flowvis,
              floral=floral,
              nest=nest,
              pollres=pollres))
}

