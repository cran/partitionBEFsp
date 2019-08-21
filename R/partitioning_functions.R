#' Calculate change in relative yield
#'
#' calculates change in relative yield, DRY, comparing observed relative yield to the expected yield 1/Q
#' @param P biomass of species grown in polyculture
#' @param M biomass of species grow in monoculture - note, must include the same species as P, listed in the same order
#' @param Q number of species in the community -defaults to length(M), but note that if you are calculating DRY for a large community of Q species of which only N are observed, you should set Q=Q, rather than Q=N.
#' @return a list of changes in relative yields
#' @export
#' @examples
#' # Please see package help file (?partitionBEFsp) for examples.

calculate_DRY<-function(P, M, Q=length(M)) {
  DRY<-P/M - 1/Q
  return(DRY)
}

#' Calculate classic partition
#'
#' calculates the classic selection and complementarity effects, sensu Loreau and Hector 2001
#' @param DRY change in relative yield, as calculated by the calculate_DRY function
#' @param M monoculture biomass
#' @param N number of species in the sample - defaults to length(M)
#' @param Q number of species in the full population - defaults to N - only required if uncorrected_cov="COMP"
#' @param uncorrected_cov A character, which can be TRUE, FALSE, or COMP. Tells whether to use the standard sample-size corrected covariance function (FALSE), or
#' a covariance function that is not corrected for sample size (TRUE), or
#' a "compromise" function that resembles the standard function for N < Q, and that resembles the non-corrected function for N ~ Q
#' If TRUE, then SS + CS = YO - YE, sensu Loreau and Hector 2001
#' defaults to FALSE
#' note - we do not recommend setting this to TRUE or "COMP", unless you require SS+CS=YO-YE
#' @return a list with elements S (the selection effect) and C (the complementarity effect)
#' @import stats
#' @export
#' @examples
# Plese see package help file (?partitionBEFsp) for examples.

classic_partition<-function(DRY, M, N=length(M), Q=N, uncorrected_cov=FALSE) {
  if(uncorrected_cov==TRUE) {
    S<-sum((mean(DRY)-DRY)*(mean(M)-M))
  } else if(uncorrected_cov==FALSE) {
    S<-N*cov(DRY, M)
  } else if(uncorrected_cov=="COMP") {
    S<-sum((mean(DRY)-DRY)*(mean(M)-M))*N/(N-(Q-N)/Q)
  }
  C<-N*mean(DRY)*mean(M)

  return(list(S=S, C=C))
}

#' Calculate population-level partition
#'
#' takes a random but incomplete sample of species of size N from a larger community Q,
#' and estiamtes population-level selection and complementarity effects
#' @param DRY change in relative yield, as calculated by the calculate_DRY function
#' @param M monoculture biomass
#' @param N number of species in the sample of the full community (i.e. the "sample") - defaults to length(M)
#' @param Q total number of species in the full community (i.e. the "population")
#' @param smallQ_correction tells whether to apply the correction for small Q, as shown in Eq. 3c in the main text - defaults to TRUE
#' @param uncorrected_cov A character, which can be TRUE, FALSE, or COMP. Tells whether to use the standard sample-size corrected covariance function (FALSE), or
#' @param nboot Number of bootstrap iterations to run for estimating confidence intervals for selection and complementarity effects. Defaults to NA - i.e. no bootstrapping.
#' a covariance function that is not corrected for sample size (TRUE), or
#' a "compromise" function that resembles the standard function for N < Q, and that resembles the non-corrected function for N ~ Q
#' If TRUE, then SS + CS = YO - YE, sensu Loreau and Hector 2001
#' defaults to FALSE
#' note - we do not recommend setting this to TRUE or "COMP", unless you require SS+CS=YO-YE
#' @return a list with elements SS (the sample-level selection effect), CS (the sample-level complementarity effect),
#' SP (the population-level selection effect), CP (the population-level complementarity effect),
#' and confint, which is a list that includes summary data and the full bootstrapped for estimates of the confidence intervals (if nboot != NA)
#' @import stats
#' @export
#' @examples
#' # Please see package help file (?partitionBEFsp) for examples.

sample_to_population_partition<-function(DRY, M, N=length(M), Q, smallQ_correction=TRUE, uncorrected_cov=FALSE, nboot=NA) {
  #calculate sample-level statistics
  if(uncorrected_cov==TRUE) {
    SS<-sum((mean(DRY)-DRY)*(mean(M)-M))
  } else if(uncorrected_cov==FALSE) {
    SS<-N*cov(DRY, M)
  } else if(uncorrected_cov=="COMP") {
    SS<-sum((mean(DRY)-DRY)*(mean(M)-M))*N/(N-(Q-N)/Q)
  }
  CS<-N*mean(DRY)*mean(M)

  #calculate population-level statistics
  SP<-Q/N*SS
  if(smallQ_correction) {
    CP<-Q/N*(CS-(Q-N)/Q*(1/N*SS))
  } else {
    CP<-Q/N*(CS-(1/N*SS))
  }

  #bootstrap confidence intervals
  if(!is.na(nboot)) {
    bootdat_full<-matrix(nrow=nboot, ncol=4)
    colnames(bootdat_full)<-c("SS", "CS", "SP", "CP")
    
    for(j in 1:nboot) {
      smp<-sample(N, replace=TRUE)
      cov(DRY[smp], M[smp])
      
      if(uncorrected_cov==TRUE) {
        bootdat_full[j,"SS"]<-sum((mean(DRY[smp])-DRY[smp])*(mean(M[smp])-M[smp]))
      } else if(uncorrected_cov==FALSE) {
        bootdat_full[j,"SS"]<-N*cov(DRY[smp], M[smp])
      } else if(uncorrected_cov=="COMP") {
        bootdat_full[j,"SS"]<-sum((mean(DRY[smp])-DRY[smp])*(mean(M[smp])-M[smp]))*N/(N-(Q-N)/Q)
      }
      bootdat_full[j,"CS"]<-N*mean(DRY[smp])*mean(M[smp])
      
      
      bootdat_full[j,"SP"]<-Q/N*bootdat_full[j,"SS"]
      if(smallQ_correction) {
        bootdat_full[j,"CP"]<-Q/N*(bootdat_full[j,"CS"]-(Q-N)/Q*(1/N*bootdat_full[j,"SS"]))
      } else {
        bootdat_full[j,"CP"]<-Q/N*(bootdat_full[j,"CS"]-(1/N*bootdat_full[j,"SS"]))
      }
    }
    
    bootdat_summary<-apply(bootdat_full, 2, function(x) cbind(mean(x), sd(x), mean(x>0)))
    bootdat_summary[3,][bootdat_summary[3,]>0.5]<-1-bootdat_summary[3,][bootdat_summary[3,]>0.5]
    bootdat_summary[3,][bootdat_summary[3,]==0]<-1/nboot
    rownames(bootdat_summary)<-c("mean", "sd", "p-value")
    
  } else {
    bootdat_full<-NA; bootdat_summary<-NA
  }
  
  return(list(SS=SS, CS=CS, SP=SP, CP=CP, confint=list(bootdat_summary=bootdat_summary, bootdat_full=bootdat_full)))
}

#' Calculate sample-level partition
#'
#' takes a complete sample of all Q species in a community,
#' and estimates sample-level selection and complementarity effects expected
#' from a subset of N species drawn randomly from that community
#' @param DRY change in relative yield, as calculated by the calculate_DRY function
#' @param M monoculture biomass
#' @param N number of species in the sample of the full community (i.e. the "sample") - defaults to length(M)
#' @param Q total number of species in the full community (i.e. the "population")
#' @param smallQ_correction tells whether to apply the correction for small Q, as shown in Eq. 3c in the main text - defaults to TRUE
#' @param uncorrected_cov A character, which can be TRUE, FALSE, or COMP. Tells whether to use the standard sample-size corrected covariance function (FALSE), or
#' a covariance function that is not corrected for sample size (TRUE), or
#' a "compromise" function that resembles the standard function for N < Q, and that resembles the non-corrected function for N ~ Q
#' If TRUE, then SS + CS = YO - YE, sensu Loreau and Hector 2001
#' defaults to FALSE
#' note - we do not recommend setting this to TRUE or "COMP", unless you require SS+CS=YO-YE
#' @return a list with elements SS (the sample-level selection effect), CS (the sample-level complementarity effect),
#' SP (the population-level selection effect), and CP (the population-level complementarity effect),
#' @import stats
#' @export
#' @examples
#' # Please see package help file (?partitionBEFsp) for examples.

population_to_sample_partition<-function(DRY, M, N, Q=length(M), smallQ_correction=TRUE, uncorrected_cov=FALSE) {
  #calculate sample-level statistics
  if(uncorrected_cov==TRUE) {
    SP<-sum((mean(DRY)-DRY)*(mean(M)-M))
  } else if(uncorrected_cov==FALSE) {
    SP<-Q*cov(DRY, M)
  }
  CP<-Q*mean(DRY)*mean(M)

  #calculate population-level statistics
  SS<-N/Q*SP
  if(smallQ_correction) {
    CS<-N/Q*CP+(Q-N)/Q*(1/N*SS)
  } else {
    CS<-N/Q*CP+(1/N*SS)
  }

  return(list(SS=SS, CS=CS, SP=SP, CP=CP))
  #SS and SP are sample and population-level selection effects
  #CS and CP are sample and population-level complementarity effects
}

