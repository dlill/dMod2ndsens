#' Difference between two profiles
#'
#' @param param The parameter to be examined
#' @param profile1,profile2 The two different profiles
#'
#'@export
profilediff <- function(param, profile1, profile2) {
  prof1sub <- subset.data.frame(profile1, whichPar == c("a","b")[param])
  prof2sub <- subset.data.frame(profile2, whichPar == c("a","b")[param])
  
  indices <- match(unlist(round(prof1sub[,c("a","b")[param]], 2)), unlist(round(prof2sub[,c("a","b")[param]], 2)))
  
  prof1vals <- prof2sub[!is.na(indices),]
  prof2vals <- prof2sub[na.omit(indices),]
  
  profvals[,"value"] <- prof1vals[,"value"]- prof2vals[,"value"]
  return(profvals)
}



#' LL along an approximated geodesic curve
#' This function evaluates the objective function along the geodesics that are returned by geodesic(direction, tau)
#' The function geodesic(direction, tau) has to exist already.
#' @param direction Vector of the initial velocity of the geodesic
#' @param tau Vector with stepwidths of the geodesic.
#' @export
LLgeod <- function(direction, tau) {
  
  mygeodesic <- geodesic(direction, tau)
  
  # whichdir <- which(directions == 1) #not yet suitable for arbitrary directions!
  # mydirparvalues <- rbind(mygeodesic)[,whichdir]
  mygeodesicpars <- lapply(1:ncol(mygeodesic), function(i) mygeodesic[,i])
  mygeodesic <- lapply(1:nrow(mygeodesic), function(i) mygeodesic[i,])
  
  myprofilevalues <- lapply(mygeodesic, function(p) {
    myobj <- obj(p, deriv = 0)
    out <- myobj$value
    attr(out, "valueData") <- attr(myobj, "valueData")
    attr(out, "valuePrior") <- attr(myobj, "valuePrior")
    return(out)
  })
  mydatavalues <- unlist(lapply(myprofilevalues, function(i) attr(i,"valueData")))
  mypriorvalues <- unlist(lapply(myprofilevalues, function(i) attr(i,"valuePrior")))
  return(list(as.vector(tau), unlist(myprofilevalues), mydatavalues, mypriorvalues, mygeodesicpars))
}


#' Brute force straight "profiles"
#' This function evaluates the objective function along straight lines that are specified by the initial vecotr bestfitpars and direction
#' 
#' @param direction Vector of the initial velocity of the line
#' @param tau Vector with stepwidths of the line.
#' @param bestfitpars Vector with the parameter values of the best fit.
#'@export
LLgeod <- function(direction, tau, bestfitpars) {
  
  mypath <-outer(tau, direction)
  mypath <- t(apply(mypath,1, function(t) {t+bestfitpars}))
  
  mypathpars <- lapply(1:ncol(mypath), function(i) mypath[,i])
  
  mypath <- lapply(1:nrow(mypath), function(i) mypath[i,])
  myprofilevalues <- lapply(mypath, function(p) {
    myobj <- obj(p, deriv = 0)
    out <- myobj$value
    attr(out, "valueData") <- attr(myobj, "valueData")
    attr(out, "valuePrior") <- attr(myobj, "valuePrior")
    return(out)
  })
  mydatavalues <- unlist(lapply(myprofilevalues, function(i) attr(i,"valueData")))
  mypriorvalues <- unlist(lapply(myprofilevalues, function(i) attr(i,"valuePrior")))
  return(list(tau, unlist(myprofilevalues), mydatavalues, mypriorvalues, mypathpars))  
}


#'Plot "profiles" and the behaviour of the parameters
#'Pick a direction of your directionslist by direcitonindex and have a look on the value of the objective function and the behaviour of
#'the parameters wrt tau, the parameter of your curve (geodesic or line)
#'@param directionindex In dMod description.R the profilelist are called via a list of directionvectors. Directionindex is the index of the direction you want to look at.
#'@export
plotprofpars <- function(profilelist, directionindex) {
  par(mfrow = c(3,3))
  plot(profilelist[[directionindex]][[1]], profilelist[[directionindex]][[2]])
  lapply(1:length(pouter), function(i) plot(profilelist[[directionindex]][[1]], profilelist[[directionindex]][[5]][[i]], ylab = names(pouter[i])))
}

#' Compute the difference of the objective functions of proflist1 and proflist2
#' proflist1-Profile haben kleineres chi^2 als die die von Proflist2, wenn der Output positiv ist.
#' It can also be used to compare the profiles parameterwise, but then one needs to change the indices of the list in the function
#' ie in indices <- ... proflist2[[i]][[1]]... to ...proflist2[[i]][[5]][[1]]... etc.
#' @param proflist1/proflist2 List of "profiles" of the functions above. The first list index comes from calling these functions with lapply(directionlist)
#' @export
simpleLLprofiledifference <- function(proflist1, proflist2){
  lapply(1:length(pouter), function(i) { #simple because it doesn't involve the actual profiles, but just how the cost function behaves when going into the direction without repotimizing
  indices <- match(round(proflist2[[i]][[1]], 2), round(proflist1[[i]][[1]], 2))
  geodesicindices <- na.omit(indices)
  nongeodesicindices <- !is.na(indices)
  xval <- proflist1[[i]][[1]][geodesicindices]
  difference <- proflist2[[i]][[2]][nongeodesicindices] - proflist1[[i]][[2]][geodesicindices] #geodäten sind besser, wenn dieser wert positiv ist.
  return(list(xval, difference))
})
}

#' Analyse the Profiles statistically
#' Find the values for which the statistical significance level is reached. I also assumed a Likelihood ratio test
#' with one parameter less which corresponds to walking in one direction...
#' @param profilelist A profilelist that is generated by calling the functions above (bruteforceprofiles or LLgeod or elegantprofiles) in an lapply(directionlist, function)
#' @param directionindex the index of the direction in the direcitonlist.
#' @export
analyseLLProfiles <- function(profilelist, directionindex) {
  objbf <- profilelist[[directionindex]][[2]][which(tau==0)]
  
  parameteranalysis <- lapply(1:length(directions), function(par) {
    confidenceindices <- which(profilelist[[directionindex]][[2]] <= (objbf+1.82))
    parname <- names(bestfit)[par]
    par0 <- bestfit[par]
    maxpar <- max(profilelist[[directionindex]][[5]][[par]][confidenceindices])
    minpar <- min(profilelist[[directionindex]][[5]][[par]][confidenceindices])
    return(data.frame(directionindex=directionindex, parname=parname, bestfit = par0, maxpar =maxpar, minpar =minpar, range = (maxpar-minpar)))
  })
  return(do.call(rbind,parameteranalysis))
}


#' Compute and compare the different ranges after which the objective function has a value higher than the level of significance
#' 
#'@export
rangedifference <- function(profilelist1, profilelist2, directionindex) {
  ana1 <- analyseLLProfiles(profilelist1, directionindex)
  range1 <- ana1[,"range"]
  range2 <- analyseLLProfiles(profilelist2, directionindex)[,"range"]
  range <- range1-range2
  return(data.frame(ana1[,c(1,2,3)], range1=range1, range2=range2, rangedifference =range))
}


#' #' Have a new objective function designed, which behaves quadratic. To be compared to the real objective function.
#' #' If the theory is applicable the output of this function should be the same as the one of geodesicprofiles.
#' #' This would be the behaviour if the nonlinear parts were negligible. 
#' #' 
#' #' 
#' #' A first test shows that it is not enough to test the nonlinear parts only at the point of the best fit. 
#' #' I assumed that they vary slowly but of course this doesn't have to be the case. 
#' #' An idea would be to test the nonlinear parts at least along the stiff directions of the metric.
#' #' Another idea to ponder is playing with the data, reducing the space normal to the tangent space.
#' #' But in my model there is one direction which is non identifiable anyways, so maybe the model isn't suitable for a first test if it works with general models?
#' #' 
#' #' Another thought that arises is that the profiles don't look perfectly parabolic so there is the problem 
#' #' that the geodesics 1. aren't straight lines in residual space and 2. their direction is not necessarily 
#' #' perpendicular to the hyperspheres centered around the origin of residual space.
#' #' 
#' #' 
#' #' @param residual as usual
#' #' @param direction as usual
#' #' @param tau as usual
#' #' @export
#' elegantprofiles <- function(bestfit, residual, direction, tau) {
#'   
#'   #residuals and objective function at bestfit
#'   myresidual <- residual[,"residual"]
#'   myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
#'   mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
#'   
#'   bestfitobj <- obj(bestfit, deriv = 0)$value
#'   
#'   # construct coordinates of paths in RNC (ie straight lines)
#'   mypath <-outer(tau, direction)
#'   mypath <- t(apply(mypath,1, function(t) {t}))
#'   mypath <- lapply(1:nrow(mypath), function(i) mypath[i,])
#'   
#'   # compute the coefficients of the metric going into a certain direcction
#'   metric <- t(myderiv)%*%myderiv
#'   metric <- metric %*% direction
#'   metric <- direction %*% metric
#'   
#'   # normalize it so that at tau==1 one has walked a unit distance
#'   vabs <- direction %*% direction
#'   metric <- metric/vabs
#'   
#'   # compute the part of the second direcitonal derivative that is not constant
#'   u <- svd(myderiv)$u
#'   PN<-diag(nrow = nrow(u)) - u%*%t(u)
#'   mysderivarray <- array(mysderiv, dim=c(nrow(mysderiv), dim(myderiv)[2], dim(myderiv)[2])) #first index switches between the rm the second indices, the second and third go through the parameter indices
#'   mysderivarray <- apply(mysderivarray, 2, "%*%", direction) # multiply third index with initial velocity
#'   mysderivarray <- PN%*%mysderivarray%*%direction # multiply second index with velocity and project it with PN
#'   nonquadcoeffs <- myresidual%*%mysderivarray
#'   
#'   # add the metric part and the nonquadratic part to get the value of the second directional derivative at the point where residual was comnputed
#'   metric <- metric+nonquadcoeffs #put it in "metric" to save variables
#'   
#'   # compute the values of the profiles
#'   myprofilevalues <- lapply(mypath, function(p) {
#'     out <- bestfitobj+metric*((p)%*%(p))
#'   })
#'   
#'   return(list(as.vector(tau), unlist(myprofilevalues), mydatavalues=0, mypriorvalues=0, mypath))
#' }
#' 
#' 









#' As Proflist (funktioniert nicht)

asproflist <- function(profilelist, bestfit, tau) {
  value <- c(unlist(lapply(1:length(profilelist), function(x) profilelist[[x]][[2]])))
  myparcols <- lapply(1:length(profilelist), function(x) rep(profilelist[[x]][[1]], length(profilelist)))
  myparcols <- matrix(unlist(myparcols), ncol=length(myparcols))
  colnames(myparcols) <- names(bestfit)
  valueData <- c(unlist(lapply(1:length(profilelist), function(x) profilelist[[x]][[3]])))
  constraint <- c(unlist(lapply(1:length(profilelist), function(x) profilelist[[x]][[4]])))
  valuePrior <- abs(constraint)
  whichPar <- c(unlist(lapply(names(bestfit), function(x) rep(x, length(tau)))))
  stepsize <- rep(0,  length(tau)) #zu faul
  gamma <-rep(1,  length(tau)) #weiß nicht was das ist
  mydataframe <- parframe(data.frame(value, constraint, stepsize, gamma, whichPar, valueData, valuePrior, myparcols, stringsAsFactors = FALSE))
  attr(mydataframe, "parameters") <- names(bestfit)
  attr(mydataframe, "metanames") <- c("value", "constraint", "stepsize", "gamma", "whichPar")
  attr(mydataframe, "obj.attributes")= c("valueData","valuePrior")
  return(mydataframe)
}