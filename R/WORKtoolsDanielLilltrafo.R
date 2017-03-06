#' Theta(p)

#' Coordinates Theta(p)
#' @details Returns a function of p, the riemann normal coordinates which computes the coordinates of theta in dependence of p
#' 
#' @param bestfit The point "theta_0" around which the rnc coordinates p are constructed
#' @param residual The residual at the point "bestfit"
#' @param christoffelsymbols The christoffel symbols at the point "bestfit"
#' 
#' @param p Argument of the returned function. A matrix with coordinate parameters, each row corresponds to one set of rnc parameters p
#' 
#' @export
Thetafun <- function(bestfit, residual, christoffelsymbols = christoffels(residual)) {
function(p) {
      th <- apply(christoffelsymbols, 3, "%*%", p) #first step in the evaluation of the quadratic term, called th to save variables
      th <- as.vector(bestfit+p-0.5*th%*%p)
  }
}


#' This function returns a function for the coordinates of a geodesic in Theta (parameter-) space.
#' @param residual The residual at the bestfit
#' @param bestfit The parameter vector with the best bestfit parameters
#' @param thirdpower Include the third order term or not. (Currently commented out.)
#' 
#' @param direction Argument of the returned function.The vector of the initial direction in both Theta and p space
#' @param tau Argument of the returned function. The stepwidth along the geodesic.
#' @export
Thetageofun <- function(bestfit, residual, christoffelsymbols = christoffels(residual)){
  
  # get 1st and 2nd derivatives 
  myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
  mysderiv <- data.matrix(attr(residual, "sderiv")[,-c(1,2)])
  
  # get parameter names to whose respect the derivatives were taken (for debugging reasons while coding)
  mypars <- colnames(myderiv)
  # myspars <- colnames(mysderiv)
  
  #get best bestfit parameters
  mybestpars <- bestfit#$argument #macht nen error, obwohls früher immer funktioniert hat
  

  geodesic <- NULL
  # if(thirdpower==FALSE ) {
    
    geodesic <- function(direction, tau) {
      mycoordinates <- t(matrix(unlist(lapply(tau, function(mytau) {
        #evaluate the quadratic coefficients.
        mycoeffs <- apply(christoffelsymbols, 3, "%*%", direction)
        mycoeffs <- mycoeffs %*% direction
        # sum up all the coordinates
        coordinates <- as.vector(mybestpars) + as.vector(direction * mytau) - as.vector(0.5 * mycoeffs * (mytau ^ 2))
      })), ncol = length(tau)))

      colnames(mycoordinates) <- mypars
      rownames(mycoordinates) <-  as.character(tau)
      attr(mycoordinates, "direction") <- direction

      return(mycoordinates)
    }
    
  # } else {
  #   mychristoffelderiv <- christoffelderiv(residual)
  #   
  #   geodesic <- function(direction, tau, thirdpower = FALSE) {
  #     mycoordinates <- t(matrix(unlist(lapply(tau, function(mytau) {
  #       #evaluate the quadratic coefficients.
  #       mycoeffs <- apply(mychristoffels, 3, "%*%", direction)
  #       mycoeffs <- mycoeffs %*% direction
  #       #evaluate the cubic coefficients.
  #       mycubics <- array(apply(mychristoffels, c(2,3), function(chrisupper) {
  #         apply(mychristoffels, c(1,2), function(chrislower) sum(chrislower*chrisupper))
  #       }), dim = rep(ncol(myderiv), 4))
  #       mycubics <- 2*mycubics - mychristoffelderiv
  #       mycubics <- apply(mycubics, c(2,3), "%*%", direction)
  #       mycubics <- apply(mycubics, 2, "%*%", direction)
  #       mycubics <- mycubics %*% direction
  #       
  #       # sum up all the coordinates
  #       coordinates <- as.vector(mybestpars) + as.vector(direction * mytau) - as.vector(0.5 * mycoeffs * (mytau ^ 2) + as.vector(1 /
  #                                                                                                                                  6 * (mycubics) * (mytau ^ 3)))
  #     })), ncol = length(tau)))
  #     mycoordinates <- rbind(direction, mycoordinates)
  #     mycoordinates <- cbind(c(0, tau), mycoordinates)
  #     colnames(mycoordinates) <- c("stepwidth", mypars)
  #     rownames(mycoordinates) <- c("direction", as.character(tau))
  #     return(mycoordinates)
  #   }  
    return(geodesic)
}




#' Jacobian dTheta/dp
#' @details Returns the Jacobian (dTheta^mu/dp^nu)^mu_nu at the point p
#'
#' @param bestfit The point "Theta_0" around which the rnc coordinates p are constructed
#' @param residual The residual at the point "bestfit"
#' @param christoffelsymbols The christoffel symbols at the point "bestfit"
#' 
#' @param p Argument of the returned function. A matrix with coordinate parameters, each row corresponds to one set of rnc parameters p
#' 
#' @export
JTpfun  <- function(bestfit, residual, christoffelsymbols = christoffels(residual)) {
 function(p){
    diag(1, length(bestfit)) - apply(christoffelsymbols, 3, "%*%", p )
 }
}


#' Metric in parameter (Theta)-space
#' @param residual the residual evaluated at the point of bestfit
#' 
#' @export
Tmetric <- function(residual) {
    myderiv <- data.matrix(attr(residual, "deriv")[,-c(1,2)])
    myderiv <- t(myderiv)%*%myderiv # actually the metric, called myderiv to save variables
}











#' p(Theta)



#' Coordinates p(Theta)
#' @details Returns a function of Theta, to give the value of the riemann normal coordinates, which are constructed around the point Theta_0
#' 
#' @param bestfit The point "Theta_0" around which the rnc coordinates p are constructed
#' @param residual The residual at the point "bestfit"
#' @param christoffelsymbols The christoffel symbols at the point "bestfit"
#' 
#' @param Theta Parameter of the returned function. A matrix with coordinate parameters, each row corresponds to one set of rnc parameters p
#' 
#' @export
pfun <- function(bestfit, residual, christoffelsymbols = christoffels(residual)) {
  
  # cubiccoeffs <- array(apply(christoffelsymbols, c(2,3), function(chrisupper) {
    # apply(christoffelsymbols, c(1,2), function(chrislower) sum(chrislower*chrisupper))
  # }), dim = rep(length(bestfit), 4))
  # cubiccoeffs <- cubiccoeffs + christoffelderiv(residual)
  
  p <- function(Theta) {

    # Theta <- matrix(Theta, ncol=1)
    # bestfit <- matrix(bestfit, ncol= 1)

      # t((Theta[i,]-bestfit)+0.5*apply(christoffelsymbols, 3, "%*%", (Theta[i,]-bestfit))%*%(Theta[i,]-bestfit))
      as.vector((Theta-bestfit))+as.vector(0.5*apply(christoffelsymbols, 3, "%*%", (Theta-bestfit))%*%(Theta-bestfit) ) #+1/6*apply(apply(cubiccoeffs, c(2,3), "%*%", (Theta[i,]-bestfit)), 2, "%*%", (Theta[i,]-bestfit))%*%(Theta[i,]-bestfit))

  }
  return(p)
}

#' Coordinates p(Theta)
#' @details Returns a function of Theta, that computes numerically the riemann normal coordinates, which are constructed around the point Theta_0
#' It makes use of the formula for cubic coefficients of the expansion of the geodesic, since if you only have it up to second order, there could be points Theta for which the quadratic formula cannot be inverted.
#' @param bestfit The point "Theta_0" around which the rnc coordinates p are constructed
#' @param residual The residual at the point "bestfit"
#' @param christoffelsymbols The christoffel symbols at the point "bestfit"
#' 
#' @param Theta Parameter of the returned function. A matrix with coordinate parameters, each row corresponds to one set of rnc parameters p
#' 
#' @export
pfunnum<- function(bestfit, residual, christoffelsymbols = christoffels(residual)) {

        cubiccoeffs <- array(apply(christoffelsymbols, c(2,3), function(chrisupper) {
          apply(christoffelsymbols, c(1,2), function(chrislower) sum(chrislower*chrisupper))
        }), dim = rep(length(bestfit), 4))
        cubiccoeffs <- 2*cubiccoeffs - christoffelderiv(residual)

  
  pnum <- function(Theta) { 
    Theta <- as.vector(Theta)
    
    solvefun <- function(p) {Theta - bestfit+p-as.vector(0.5*apply(christoffelsymbols, 3, "%*%", p)%*%p)+as.vector(1/6*apply(apply(cubiccoeffs, c(2,3), "%*%", p), 2, "%*%", p)%*%p)}
    # solvefun <- function(p) {Theta - as.vector(bestfit+p-0.5*apply(christoffelsymbols, 3, "%*%", p)%*%p)}
    
    myroot <- as.numeric(multiroot(solvefun, bestfit)$root)
    names(myroot)<- names(bestfit)
    return(myroot)
  }
  return(pnum)
}

#' This function returns a function for the coordinates of a geodesic in Theta (parameter-) space.
#' @param direction Argument of the returned function.The vector of the initial direction in both Theta and p space
#' @param tau Argument of the returned function. The stepwidth along the geodesic.
#' @export
pgeofun <- function() {
  function(direction, tau) {
    do.call(rbind,lapply(tau, function(tau){direction*tau}))
  }
}

#' Jacobian dp/dTheta
#' @details Returns the Jacobian (dp^mu/dTheta^nu)^mu_nu at the point Theta
#' 
#' @param bestfit The point "Theta_0" around which the rnc coordinates p are constructed
#' @param residual The residual at the point "bestfit"
#' @param christoffelsymbols The christoffel symbols at the point "bestfit"
#' 
#' @param Theta Argument of the returned function. Theta A matrix with coordinate parameters, each row corresponds to one set of parameters Theta
#' 
#' @export
JpTfun <- function(bestfit, residual, christoffelsymbols = christoffels(residual)) {
  function(Theta){
        return(diag(1, length(bestfit)) + apply(christoffelsymbols, 3, "%*%", (Theta-bestfit)))
  }
}


#' Metric in parameter (p)-space
#' @description The Tmetric is being transformed into the metric of p.
#' For this the coordinate transformation has to be defined already with the respective functions.
#' The underlying formula is gmunu(p)=t(JTp(p))mualpha%*%galphabeta(Theta(p))%*%JTpbetanu
#' @details Prerequisites of this function to work is to have
#' a) a function called resi(Theta) which gives out the residuals at point Theta
#' b) a function called Theta(p)
#' c) a function called JTP(p)
#' @param p The coordinate value of p
#' 
#' 
#' @export
pmetric <- function(p) {
    t(JTp(p))%*%Tmetric(resi(Theta(p)))%*%JTp(p)
}







#' Wrapper function for the coordinate Trafos
#' @param bestfit The point "Theta_0" around which the rnc coordinates p are constructed
#' @param residual The residual at the point "bestfit"
#' @param christoffelsymbols The christoffel symbols at the point "bestfit"
#' 
#' @export
alltrafos <- function(bestfit, residual, christoffelsymbols = christoffels(residual)){
  # Theta(p)
  Theta <<- Thetafun(bestfit, residual, christoffelsymbols = christoffels(residual))
  Thetageo <<- Thetageofun(bestfit, residual, christoffelsymbols = christoffels(residual))
  JTp <<- JTpfun (bestfit, residual, christoffelsymbols = christoffels(residual))
  # Tmetric <- Tmetricfun()
  # p(Theta)
  p <<- pfun(bestfit, residual, christoffelsymbols = christoffels(residual))
  pgeo <<- pgeofun()
  JpT <<- JpTfun(bestfit, residual, christoffelsymbols = christoffels(residual))
  # pmetric <- pmetricfun(bestfit, residual, christoffelsymbols = christoffels(residual))
  # p numerically
   pnum <<- pfunnum(bestfit, residual, christoffelsymbols = christoffels(residual))
    
  # list(
  #   Theta = Theta,
  #   Thetageo = Thetageo,
  #   JTp = JTp,
  #   # Tmetric = Tmetric,
  #   p = p,
  #   pgeo = pgeo,
  #   JpT = JpT#,
  #   # pmetric = pmetric,
    # pnum =pnum
  # )
}






#' Forward and backward coordinate transformation
#' 
#' @export
ThetapTheta <- function(Theta) {
  out<-Theta(p(Theta))
  attr(out, "difference") <- out - Theta
  return(out)
}

#' Backward and forward coordinate transformation
#' 
#' @export
pThetap <- function(pmu) {
  out <- p(Theta(pmu))
  attr(out, "difference") <- out - pmu
  return(out)
}






