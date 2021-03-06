
#' @export
coordTransform <- function(data, transformations) {
  
  mynames <- unique(as.character(data$name))
  
  # Replicate transformation if not a list
  if (!is.list(transformations))
    transformations <- as.list(structure(rep(transformations, length(mynames)), names = mynames))
  
  out <- do.call(rbind, lapply(mynames, function(n) {
    
    subdata <- subset(data, name == n)
    
    if (n %in% names(transformations)) {
      
      mysymbol <- getSymbols(transformations[[n]])[1]
      mytrafo <- replaceSymbols(mysymbol, "value", transformations[[n]])
      mytrafo <- parse(text = mytrafo)
      
      if ("sigma" %in% colnames(subdata))
        subdata$sigma <- abs(with(subdata, eval(D(mytrafo, "value")))) * subdata$sigma
      subdata$value <- with(subdata, eval(mytrafo))
      
    }
    
    return(subdata)
    
  }))
  

  return(out)
  
  
}



#' @export
theme_dMod <- function (base_size = 12, base_family = "") {
  colors <- list(
    medium = c(gray = '#737373', red = '#F15A60', green = '#7AC36A', blue = '#5A9BD4', orange = '#FAA75B', purple = '#9E67AB', maroon = '#CE7058', magenta = '#D77FB4'),
    dark = c(black = '#010202', red = '#EE2E2F', green = '#008C48', blue = '#185AA9', orange = '#F47D23', purple = '#662C91', maroon = '#A21D21', magenta = '#B43894'),
    light = c(gray = '#CCCCCC', red = '#F2AFAD', green = '#D9E4AA', blue = '#B8D2EC', orange = '#F3D1B0', purple = '#D5B2D4', maroon = '#DDB9A9', magenta = '#EBC0DA')
  )
  gray <- colors$medium["gray"]
  black <- colors$dark["black"]
  
  theme_bw(base_size = base_size, base_family = base_family) + 
    theme(line = element_line(colour = gray), 
          rect = element_rect(fill = "white", colour = NA), 
          text = element_text(colour = black), 
          axis.ticks = element_line(colour = gray), 
          legend.key = element_rect(colour = NA), 
          panel.border = element_rect(colour = gray), 
          panel.grid = element_line(colour = gray, size = 0.2), 
          strip.background = element_rect(fill = "white", colour = NA))
}


ggplot <- function(...) ggplot2::ggplot(...) + theme_dMod()


#' Plot a list of model predictions
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by \link{Xs}.
#' @param ... Further arguments going to \code{subset}. 
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet Either \code{"wrap"} or \code{"grid"}
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @export
#' @import ggplot2
plotPrediction <- function(prediction, ..., scales = "free", facet = "wrap", transform = NULL) {
  
  prediction <- subset(wide2long.list(prediction), ...)
  
  if (!is.null(transform)) prediction <- coordTransform(prediction, transform)
  
  if (facet == "wrap")
    p <- ggplot(prediction, aes(x = time, y = value, group = condition, color = condition)) + facet_wrap(~name, scales = scales)
  if (facet == "grid")
    p <- ggplot(prediction, aes(x = time, y = value)) + facet_grid(name ~ condition, scales = scales)
 
  p <- p + geom_line() 
  
  attr(p, "data") <- prediction
  return(p)
   
}


#' Plot a list of model predictions and a list of data points in a combined plot
#' 
#' @param prediction Named list of matrices or data.frames, usually the output of a prediction function
#' as generated by \link{Xs}.
#' @param data Named list of data.frames as being used in \link{res}, i.e. with columns \code{name}, \code{time}, 
#' \code{value} and \code{sigma}.
#' @param ... Further arguments going to \code{subset}. 
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet Either \code{"wrap"} or \code{"grid"}
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{sigma},
#' \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @export
plotCombined <- function (prediction, data = NULL, ..., scales = "free", facet = "wrap", transform = NULL) {
  
  mynames <- c("time", "name", "value", "sigma", "condition")
  
  if (!is.null(prediction)) {
    prediction <- cbind(wide2long(prediction), sigma = NA)
    prediction <- subset(prediction, ...)
    
    if (!is.null(transform)) prediction <- coordTransform(prediction, transform)
    
  }
  
  if (!is.null(data)) {
    data <- lbind(data)
    data <- subset(data, ...)
    
    if (!is.null(transform)) data <- coordTransform(data, transform)
    
  }
  
  total <- rbind(prediction[, mynames], data[, mynames])
  
  if (facet == "wrap")
    p <- ggplot(total, aes(x = time, y = value, ymin = value - sigma, ymax = value + sigma, 
                           group = condition, color = condition)) + facet_wrap(~name, scales = scales)
  if (facet == "grid")
    p <- ggplot(total, aes(x = time, y = value, ymin = value - sigma, ymax = value + sigma)) + facet_grid(name ~ condition, scales = scales)
  
  if (!is.null(prediction))
    p <- p +  geom_line(data = prediction)
  
  if (!is.null(data))
    p <- p + geom_point(data = data) + geom_errorbar(data = data, width = 0)
  
  attr(p, "data") <- list(data = data, prediction = prediction)
  return(p)
  
}

#' Plot predictions that are already in long format
#'
#' @param prediction 
#' @param data 
#' @param ... 
#' @param scales 
#' @param facet 
#' @param transform 
#'
#' @return
#' @export
#'
#' @examples
plotCombined_long <- function (prediction, data = NULL, ..., scales = "free", facet = "wrap", transform = NULL, whichPar = 1) {
  
  mynames <- c("time", "name", "value", "condition", outer_pars[whichPar])
  mygrouping <- as.name(outer_pars[whichPar])
  if (!is.null(data)) {
    data <- lbind(data)
    data <- subset(data, ...)
    
    if (!is.null(transform)) data <- coordTransform(data, transform)
    
  }
  
  total <- rbind(prediction[, mynames])
  
  if (facet == "wrap")
    p <- ggplot(total, aes_(x = ~time, y = ~value, ymin = ~(value - sigma), ymax = ~(value + sigma), 
                           group = mygrouping, color = mygrouping)) + facet_wrap(~name, scales = scales)
  if (facet == "grid")
    p <- ggplot(total, aes(x = time, y = value, ymin = value - sigma, ymax = value + sigma)) + facet_grid(name ~ condition, scales = scales)
  
  if (!is.null(prediction))
    p <- p +  geom_line(data = prediction)
  
  if (!is.null(data))
    p <- p + geom_point(data = data) + geom_errorbar(data = data, width = 0)
  
  attr(p, "data") <- list(data = data, prediction = prediction)
  return(p)
  
}




#' Plot a list data points
#' 
#' @param data Named list of data.frames as being used in \link{res}, i.e. with columns \code{name}, \code{time}, 
#' \code{value} and \code{sigma}.
#' @param ... Further arguments going to \code{subset}. 
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet Either \code{"wrap"} or \code{"grid"}
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{sigma},
#' \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @export
plotData <- function (data, ..., scales = "free", facet = "wrap", transform = NULL) {
  
  data <- subset(lbind(data), ...)
  
  if (!is.null(transform)) data <- coordTransform(data, transform)

  if (facet == "wrap")
    p <-  ggplot(data, aes(x = time, y = value, ymin = value - sigma, 
                           ymax = value + sigma, group = condition, color = condition)) + facet_wrap(~name, scales = scales)
  if (facet == "grid")
    p <- ggplot(data, aes(x = time, y = value, ymin = value - sigma, 
                          ymax = value + sigma)) +  facet_grid(name ~ condition, scales = scales)
  
  p <- p + geom_point() + geom_errorbar(width = 0)
  
  
  attr(p, "data") <- data
  return(p)
  
}

#' Profile likelihood plot
#' 
#' @param ... Lists of profiles as being returned by \link{profile}.
#' @param maxvalue Numeric, the value where profiles are cut off.
#' @param parlist Matrix or data.frame with columns for the parameters to be added to the plot as points.
#' If a "value" column is contained, deltas are calculated with respect to lowest chisquare of profiles.
#' @return A plot object of class \code{ggplot}.
#' @export
plotProfile <- function(..., maxvalue = 5, parlist = NULL) {
  
  
  arglist <- list(...)
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    proflist <- arglist[[i]]
    
    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    do.valueData <- "valueData" %in% colnames(proflist[[1]])
    do.valuePrior <- "valuePrior" %in% colnames(proflist[[1]])
    
    
    # Discard faulty profiles
    proflistidx <- sapply(proflist, function(prf) any(class(prf) == "data.frame"))
    proflist <- proflist[proflistidx]
    if (sum(!proflistidx) > 0) {
      warning(sum(!proflistidx), " profiles discarded.", call. = FALSE)
    }
    
    subdata <- do.call(rbind, lapply(names(proflist), function(n) {
      
      values <- proflist[[n]][, "value"]
      origin <- which.min(abs(proflist[[n]][, "constraint"]))
      zerovalue <- proflist[[n]][origin, "value"]
      parvalues <- proflist[[n]][, n]
      deltavalues <- values - zerovalue

      sub <- subset(data.frame(name = n, delta = deltavalues, par = parvalues, proflist = i, mode="total", is.zero = 1:nrow(proflist[[n]]) == origin), delta <= maxvalue)
      
      obj.attributes <- attr(proflist[[n]], "obj.attributes")
      if(!is.null(obj.attributes)) {
        for(mode in obj.attributes) {
          valuesO <- proflist[[n]][, mode]
          originO <- which.min(abs(proflist[[n]][, "constraint"]))
          zerovalueO <- proflist[[n]][originO, mode]
          deltavaluesO <- valuesO - zerovalueO
          sub <- rbind(sub,subset(data.frame(name = n, delta = deltavaluesO, par = parvalues, proflist = i, mode=mode, is.zero = 1:nrow(proflist[[n]]) == originO), delta <= maxvalue))
        }
      }
      
      return(sub)
    }))
    return(subdata)
  }))
  
  data$proflist <- as.factor(data$proflist)
  data.zero <- subset(data, is.zero)

  threshold <- c(1, 2.7, 3.84)
  
  p <- ggplot(data, aes(x=par, y=delta, group=interaction(proflist,mode), color=proflist, linetype=mode)) + facet_wrap(~name, scales="free_x") + 
    geom_line() + #geom_point(aes=aes(size=1), alpha=1/3) +
    geom_point(data = data.zero) +
    geom_hline(yintercept=threshold, lty=2, color="gray") + 
    ylab(expression(paste("CL /", Delta*chi^2))) +
    scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84")) +
    xlab("parameter value")
  
  if(!is.null(parlist)){
    delta <- 0
    if("value" %in% colnames(parlist)){
      minval <- min(unlist(lapply(1:length(arglist), function(i){ 
        origin <- which.min(arglist[[i]][[1]][, "constraint"])
        zerovalue <- arglist[[i]][[1]][origin, 1]  
      })))
      values <- parlist[,"value"]
      parlist <- parlist[,!(colnames(parlist) %in% "value")]
      delta <- as.numeric(values - minval)
    }
    points <- data.frame(par = as.numeric(as.matrix(parlist)), name = rep(colnames(parlist), each = nrow(parlist)), delta = delta)

    #points <- data.frame(name = colnames(parlist), par = as.numeric(parlist), delta=0)
    p <- p + geom_point(data=points, aes(x=par, y=delta, group=NULL, linetype = NULL), color = "black")
  }
  attr(p, "data") <- data
  return(p)
  
}

#' Profile likelihood: plot of the parameter paths.
#' 
#' @param ... Lists of profiles as being returned by \link{profile}.
#' @param whichPar Character or index vector, indicating the parameters that are taken as possible reference (x-axis)
#' @param sort Logical. If paths from different parameter profiles are plotted together, possible
#' combinations are either sorted or all combinations are taken as they are.
#' @return A plot object of class \code{ggplot}.
#' @export
plotPaths <- function(..., whichPar = NULL, sort = FALSE, relative = TRUE, scales = "free") {
  
  arglist <- list(...)
  
  
  data <- do.call(rbind, lapply(1:length(arglist), function(i) {
    # choose a proflist
    proflist <- arglist[[i]]
    
    if(is.data.frame(proflist)) {
      whichPars <- unique(proflist$whichPar)
      proflist <- lapply(whichPars, function(n) {
        with(proflist, proflist[whichPar == n, ])
      })
      names(proflist) <- whichPars
    }
    
    if(is.null(whichPar)) whichPar <- names(proflist)
    if(is.numeric(whichPar)) whichPar <- names(proflist)[whichPar]
    subdata <- do.call(rbind, lapply(whichPar, function(n) {
      # chose a profile
      parameters <- attr(proflist[[n]], "parameters")
      # matirx
      paths <- as.matrix(proflist[[n]][, parameters])
      values <- proflist[[n]][,1]
      if(relative) 
        for(j in 1:ncol(paths)) paths[, j] <- paths[, j] - paths[which.min(values), j]
      combinations <- expand.grid.alt(whichPar, colnames(paths))
      if(sort) combinations <- apply(combinations, 1, sort) else combinations <- apply(combinations, 1, identity)
      combinations <- submatrix(combinations, cols = -which(combinations[1,] == combinations[2,]))
      combinations <- submatrix(combinations, cols = !duplicated(paste(combinations[1,], combinations[2,])))
      
      
      
      
      path.data <- do.call(rbind, lapply(1:dim(combinations)[2], function(j) {
        data.frame(chisquare = values, 
                   name = n,
                   proflist = i,
                   combination = paste(combinations[,j], collapse = " - \n"),
                   x = paths[, combinations[1,j]],
                   y = paths[, combinations[2,j]])
      }))
      
      return(path.data)
      
    }))
    
    return(subdata)
    
  }))
  
  data$proflist <- as.factor(data$proflist)
  
  
  if(relative)
    axis.labels <- c(expression(paste(Delta, "parameter 1")), expression(paste(Delta, "parameter 2")))  
  else
    axis.labels <- c("parameter 1", "parameter 2")
  
  
  p <- ggplot(data, aes(x=x, y=y, group=interaction(name, proflist), color=name, lty=proflist)) + 
    facet_wrap(~combination, scales = scales) + 
    geom_path() + geom_point(aes=aes(size=1), alpha=1/3) +
    xlab(axis.labels[1]) + ylab(axis.labels[2]) +
    scale_linetype_discrete(name = "profile\nlist") +
    scale_color_discrete(name = "profiled\nparameter")
  
  attr(p, "data") <- data
  return(p)
  
}





plotPredictionCont <- function(out, ...) {
  
  require(ggplot2)
  
  mynames <- c("time", "name", "value", "sigma", "condition")
  
  
  forclist <- lapply(out, function(o) attr(o, "forc"))
  out <- cbind(wide2long.list(out), sigma=NA)
  out <- subset(out, ...)
  targets <- as.character(unique(out$name))
  
  print(targets)
  
  forc <- lbind(lapply(forclist, function(fo) {
    names(fo) <- substr(names(fo), 1, nchar(names(fo))-1)
    print(names(fo))
    data <- do.call(rbind, lapply(targets[targets%in%names(fo)], function(t) {
      print(t)
      data.frame(time = fo[[t]][,1], name = t, value = fo[[t]][,2], sigma=1/sqrt(fo[[paste0("weight", t)]][,2]))
    }))
    
    return(data)
        
  }))
  
  ggplot(rbind(out[, mynames], forc[, mynames]), 
         aes(x = time, 
             y = value, ymin = value - sigma, ymax = value + sigma, 
             group = condition, color = condition, fill=condition)) + 
    facet_wrap(~name, scales = "free") + 
    geom_line(data = out) + 
    geom_line(data = forc, lty=2) + 
    geom_ribbon(data = forc, alpha=0.3, lty=0)
  
  
}



#' Plot an array of model predictions for a list of parameters
#' 
#' @param parframe Object of class \code{parframe}, e.g. returned by \link{mstrust} or \link{profile}
#' @param x The model prediction function \code{x(times, pars, fixed, ...)}
#' @param times Numeric vector of time points for the model prediction
#' @param data Named list of data.frames as being used in \link{res}, i.e. with columns \code{name}, \code{time}, 
#' \code{value} and \code{sigma}.
#' @param ... Further arguments going to \code{subset}.
#' @param fixed Named numeric vector with fixed parameters
#' @param deriv Logical. If \code{x} supports the argument \code{deriv}, it is used.
#' @param scales The scales argument of \code{facet_wrap} or \code{facet_grid}, i.e. \code{"free"}, \code{"fixed"}, 
#' \code{"free_x"} or \code{"free_y"}
#' @param facet Either \code{"wrap"} or \code{"grid"}
#' @details The data.frame being plotted has columns \code{time}, \code{value}, \code{sigma},
#' \code{name} and \code{condition}.
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @export
plotArray <- function(parlist, x, times, data = NULL, ..., fixed = NULL, deriv = FALSE, scales = "free", facet = "wrap") {

  parameters <- attr(parlist, "parameters")
  
  pars<- lapply(1:nrow(parlist), function(i) unlist(parlist[i, c("value", parameters)]))
  
  
  prediction <- lapply(pars, function(p) {
    pred <- x(times, p, fixed, deriv = deriv)
    newnames <- sapply(names(pred), function(cond) paste(colnames(pred[[cond]])[-1], cond, sep = ",\n "))
    pred <- do.call(cbind, pred)
    pred <- pred[, -which(colnames(pred) == "time")]
    pred <- cbind(times, pred)
    colnames(pred) <- c("time", newnames)
    return(pred)
  }); names(prediction) <- parlist[ ,"value"]
  
  
  if(!is.null(data)) {
    for(n in names(data)) data[[n]]$name <- paste(data[[n]]$name, n, sep = ",\n ")
    data <- do.call(rbind, data)
    data <- list(data)
    names(data) <- names(prediction)[1]
  }
  
  plotCombined(prediction, data, ..., scales = scales, facet = facet) 
  
  
}

#' Plot Fluxes given a list of flux Equations
#' 
#' @param pouter parameters
#' @param x The model prediction function \code{x(times, pouter, fixed, ...)} needs to return pinner as attribute,
#' e.g.:\cr
#'  \code{x <- function(times, pouter, fixed=NULL, attach=TRUE, ...) {\cr
#'  out <- lapply(conditions, function(cond) { \cr
#'  pinner <- pL[[cond]](pouter, fixed) \cr
#'  prediction <- xL[[cond]](times, pinner, ...)\cr
#'  observation <- g(prediction, pinner, attach.input = attach)\cr
#'  attr(observation, "pinner") <- pinner\cr
#'  return(observation)\cr
#' }); names(out) <- conditions\cr
#' return(out)}
#' }
#' @param fluxEquations list of chars containing expressions for the fluxes, 
#' if names are given, they are shown in the legend. Easy to obtain via \link{subset.eqnList}, see Examples.
#' @param times Numeric vector of time points for the model prediction
#' @param fixed Named numeric vector with fixed parameters
#'  
#' 
#' @return A plot object of class \code{ggplot}.
#' @examples 
#' \dontrun{
#' plotFluxes(bestfit,times,attr(subset(f,"B"%in%Product), "rates"),nameFlux="B production")
#' }
#' @export
plotFluxes <- function(pouter, x, times, fluxEquations, nameFlux = "Fluxes:", fixed = NULL){
  if(is.null(names(fluxEquations))) names(fluxEquations) <- fluxEquations
  flux <- funC0(fluxEquations)
  prediction.all <- x(times, pouter, fixed, deriv = FALSE)
  
  out <- lapply(names(prediction.all), function(cond) {
    prediction <- prediction.all[[cond]]
    pinner <- attr(prediction,"pinner")
    pinner.list <- as.list(pinner)
    prediction.list <- as.list(as.data.frame(prediction))
    fluxes <- cbind(time=prediction[,"time"],flux(c(prediction.list, pinner.list)))
    return(fluxes)
  }); names(out) <- names(prediction.all)
  out <- wide2long(out)
  
  cbPalette <- c("#999999", "#E69F00", "#F0E442", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7","#CC6666", "#9999CC", "#66CC99","red", "blue", "green","black")
  P <- ggplot(out, aes(x=time, y=value, group=name, fill=name, log="y"))+ ylab("flux") + 
    facet_wrap(~condition) + scale_fill_manual(values=cbPalette, name=nameFlux) +
    geom_density(stat="identity", position="stack", alpha=0.3, color="darkgrey", size=0.4) +
    xlab("time [min]") + ylab("contribution [a.u.]")
  
  return(P)
  
}

#plotArray <- function(x, times, fitlist, data = NULL, fixed=NULL, binwidth=10) {
#  
#  require("ggplot2")
#  require("wq")
#  
#  n <- dim(fitlist)[1]
#  chisquare <- as.data.frame(log10(fitlist[,1]))
#  colnames(chisquare) <- "logchisquare"
#  
#  
#  out <- do.call(rbind, lapply(1:n, function(i) {
#    myout <- x(times, c(fitlist[i,-1], fixed), deriv=FALSE)
#    myout <- wide2long(myout)
#    myout <- cbind(myout, logchisquare = log10(fitlist[i,1]))
#    return(myout)
#  }))
#  
#  
#  P1 <- ggplot(out, aes(x=time, y=value, group=logchisquare, color=logchisquare)) +
#    facet_wrap(~name, scales="free") 
#  if(!is.null(data)) {
#    data <- cbind(data, logchisquare = NA)
#    P1 <- P1 + 
#      geom_point(aes(x=time, y=value), data = data, color="gray") + 
#      geom_errorbar(aes(x=time, y=value, ymin=value-sigma, ymax=value+sigma), data=data, color="gray", width=0)
#  }
#  P1 <- P1 +
#    geom_line(alpha=.3) +
#    scale_color_gradientn(colours = rainbow(7)) +
#    theme(legend.position="none")
#  
#  
#  #   P3 <- ggplot(out, aes(x=logchisquare, group=logchisquare, fill=logchisquare)) + 
#  #     geom_histogram(binwidth=binwidth) + 
#  #     scale_fill_gradientn(colours = rainbow(7)) +
#  #     theme(legend.position="none") +
#  #     xlab("") +
#  #     coord_flip()
#  
#  
#  
#  P2 <- ggplot(chisquare, aes(x = 1:length(logchisquare), y=logchisquare, color=logchisquare)) +
#    geom_point() +
#    scale_color_gradientn(colours = rainbow(7)) + 
#    theme(legend.position="none") +
#    xlab("sorted index") +
#    ylab(expression(log[10](chi^2)))
#  
#  layOut(list(P2, 1, 1),
#         list(P1, 1, 2:5))
#  
#  
#}


plotObjective <- function(out) {
  require("ggplot2")
  require("wq")
  
  value <- out$value
  
  gradient <- out$gradient
  npar <- length(gradient)
  names <- factor(paste(names(gradient), 1:npar, sep=", "), levels = paste(names(gradient), 1:npar, sep=", "))
  gradient.data <- data.frame(name = names, value = gradient)
  
  
  
  hessian <- out$hessian
  hessian.data <- data.frame(x = as.factor(1:npar), y=as.factor(rep(1:npar, each=npar)), hessian = as.vector(hessian))
  
  
  P1 <- ggplot(gradient.data, aes(x=name, y= value)) + 
    geom_bar(stat="identity") + 
    coord_flip() + ylab("gradient value") + xlab("parameter name, i") + 
    ggtitle(paste("obj. value:", value)) 
  
  
  P2 <- ggplot(hessian.data, aes(x=x, y=y, z=hessian, fill=hessian)) + 
    geom_tile(color="gray") + scale_fill_gradient2() + xlab("i") + ylab("j") + 
    theme(legend.position=c(0, 1), legend.justification=c(0,1))
  
  layOut(list(P2, 1, 1:2),
         list(P1, 1, 3))
  
  
  
}

#' @export
plotValues <- function(pars, values = "value") {
  
  mycolnames <- colnames(pars)
  mycolnames[mycolnames == values] <- "value"
  colnames(pars) <- mycolnames
 
  pars <- cbind(index = 1:nrow(pars), pars)
   
  ggplot(pars, aes(x = index, y = value, pch = converged, color = iterations)) + geom_point() + 
    xlab("index") + ylab("value")
  
}

#' Plot parameter values for a fitlist
#' 
#' @param myparframe fitlist as obtained by as.parframe(mstrust)
#' @param whichFits indexlist e.g. 1:10
#' @export
plotPars <- function(myparframe, whichFits = 1:length(myparframe)){
  parNames <- attr(myparframe,"parameters")
  parOut <- wide2long.data.frame(out = ((myparframe[whichFits,c("value",parNames)])) , keep = 1)
  names(parOut) <- c("value","name","parvalue")
  plot <- ggplot(parOut, aes(x = name, y = parvalue, color = value)) + geom_point() + theme(axis.text.x = element_text(angle = 270, hjust = 0))
  return(plot)
}

plotFitList <- function(fitlist) {
  require(ggplot2)
    
  ggplot(fitlist, aes(x=chisquare, y=value)) + 
    facet_wrap(~name, scales="free") + 
    geom_point()  
  
}


plotFluxesOld <- function(out, fluxEquations, pars) {
  
  require(scales)
  
  if(any(class(out)=="list")) out <- wide2long(out)
  
  nFluxes <- length(fluxEquations)
  if(is.null(names(fluxEquations))) names(fluxEquations) <- paste0("reaction", 1:nFluxes)
  fluxEquations <- c(fluxEquations, sum = paste(fluxEquations, collapse="+"))
  
  # Evaluate fluxes
  fluxes <- with(c(as.list(out), as.list(pars)), {
    flux <- do.call(rbind, lapply(1:(nFluxes+1), function(i) {
      ev <- eval(parse(text=fluxEquations[i]))
      nout <- data.frame(time = out$time, 
                         name = names(fluxEquations)[i], 
                         value = ev, 
                         condition = out$condition)
    }))
    return(flux)
  })
  
  fluxes1 <- subset(fluxes, name != "sum")
  fluxes2 <- subset(fluxes, name == "sum")
  
  
  P <- ggplot(out, aes(x=time, y=value, group=name, fill=name)) + facet_wrap(~condition) +
    geom_density(stat="identity", position="stack", alpha=0.3, color="darkgrey", size=0.4, data=fluxes1) +
    geom_line(aes(x=time, y=value, group=NULL, fill=NULL), color="black", data=fluxes2)
  
  
  return(P)
  
}
