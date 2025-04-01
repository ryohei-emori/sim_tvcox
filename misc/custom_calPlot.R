custom_calPlot <- function (object, time, formula, data, splitMethod = "none", 
          B = 1, M, pseudo, type, showPseudo, pseudo.col = NULL, pseudo.pch = NULL, 
          method = "nne", round = TRUE, bandwidth = NULL, q = 10, bars = FALSE, 
          hanging = FALSE, names = "quantiles", showFrequencies = FALSE, 
          jack.density = 55, plot = TRUE, add = FALSE, diag = !add, 
          legend = !add, axes = !add, xlim = c(0, 1), ylim = c(0, 1), 
          xlab, ylab, col, lwd, lty, pch, cause = 1, percent = TRUE, 
          giveToModel = NULL, na.action = na.fail, cores = 1, verbose = FALSE, 
          cex = 1, ...) 
{
  if (missing(pseudo)) {
    if (method == "quantiles" || bars == TRUE) 
      pseudo <- FALSE
    else pseudo <- TRUE
  }
  if (missing(showPseudo)) 
    showPseudo <- ifelse(add || (pseudo != FALSE), FALSE, 
                         TRUE)
  if (!inherits(x = object, what = "list")) {
    object <- list(object)
  }
  if (is.null(names(object))) 
    names(object) <- paste0("Model.", 1:length(object))
  if (bars) {
    method = "quantile"
    if (!(length(object) == 1)) 
      stop(paste0("Barplots work only for one prediction at a time. Provided are ", 
                  length(object), "predictions"))
  }
  if (is.null(names(object))) {
    names(object) <- sapply(object, function(o) class(o)[1])
    names(object) <- make.names(names(object), unique = TRUE)
  }
  else {
    names(object)[(names(object) == "")] <- sapply(object[(names(object) == 
                                                             "")], function(o) class(o)[1])
  }
  NF <- length(object)
  if (missing(lwd)) 
    lwd <- rep(3, NF)
  if (missing(col)) {
    if (bars) 
      col <- c("grey90", "grey30")
    else col <- 1:NF
  }
  if (missing(lty)) 
    lty <- rep(1, NF)
  if (missing(pch)) 
    pch <- rep(1, NF)
  if (length(lwd) < NF) 
    lwd <- rep(lwd, NF)
  if (length(lty) < NF) 
    lty <- rep(lty, NF)
  if (length(col) < NF) 
    col <- rep(col, NF)
  if (length(pch) < NF) 
    pch <- rep(pch, NF)
  if (missing(data)) {
    data <- eval(object[[1]]$call$data)
    if (!inherits(x = data, what = "data.frame")) 
      stop("Argument data is missing.")
    else if (verbose) 
      warning("Argument data is missing. I use the data from the call to the first model instead.")
  }
  if (missing(formula)) {
    if (length(grep("~", as.character(object[[1]]$call$formula))) == 
        0) {
      stop(paste("Argument formula is missing and first model has no usable formula:", 
                 as.character(object[[1]]$call$formula)))
    }
    else {
      ftry <- try(formula <- eval(object[[1]]$call$formula), 
                  silent = TRUE)
      if (inherits(x = ftry, what = "try-error") || !inherits(x = formula, 
                                                              what = "formula")) 
        stop("Argument formula is missing and first model has no usable formula.")
      else if (verbose) 
        warning("Formula missing. Using formula from first model")
    }
  }
  m <- model.frame(formula, data, na.action = na.action)
  response <- model.response(m)
  if (inherits(x = response, what = "Surv")) 
    model.type <- "survival"
  else model.type <- attr(response, "model")
  if (is.null(model.type) & length(unique(response)) == 2) 
    stop("This function works only for survival and competing risks models.")
  if (missing(type)) 
    type <- ifelse(model.type == "survival", "survival", 
                   "risk")
  if (missing(ylab)) 
    if (bars) 
      ylab = ""
  else ylab <- ifelse(type == "survival", "Observed survival frequencies", 
                      "Observed event frequencies")
  if (type == "survival" && !(model.type %in% c("survival", 
                                                "binary"))) 
    stop(paste0("Type survival works only in survival or binary outcome models. This is a ", 
                model.type, " model"))
  if (!(model.type == "binary")) {
    neworder <- order(response[, "time"], -response[, "status"])
    response <- response[neworder, , drop = FALSE]
    Y <- response[, "time"]
    data <- data[neworder, ]
    if (missing(time)) 
      time <- median(Y)
    else if (length(time) > 1) 
      stop("Please specify only one time point.")
  }
  predictHandlerFun <- switch(model.type, binary = "predictStatusProb", 
                              competing.risks = "predictEventProb", survival = "predictSurvProb")
  if (pseudo == FALSE && splitMethod != "none") 
    stop(paste0("Split method ", splitMethod, " is only implemented for : 'pseudo=TRUE'."))
  if (model.type == "binary") 
    if (is.factor(response)) 
      jack <- as.numeric(response == levels(response)[2])
  else jack <- as.numeric(response)
  else {
    if (pseudo == TRUE) {
      margForm <- update(formula, paste(".~1"))
      margFit <- prodlim::prodlim(margForm, data = data)
      jack <- prodlim::jackknife(margFit, cause = cause, 
                                 times = time)
    }
    else {
      jack <- NULL
    }
  }
  axis1.DefaultArgs <- list(side = 1, las = 1, at = seq(0, 
                                                        ylim[2], ylim[2]/4))
  axis2.DefaultArgs <- list(side = 2, las = 2, at = seq(0, 
                                                        ylim[2], ylim[2]/4), mgp = c(4, 1, 0))
  if (bars) {
    legend.DefaultArgs <- list(legend = names(object), col = col, 
                               cex = cex, bty = "n", x = "topleft")
    names.DefaultArgs <- list(cex = 0.7 * par()$cex, y = c(-abs(diff(ylim))/15, 
                                                           -abs(diff(ylim))/25))
    frequencies.DefaultArgs <- list(cex = 0.7 * par()$cex, 
                                    percent = FALSE, offset = 0)
  }
  else {
    legend.DefaultArgs <- list(legend = names(object), lwd = lwd, 
                               col = col, lty = lty, cex = cex, bty = "n", y.intersp = 1.3, 
                               x = "topleft")
  }
  if (bars) {
    if (type == "survival") 
      legend.DefaultArgs$legend <- c("Predicted survival", 
                                     "Observed frequencies")
    else legend.DefaultArgs$legend <- c("Predicted risks", 
                                        "Observed frequencies")
  }
  lines.DefaultArgs <- list(type = "l")
  abline.DefaultArgs <- list(lwd = 1, col = "red")
  if (missing(ylim)) {
    if (showPseudo && !bars) {
      ylim <- c(min(jack), max(jack))
    }
    else ylim <- c(0, 1)
  }
  if (missing(xlim)) {
    xlim <- c(0, 1)
  }
  if (missing(xlab)) 
    if (bars) 
      xlab <- ifelse(type == "survival", "Survival groups", 
                     "Risk groups")
  else xlab <- ifelse(type == "survival", "Predicted survival probability", 
                      "Predicted event probability")
  plot.DefaultArgs <- list(x = 0, y = 0, type = "n", ylim = ylim, 
                           xlim = xlim, ylab = ylab, xlab = xlab)
  barplot.DefaultArgs <- list(ylim = ylim, col = col, axes = FALSE, 
                              ylab = ylab, xlab = xlab, beside = TRUE, legend.text = NULL, 
                              cex.axis = cex, cex.lab = par()$cex.lab, cex.names = cex)
  if (bars) 
    control <- prodlim::SmartControl(call = list(...), keys = c("barplot", 
                                                                "legend", "axis2", "abline", "names", "frequencies"), 
                                     ignore = NULL, ignore.case = TRUE, defaults = list(barplot = barplot.DefaultArgs, 
                                                                                        abline = abline.DefaultArgs, legend = legend.DefaultArgs, 
                                                                                        names = names.DefaultArgs, frequencies = frequencies.DefaultArgs, 
                                                                                        axis2 = axis2.DefaultArgs), forced = list(abline = list(h = 0)), 
                                     verbose = TRUE)
  else control <- prodlim::SmartControl(call = list(...), keys = c("plot", 
                                                                   "lines", "legend", "axis1", "axis2"), ignore = NULL, 
                                        ignore.case = TRUE, defaults = list(plot = plot.DefaultArgs, 
                                                                            lines = lines.DefaultArgs, legend = legend.DefaultArgs, 
                                                                            axis1 = axis1.DefaultArgs, axis2 = axis2.DefaultArgs), 
                                        forced = list(plot = list(axes = FALSE), axis1 = list(side = 1)), 
                                        verbose = TRUE)
  splitMethod <- resolvesplitMethod(splitMethod = splitMethod, 
                                    B = B, N = NROW(data), M = M)
  k <- splitMethod$k
  B <- splitMethod$B
  N <- splitMethod$N
  NF <- length(object)
  apppred <- do.call("cbind", lapply(1:NF, function(f) {
    fit <- object[[f]]
    if (inherits(x = fit, what = "numeric") || inherits(x = fit, 
                                                        what = "double")) 
      fit <- matrix(fit, ncol = 1)
    apppred <- switch(model.type, competing.risks = {
      p <- as.vector(do.call(predictHandlerFun, list(fit, 
                                                     newdata = data, times = time, cause = cause)))
      if (inherits(x = fit, what = "numeric") || inherits(x = fit, 
                                                          what = "matrix")) p <- p[neworder]
      p
    }, survival = {
      p <- as.vector(do.call(predictHandlerFun, list(fit, 
                                                     newdata = data, times = time)))
      if (inherits(x = fit, what = "numeric") || inherits(x = fit, 
                                                          what = "matrix")) p <- p[neworder]
      p
    }, binary = {
      p <- do.call(predictHandlerFun, list(fit, newdata = data))
      if (inherits(x = fit, what = "numeric") || inherits(x = fit, 
                                                          what = "matrix")) p <- p[neworder]
      p
    })
  }))
  colnames(apppred) <- names(object)
  if (pseudo == TRUE) 
    apppred <- data.frame(jack = jack, apppred)
  else apppred <- data.frame(apppred)
  if (splitMethod$internal.name %in% c("noPlan")) {
    predframe <- apppred
  }
  if (splitMethod$internal.name %in% c("crossval", "loocv")) {
    groups <- splitMethod$index[, 1, drop = TRUE]
    cv.list <- lapply(1:k, function(g) {
      if (verbose == TRUE) 
        internalTalk(g, k)
      id <- groups == g
      train.k <- data[!id, , drop = FALSE]
      val.k <- data[id, , drop = FALSE]
      model.pred <- lapply(1:NF, function(f) {
        extraArgs <- giveToModel[[f]]
        fit <- object[[f]]
        fit.k <- internalReevalFit(object = fit, data = train.k, 
                                   step = paste("CV group=", k), silent = FALSE, 
                                   verbose = verbose)
        switch(model.type, competing.risks = {
          do.call(predictHandlerFun, list(object = fit.k, 
                                          newdata = val.k, times = time, cause = cause))
        }, survival = {
          p <- do.call(predictHandlerFun, c(list(object = fit.k, 
                                                 newdata = val.k, times = time), extraArgs))
          p
        }, binary = {
          p <- do.call(predictHandlerFun, list(object = fit.k, 
                                               newdata = val.k))
          p
        })
      })
      model.pred
    })
    predframe <- do.call("cbind", lapply(1:NF, function(f) {
      pred <- do.call("rbind", lapply(cv.list, function(x) x[[f]]))
      if (splitMethod$internal.name != "loocv") {
        pred <- pred[order(order(groups)), ]
      }
      pred
    }))
    colnames(predframe) <- names(object)
    if (pseudo == TRUE) 
      predframe <- cbind(data.frame(jack = jack), predframe)
  }
  if (splitMethod$internal.name %in% c("Boot632plus", "BootCv", 
                                       "Boot632")) {
    if (splitMethod$internal.name %in% c("Boot632plus", "Boot632")) {
      stop("Don't know how to do the 632(+) for the calibration curve.")
    }
    ResampleIndex <- splitMethod$index
    pred.list <- parallel::mclapply(1:B, function(b) {
      if (verbose == TRUE) 
        internalTalk(b, B)
      jackRefit <- FALSE
      vindex.b <- match(1:N, unique(ResampleIndex[, b]), 
                        nomatch = 0) == 0
      val.b <- data[vindex.b, , drop = FALSE]
      if (jackRefit) {
        margFit.b <- prodlim::prodlim(margForm, data = val.b)
        jack.b <- prodlim::jackknife(margFit.b, cause = cause, 
                                     times = time)
      }
      else {
        jack.b <- jack[match(1:N, unique(ResampleIndex[, 
                                                       b]), nomatch = 0) == 0]
      }
      train.b <- data[ResampleIndex[, b], , drop = FALSE]
      frame.b <- data.frame(jack = jack.b)
      bootpred <- do.call("cbind", lapply(1:NF, function(f) {
        fit <- object[[f]]
        fit.b <- internalReevalFit(object = fit, data = train.b, 
                                   step = b, silent = FALSE, verbose = verbose)
        extraArgs <- giveToModel[[f]]
        try2predict <- try(pred.b <- switch(model.type, 
                                            competing.risks = {
                                              do.call(predictHandlerFun, list(object = fit.b, 
                                                                              newdata = val.b, times = time, cause = cause))
                                            }, survival = {
                                              p <- do.call(predictHandlerFun, c(list(object = fit.b, 
                                                                                     newdata = val.b, times = time), extraArgs))
                                              p
                                            }, binary = {
                                              p <- do.call(predictHandlerFun, list(object = fit.b, 
                                                                                   newdata = val.b))
                                              p
                                            }), silent = TRUE)
        if (inherits(try2predict, "try-error") == TRUE) {
          rep(NA, NROW(val.b))
        }
        else {
          pred.b
        }
      }))
      colnames(bootpred) <- names(object)
      cbind(frame.b, bootpred)
    }, mc.cores = cores)
    predframe <- do.call("rbind", pred.list)
    rm(pred.list)
  }
  method <- match.arg(method, c("quantile", "nne"))
  getXY <- function(f) {
    if (pseudo == TRUE) {
      p <- predframe[, f + 1]
      jackF <- predframe[, 1]
    }
    else {
      p <- predframe[, f]
    }
    switch(method, quantile = {
      if (length(q) == 1) groups <- quantile(p, seq(0, 
                                                    1, 1/q)) else {
                                                      groups <- q
                                                    }
      xgroups <- (groups[-(length(groups))] + groups[-1])/2
      pcut <- cut(p, groups, include.lowest = TRUE)
      if (pseudo == TRUE) {
        plotFrame = data.frame(Pred = tapply(p, pcut, 
                                             mean), Obs = pmin(1, pmax(0, tapply(jackF, 
                                                                                 pcut, mean))))
        attr(plotFrame, "quantiles") <- groups
        plotFrame
      } else {
        form.pcut <- update(formula, paste(".~pcut"))
        if (inherits(x = data, what = "data.table")) pdata <- cbind(data[, 
                                                                         all.vars(update(formula, ".~1")), drop = FALSE, 
                                                                         with = FALSE], pcut = pcut) else pdata <- cbind(data[, 
                                                                                                                              all.vars(update(formula, ".~1")), drop = FALSE], 
                                                                                                                         pcut = pcut)
                                                                         y <- unlist(predict(f <- prodlim::prodlim(form.pcut, 
                                                                                                                   data = pdata), cause = cause, newdata = data.frame(pcut = levels(pcut)), 
                                                                                             times = time, type = ifelse(model.type == "survival", 
                                                                                                                         "surv", "cuminc")))
                                                                         if (model.type == "survival") y[is.na(y)] <- min(y, 
                                                                                                                          na.rm = TRUE) else y[is.na(y)] <- max(y, na.rm = TRUE)
                                                                         plotFrame = data.frame(Pred = tapply(p, pcut, 
                                                                                                              mean), Obs = y)
                                                                         attr(plotFrame, "quantiles") <- groups
                                                                         plotFrame
      }
    }, nne = {
      if (pseudo == TRUE) {
        if (round == TRUE) {
          if (!is.null(bandwidth) && bandwidth >= 1) {
          } else {
            p <- round(p, 2)
          }
        }
        p <- na.omit(p)
        if (no <- length(attr(p, "na.action"))) warning("calPlot: removed ", 
                                                        no, " missing values in risk prediction.", 
                                                        call. = FALSE, immediate. = TRUE)
        if (is.null(bandwidth)) {
          if (length(p) > length(apppred[, f + 1])) {
            bw <- prodlim::neighborhood(apppred[, f + 
                                                  1])$bandwidth
          } else {
            bw <- prodlim::neighborhood(p)$bandwidth
          }
        } else {
          bw <- bandwidth
        }
        if (bw >= 1) {
          plotFrame <- data.frame(Pred = mean(p), Obs = mean(jackF))
        } else {
          nbh <- prodlim::meanNeighbors(x = p, y = jackF, 
                                        bandwidth = bw)
          plotFrame <- data.frame(Pred = nbh$uniqueX, 
                                  Obs = nbh$averageY)
        }
        attr(plotFrame, "bandwidth") <- bw
        plotFrame
      } else {
        form.p <- update(formula, paste(".~p"))
        if (inherits(x = data, what = "data.table")) pdata <- cbind(data[, 
                                                                         all.vars(update(formula, ".~1")), drop = FALSE, 
                                                                         with = FALSE], p = p) else pdata <- cbind(data[, 
                                                                                                                        all.vars(update(formula, ".~1")), drop = FALSE], 
                                                                                                                   p = p)
                                                                         y <- unlist(predict(prodlim::prodlim(form.p, 
                                                                                                              data = pdata), cause = cause, newdata = data.frame(p = sort(p)), 
                                                                                             times = time, type = ifelse(type == "survival", 
                                                                                                                         "surv", "cuminc")))
                                                                         plotFrame <- data.frame(Pred = sort(p), Obs = y)
                                                                         plotFrame
      }
    })
  }
  plotFrames <- lapply(1:NF, function(f) {
    getXY(f)
  })
  names(plotFrames) <- names(object)
  if (bars) {
    if (model.type == "survival" && type == "risk") 
      plotFrames[[1]] <- plotFrames[[1]][NROW(plotFrames[[1]]):1, 
      ]
    if ((is.logical(names[1]) && names[1] == TRUE) || names[1] %in% 
        c("quantiles.labels", "quantiles")) {
      qq <- attr(plotFrames[[1]], "quantiles")
      if (model.type == "survival" && type == "risk") 
        qq <- rev(1 - qq)
      if (names[1] == "quantiles.labels") {
        pp <- seq(0, 1, 1/q)
        names <- paste0("(", sprintf("%1.0f", 100 * pp[-length(pp)]), 
                        ",", sprintf("%1.0f", 100 * pp[-1]), ")\n", 
                        sprintf("%1.1f", 100 * qq[-length(qq)]), " - ", 
                        sprintf("%1.1f", 100 * qq[-1]))
      }
      else names <- paste0(sprintf("%1.1f", 100 * qq[-length(qq)]), 
                           " - ", sprintf("%1.1f", 100 * qq[-1]))
    }
  }
  summary <- list(n = NROW(data))
  if (model.type %in% c("survival", "competing.risks")) 
    summary <- c(summary, list(Event = table(response[response[, 
                                                               "status"] != 0 & response[, "time"] <= time, ifelse(model.type == 
                                                                                                                     "survival", "status", "event")]), Lost = sum(response[, 
                                                                                                                                                                           "status"] == 0 & response[, "time"] <= time), Event.free = NROW(response[response[, 
                                                                                                                                                                                                                                                             "time"] > time, ])))
  out <- list(plotFrames = plotFrames, predictions = apppred, 
              time = time, cause = cause, pseudo = pseudo, summary = summary, 
              control = control, legend = legend, bars = bars, diag = diag, 
              add = add, legend = legend, names = names, method = method, 
              model.type = model.type, type = type, axes = axes, percent = percent, 
              hanging = hanging, showFrequencies = showFrequencies, 
              col = col, ylim = ylim, xlim = xlim, ylab = ylab, xlab = xlab, 
              lwd = lwd, lty = lty, pch = pch, lty = lty, NF = NF, 
              pseudo.col = pseudo.col, pseudo.pch = pseudo.pch, showPseudo = showPseudo, 
              jack.density = jack.density)
  if (method == "nne") 
    out <- c(out, list(bandwidth = sapply(plotFrames, function(x) attr(x, 
                                                                       "bandwidth"))))
  if (plot) {
    coords <- plot.calibrationPlot(out)
    out <- c(out, coords)
  }
  class(out) <- "calibrationPlot"
  invisible(out)
}
<bytecode: 0x1272bd9d0>
  <environment: namespace:pec>