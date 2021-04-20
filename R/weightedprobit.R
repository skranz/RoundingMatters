#' Function to estimate clustered weighted probit regression with clustered standard errors
# and to compute marginal effects, borrowed from function probitmfxest from mfx package. The only difference
# is that the new function defined here allows for weights
probitweighted <- function (formula, data, atmean = TRUE, robust = FALSE, clustervar1 = NULL,
            clustervar2 = NULL, start = NULL, control = list(), weights=NULL) {
    if (is.null(formula)) {
      stop("formula is missing")
    }
    if (!is.data.frame(data)) {
      stop("data arguement must contain data.frame object")
    }
    if (is.null(clustervar1) & !is.null(clustervar2)) {
      stop("use clustervar1 arguement before clustervar2 arguement")
    }
    if (!is.null(clustervar1)) {
      if (is.null(clustervar2)) {
        if (!(clustervar1 %in% names(data))) {
          stop("clustervar1 not in data.frame object")
        }
        # add weights to data frame if appropriate
        if (is.null(weights))
        {
          data = data.frame(model.frame(formula, data, na.action = NULL),
                            data[, clustervar1])
          names(data)[dim(data)[2]] = clustervar1
        } else {
        data = data.frame(model.frame(formula, data, na.action = NULL),
                          data[, clustervar1],data[, weights])
        names(data)[dim(data)[2]-1] = clustervar1
        names(data)[dim(data)[2]] = weights
        }
        data = na.omit(data)
      }
      if (!is.null(clustervar2)) {
        if (!(clustervar1 %in% names(data))) {
          stop("clustervar1 not in data.frame object")
        }
        if (!(clustervar2 %in% names(data))) {
          stop("clustervar2 not in data.frame object")
        }
        data = data.frame(model.frame(formula, data, na.action = NULL),
                          data[, c(clustervar1, clustervar2)])
        names(data)[c(dim(data)[2] - 1):dim(data)[2]] = c(clustervar1,
                                                          clustervar2)
        data = na.omit(data)
      }
    }
    # weights are allowed
    if (is.null(weights))
      data$w <- rep(1,nrow(data)) else
        data$w <- data[,weights]
    fit = glm(formula, data = data, family = binomial(link = "probit"),
              x = T, start = start, control = control, weights=w)
    x1 = model.matrix(fit)
    if (any(alias <- is.na(coef(fit)))) {
      x1 <- x1[, !alias, drop = FALSE]
    }
    xm = as.matrix(colMeans(x1))
    be = as.matrix(na.omit(coef(fit)))
    k1 = length(na.omit(coef(fit)))
    xb = t(xm) %*% be
    fxb = ifelse(atmean == TRUE, dnorm(xb), mean(dnorm(x1 %*%
                                                         be)))
    vcv = vcov(fit)
    if (robust) {
      if (is.null(clustervar1)) {
        vcv = vcovHC(fit, type = "HC0")
      }
      else {
        if (is.null(clustervar2)) {
          vcv = clusterVCV(data = data, fm = fit, cluster1 = clustervar1,
                           cluster2 = NULL)
        }
        else {
          vcv = clusterVCV(data = data, fm = fit, cluster1 = clustervar1,
                           cluster2 = clustervar2)
        }
      }
    }
    if (robust == FALSE & is.null(clustervar1) == FALSE) {
      if (is.null(clustervar2)) {
        vcv = clusterVCV(data = data, fm = fit, cluster1 = clustervar1,
                         cluster2 = NULL)
      }
      else {
        vcv = clusterVCV(data = data, fm = fit, cluster1 = clustervar1,
                         cluster2 = clustervar2)
      }
    }
    mfx = data.frame(mfx = fxb * be, se = NA)
    if (atmean) {
      gr = as.numeric(fxb) * (diag(k1) - as.numeric(xb) * (be %*%
                                                             t(xm)))
      mfx$se = sqrt(diag(gr %*% vcv %*% t(gr)))
    }
    else {
      gr = apply(x1, 1, function(x) {
        as.numeric(as.numeric(dnorm(x %*% be)) * (diag(k1) -
                                                    as.numeric(x %*% be) * (be %*% t(x))))
      })
      gr = matrix(apply(gr, 1, mean), nrow = k1)
      mfx$se = sqrt(diag(gr %*% vcv %*% t(gr)))
    }
    temp1 = apply(x1, 2, function(x) length(table(x)) == 1)
    const = names(temp1[temp1 == TRUE])
    mfx = mfx[row.names(mfx) != const, ]
    temp1 = apply(x1, 2, function(x) length(table(x)) == 2)
    disch = names(temp1[temp1 == TRUE])
    if (length(disch) != 0) {
      for (i in 1:length(disch)) {
        if (atmean) {
          disx0 = disx1 = xm
          disx1[disch[i], ] = max(x1[, disch[i]])
          disx0[disch[i], ] = min(x1[, disch[i]])
          mfx[disch[i], 1] = pnorm(t(be) %*% disx1) - pnorm(t(be) %*%
                                                              disx0)
          gr = dnorm(t(be) %*% disx1) %*% t(disx1) - dnorm(t(be) %*%
                                                             disx0) %*% t(disx0)
          mfx[disch[i], 2] = sqrt(gr %*% vcv %*% t(gr))
        }
        else {
          disx0 = disx1 = x1
          disx1[, disch[i]] = max(x1[, disch[i]])
          disx0[, disch[i]] = min(x1[, disch[i]])
          mfx[disch[i], 1] = mean(pnorm(disx1 %*% be) -
                                    pnorm(disx0 %*% be))
          gr = as.numeric(dnorm(disx1 %*% be)) * disx1 -
            as.numeric(dnorm(disx0 %*% be)) * disx0
          avegr = as.matrix(colMeans(gr))
          mfx[disch[i], 2] = sqrt(t(avegr) %*% vcv %*%
                                    avegr)
        }
      }
    }
    mfx$discretechgvar = ifelse(rownames(mfx) %in% disch, 1,
                                0)
    output = list(fit = fit, mfx = mfx)
    return(output)
  }

#' Weighted probit with marginal effects and cluster robust se
probitmfxweights <- function (formula, data, atmean = TRUE, robust = FALSE, clustervar1 = NULL, clustervar2 = NULL, start = NULL, control = list(), weights=NULL) {
  res = probitweighted(formula, data, atmean, robust, clustervar1,
                     clustervar2, start, control, weights)
  est = NULL
  est$mfxest = cbind(dFdx = res$mfx$mfx, StdErr = res$mfx$se,
                     z.value = res$mfx$mfx/res$mfx$se, p.value = 2 * pt(-abs(res$mfx$mfx/res$mfx$se),
                                                                        df = Inf))
  colnames(est$mfxest) = c("dF/dx", "Std. Err.",
                           "z", "P>|z|")
  rownames(est$mfxest) = rownames(res$mfx)
  est$mfxest %<>% as_tibble(rownames="Var")
  est$fit = res$fit
  varnames_full <- data.frame(Var=names(res$fit$coefficients)[-1])  %>% as_tibble()
  est$mfxest <- varnames_full %>% left_join(est$mfxest) %>% column_to_rownames(var = "Var")
  est$dcvar = rownames(res$mfx[res$mfx$discretechgvar == 1,
  ])
  est$call = match.call()
  class(est) = "probitmfx"
  est
}

# assign environments of the two functions above to the environment of the probitmfx function from the mfx package
#environment(probitmfxweights) <- environment(probitmfx)
#environment(probitweighted) <- environment(probitmfx)




