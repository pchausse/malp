## The malp object and its methods
####################################

malp <- function(formula, data)
{
    cl <- match.call()
    mf <- model.frame(formula, data)
    if (attr(terms(mf), "intercept")!=1)
        stop("The model must include an intercept")
    Y <- model.response(mf)
    muY <- mean(Y)
    fit <- lm(formula, data)
    gamma <- cor(Y, fitted(fit))
    alpha1 <- (1-1/gamma)*muY+coef(fit)[1]/gamma
    alpha2 <- coef(fit)[-1]/gamma
    coef <- c(alpha1, alpha2)
    names(coef) <- names(coef(fit))
    obj <- list(coefficients=coef, gamma=gamma, lm=fit, varY=var(Y),
                muY=muY, call=cl, data=data)
    class(obj) <- "malp"
    obj
}

coef.malp <- function(object, ...)
    object$coefficients

print.malp <- function (x, ...) 
{
    class(x) <- "lm"
    print(x, ...)
}

predict.malp <- function(object, newdata=NULL, se.fit = FALSE,
                         interval = c("none","confidence"), level=0.95, 
                         includeLS=FALSE, LSdfCorr=FALSE,
                         vcovMet=c("Asymptotic", "Boot", "Jacknife"), ...)
{
    vcovMet <- match.arg(vcovMet)
    interval <- match.arg(interval)
    pr <- predict(object$lm, newdata, se.fit=TRUE)
    pr2 <- coef(object)[1]+pr$fit/object$gamma-coef(object$lm)[1]/object$gamma
    n <- nobs(object$lm)
    df <- object$lm$df.residual
    if ((interval=="none") & !se.fit)
        return(list(MALP=pr2, LSLP=pr$fit))
    sigLS <- pr$se.fit^2*n 
    if (vcovMet == "Asymptotic")
        {
            if (!LSdfCorr)
                sigLS <- sigLS*df/n +
                    object$varY*(1-object$gamma^2)/n
            g2 <- object$gamma^2
            tmp <- 2*g2/(1+object$gamma)-1 -
                (1-g2)/(object$varY*g2)*(pr$fit-object$muY)^2
            sigMA <- (sigLS+object$varY*(1-g2)*tmp)/g2
        } else {
            tt <- terms(object$lm)
            Terms <- delete.response(tt)
            if (is.null(newdata)) {
                model.frame(Terms, object$data, xlev = object$lm$xlevels)
            } else {
                m <- model.frame(Terms, newdata, xlev = object$lm$xlevels)
            }
            X <- model.matrix(Terms, m)
            V <- vcov(object, method=vcovMet, ...)
            sigMA <- apply(X, 1, function(x)
                c(t(x)%*%V%*%x)*nobs(object$lm))
        }
    if (interval == "confidence")
    {
        crit <- qnorm(0.5+level/2)
        pr <- cbind(fit=pr$fit, lwr=pr$fit-crit*sigLS, upr=pr$fit+crit*sigLS)
        pr2 <- cbind(fit=pr2, lwr=pr2-crit*sigMA, upr=pr2+crit*sigMA)
    } else {
        pr <- pr$fit
    }
    if (!se.fit)        
        return(list(MALP=pr2, LSLP=pr))
    list(MALP=list(fit=pr2, se.fit=sigMA),
         LSLP=list(fit=pr, se.fit=sigLS))
}

plot.malp <- function (x, y=NULL, includeLS=FALSE,
                       pch=21:22, col=2:3, bg=2:3, ...)
{
    yhat <- predict(x)$MALP
    y <- model.response(model.frame(x$lm))
    plot(yhat, y, pch=pch[1], col=col[1], bg=bg[1], ...)
    abline(1,1)
    if (includeLS)
    {
        yhat2 <- predict(x)$LSLP
        points(yhat2, y, pch=pch[2], col=col[2], bg=bg[2], ...)
        legend("topleft", c("MALP","LSLP"), pch=pch, col=col,
               pt.bg=bg, bty='n', lty=NULL)
    }
    grid()
    invisible()
}

vcov.malp <- function(object, method=c("Boot", "Jackknife"), B=400, ...)
{
    method <- match.arg(method)
    n <- nrow(object$data)
    B <- ifelse(method == "Boot", B, n)
    alpha <- sapply(1:B, function(i) {
        ind <- if (method=="Boot") {
                   sample(nobs(object$lm), replace=TRUE)
               } else { -i}
        coef(update(object, data=object$data[ind,]))
    })
    V <- var(t(alpha))
    if (method == "Jackknife")
        V <- V*(n-1)^2/n
    V
}

## The summary.malp object and its methods

summary.malp <- function(object, vcovMet=c("Boot", "Jackknife"), ...)
{
    vcovMet <- match.arg(vcovMet)
    V <- vcov(object, method=vcovMet, ...)
    se <- sqrt(diag(V))
    b <- coef(object)
    t <- b/se
    pv <- 2*pnorm(-abs(t))
    coefs <- cbind(b, se, t, pv)
    dimnames(coefs) <- list(names(b),
                            c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    yhat <- predict(object)$MALP
    y <- model.response(model.frame(object$lm))
    CCC <- 2*cov(y, yhat)/(var(yhat)+object$varY+(mean(yhat)-object$muY)^2)
    PCC <- cor(y,yhat)
    ans <- list(CCC=CCC, PCC=PCC, coefficients=coefs,
                call=object$call)
    class(ans) <- "summary.malp"
    ans
}

print.summary.malp <- function(x, digits=5,
                               signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    cat("\nCCC: ", formatC(x$CCC, digits = digits), "\n", sep="")
    cat("PCC: ", formatC(x$PCC, digits = digits), "\n", sep="")
}
