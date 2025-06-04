## The malp object and its methods
####################################

malp <- function(formula, data, obj=TRUE)
{
    cl <- match.call()
    mf <- model.frame(formula, data)
    if (attr(terms(mf), "intercept")!=1)
        stop("The model must include an intercept")
    Y <- model.response(mf)
    X <- model.matrix(mf, data)
    omit <- attr(na.omit(cbind(Y,X)), "na.action")
    if (!is.null(omit))
    {
        Y <- Y[-omit]
        X <- X[-omit,,drop=FALSE]
        data <- data[-omit,,drop=FALSE]
        if (length(Y)==0)
            stop("No observations left after removing missing values")
    }
    muY <- mean(Y)    
    if (!obj)
    {
        fit <- lm.fit(X, Y)
        b <- fit$coefficients
        Yhat <- fit$fitted.values
        gamma <- cor(Y, Yhat)
        b2 <- c((1-1/gamma)*muY+b[1]/gamma,
                b[-1]/gamma)
        names(b2) <- names(b)
        return(b2)
    }
    fit <- lm(formula, data)
    gamma <- cor(Y, fitted(fit))
    alpha1 <- (1-1/gamma)*muY+coef(fit)[1]/gamma
    alpha2 <- coef(fit)[-1]/gamma
    coef <- c(alpha1, alpha2)
    names(coef) <- names(coef(fit))
    obj <- list(coefficients=coef, gamma=gamma, lm=fit, varY=var(Y),
                muY=muY, call=cl, data=data, na.action=omit)
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

ccc <- function(x,y, type="Lin")
{
    type <- match.arg(type)
    sx <- sd(x)
    sy <- sd(y)
    if (type=="Lin")
        2*cor(x,y)*sx*sy/(sx^2+sy^2+(mean(x)-mean(y))^2)
}

predict.malp <- function (object, newdata = NULL, se.fit = FALSE,
                          interval = c("none", "confidence", "prediction"),
                          level = 0.95,
                          includeLS = FALSE, LSdfCorr = FALSE, 
                          vcovMet = c("Asymptotic", "Normal", "Boot", "Jackknife"),
                          bootInterval=FALSE,
                          bootIntType=c("all", "norm", "basic",
                                        "stud", "perc", "bca"),
                          Bse.=100, B.=300,
                          parallel = c("no", "multicore", "snow"),
                          ncpus = getOption("boot.ncpus", 1L), cl = NULL, ...) 
{
    vcovMet <- match.arg(vcovMet)
    interval <- match.arg(interval)
    parallel <- match.arg(parallel)
    bootIntType <- match.arg(bootIntType)
    if (bootInterval)
    {
        if (interval=="prediction")
            warning("interval='prediction' is ignored for bootstrap intervals.")
        res <- bootMALP(object, newdata, B., Bse., se.fit, vcovMet,
                        parallel, ncpus, cl)
        ci <- confint(res, level=level, type.=bootIntType)
        return(ci)
    }
    pr <- predict(object$lm, newdata, se.fit = TRUE)
    pr2 <- coef(object)[1] + pr$fit/object$gamma - coef(object$lm)[1]/object$gamma
    n <- nobs(object$lm)
    df <- object$lm$df.residual
    if ((interval == "none") & !se.fit)
    {
        if (includeLS)
            return(list(MALP = pr2, LSLP = pr$fit))
        else
            return(pr2)
    }
    sigLS <- pr$se.fit^2 * n
    if (vcovMet == "Normal") {
        if (!LSdfCorr) 
            sigLS <- sigLS * df/n + object$varY * (1 - object$gamma^2)/n
        g2 <- object$gamma^2
        tmp <- 2 * g2/(1 + object$gamma) - 1 - (1 - g2)/(object$varY * 
            g2) * (pr$fit - object$muY)^2
        sigMA <- (sigLS + object$varY * (1 - g2) * tmp)/g2
    } else {
        tt <- terms(object$lm)
        Terms <- delete.response(tt)
        if (is.null(newdata)) {
            m <- model.frame(Terms, object$data, xlev = object$lm$xlevels)
        }
        else {
            m <- model.frame(Terms, newdata, xlev = object$lm$xlevels)
        }
        X <- model.matrix(Terms, m)
        V <- vcov(object, method = vcovMet, B=Bse.)
        sigMA <- apply(X, 1, function(x) c(t(x) %*% V %*% x) * 
            nobs(object$lm))
    }
    if (interval == "prediction")
    {
        sige <- object$varY*(1-object$gamma^2)
        if (LSdfCorr)
            sige <- (n-1)/df*sige
    } else {
        sige <- 0
    }
    sigMA <- sqrt(sigMA/n+sige)
    sigLS <- sqrt(sigLS/n+sige)
    if (interval != "none") {
        crit <- qnorm(0.5 + level/2)
        pr <- cbind(fit = pr$fit, lwr = pr$fit - crit * sigLS, 
                    upr = pr$fit + crit * sigLS)
        if (interval == "prediction")
            pr3 <- pr[,1]
        else
            pr3 <- pr2
        pr2 <- cbind(fit = pr2, lwr = pr3 - crit * sigMA, upr = pr3 + 
            crit * sigMA)
    } else {
        pr <- pr$fit
    }
    if (!se.fit)
    {
        if (includeLS)
            return(list(MALP = pr2, LSLP = pr))
        else
            return(pr2)
    }
    if (includeLS)
        list(MALP = list(fit = pr2, se.fit = sigMA),
             LSLP = list(fit = pr, se.fit = sigLS))
    else 
        list(fit = pr2, se.fit = sigMA)
}

plot.malp <- function (x, y=NULL, which=c("MALP", "LSLP", "Both"),
                       pch=21:22, col=2:3, bg=2:3, ...)
{
    which <- match.arg(which)
    yhat <- predict(x)
    yhat2 <- predict(x$lm)
    y <- model.response(model.frame(x$lm))
    if (which %in% c("MALP", "Both"))
        plot(y, yhat, pch=pch[1], col=col[1], bg=bg[1], ...)
    else
        plot(y, yhat2, pch=pch[2], col=col[2], bg=bg[2], ylab="yhat", ...)
    abline(0,1, lty=2, lwd=2)
    if (which == "Both")
    {
        points(y, yhat2, pch=pch[2], col=col[2], bg=bg[2], ...)
        legend("topleft", c("MALP","LSLP"), pch=pch, col=col,
               pt.bg=bg, bty='n', lty=NULL)
    }
    grid()
    invisible()
}

.Ximat <- function(obj)
{
    X <- model.matrix(obj$lm)[,-1, drop=FALSE]
    Y <- model.response(model.frame(obj$lm))
    T <- cbind(Y, X)
    X <- scale(X, scale=FALSE)
    Y <- Y-mean(Y)
    T <- cbind(T, Y^2, X*Y)
    if (ncol(X)==1)
        T2 <- c(X)^2
    else
        T2 <- t(sapply(1:nrow(X), function(i) c(X[i,]%*%t(X[i,]))))
    T <- cbind(T, T2)
    dimnames(T) <- NULL
    var(T)
}

vcov.malp <- function(object, method=c("Asymptotic", "Normal", "Boot", "Jackknife"), B=400,
                      LSdfCorr = FALSE, ...)
{
    method <- match.arg(method)
    if (method=="Normal")
    {
        n <- nobs(object$lm)
        dfC1 <- object$lm$df.residual/(n-1)
        dfC2 <- object$lm$df.residual/n
        if (LSdfCorr)
            dfC1 <- dfC2 <- 1
        Omega <- vcov(object$lm)[-1,-1]*dfC2
        g <- object$gamma
        b <- coef(object$lm)[-1]
        Ybar0 <- object$muY-coef(object$lm)[1]
        Xbar <- colMeans(model.matrix(object$lm))[-1]
        sig2 <- sum(residuals(object$lm)^2)/object$lm$df.residual
        tmp <- (1-g^2)^2/(n*g^4)
        tmp2 <- c(crossprod(Xbar, Omega))
        V11 <- 2*sig2/(n*(1+g))*dfC1 - tmp*Ybar0^2 + sum(tmp2*Xbar)/g^2
        V12 <- -tmp2/g^2+tmp*Ybar0*b
        V22 <- Omega/g^2-tmp*(b%*%t(b))
        V <- cbind(c(V11,V12), rbind(c(V12), V22))
        dimnames(V) <- list(names(coef(object$lm)), names(coef(object$lm)))
        return(V)
    }
    if (method=="Asymptotic")
    {
        Xi <- .Ximat(object)
        p <- length(object$coef)-1
        Xi11 <- Xi[1:(p+1), 1:(p+1)]
        Xi22 <- Xi[-(1:(p+1)), -(1:(p+1))]
        Xi21 <- Xi[-(1:(p+1)), 1:(p+1)]
        n <- nobs(object$lm)
        b1 <- coef(object$lm)[-1]
        g <- object$gamma
        muX <- colMeans(model.matrix(object$lm)[,-1, drop=FALSE])
        sigY <- object$varY
        sigXinv <- chol2inv(object$lm$qr$qr)[-1,-1, drop=FALSE]*n
        A <- rbind(t(b1)/(2*g*sigY),
                   sigXinv/g-b1%*%t(b1)/(g^3*sigY),
                   -kronecker(b1, sigXinv)/g+c(b1%*%t(b1))%*%t(b1)/(2*g^3*sigY))
        B <- c(1, -b1/g)
        V <- matrix(0, p+1, p+1)
        dimnames(V) <- list(names(coef(object$lm)), names(coef(object$lm)))
        AA <- crossprod(A,Xi22)%*%A
        BB <- 2*sigY*(1-g)
        AB <- crossprod(A, Xi21)%*%B
        V[2:(p+1), 2:(p+1)] <- AA
        V[2:(p+1),1] <- AB - AA%*%muX
        V[1, 2:(p+1)] <- t(V[2:(p+1),1])
        V[1, 1] <- c(crossprod(muX, AA)%*%muX) -2*c(t(muX)%*%AB) + BB
        return(V/n)
    }
    
    n <- nrow(object$data)            
    B <- ifelse(method == "Boot", B, n)
    alpha <- sapply(1:B, function(i) {
        ind <- if (method=="Boot") {
                   sample(nobs(object$lm), replace=TRUE)
               } else { -i}
        form <- formula(object$lm)
        environment(form) <- environment()
        malp(form, data = object$data[ind, ], obj=FALSE)        
    })
    V <- var(t(alpha))
    if (method == "Jackknife")
        V <- V*(n-1)^2/n
    V
}

## The summary.malp object and its methods

.goodFit <- function(object, which=c("malp","lslp"))
{
    which <- match.arg(which)
    y <- model.response(model.frame(object$lm))
    yhat <- if (which == "malp") predict(object) else predict(object$lm)
    CCC <- ccc(y, yhat)
    PCC <- cor(y,yhat)
    MSE <- mean((y-yhat)^2)
    c(PCC=PCC, CCC=CCC, MSE=MSE)
}


summary.malp <- function(object, vcovMet=c("Asymptotic", "Normal", "Boot", "Jackknife"),
                         se=TRUE, LSdfCorr=FALSE, ...)
{
    vcovMet <- match.arg(vcovMet)
    b <- coef(object)    
    if (se)
    {
        V <- vcov(object, method=vcovMet, LSdfCorr=LSdfCorr, ...)
        se <- sqrt(diag(V))
        t <- b/se
        pv <- 2*pnorm(-abs(t))
    } else {
        se <- t <- pv <- rep(NA, length(b))
    }
    coefs <- cbind(b, se, t, pv)
    dimnames(coefs) <- list(names(b),
                            c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    fm <- .goodFit(object, "malp")
    fl <- .goodFit(object, "lslp")
    ans <- list(fitMALP=fm, fitLSLP=fl, coefficients=coefs,
                call=object$call)
    class(ans) <- "summary.malp"
    ans
}

print.summary.malp <- function(x, digits=5,
                               signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (all(is.na(x$coefficients[,2])))
        cat("(Summary without standard errors)\n")
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    cat("\nPCC: ", formatC(x$fitMALP["PCC"], digits = digits), "\n", sep="")
    cat("CCC: ", formatC(x$fitMALP["CCC"], digits = digits), "\n", sep="")
    cat("MSE: ", formatC(x$fitMALP["MSE"], digits = digits), "\n", sep="")
}

bootMALP <- function (object, newdata = NULL, B=300, Bse=100, se.fit=FALSE,
                      vcovMet = c("Asymptotic", "Boot", "Jackknife"),
                      parallel = c("no", "multicore", "snow"),
                      ncpus = getOption("boot.ncpus", 1L), cl = NULL) 
{
    vcovMet <- match.arg(vcovMet)
    parallel <- match.arg(parallel)
    .bootPr <- function(data, i, obj, newdata.=NULL, vcovMet., Bse., se.fit.)
    {
        fitB <- update(obj, data=data[i,])
        pr <- predict(fitB, newdata=newdata., se.fit=se.fit.,
                      vcovMet=vcovMet., Bse.=Bse.)
        if (se.fit.)
            c(pr$fit, pr$se.fit^2)
        else pr
    }
    res <- boot(object$data, .bootPr, R=B, newdata.=newdata, vcovMet.=vcovMet,
                Bse.=Bse, obj=object, se.fit.=se.fit,
                parallel=parallel, ncpus=ncpus, cl=cl)
    res$se.fit <- se.fit
    class(res) <- "bootMALP"
    res
}

confint.bootMALP <- function(object, parm, level=0.95,
                             type.=c("norm", "basic", "stud",
                                     "perc", "bca", "all"), ...)
{
    type. <- match.arg(type.)
    if (type.=="stud" & !object$se.fit)
        stop("The bootMALP object must be computed with se.fit=TRUE for studentized intervals")
    if (type. == "all")
    {
        type. <- c("norm", "basic", "stud", "perc", "bca")
        if (!object$se.fit)
            type. <- type.[-3]
    }
    n <- ifelse(object$se.fit, ncol(object$t)/2, ncol(object$t))
    all <- lapply(type., function(ti) {
        conf <- sapply(1:n, function(i)
        {
            ind <- if (object$se.fit) c(i,i+n)
                   else i
            res <- try(boot.ci(object, type=ti, conf=level, index=ind),
                       silent=TRUE)
            if (inherits(res, "try-error"))
                c(object$t0[i], rep(NA, 2))
            else
                c(res[[2L]], tail(res[[4L]][1,],2))
        })
        rownames(conf) <- c("fit","lower","upper")
        t(conf)})
    names(all) <- type.
    all
}
