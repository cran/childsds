##' Calculate raw values for percentile curve
##'
##' calculates quantile values for given RefGroup and given
##' percentiles
##' @title calculate raw values 
##' @param ref Refgroup object
##' @param item name of the measurement item
##' @param perc vector of percentiles to be calculated
##' @param stack wether or not the data should be stacked, stacked data
##' would most possibly be used in ggplot2
##' @param age desired values of age 
##' @param include.pars indicator whether or not parameters should be included
##' @param digits specification of number of decimal places
##' @param sex name of the sex variable (character) if different from sex, not
##' functional in this version and therefore ignored
##' @return data frame either with the different percentiles as columns
##' or, if stacked, as data frame with four columns: age, sex, variable, value
##' @author Mandy Vogel
##' @examples
##' ptab <- make_percentile_tab(ref = kro.ref,
##'                             item = "height",
##'                            perc = c(2.5,10,50,90,97.5),
##'                            stack = TRUE)
##'
##' ggplot2::ggplot(ptab, ggplot2::aes(x = age, y = value, colour = variable)) +
##'    ggplot2::geom_line() +
##'    ggplot2::facet_wrap(~ sex, nrow = 2)
##' @export
make_percentile_tab <- function (ref, item, perc = c(2.5, 5, 50, 95, 97.5), stack = F,
                                 age = NULL, include.pars = T, digits = 4, sex ) {
    if(stack) include.pars <- F
    dists <- unlist(ref@refs[[item]]@dist)    
    refs <- ref@refs[[item]]@params
    sexes <- names(refs)
    if(!is.null(age) & !inherits(age,"factor" )) {
        refs <- purrr::map(refs, function(ref){
            age <- age[dplyr::between(age, min(ref$age), max(ref$age))]
            purrr::map(ref,
                       function(param) stats::approx(x = ref$age,
                                                     y = param,
                                                     xout = age,
                                                     rule = 1)$y)
        })
    } 
    nam <- paste(sprintf("perc_%02d", floor(perc)),
                 gsub("0.", "", perc - floor(perc)), sep = "_")
    perc <- perc/100
    pertab <- purrr::map2(refs, sexes, function(ref,sex) {
        params <- base::as.list(ref)
        params$age <- NULL        
        res <- purrr::map(perc,function(p) {
            params$p <- p 
            do.call(get(paste0("q",dists[sex]), asNamespace("gamlss.dist")), params)
        })
        names(res) <- nam
        res <- as.data.frame(res)
        if(include.pars)
            res <- dplyr::bind_cols(res, round(as.data.frame(ref[-1]), digits))
        res$sex <- sex
        res$age <- ref$age        
        res 
    })        
    res <- Reduce(dplyr::bind_rows, pertab)
    if (requireNamespace("reshape2") & stack){
        return(reshape2::melt(res, id.vars = c("age", "sex")))
    } 
        
    if (!requireNamespace("reshape2") & stack) 
        print("For stacking the package reshape2 is required")
    return(dplyr::select(res, sex, age, dplyr::everything()))
}

##' Calculate raw values for percentile curve
##'
##' calculates quantile values for given RefGroup and given
##' percentiles
##' @title calculate raw values 
##' @param ref Refgroup object
##' @param item name of the measurement item
##' @param perc vector of percentiles to be calculated
##' @param stack wether or not the data should be stacked, stacked data
##' would most possibly be used in ggplot2
##' @param age desired values of age
##' @param x2 desired values of second covariate
##' @param include.pars indicator whether or not parameters should be included
##' @param digits specification of number of decimal places
##' @param sex name of the sex variable (character) if different from sex, not
##' functional in this version and therefore ignored
##' @return data frame either with the different percentiles as columns
##' or, if stacked, as data frame with four columns: age, sex, variable, value
##' @author Mandy Vogel
##' @examples
##' ptab <- make_percentile_tab(ref = kro.ref,
##'                             item = "height",
##'                            perc = c(2.5,10,50,90,97.5),
##'                            stack = TRUE)
##'
##' ggplot2::ggplot(ptab, ggplot2::aes(x = age, y = value, colour = variable)) +
##'    ggplot2::geom_line() +
##'    ggplot2::facet_wrap(~ sex, nrow = 2)
##' @export
make_percentile_tab2d <- function (ref, item,
                                   perc = c(2.5, 5, 50, 95, 97.5), stack = F,
                                   age = NULL, x2 = NULL,
                                   include.pars = T, digits = 4, sex ) {
    if(is.null(age ) | is.null(x2 ) ) {
        stop("age and x2 are not allowed to be NULL" )
    }
    if(stack) include.pars <- F
    dists <- unlist(ref@refs[[item]]@dist)    
    refs <- ref@refs[[item]]@params
    sexes <- names(refs)
    if(!is.null(age)) {
        refs <- purrr::map(refs, function(ref){
            age <- age[dplyr::between(age, min(ref$age), max(ref$age))]
            tmpdf <- expand.grid(age = age,
                                 x2 = x2)
            purrr::map(ref,
                       function(param) interp::interp(x = ref$age,
                                                     y = ref$x2,
                                                     z = param,
                                                     xo = tmpdf$age,
                                                     yo = tmpdf$x2,
                                                     output = "points")$z
                       )
        })
    } 
    nam <- paste(sprintf("perc_%02d", floor(perc)),
                 gsub("0.", "", perc - floor(perc)), sep = "_")
    perc <- perc/100
    pertab <- purrr::map2(refs, sexes, function(ref,sex) {
        params <- base::as.list(ref)
        params$age <- NULL
        params$x2 <- NULL                
        res <- purrr::map(perc,function(p) {
            params$p <- p 
            do.call(get(paste0("q",dists[sex]), asNamespace("gamlss.dist")), params)
        })
        names(res) <- nam
        res <- as.data.frame(res)
        if(include.pars)
            res <- dplyr::bind_cols(res, round(as.data.frame(ref[-1]), digits))
        res$sex <- sex
        res$age <- ref$age
        res$x2 <- ref$x2
        res 
    })        
    res <- Reduce(rbind, pertab)
    if (requireNamespace("reshape2") & stack){
        return(reshape2::melt(res, id.vars = c("age", "sex","x2")))
    } 
        
    if (!requireNamespace("reshape2") & stack) 
        print("For stacking the package reshape2 is required")
    return(dplyr::select(res, sex, age, x2, dplyr::everything()))
}




##' Worm plot ggplot version
##'
##' creates a wormplot for a gamlss model or a given vector of
##' normalized quantile residuals, either for all residuals or
##' grouped by age intervals
##' @title Worm Plot ggplot version
##' @param m a gamlss model
##' @param residuals nlormalized quantile residuals
##' @param age numeric vector of ages
##' @param n.inter number of age intervals or cut points
##' @param y.limits limits of the y-axis
##' @return ggplot object
##' @export 
wormplot_gg <- function(m=NULL, residuals=NULL, age=NA, n.inter=1, y.limits = c(-1,1)){
    if(inherits(m,"gamlss"))
        residuals <- residuals(m)
    mm <- tibble::tibble(x = seq(-4,4,l = 1000),
                         yu = 1.96 * sqrt(stats::pnorm(.data$x)*(1-stats::pnorm(.data$x))/length(residuals))/stats::dnorm(.data$x),
                         yl = -.data$yu)
    mm <- reshape2::melt(mm, id.var = "x")
    tmp <- data.frame(residuals = residuals,
                      age = age)
    if(n.inter > 1) {
        if(all(is.na(age))) stop("intervals only possible if a vector of ages  is given")
        tmp$ag <- cut(tmp$age, breaks = stats::quantile(tmp$age,
                                                        probs = (0:n.inter)/n.inter),
                      include.lowest = T,
                      na.rm = T)
    } else {
        tmp$ag <- "all ages"

    }
    tmp <- dplyr::group_by(tmp, .data$ag) %>% tidyr::nest()
    tmp <- dplyr::mutate(tmp, qq = purrr::map(.data$data, function(x){
        qq <- as.data.frame(stats::qqnorm(x$residuals, plot.it = F))
        qq$y <- qq$y - qq$x
        qq
    }))
    tmp <- tidyr::unnest(tmp, .data$qq)
    ggplot2::ggplot(tmp, ggplot2::aes_string(x = "x", y = "y")) +
        ggplot2::geom_point(shape = 21, size = 1, colour = "#2171B5") +
        ggplot2::geom_line(data = mm, inherit.aes = F,
                  ggplot2::aes_string(x = "x", y = "value", group = "variable"),
                  linetype = 2, colour = "forestgreen") +
        ggplot2::geom_vline(xintercept = 0, colour = "firebrick3", linetype = 4) +
        ggplot2::scale_y_continuous(limits = y.limits) +
        ggplot2::facet_wrap( ~ag) +
        ggplot2::labs(x = "Unit normal quantiles",
             y = "Deviation") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
              axis.title = ggplot2::element_text(size = 13, colour ="black"),
              axis.text = ggplot2::element_text(size = 13, colour ="black"),
              title = ggplot2::element_text(colour = "black"))
}



##' Worm plot ggplot version, 2 covariates
##'
##' creates a wormplot for a gamlss model or a given vector of
##' normalized quantile residuals, either for all residuals or
##' grouped by age intervals
##' @title Worm Plot ggplot version
##' @param m a gamlss model
##' @param residuals nlormalized quantile residuals
##' @param age numeric vector of ages
##' @param x2 numeric vector of second covariate
##' @param n.inter.age number of age intervals or cut points
##' @param n.inter.x2 number of x2 intervals or cut points
##' @param name.x2 name of x2 for use in graphics
##' @param y.limits limits of the y-axis
##' @return ggplot object
##' @export 

wormplot2d_gg <- function(m=NULL, residuals=NULL, age=NA, x2 = NA, name.x2, n.inter.age=1, n.inter.x2 = 1, y.limits = c(-1,1)){
    if(inherits(m,"gamlss"))
        residuals <- residuals(m)
    mm <- tibble::tibble(x = seq(-4,4,l = 1000),
                         yu = 1.96 * sqrt(stats::pnorm(.data$x)*(1-stats::pnorm(.data$x))/length(residuals))/stats::dnorm(.data$x),
                         yl = -.data$yu)
    mm <- reshape2::melt(mm, id.var = "x")
    tmp <- data.frame(residuals = residuals,
                      age = age,
                      x2 = x2)
    if(n.inter.age > 1) {
        if(all(is.na(age))) stop("intervals only possible if a vector of ages is given")
        tmp$ag <- cut(tmp$age, breaks = stats::quantile(tmp$age,
                                                 probs = (0:n.inter.age)/n.inter.age,
                                                 include.lowest = T,
                                                 na.rm = T))
    } else {
        tmp$ag <- "all ages"
    }
    if(n.inter.x2 > 1) {
        if(all(is.na(age))) stop("intervals only possible if a vector x2 is given")
        tmp$x2g <- cut(tmp$x2, breaks = stats::quantile(tmp$x2,
                                                 probs = (0:n.inter.x2)/n.inter.x2,
                                                 include.lowest = T,
                                                 na.rm = T))
    } else {
        tmp$x2g <- "all x2"
    }
    tmp <- dplyr::filter(tmp, !is.na(.data$age ) & !is.na(.data$x2g))
    tmp <- dplyr::group_by(tmp, .data$ag, .data$x2g) %>% tidyr::nest()
    tmp <- dplyr::mutate(tmp, qq = purrr::map(.data$data, function(x){
        qq <- as.data.frame(stats::qqnorm(x$residuals, plot.it = F))
        qq$y <- qq$y - qq$x
        qq
    }))
    tmp <- tidyr::unnest(tmp, .data$qq)
    ggplot2::ggplot(tmp, ggplot2::aes_string(x = "x", y = "y")) +
        ggplot2::geom_point(shape = 21, size = 1, colour = "#2171B5") +
        ggplot2::geom_line(data = mm, inherit.aes = F,
                  ggplot2::aes_string(x = "x", y = "value", group = "variable"),
                  linetype = 2, colour = "forestgreen") +
        ggplot2::geom_vline(xintercept = 0, colour = "firebrick3", linetype = 4) +
        ggplot2::scale_y_continuous(limits = y.limits) +
        ggplot2::facet_wrap( ~ag + x2g) +
        ggplot2::labs(x = "Unit normal quantiles",
             y = "Deviation") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none",
              axis.title = ggplot2::element_text(size = 13, colour ="black"),
              axis.text = ggplot2::element_text(size = 13, colour ="black"),
              title = ggplot2::element_text(colour = "black"))
}


##' dummy function to call a function of the gamlss.dist package
##'
##' dummy function to call a function of the gamlss.dist package
##' @title dummy
##' @param ref reference object
##' @param item item name specifiying the reference
##' @return dummy output
##' @author Mandy Vogel
dummy.refs <- function(ref, item){
    minmax <- purrr::map(ref@refs[[item]]@params, ~.[c(1,nrow(.)),])
    dists <- purrr::map(ref@refs[[item]]@dist, ~.)
    purrr::map2(dists, minmax, function(d,dt){
        fct <- ifelse(base::grepl("BCT", d), gamlss.dist::pBCT,
                      ifelse(base::grepl("BCPE", d), gamlss.dist::pBCPE))
        try(fct(q = dt$mu, mu = dt$mu, sigma = dt$sigma, tau = dt$tau))
    })
}



