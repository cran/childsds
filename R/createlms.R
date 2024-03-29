##' prepare data for repeated iteration process
##'
##' given a dataframe, the column name of the subject identifier, sex, age,
##' value and group colums, the function creates a dataframe containing only
##' these five columns with the standard column names group, subject, sex, age, value.
##' lines containing missing values are removed. 
##' @title prepare data for iteration process
##' @param data dataframe containing measurement values, age, sex, and subject identifier
##' @param group optional variable indicating groups of subjects within the data frame in most cases (families)
##' @param subject subject identifier
##' @param sex column containing the sex (or any other stratum), ideally of type character, iteration process will run on each of the levels separately
##' @param value numeric column containing the measurement values
##' @param age numeric column containing the age
##' @param lb optional - lower bound for age
##' @param ub optional - upper bound for age
##' @return list of dataframes containing the columns group, subject, sex, age, value; one dataframe for every level of sex
##' @author Mandy Vogel
##' @importFrom magrittr %<>%
##' @importFrom magrittr  %>%
##' @importFrom dplyr n
##' @export
prepare_data <- function(data, group = NULL, subject = "SIC", sex = NULL, value = "value", age = "age", lb = -Inf, ub = Inf){
    if(is.null(group)){
        data$group <- NA
        group <- "group"
    }
    if(is.null(sex)){
        print("No sex variable is given. Data will not be grouped.")
        data$sex <- "all"
        sex <- "sex"
    }    
    data %<>% dplyr::select(group, subject, value, age, sex)
    data %<>% dplyr::rename("subject" = subject,
                      "group" = group,
                      "sex" = sex,
                      "value" = value,
                      "age" = age) 
    data %<>% tidyr::drop_na(subject, sex, value, age) %>% dplyr::ungroup()
    data %<>% dplyr::filter(dplyr::between(age, lb, ub))
    split(data, data$sex)
}

##' Select groups (families)
##'
##' function selects a given proportion of groups/families from the data
##' if no grouping variable is given the original data set is returned
##' function is called inside \code{\link{do_iterations}} and may not called directly
##' @title select families
##' @param data dataframe as returned by prepare data
##' @param prop proportion of families to be sampled
##' @param group name of the group variable (character) if not "group", ignored
##' @param verbose if TRUE information about sample size is printed out
##' @return dataframe containing only prop.fam percent the families in data
##' @author Mandy Vogel
##' @export
select_fams <- function(data, prop = 0.75, group, verbose = F){
    if(sum(is.na(data$group)) > 0) {
        if(verbose) print("select_fams: Missings in group variable. Returning original data set.")
        return(data)}
    weights <- dplyr::group_by(data, group) %>% dplyr::summarise(n=n(), wgt = 1-1/n()+1)
    weights <- weights$group[sample(1:nrow(weights),size = (nrow(weights) * prop), prob = weights$wgt)]
    data[data$group %in% weights,]
}


##' Choose one measurement per subject
##'
##' function samples one measurement per subject, if prop < 1 additional
##' a prop*100 percent will be sampled from the measurements
##' the function is called inside \code{\link{do_iterations}} and may not called directly
##' @title choose one measurement per subject
##' @param data dataframe as returned by prepare data
##' @param subject name of the column containing the subject identifier
##' @param prop optional - proportion of measurements to sample
##' @param verbose if TRUE information about sample size is printed out
##' @return dataframe containing the sampled rows
##' @author Mandy Vogel
##' @export
select_meas <- function(data, subject = "subject", prop = 1, verbose = F){
    n1 <- nrow(data)
    if(sum(duplicated(data$subject)) > 0 ) data %<>% dplyr::group_by(subject) %>% dplyr::sample_n(1)
    n2 <- nrow(data)
    if(prop < 1) data %<>% dplyr::ungroup() %>% dplyr::sample_frac(prop)
    n3 <- nrow(data)
    if(verbose) print(paste("sampled", n2, "lines from", n1, "lines. Choose final fraction of", prop,"makes" ,n3, "lines"))
    data
}

##' fit gamlss
##'
##' wrapper around the \code{\link[gamlss]{lms}} function in the gamlss package
##' returns the fitted lms-parameter at given age points
##' the function is called inside \code{\link{do_iterations}} and may not called directly
##' @title fit lms
##' @param data dataframe as return by select_meas()
##' @param age.min lower bound of age
##' @param age.max upper bound of age
##' @param age.int stepwidth of the age variable
##' @param keep.models indicator whether or not models in each iteration should be kept
##' @param dist distribution used for the fitting process, has to be one of BCCGo, BCPEo, BCTo as they are accepted by lms()
##' @param mu.df degree of freedem location parameter
##' @param sigma.df degree of freedem spread parameter
##' @param nu.df degree of freedem skewness parameter
##' @param tau.df degree of freedem kurtosis parameter
##' @param trans.x indicator wether age should be transformed or not
##' @param lim.trans limits for the exponent of transformation of age
##' @param value names of the value variable (character) if different from value, ignored
##' @param tmpdata ignored
##' @return list containing a dataframe of the fitted lms parameter at the given age points and the fitted model
##' @author Mandy Vogel
fit_gamlss <- function(data, age.min = 0.25, age.max = 18, age.int = 1/12, keep.models = F,
                       dist = "BCCGo", mu.df = 4,sigma.df = 3, nu.df = 2, tau.df = 2,
                       trans.x = F, lim.trans = c(0,1.5), value, tmpdata){
    tmpdata <<- dplyr::select(data, -group)
    tmpdata <- dplyr::select(data, -group)
    tr.obj <- try(mm <- gamlss::lms(value, age, data = tmpdata,
                            families = dist,method.pb = "ML", k = 2,trace = F,
                            mu.df = mu.df, sigma.df = sigma.df,
                            nu.df = nu.df, tau.df = tau.df,
                            trans.x = trans.x, lim.trans = lim.trans))
    if("try-error" %in% class(tr.obj) ){
        tr.obj <- try(mm <- gamlss::lms(value, age, data = tmpdata,
                                        families = dist,method.pb = "ML", k = 2,trace = F,
                                        mu.df = mu.df, sigma.df = sigma.df,
                                        nu.df = nu.df, tau.df = 1,
                                        trans.x = trans.x, lim.trans = lim.trans)) 
    }
    if(!exists("mm") || is.null(mm)) invisible(return(NULL)) 
    age <- seq(age.min, age.max, by = age.int)
    if (mm$family[1] == dist & !("try-error" %in% class(tr.obj))) {
        lms <- as.data.frame(gamlss::predictAll(mm,
                                        newdata = data.frame(age = age)))
        lms$age <- age
        lms %<>% dplyr::select(age, dplyr::everything())
        if(!keep.models) mm <- NULL
        return(list(lms = lms, model = mm))
    }
    invisible(return(NULL))
}

##' fit gamlss
##'
##' wrapper around the \code{\link[VGAM]{vgam}} function in the VGAM package
##' returns the fitted lms-parameter at given age points
##' the function is called inside \code{\link{do_iterations}} and may not called directly
##' @title fit lms parameters via VGAM
##' @param data dataframe as return by select_meas()
##' @param age.min lower bound of age
##' @param age.max upper bound of age
##' @param age.int stepwidth of the age variable
##' @param keep.models indicator whether or not models in each iteration should be kept
##' @param dist distribution used for the fitting process, has to be one of BCCGo, BCPEo, BCTo as they are accepted by lms()
##' @param mu.df degree of freedem location parameter
##' @param sigma.df degree of freedem spread parameter
##' @param nu.df degree of freedem skewness parameter
##' @param value names of the value variable (character) if different from value, ignored
##' @return list containing a dataframe of the fitted lms parameter at the given age points and the fitted model
##' @author mandy
fit_vgam <- function(data, age.min = 0.25, age.max = 18, age.int = 1/12, keep.models = F,
		       dist = "BCN", mu.df = 4,sigma.df = 3, nu.df = 2, value){
    tr.obj <- try(mm <-VGAM::vgam(value ~  VGAM::s(age, df = c(mu.df, sigma.df, nu.df)),
                                  data = dplyr::select(data, -group),
                                  family = switch(dist,
                                                  BCN = VGAM::lms.bcn(zero = NULL),
                                                  YJN = VGAM::lms.yjn(zero = NULL)), trace = F))
    if("try-error" %in% class(tr.obj) ){
        tr.obj <- try(mm <-VGAM::vgam(value ~  VGAM::s(age, df = c(mu.df, sigma.df, 1)), data = dplyr::select(data, -group),
                                      family = VGAM::lms.bcn(zero = NULL), trace = F))
    }
    if(!exists("mm") || is.null(mm)) invisible(return(NULL)) 
    age <- seq(age.min, age.max, by = age.int)
    if (!("try-error" %in% class(tr.obj))) {
	lms <- as.data.frame(VGAM::predict(mm,
					newdata = data.frame(age = age)))
        names(lms) <- c("nu","mu","sigma")
        lms$age <- age
        lms$sigma <- exp(lms$sigma)
	lms %<>% dplyr::select(age, dplyr::everything())
	if(!keep.models) mm <- NULL
	return(list(lms = lms, model = mm))
    }
    invisible(return(NULL))
}


##' one iteration
##'
##' function samples families then measurements and fits the model
##' the function is called inside \code{\link{do_iterations}} and may not called directly
##' @title one iteration
##' @param data.list list of dataframes as returned by prepare_data
##' @param method use vgam or gamlss
##' @param prop.fam proportion of families to be sampled
##' @param prop.subject proportion of subject to be sampled
##' @param age.min lower bound of age
##' @param age.max upper bound of age
##' @param age.int stepwidth of the age variable
##' @param keep.models indicator whether or not models in each iteration should be kept
##' @param dist distribution used for the fitting process, has to be one of BCCGo, BCPEo, BCTo as they are accepted by lms()
##' @param formula formula for the location parameter
##' @param sigma.formula formula for the sigma parameter
##' @param nu.formula formula for the nu parameter
##' @param tau.formula formula for the tau parameter
##' @param sigma.df degree of freedem spread parameter
##' @param nu.df degree of freedem skewness parameter
##' @param mu.df degree of freedem location parameter
##' @param tau.df degree of freedem kurtosis parameter
##' @param verbose whether or not information about sampling will be printed during while iterate
##' @param trans.x indicator wether age should be transformed or not
##' @param lim.trans limits for the exponent of transformation of age
##' @param method.pb GAIC or ML
##' @return list of lists each containing a dataframe of the fitted lms parameter at the given age points and the fitted model
##' @author Mandy Vogel
##' @export
one_iteration <- function(data.list, method, prop.fam = 0.75, prop.subject = 1,
                          age.min = 0, age.max = 18, age.int = 1/12,
                          keep.models = F,
                          dist = "BCCGo",
                          formula = NULL,
                          sigma.df = 3, nu.df = 2, mu.df = 4, tau.df = 2,
                          sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
                          verbose = F,trans.x = F, lim.trans = c(0,1.5),
                          method.pb = "ML"

                          ){
    tmp.l <- lapply(data.list, select_fams, prop = prop.fam, verbose = verbose)
    tmp.l <- lapply(tmp.l, select_meas, prop = prop.subject, verbose = verbose)
    switch(method,
           gamlss = lapply(tmp.l, fit_gamlss, dist = dist, keep.models = keep.models,
                           sigma.df = sigma.df, nu.df = nu.df, mu.df = mu.df, tau.df = tau.df,
                           age.min = age.min, age.max = age.max, trans.x = trans.x, lim.trans = lim.trans),
           vgam = lapply(tmp.l, fit_vgam, dist = dist, keep.models = keep.models,
                         sigma.df = sigma.df, nu.df = nu.df, mu.df = mu.df,
                         age.min = age.min, age.max = age.max),
           gamlss1 = lapply(tmp.l, fit_gamlss1, dist = dist, keep.models = keep.models,
                            sigma.formula = sigma.formula,
                            nu.formula = nu.formula,
                            tau.formula = tau.formula,
                            formula = formula, method.pb = method.pb,
                            age.min = age.min, age.max = age.max)
           )
}

##' Do lms iterations 
##'
##' function samples families, samples measurements (and subjects), fits the model for a
##' given number of iterations
##' @title do lms iterations    
##' @param n number of desired fits
##' @param max.it maximum number of iterations that will be run
##' @param verbose whether or not information about sampling will be printed during while iterate
##' @inheritParams one_iteration
##' @return list of lists for models and fitted parameters
##' @author Mandy Vogel
##' @export
do_iterations <- function(data.list, n = 10, max.it = 1000,
                          method = "gamlss", prop.fam = 0.75, prop.subject = 1,
                          age.min = 0, age.max = 18, age.int = 1/12,keep.models = F,
                          dist = "BCCGo",
                          mu.df = 4, sigma.df = 3, nu.df = 2, tau.df = 2, verbose = F,
                          formula = NULL, 
                          sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
                          method.pb = "ML",
                          trans.x = F, lim.trans = c(0,1.5)){
    range.age <- range(unlist(lapply(data.list, function(df) df$age)))
    if(age.min < (range.age[1] - 0.01*range.age[1])) {
        age.min <- find.nearest.month(range.age[1])
        print(paste("requested age.min is smaller than max age in the data. set age.min to min(age):", round(find.nearest.month(range.age[1]),2)))}
    if(age.max > (range.age[2] + 0.01*range.age[2])){
        age.max <- find.nearest.month(range.age[2])
        print(paste("requested age.max is greater than age range in data. set age.max to max(age):", round(find.nearest.month(range.age[2]),2)))
    } 
    if(sum(is.na(data.list[[1]]$group)) > 0) print("no grouping variable is given. Therefore, no grouping will be done.")
    sexes <- names(data.list)
    res <- list()
    i <- 1
    counter <- as.data.frame(t(numeric(length(sexes))))
    names(counter) <- sexes
    while(i <= max.it & min(counter[1,]) < n ){
        tmp.res <- one_iteration(data.list = data.list,
                                 method = method,
                                 dist = dist,
                                 keep.models = keep.models,
                                 mu.df = mu.df,
                                 sigma.df = sigma.df,
                                 nu.df = nu.df,
                                 tau.df = tau.df,
                                 prop.fam = prop.fam,
                                 prop.subject = prop.subject,
                                 age.min = age.min,
                                 age.max = age.max,
                                 verbose = verbose,
                                 trans.x = trans.x,
                                 lim.trans = lim.trans,
                                 formula = formula,
                                 sigma.formula = sigma.formula,
                                 tau.formula = tau.formula,
                                 nu.formula = nu.formula
                                 )
        fit_succ <- as.data.frame(lapply(tmp.res, function(x) as.numeric(!is.null(x))))
        counter <- as.data.frame(t(colSums(dplyr::bind_rows(counter,fit_succ))))
        print(fit_succ)
        print(counter)         
        res[[length(res) + 1]] <- tmp.res
        print(paste(i,"iterations done."))
        i <- i + 1 
    }
    lms <- lapply(sexes, function(sex) {
        lms <- lapply(res, function(x) x[[sex]]$lms)
        lapply(lms, function(ls) ls[sapply(ls, function(x) !is.null(x))])
    })
    names(lms) <- sexes
    print(paste("fitted out of", n, "iterations:"))
    print(sapply(lms, length))
    if( keep.models ){
        models <- lapply(sexes, function(sex) {
            models <- lapply(res, function(x) x[[sex]]$model)
            lapply(models, function(ls) ls[sapply(ls, function(x) !is.null(x))])
        })
        names(models) <- sexes
    } else {
        models <- NULL
    }
    attr(lms, "distribution") <- dist
    list(lms = lms, models = models)

}

##' aggregate lms parameters
##'
##' function takes the lms part of the result from the do_iterations() function and returns
##' the mean parameters
##' @title aggregate lms parameters
##' @param lms.list list of parameter tables as returned by do_iterations()
##' @return list of dataframes containing the aggregated parameters, each for every level of sex
##' @author Mandy Vogel
##' @export 
aggregate_lms <- function(lms.list){
    dist <- attr(lms.list, "distribution")
    columns <- names(lms.list[[1]][[1]])
    lms.agg <- lapply(lms.list, function(lms) {
        as.data.frame(sapply(columns, function(col) rowMeans(Reduce(cbind,lapply(lms, function(x) x[,col]))),
                             USE.NAMES = T))})
    attr(lms.agg, "distribution") <- dist
    lms.agg
}


##' Calculate confidence intervals
##'
##' The function takes a lms list as returned by \code{\link{do_iterations}} and
##' calculates the confidence bands for a given set of percentiles using
##' \code{\link[boot]{envelope}} from the boot package 
##' @title Calculate confidence intervals
##' @param lms.list lms part of the returned list of \code{\link{do_iterations}}
##' @param perc percentiles for which the confidence bands are calculated
##' @param level confidence level
##' @param type for now only point is a valid value
##' @return list containing the respective confidence envelopes
##' @author mandy
##' @export
calc_confints <- function(lms.list, perc = c(2.5,5,50,95,97.5), level = 0.95, type = c("point")){
    dist <- attributes(lms.list)$distribution
    nam <- paste(sprintf("perc_%02d",floor(perc)),
                 gsub("0.","", perc-floor(perc)), sep = "_") 
    perc <- perc/100
    sexes <- names(lms.list)
    lapply(sexes, function(sex){
        res.l <- list()
        for(i in 1:length(lms.list[[sex]])){
            sex.df <- lms.list[[sex]][[i]]
            perc.values <-  lapply(perc, function(p, sex.df){
                eval(parse(
                    text = paste0("gamlss.dist::q",dist,"(",p,",",
                                  paste(
                                      paste0(names(sex.df)[-which(names(sex.df) == "age")],
                                             "=sex.df$",
                                             names(sex.df)[-which(names(sex.df) == "age")]) ,
                                      collapse = ","),
                                  ")")))}, sex.df = sex.df)
            names(perc.values) <- nam
            perc.values$age <- sex.df$age
            perc.values$sex <- sex
            res.l[[length(res.l) + 1]] <- as.data.frame(perc.values, stringsAsFactors = F)                
        }
        confenv <- lapply(nam, function(perc){
            pm <- sapply(res.l, function(xx, perc){
                xx[,perc]
            }, perc = perc)
            pm <- as.data.frame(t(boot::envelope(mat = t(pm), level = level)[[type]]))
            names(pm) <- c("upper","lower")
            pm
        })
        names(confenv) <- nam
        confenv
    })
}

 

## find nearest month
##
## in the iteration process I assume age given in years
## and try to set the min and max age to full month 
## @param x - original age
## @return numeric - the nearest full month
## @author mandy
find.nearest.month <- function(x){
    full <- x %/% 1
    rem <- x %% 1
    ind <- which.min(abs(rem - 0:11/12))
    rem <- (0:11/12)[ind]
    full + rem
}


##' fit_gamlss
##' 
##' wrapper around the \code{\link[gamlss]{gamlss}} function from the gamlss package
##' returns the fitted lms-parameter at given age points
##' the function is called inside \code{\link{do_iterations}} and may not be called directly
##' @title fit_gamlss1
##' @param data dataframe as return by select_meas()
##' @param age.min lower bound of age
##' @param age.max upper bound of age
##' @param age.int stepwidth of the age variable
##' @param keep.models indicator whether or not models in each iteration should be kept
##' @param dist distribution used for the fitting process, has to be one of BCCGo, BCPEo, BCTo as they are accepted by lms()
##' @param formula formula for the location parameter
##' @param sigma.formula formula for the sigma parameter
##' @param nu.formula formula for the nu parameter
##' @param tau.formula formula for the tau parameter
##' @param method.pb GAIC or ML
##' @return list containing a dataframe of the fitted lms parameter at the given age points and the fitted model
##' @author Mandy Vogel
fit_gamlss1 <- function(data, age.min = 0, age.max = 80, age.int = 1/12, keep.models = F,
                        dist = "BCCGo", formula = NULL, 
                        sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
                        method.pb = "ML"){
    formula <- formula(formula)
    tr.obj <- try(mm <- gamlss::gamlss(formula,
                                       sigma.formula = sigma.formula,
                                       nu.formula = nu.formula,
                                       tau.formula = tau.formula,
                                       family = dist,
                                       method.pb = method.pb,
                                       k = 2,trace = F,
                                       data = data[,-grep("group",names(data))],
                                       ))
    if("try-error" %in% class(tr.obj) ){
            tr.obj <- try(mm <- gamlss::gamlss(formula,
                                               sigma.formula = sigma.formula,
                                               nu.formula = nu.formula,
                                               tau.formula = ~ 1,
                                               family = dist,
                                               method.pb = "GAIC",
                                               k = 2,trace = F,
                                               data = data[,-grep("group",names(data))],
                                               ))
    }
    if(!exists("mm") || is.null(mm)) invisible(return(NULL)) 
    age <- seq(age.min, age.max, by = age.int)
    if (mm$family[1] == dist & !("try-error" %in% class(tr.obj))) {
        lms <- as.data.frame(gamlss::predictAll(mm,
                                                newdata = data.frame(age = age),
                                                data = data
                                                ))
        lms$age <- age
        lms <- lms %>% dplyr::select(age, dplyr::everything())
        if(!keep.models) mm <- NULL
        return(list(lms = lms, model = mm))
    }
    invisible(return(NULL))
}
