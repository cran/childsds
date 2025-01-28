##' Calculate SDS values
##'
##' The function takes a vector of measurement values, and of age and of sex
##' and a RefGroup object as arguments. It calculates the sds or percentile
##' values.
##' @title Calculate SDS Values
##' @param value vector of measurement values
##' @param age vector of age values
##' @param sex vector of sex 
##' @param item name of the item e.g. "height"
##' @param ref RefGroup object
##' @param type "SDS" or "perc"
##' @param male coding of sex for male
##' @param female coding of sex for female
##' @return vector containing SDS or percentile values
##' @examples
##' anthro <- data.frame(age = c(11.61,12.49,9.5,10.42,8.42,10.75,9.57,10.48),
##'                      height = c(148.2,154.4,141.6,145.3,146,140.9,145.5,150),
##'                      sex = sample(c("male","female"), size = 8, replace = TRUE),
##'                      weight = c(69.5,72.65,47.3,51.6,45.6,48.9,53.5,58.5))
##' anthro$height_sds <- sds(anthro$height,
##'                          age = anthro$age,
##'                          sex = anthro$sex, male = "male", female = "female",
##'                          ref = kro.ref,
##'                          item = "height",
##'                          type = "SDS")
##' 
##' anthro$bmi <- anthro$weight/(anthro$height**2) * 10000
##' anthro$bmi_perc <- sds(anthro$bmi,
##'                        age = anthro$age,
##'                        sex = anthro$sex, male = "male", female = "female",
##'                        ref = kro.ref,
##'                        item = "bmi",
##'                        type = "perc")
##' data(who.ref)
##' x <- data.frame(height=c(50,100,60,54),
##'                 sex=c("m","f","f","m"),
##'                 age=c(0,2.9,0.6,0.2))
##' sds(value = x$height, age = x$age, sex = x$sex, male = "m", female = "f",
##'     ref = who.ref, item = "height")
##' @author Mandy Vogel
##' @export
sds <- function(value, age, sex, item, ref, type = "SDS", male = "male", female = "female"){
    if(!(length(value) == length(age) & length(value) == length(sex))){
        print("value, age, and sex must be of the same length")
        invisible(return(NULL))
    }
    if(!item %in% names(ref@refs)){
        stop(paste(item,"is not available in the given refs object.",
                   "Please choose one of the following items:",
                   paste(names(ref@refs), collapse = " - ")))
        return(invisible(NULL))
    }
    value <- tidyr::replace_na(value, Inf)
    sex <- as.character(factor(sex, levels = c(male, female), labels = c("male", "female")))
    refs <- ref@refs[[item]]@params
    dists <- ref@refs[[item]]@dist
    tmpdf <- data.frame(age = age, sex = sex, value = value)
    params <- purrr::map2(refs, names(refs), function(ref,sex){
        mm <- purrr::map(ref,
                         function(param) stats::approx(ref$age,
                                                       param, xout = tmpdf$age, rule = 1)$y)
        mm$q <- tmpdf$value
        mm$age <- NULL
        compl <- stats::complete.cases(as.data.frame(mm))
        mm2 <- purrr::map(mm, function(x) x[compl])
        rescol <- rep(NA_real_, nrow(tmpdf))
        rescol[compl] <- do.call(get(paste0("p",dists[[sex]]), asNamespace("gamlss.dist")), mm2)
        rescol
    })
    perc <- as.data.frame(params)
    perc$res <- ifelse(sex == "male", perc$male, perc$female)
    perc$res[is.infinite(value)] <- NA
    if(type == "SDS") return(stats::qnorm(perc$res))
    perc$res
}


##' Calculate SDS values - old version for comparison
##'
##' The function takes a vector of measurement values, and of age and of sex
##' and a RefGroup object as arguments. It calculates the sds or percentile
##' values.
##' @title Calculate SDS Values
##' @param value vector of measurement values
##' @param age vector of age values
##' @param sex vector of sex 
##' @param item name of the item e.g. "height"
##' @param ref RefGroup object
##' @param type "SDS" or "perc"
##' @param male coding of sex for male
##' @param female coding of sex for female
##' @return vector containing SDS or percentile values
##' @examples
##' anthro <- data.frame(age = c(11.61,12.49,9.5,10.42,8.42,10.75,9.57,10.48),
##'                      height = c(148.2,154.4,141.6,145.3,146,140.9,145.5,150),
##'                      sex = sample(c("male","female"), size = 8, replace = TRUE),
##'                      weight = c(69.5,72.65,47.3,51.6,45.6,48.9,53.5,58.5))
##' anthro$height_sds <- sds(anthro$height,
##'                          age = anthro$age,
##'                          sex = anthro$sex, male = "male", female = "female",
##'                          ref = kro.ref,
##'                          item = "height",
##'                          type = "SDS")
##' 
##' anthro$bmi <- anthro$weight/(anthro$height**2) * 10000
##' anthro$bmi_perc <- sds(anthro$bmi,
##'                        age = anthro$age,
##'                        sex = anthro$sex, male = "male", female = "female",
##'                        ref = kro.ref,
##'                        item = "bmi",
##'                        type = "perc")
##' data(who.ref)
##' x <- data.frame(height=c(50,100,60,54),
##'                 sex=c("m","f","f","m"),
##'                 age=c(0,2.9,0.6,0.2))
##' sds(value = x$height, age = x$age, sex = x$sex, male = "m", female = "f",
##'     ref = who.ref, item = "height")
##' @author Mandy Vogel
sdsold <- function(value, age, sex, item, ref, type = "SDS", male = "male", female = "female"){
    if(!(length(value) == length(age) & length(value) == length(sex))){
        print("value, age, and sex must be of the same length")
        invisible(return(NULL))
    }
    if(!item %in% names(ref@refs)){
        stop(paste(item,"is not available in the given refs object.",
                    "Please choose one of the following items:",
                    paste(names(ref@refs), collapse = " - ")))
        return(invisible(NULL))
    }
    sex <- as.character(factor(sex, levels = c(male, female), labels = c("male", "female")))
    refs <- ref@refs[[item]]@params
    dists <- ref@refs[[item]]@dist
    par.appr <- lapply(refs, function(df){
        as.data.frame(lapply(df, function(param) stats::approx(df$age, param, xout = age, rule = 1)$y))})
    res <- numeric()
    for(i in 1:length(value)) {
        if(is.na(value[i]) | any(is.na(unlist(par.appr[[sex[i]]][i,-1]))) ){
            res[i] <- NA
        } else{
        res[i] <- eval(parse(text = paste0("gamlss.dist::p",dists[[sex[i]]],"(",value[i],",",
                                           paste(paste(names(par.appr[[sex[i]]])[-1],"=", par.appr[[sex[i]]][i,-1]), collapse = ","),
                                           ")")))
        }
    }
    if(type == "SDS") return(stats::qnorm(res))
    return(round(res * 100, 2))
}



##' Calculate SDS values for 2-dimensional matrix of covariates -- old version
##'
##' The function takes a vector of measurement values, and of age and a
##' second covariate (like age and height for blood pressure) of sex
##' and a RefGroup object as arguments. It calculates the sds or percentile
##' values. This function is beta.
##' @title Calculate SDS Values for 2-dimensional matrix of covariates
##' @param value vector of measurement values
##' @param age vector of age values
##' @param x2 second vector of covariates
##' @param sex vector of sex 
##' @param item name of the item e.g. "height"
##' @param ref RefGroup object
##' @param type "SDS" or "perc"
##' @param male coding of sex for male
##' @param female coding of sex for male
##' @details the function searches for the nearest given point in the reference grid.
##' From there, the SDS/percentile value will be calculated. Different from \code{\link{sds}},
##' no interpolation will be applied. The procedure is according to Neuhauser et al. Blood
##' Pressure Percentiles by Age and Height from Nonoverweight Children and Adolescents
##' in Germany. 2011.
##' @return vector containing SDS or percentile values
##' @author Mandy Vogel
##' @export
sds_2d <- function(value, age, x2, sex, item, ref, type = "SDS", male = "male", female = "female"){
    if(!(length(value) == length(age) &
         length(value) == length(sex) &
         length(value) == length(x2))){
        print("value, age, x2, and sex must be of the same length")
        invisible(return(NULL))
    }
    if(!item %in% names(ref@refs)){
        stop(paste(item,"is not available in the given refs object.",
                    "Please choose one of the following items:",
                    paste(names(ref@refs), collapse = " - ")))
        return(invisible(NULL))
    }
    sex <- as.character(factor(sex,
                               levels = c(male, female),
                               labels = c("male", "female")))
    refs <- ref@refs[[item]]@params
    dists <- ref@refs[[item]]@dist
    par.appr <- list()
    par.appr <- lapply(refs, function(df){
        tmpdf <- data.frame(
                   mu = df$mu[class::knn(df[,c("age","x2")],
                                         data.frame(age = age,
                                                    x2 = x2),
                                         k = 1,cl = 1:nrow(df))],
                   nu = df$nu[class::knn(df[,c("age","x2")],
                                         data.frame(age = age,
                                                    x2 = x2),
                                         k = 1,cl = 1:nrow(df))],
                   sigma = df$sigma[class::knn(df[,c("age","x2")],
                                               data.frame(age = age,
                                                          x2 = x2),
                                               k = 1,cl = 1:nrow(df))]
                   ) 
        tmpdf
    })
    res <- numeric()
    for(i in 1:length(value)) {
        if(is.na(value[i]) | any(is.na(unlist(par.appr[[sex[i]]][i,]))) ){
            res[i] <- NA
        } else{
            print(paste0("gamlss.dist::p",dists[[sex[i]]],"(",value[i],",",
                                           paste(paste(names(par.appr[[sex[i]]]),"=", par.appr[[sex[i]]][i,]), collapse = ","),
                                           ")"))
        res[i] <- eval(parse(text = paste0("gamlss.dist::p",dists[[sex[i]]],"(",value[i],",",
                                           paste(paste(names(par.appr[[sex[i]]]),"=", par.appr[[sex[i]]][i,]), collapse = ","),
                                           ")")))
        }
    }
    if(type == "SDS") return(stats::qnorm(res))
    return(round(res * 100, 2)) 
}

##' Calculate SDS values for 2-dimensional matrix of covariates 
##'
##' The function takes a vector of measurement values, and of age and a
##' second covariate (like age and height for blood pressure) of sex
##' and a RefGroup object as arguments. It calculates the sds or percentile
##' values. This function is beta.
##' @title Calculate SDS Values for 2-dimensional matrix of covariates
##' @param value vector of measurement values
##' @param age vector of age values
##' @param x2 second vector of covariates
##' @param sex vector of sex 
##' @param item name of the item e.g. "height"
##' @param ref RefGroup object
##' @param type "SDS" or "perc"
##' @param male coding of sex for male
##' @param female coding of sex for male
##' @details the function searches for the nearest given point in the reference grid.
##' From there, the SDS/percentile value will be calculated. Different from \code{\link{sds}},
##' no interpolation will be applied. The procedure is according to Neuhauser et al. Blood
##' Pressure Percentiles by Age and Height from Nonoverweight Children and Adolescents
##' in Germany. 2011.
##' @return vector containing SDS or percentile values
##' @author Mandy Vogel
##' @importFrom interp interp
##' @export
sds2d <- function(value, age, x2, sex, item, ref, type = "SDS", male = "male", female = "female"){
    if(!(length(value) == length(age) &
         length(value) == length(sex) &
         length(value) == length(x2))){
        print("value, age, x2, and sex must be of the same length")
        invisible(return(NULL))
    }
    if(!item %in% names(ref@refs)){
        stop(paste(item,"is not available in the given refs object.",
                   "Please choose one of the following items:",
                   paste(names(ref@refs), collapse = " - ")))
        return(invisible(NULL))
    }
    value <- tidyr::replace_na(value, Inf)
    sex <- as.character(factor(sex, levels = c(male, female), labels = c("male", "female")))
    refs <- ref@refs[[item]]@params
    dists <- ref@refs[[item]]@dist
    tmpdf <- data.frame(age = age, x2 = x2, sex = sex, value = value)
    params <- purrr::map2(refs, names(refs), function(ref,sex){
        mm <- purrr::map(dplyr::select(ref, -age, -x2),
                         function(param) {
                             tmp <- interp::interp(x = ref$age, y = ref$x2,
                                            z = param,
                                            xo = tmpdf$age, yo = tmpdf$x2,
                                            ## method = "akima",
                                            output = "points")$z
                             tmp
                         })
        mm$age <- tmpdf$age
        mm$x2 <- tmpdf$x2
        mm$q <- tmpdf$value
        mm$age <- NULL
        mm$x2 <- NULL
        mm$mu[mm$mu < 0] <- NA
        compl <- stats::complete.cases(as.data.frame(mm))
        mm2 <- purrr::map(mm, function(x) x[compl])
        rescol <- rep(NA_real_, nrow(tmpdf))
        rescol[compl] <- do.call(get(paste0("p",dists[[sex]]), asNamespace("gamlss.dist")), mm2)
        rescol
    })
    perc <- as.data.frame(params)
    perc$res <- ifelse(sex == "male", perc$male, perc$female)
    perc$res[is.infinite(value)] <- NA
    if(type == "SDS") return(stats::qnorm(perc$res))
    perc$res
}


##' Calculate SDS values depending on the Tanner stage
##'
##' The function takes a vector of measurement values, and of tanner stage and of sex
##' and a RefGroup object as arguments. It calculates the sds or percentile
##' values.
##' @title Calculate SDS Values
##' @param value vector of measurement values
##' @param pubstatus vector of Tanner stages coded 1 to 5
##' @param sex vector of sex 
##' @param item name of the item e.g. "height"
##' @param ref RefGroup object
##' @param type "SDS" or "perc"
##' @param male coding of sex for male
##' @param female coding of sex for female
##' @return vector containing SDS or percentile values
##' @author Mandy Vogel
##' @export
sds_pub <- function(value, pubstatus, sex, item, ref, type = "SDS",
                    male = "male", female = "female" ) {
    if(!(length(value) == length(pubstatus) &
         length(value) == length(sex))){
        print("value, pubstatus and sex must be of the same length")
        invisible(return(NULL))
    }
    if(!item %in% names(ref@refs)){
        stop(paste(item,"is not available in the given refs object.",
                   "Please choose one of the following items:",
                   paste(names(ref@refs), collapse = " - ")))
        return(invisible(NULL))
    }
    value <- tidyr::replace_na(value, Inf)
    sex <- as.character(factor(sex,
                               levels = c(male, female),
                               labels = c("male", "female")))
    if(!inherits(pubstatus,"factor") ) pubstatus <- factor(pubstatus,
                                                           levels = 1:5)
    refs <- ref@refs[[item]]@params
    dists <- ref@refs[[item]]@dist
    tmpdf <- data.frame(age = pubstatus, sex = sex, value = value)
    params <- purrr::map2(refs, names(refs), function(ref,sex){
        mm <- dplyr::left_join(
                         dplyr::select(tmpdf,-sex,-value), ref,
                         by = "age" )
        mm$q <- tmpdf$value
        mm$age <- NULL
        compl <- stats::complete.cases(as.data.frame(mm))
        mm2 <- purrr::map(mm, function(x) x[compl])
        rescol <- rep(NA_real_, nrow(tmpdf))
        rescol[compl] <- do.call(get(paste0("p",dists[[sex]]), asNamespace("gamlss.dist")), mm2)
        rescol
    })
    perc <- as.data.frame(params)
    perc$res <- ifelse(sex == "male", perc$male, perc$female)
    perc$res[is.infinite(value)] <- NA
    if(type == "SDS") return(stats::qnorm(perc$res))
    perc$res
}


##' Calculate percentage relative to a given base percentile
##'
##' The function calculates the percentage of a given bmi value relative to
##' a specific percentile
##' @title Calculate percentage relative to a given base percentile
##' @param bmi vector of bmi
##' @param age vector of age
##' @param sex vector of sex (coding "male" and "female") is assumend
##' @param ref.perc single value: reference percentile (0,100)
##' @param ref RefGroup object
##' @param item item within ref
##' @param rownr indicator of order
##' @return vector containing values between 0 and 1
##' @author Mandy Vogel
##' @export
calc_percent_excess <- function (bmi = NULL, age = NULL, sex = NULL,
                                 ref.perc = 50,
                                 ref, item = "bmi", rownr = NULL) {
    tmp <- data.frame(age = age, sex = sex)
    tmp$rownr <- 1:nrow(tmp)
    dists <- unlist(ref@refs[[item]]@dist)    
    refs <- ref@refs[[item]]@params
    sexes <- names(refs)
    refs <- purrr::map(refs, function(ref){
        age <- age[dplyr::between(age, min(ref$age), max(ref$age))]
        purrr::map(ref,
                   function(param) stats::approx(x = ref$age,
                                                     y = param,
                                                 xout = age,
                                                 rule = 1)$y)
        })
    ref.perc <- ref.perc/100
    pertab <- purrr::map2(refs, sexes, function(ref,sex) {
        params <- base::as.list(ref)
        params$age <- NULL        
        res <- purrr::map(ref.perc, function(p) {
            params$p <- p 
            do.call(get(paste0("q", dists[sex]), asNamespace("gamlss.dist")), params)
        })
        names(res) <- "basis"
        res <- as.data.frame(res)
        if(base::length(res$basis) != base::length(bmi)) stop("bmi and base vector are not of the same length")
         res$percent <- bmi/res$basis
        res$sex <- sex
        res$age <- ref$age        
        res 
    })
    pertab <- purrr::map(pertab, ~ {
        .$rownr <- 1:nrow(tmp)
        return(.)})
    res <- base::Reduce(dplyr::bind_rows, pertab)
    res <- dplyr::left_join(tmp, res, by = c("rownr","age","sex"))
    res <- dplyr::arrange(res, rownr)
    res$percent
}



##' Calculate SDS values depending on the Tanner stage and a second variable
##'
##' The function takes a vector of measurement values, and of tanner stage
##' of a second variable (x2) and of sex
##' and a RefGroup object as arguments. It calculates the sds or percentile
##' values.
##' @title Calculate SDS Values
##' @param value vector of measurement values
##' @param x2 2nd predictor (vector), e.g. bmisds must be contained in reference
##' @param sex vector of sex 
##' @param item name of the item e.g. "height"
##' @param ref RefGroup object
##' @param type "SDS" or "perc"
##' @param male coding of sex for male
##' @param female coding of sex for female
##' @param pubstat vector of Tanner stages coded 1 to 5
##' @param age not used yet
##' @param id order of values
##' @return vector containing SDS or percentile values
##' @author Mandy Vogel
##' @export
sds_pub2d <- function(value, pubstat, x2, sex, item, ref,
                      type = "SDS", male = "male", female = "female", age = NULL, id = 1:length(value)){
    if(!(length(value) == length(pubstat) & length(value) == length(sex))){
        print("value, pubstat, and sex must be of the same length")
        invisible(return(NULL))
    }
    if(!item %in% names(ref@refs)){
        stop(paste(item,"is not available in the given refs object.",
                   "Please choose one of the following items:",
                   paste(names(ref@refs), collapse = " - ")))
        return(invisible(NULL))
    }
    value <- tidyr::replace_na(value, Inf)
    sex <- as.character(factor(sex,
                               levels = c(male, female),
                               labels = c("male", "female")))
    refs <- ref@refs[[item]]@params
    dists <- ref@refs[[item]]@dist
    tmpdf <- data.frame(age = pubstat,
                        sex = sex,
                        value = value,
                        x2 = x2,
                        id = 1:length(pubstat))
    pubstats <- unique(stats::na.omit(tmpdf$age))
    res <- purrr::map(pubstats, function(ts) {
      params <- purrr::map2(refs, names(refs), function(ref,sex){
        ref <- dplyr::filter(ref, age == ts)
        tmpdf <- dplyr::filter(tmpdf, age == ts)
        tmpdf$age <- NULL
        ref$age <- NULL
        mm <- purrr::map(ref,
                         function(param) stats::approx(ref$x2,
                                                       param, xout = tmpdf$x2, rule = 1)$y)
        mm$q <- tmpdf$value
        mm$age <- NULL
        mm$x2 <- NULL
        compl <- stats::complete.cases(as.data.frame(mm))
        mm2 <- purrr::map(mm, function(x) x[compl])
        rescol <- rep(NA_real_, nrow(tmpdf))
        rescol[compl] <- do.call(get(paste0("p",dists[[sex]]), asNamespace("gamlss.dist")), mm2)
        tmpdf[[sex]] <- rescol
        tmpdf$ts <- ts
        tmpdf
      })
      perc <- purrr::reduce(params, dplyr::inner_join)
      perc$res <- ifelse(perc$sex == "male", perc$male, perc$female)
      perc$res[is.infinite(perc$res)] <- NA
      dplyr::select(perc, id, res)
    }) %>% dplyr::bind_rows()
    tmpdf <- dplyr::select(tmpdf, id) %>%
      dplyr::left_join(res, by = "id")
    perc <- tmpdf$res
    if(type == "SDS") return(stats::qnorm(perc))
    perc
}


