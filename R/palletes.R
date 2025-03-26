##' change saturation and brightness of given colors
##'
##' @title change saturation and brightness of given colors
##' @param colors vector of hex colors
##' @param saturation_factor factor
##' @param brightness_factor factor
##' @return vector of adjusted colors
##' @author Mandy Vogel
##' @export
adjust_colors <- function(colors, saturation_factor = 1, brightness_factor = 1) {
  hsv_values <- DescTools::ColToHsv(colors)
  # Adjust saturation and brightness
  hsv_values["s", ] <- pmin(pmax(hsv_values["s", ] * saturation_factor, 0), 1)  # Keep within [0,1]
  hsv_values["v", ] <- pmin(pmax(hsv_values["v", ] * brightness_factor, 0), 1)  # Keep within [0,1]
  HSV <- t(hsv_values)
  colnames(HSV) <- toupper(colnames(HSV))
  colorspace::hex(colorspace::HSV(HSV))
  # Convert back to RGB and then HEX
  adjusted_colors <- apply(hsv_values, 2, function(hsv) hsv(hsv[1], hsv[2], hsv[3]))
  return(adjusted_colors)
}



##' color palette Pediatric research 
##'
##' @title color palette Pediatric research 
##' @param pal orig, sat or bright
##' @param primary first colour for two-color scales, must be one of the scale colors
##' @param other second colour for two-color scales
##' @param direction if -1 one order of colours is reversed
##' @return function of n
##' @author Mandy Vogel
##' @export
prs_pal <- function(pal = "orig", primary = "blue", other = "violet", direction = 1){
  col.names <- c("blue","bluegreen","yellow", "orangebrown", "rose","lgreen", "red","violet","green")
  prs.cols <- c("#548CC3","#7BD0C3" , "#CBC436","#D19940","#CD5EA6",
                "#6CA43E" ,"#D16256","#956CE1", "#497D4E" )
  prs.cols.sat <- adjust_colors(prs.cols, saturation_factor = 1.9, brightness_factor = 1)
  prs.cols.bright <- adjust_colors(prs.cols, saturation_factor = 1.5, brightness_factor = 1.4)
  prs.l <- list(orig = prs.cols,
            bright = prs.cols.bright,
            sat = prs.cols.sat)
  print(pal)
  (prs.l <- purrr::map(prs.l, ~set_names(., col.names)))
  pos <- grep(pal, names(prs.l))
  prs_colors <- prs.l[[pos]]
  stopifnot(primary  %in% names(prs.l[[1]]))
  function(n){
    if(n > 9) warning("only 7 colours available")
    if(n == 2) {
      other <- if(!other  %in% names(prs.l[[1]])){
                 other
               } else {
                 prs_colors[other]
               }
      color_list <- c(other, prs_colors[primary])
    } else {
      color_list <- prs.l[[pal]][1:n]
    }
    color_list <- unname(unlist(color_list))
    if(direction >= 0) color_list else rev(color_list)
  }
}


##' discrete ggplot colour scale 
##'
##' @title discrete ggplot colour scale of PRS colors
##' @param pal orig, sat, or bright
##' @param primary primary colour of two
##' @param other the other colour
##' @param direction if -1 one order of colours is reversed
##' @param ... further argument
##' @return colour scale for use with ggplot()
##' @author Mandy Vogel
##' @export
scale_colour_prs <- function(pal="orig",
                             primary = "blue", other = "violet",
                             direction = 1,...){
    ggplot2::discrete_scale(
             aesthetics = "colour",
             palette = prs_pal(pal, primary, other, direction),
             ...
           )
}

##' discrete ggplot fill scale 
##'
##' @title discrete ggplot colour scale of PRS colors
##' @param pal orig, sat, or bright
##' @param primary primary colour of two
##' @param other the other colour
##' @param direction if -1 one order of colours is reversed
##' @param ... further arguments
##' @return colour scale for use with ggplot()
##' @author Mandy Vogel
##' @export
scale_fill_prs <- function(pal="orig",
                             primary = "blue", other = "violet",
                             direction = 1,...){
    ggplot2::discrete_scale(
             aesthetics = "fill",
             palette = prs_pal(pal, primary, other, direction),
             ...
           )
}



##' Uchu colour scale
##'
##' provides the uchu colour scale as ggplot2 colour scale (https://github.com/NeverCease/uchu)
##' @title uchu colour scale
##' @name scale_uchu
##' @param pal colour palette
##' @param ... further arguments to scale
##' @return ggplot colour scale
##' @author Mandy Vogel
##' @export
scale_colour_uchu <- function(pal=c("gray","red","pink","purple","blue","green","yellow","orange","general"),...){
    ggplot2::continuous_scale(
                 aesthetics = "colour",
                 scale_name = "uchu",
                 palette = scales::colour_ramp(uchu.pal[[pal]]),
                 ...
           )
}

##' @rdname scale_uchu
##' @export
scale_fill_uchu <- function(pal=c("gray","red","pink","purple","blue","green","yellow","orange","general"),...){
    ggplot2::continuous_scale(
                 aesthetics = "fill",
                 scale_name = "uchu",
                 palette = scales::colour_ramp(uchu.pal[[pal]]),
                 ...
             )
}

