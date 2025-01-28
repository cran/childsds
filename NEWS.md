# childsds 0.9.9
* add leptin
* add sds_pub2d()
# childsds 0.9.8
* add lightness
# childsds 0.9.4
* add calc\_perc\_excess()
* add life_vegf.ref
* add life_shgb.ref
# childsds 0.9.1
## Major changes
* add life_steroide.ref
* add life_oxyandrogen.ref
## Minor changes
* correct life_thyroid.ref
* add life_hba1c.ref
# childsds 0.9.0
## Major changes
* new version sds2d
* new function make_percentile_tab2d
* new function wormplot2d_gg
* add ripka_bf.ref
* add linden_heart.ref
* add ofenheimer_bf.ref
* add cole_lobstein.ref
* add valencia_nc.ref
* add ghouili_anthro.ref
* add gomez_bmitmi.ref
* add schafmeyer_leg.ref
* add liao_igf1.ref
## Minor changes
* now, the wormplots using quantiles instead of equal spaced cut points for grouping
* remove non-valid references from nl4.ref
# childsds 0.8.0
## Major changes
* add life_fibroscan.ref
* add kawel_boehm.ref
* add life_igf.ref
* add life_heart.ref
* add momo.ref
# childsds 0.7.5
## Major changes
* add cn.refs
* add life_thyr.ref
* add metabolom.ref
* add motor.ref
# childsds 0.7.5
## Major changes
* add bone.refs
# childsds 0.7.4
## Attention
* There were some implausible values in the Wuehl referencens. We replaced these values by interpolated ones.
## Bugs
* solved bug "The percent BMI code appears to be bugged and is outputting a data.frame."
## New References
* Duran et al, bodyfat, fat mass, lean body mass index
* Doyon et al, carotid artery intima-media thickness, distensibility

# childsds 0.7.3
## Minor changes
* bug fix nl4 references (thanks to Robert Euser)
* bug fix make_percentile_tab()

# childsds 0.7.1
## Minor changes
* bug fix sds() function
* update make_percentile_tab()

# childsds 0.7.0
## Major changes
* add aga15 references
* add skinfold references
* add wuehl blood pressure references
* add whr, whtr references
* replace explicite loop in sds() function
## Minor changes

## childsds 0.6.7 ##

### Major changes ###

* add uk1990 references
* add liver enzymes references

### Minor changes ###

# childsds 0.6.6
## Major changes
* add ggplot version of the worm plot function
## Minor changes
* add include.pars to make_percentile_tab
* add transformation argument to lms() wrapper
# childsds 0.6.5
## Major changes
* add who 2007 reference tables
* add nl3 reference tables
* add preterm reference tables
* add turkish reference tables height, weight
* add igf international and japan
* add saudi arabian preschool references
## Bug fixes
* removes class labelled from who references
* add missing argument to one_it inside do_iterations
# childsds 0.6.4
## Major changes
* add colombian skinfold references
## Bug fixes
*  fix typos in italian references
# childsds 0.6.3
## Major changes
* add option whether or not models should be kept during the iteration process
## Minor changes
* add subscapular skinfold to kiggs.ref
* add function to find better upper and lower age bounds (full month)
## Bug fixes
*  remove user accessable documentation of unexported functions

# childsds 0.6.2
## Major changes
* add sds_2d for calculation of bloodpressure sds
* add kiggs bloodpressure references
* add additional new references (Belgium)

## Bug fixes
* Remove direct calls to dplyr::n()

# childsds 0.6.1

* Added a `NEWS.md` file to track changes to the package.

## Major changes
* add new references
* add function make_percentile_tab() to extract data for use in ggplot()
## Bug fixes
* fix namespace problems for gamlss distributions
* fix missing export of sds() function

