


## 0. install/load packages ----------------------

use_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

use_package("secr")



## 1. Load data  ----------------------------

captures_url <- "https://raw.githubusercontent.com/harryjobann/clouded_secr/refs/heads/main/captures.csv"
traps_url    <- "https://raw.githubusercontent.com/harryjobann/clouded_secr/refs/heads/main/traps.csv"

download.file(captures_url, destfile = "captures.csv", mode = "wb")
download.file(traps_url,    destfile = "traps.csv",    mode = "wb")

captures_df <- read.csv("captures.csv")
traps_df    <- read.csv("traps.csv")




## 2. Inspect -----------------------------
head(captures_df)
head(traps_df)


## 3. Work out number of occasions -------------------------

## Assumes the captures file has a column called 'Occasion'
## storing the sampling day as an integer.


# traps_df already loaded from traps.csv
n_occasions <- ncol(traps_df) - 3 #minus 3 because first 3 cols are Detector, X, Y
n_occasions


## 4. Build the capthist object ----------------------------

clouded_capthist <- read.capthist(
  captfile   = "captures.csv",
  trapfile   = "traps.csv",
  detector   = "proximity",
  fmt        = "trapID",
  noccasions = n_occasions
)

## If everything is formatted correctly you should see:
## "No errors found :-)"

summary(clouded_capthist)
plot(clouded_capthist, border = 3000, tracks = TRUE, varycol = TRUE)


## 5. Get a simple perimeter from trap coords -
## This uses the trap coordinates in the capthist to build a polygon

trap_coords <- as.data.frame(traps(clouded_capthist))  
head(trap_coords)


hull_idx <- chull(trap_coords$x, trap_coords$y)
perimeter_hull <- trap_coords[hull_idx, c("x", "y")]


## 6. Build the mask (state space) with 1 km buffer --------

## This is where we actually implement "coordinates + 1 km
## buffer round the edge" for the SECR state space.


buffer_width_m <- 6000   # 14 km buffer
spacing_m      <- 500    # grid spacing for the mask - using 500 here as home range slightly smaller than for tigers 

clouded_mask <- make.mask(
  traps   = traps(clouded_capthist),
  buffer  = buffer_width_m,
  spacing = spacing_m,
  type    = "trapbuffer"
)

plot(clouded_mask)
points(traps(clouded_capthist))  # check traps sit inside the mask grid


## 7. Fit SECR model -------------------------------

## Constant detection probability (g0) across individuals,
## traps and occasions.

clouded_secr <- secr.fit(
  capthist   = clouded_capthist,
  model      = g0 ~ 1,
  mask       = clouded_mask,
  CL         = TRUE,
  trace      = TRUE,
  biasLimit  = NA,
  verify     = TRUE
)



## 7a. Build a capthist with Sex as an individual covariate
clouded_capthist_sex <- read.capthist(
  captfile   = "captures.csv",
  trapfile   = "traps.csv",
  detector   = "proximity",
  fmt        = "trapID",
  noccasions = n_occasions,
  covnames   = c("Sex")
)

## Sanity check
table(covariates(clouded_capthist_sex)$Sex)

## Fit model by sex
clouded_secr_sex_group <- secr.fit(
  capthist = clouded_capthist_sex,
  model    = list(g0   = ~ g,   # g0 differs by Sex group
                  sigma = ~ 1), # same sigma for both groups (sensible with limited data)
  mask     = clouded_mask,
  CL       = FALSE,             # full likelihood needed to use groups
  groups   = "Sex",             # tell secr that Sex defines groups
  trace    = TRUE,
  verify   = TRUE,
  biasLimit = NA
)

clouded_secr_sex_group




## 7c. Model comparison
AIC(clouded_secr, clouded_secr_sex_group)



## 8. Inspect model results --------------------------------

clouded_secr          
derived(clouded_secr) # density, CV, confidence limits etc.



## 8a. Overall derived parameters
derived(clouded_secr_sex_group)

## 8b. Density by Sex
derived(clouded_secr_sex_group, group = "Sex")




## 9. Estimate population size in the region ---------------
## Here 'region' is the mask (traps + 1 km buffer).

clouded_N <- region.N(
  object = clouded_secr,
  region = clouded_mask
)

clouded_N



