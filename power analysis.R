#install.packages(c("lme4","lmerTest","dplyr","broom","broom.mixed"))
library(lme4)
library(lmerTest)
library(plot3D)
library(dplyr)
library(broom)
library(broom.mixed)

# nSamps: number of samples per site
# nSites: number of sites (or groups)
# sd.obs: standard deviation in ecological metric (e.g., size)
# sd.group: standard deviation in how the mean of the metric varies among groups

# Example: imagine we collected phosphorus readings at each of a number of sites
# phosophorus reading at eah site is our "x" variable
# we want to ask: is the performance of critters at each site influenced by phosphorus
# critter performance is our "y" variable
# site is our replicate
# individual critters are our samples

costSite <- 4 # hours of work/driving to add extra site
costSamp <- 0.5 # hours of work to add extra sample

nSamps <- 25
nSites <- 15
effect <- 2
sd.obs <- 4 # what is the sample variance
sd.group <- 1 # what is the variability across sites/groups

Nboots <- 1000 # how many bootstraps to run: the larger the better (but slow!)

# Make simulated dataframe to see what is happening in our statsig function
site <- as.factor(rep(1:nSites, each = nSamps)) # site ID
phos_bar <- rnorm(nSites, mean = 0, sd = 1) # site level phosphorous
phos_samp <- round(rep(rnorm(nSites, mean = 0, sd = 1), each = nSamps), 3) # values of phosphorous are the same for all samples because we only took one measurement at each site

# the true model is:
# y = mx + b + error, where b = 0.2, m = effect, and error is set by two processes
# group variance and sample variance

y_bar <- rep(round(0.2 + effect * phos_bar + rnorm(length(phos_bar), mean = 0, sd = sd.group), 3), each = nSamps) # create random mean fish abundance values at each site for a given effect size and phosphorous level
y_samps  = round(rnorm(length(phos_samp), mean = rep(round(0.2 + effect * phos_bar + rnorm(length(phos_bar), mean = 0, sd = sd.group), 3), each = nSamps), sd = sd.obs), 3) # measured fish abundance at each of five sampling points per site.

# View the simulated dataset
dat <- data.frame(site = site, phos_samp, y_bar, y_samps)

head(dat)


plot(dat$phos_samp,dat$y_samps)


statsig <- function(nSamps, nSites, effect, sd.obs, sd.group) {
  # imagine you have only 1 measurement of x per stream/site
  x <- rnorm(nSites, mean = 0, sd = 1)  # generate the predictor variable at the site
  x1 <- rep(x, each = nSamps)  # repeat the x variable for each sample at the site
  group <- as.factor(rep(1:nSites, each = nSamps)) # track the site ID for each samples
  y <- 0.2 + effect * x + rnorm(length(x), mean = 0, sd = sd.group) # create random ecological metric
  y1 <- rnorm(length(x1), mean = rep(y, each = nSamps), sd = sd.obs) # create random data for each group
  fit1 <- lm(y1 ~ x1)
  fit2 <- suppressMessages(lmer(y1 ~ x1 + (1 | group), REML = F))
  # get upper and lower 95% CIs
  fit1_summ <- 
    broom::tidy(fit1) %>% 
    mutate(lower_ci = estimate - 1.96 * std.error, 
           upper_ci = estimate + 1.96 * std.error)
  fit2_summ <- 
    broom.mixed::tidy(fit2) %>% 
    mutate(lower_ci = estimate - 1.96 * std.error, 
           upper_ci = estimate + 1.96 * std.error)
  # Calculate power for lm
  power1 <- effect >= fit1_summ$lower_ci[2] & 
            effect <= fit1_summ$upper_ci[2] &
            0 <= (fit1_summ$lower_ci[2] * fit1_summ$upper_ci[2])
  # Calculate power for mixed model
  power2 <- effect >= fit2_summ$lower_ci[2] & 
            effect <= fit2_summ$upper_ci[2] &
            0 <= (fit2_summ$lower_ci[2] * fit2_summ$upper_ci[2])
  
  round(c(slope1 = fit1_summ$estimate[2], R1 = summary(fit1)$adj.r.squared, 
          Type_I_error_1 = fit1_summ$p.value[2] <= 0.05, power1 = power1, 
          slope2 = fit2_summ$estimate[2], R1 = summary(fit2)$adj.r.squared, 
          Type_I_error_2 = fit2_summ$p.value[2] <= 0.05, power2 = power2), 3)
}

# set up our exploration and scenarios
effect_vec <- c(0, 1, 2) # no effect, small effect, big effect
nSites <- seq(from = 5, to = 50, by = 5) # number of sites to sample
nSamps <- c(3, 5, 10, 15, 20, 25, 30, 35, 50, 100) # number of samples per sight

results <- 
  array(NA, dim = c(length(effect_vec), length(nSites), length(nSamps), 7), 
        dimnames=list("Effect size"=effect_vec, "N Sites"=nSites, 
                      "N Samps per Site"=nSamps, 
                      "Metric"=c("slope (lm)", "R2 (lm)", "Type I error (lm)", "Power (lm)", "Slope (lme4)", "Type I error (lme4)", "Power (lme4)")))

counter <- 1
progBar <- txtProgressBar(min = 0,  max = length(effect_vec) * length(nSites) * length(nSamps), title = "More power", style = 3, initial = 0)
ptm = Sys.time()
for(i in 1:length(effect_vec)){
  for(j in 1:length(nSites)){
    for(k in 1:length(nSamps)){
      out <- replicate(Nboots, statsig(nSamps=nSamps[k], nSites=nSites[j], effect=effect_vec[i], sd.obs=sd.obs, sd.group=sd.group), simplify=T)
      results[i, j, k, ] <- rowMeans(out)
      setTxtProgressBar(progBar, counter)
      counter <- counter + 1
    }
  }
}
endtime <- Sys.time()-ptm
endtime

# analyze your cost-benefit ratios
power <- results[2, , , "Power (lme4)"]

cost <- sapply(costSamp*nSamps, function(x){x+costSite*nSites})
dimnames(cost) <- dimnames(power)

# find where you have minimal costs that maximizes your power and chance to avoid type 1 error
power2 <- power
cost2 <- cost
power2[power < 0.8] <- NA
cost2[power<0.8] <- NA

cost2 <- cost2/max(cost2,na.rm=T) # evaluate cost on the same scale as statistical power: 0-1

costBenefit <- which(cost2/power2==min(cost2/power2,na.rm=T), arr.ind=T)
paste("Sites = ", nSites[costBenefit[1]], " & ", 
        "Samples = ", nSamps[costBenefit[2]], sep="")

layout(1)
par(mar=c(5, 4, 3, 3))
filled.contour(x=nSites,  y=nSamps,  z=results[2, , , "Power (lme4)"], ylab='Number of samples', xlab='Number of sites', color.palette = colorRampPalette(c("red", 'orange', 'dodgerblue', "yellow")))

layout(1)
par(mar=c(5, 4, 3, 3))
filled.contour(x=nSites,  y=nSamps,  z=results[2, , , "Type I error (lme4)"], ylab='Number of samples', xlab='Number of sites', color.palette = colorRampPalette(c("red", 'orange', 'dodgerblue', "yellow")))
