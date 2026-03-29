#This code contains all the individual replication scripts 
# needed to replicate "The rising cost of climate change: 
# evidence from the bond market”
# by Bauer & Rudebusch (Review of Economics and Statistics, 2023). 

rm(list = ls())
#Your working directory must be organized as follows
# ROOT (your main working directory)
#  |- data (stores data)
#  |- R (stores scripts)
#  |- results (stores .RData files produced by scripts)
#        |- plots (stores plots produced by scripts)
#        |- tabs (stores tables produced by scripts)

#Set working directory to ROOT folder
setwd("ROOT")

#INSTALL AND LOAD PACKAGES
source("R/setup.R")

#ESTIMATE UC MODEL
graphics.off()

#Load functions and dependencies
source("R/uc_fns.R")
source("R/data_fns.R")
source("R/sdr_fns.R")
source("R/estimate_uc_real_rate.R")

#Load related modules
library(KFAS)
library(dplyr)

#Maturities to be estimated
rates <- c("1y","10y")

#Estimate UC
for (i in rates){
  estimate_uc_real_rate(
    rate = i,
    N = 5,
    M = 20000,
    seed = 616,
    save_results = TRUE,
    results_dir = "results",
    verbose = TRUE
  )
}

#PLOT FIGURE 1 (adjust linewidth)
source("R/figure_1.R")

#COMPUTE AND PLOT TERM STRUCTURE - UC MODEL
graphics.off()
source("R/sdr_uc.R")
rates <- c("1y","10y")

for (i in rates){
  run_uc_sdr(
    spec = i,
    seed = 616,
    plot_figure = FALSE,
    save_results = TRUE,
    results_dir = "results",
    LB = 0,
    year1 = 1990,
    year2 = 2019,
    N = 400,
    fig_width = 7,
    fig_height = 4.5
  ) 
}

#COMPUTE AND PLOT TERM STRUCTURE - AR MEANSHIFT
library(sandwich)
library(dplyr)
library(ggplot2)
graphics.off()
source("R/ar_meanshift.R") #create a function like uc

rates <- c("1y","10y")
for (i in rates){
  run_ar_meanshift_sdr(
    spec = i,
    seed = 616,
    N=400,
    byear = 1990,
    plot_figure = TRUE,
    save_results = TRUE,
    results_dir = "results",
    fig_width = 7,
    fig_height = 4.5
  )
}

#COMPUTE AND PLOT TERM STRUCTURE - AR(3) EXPONENTIAL LEARNING
library(dplyr)
library(dynlm)
library(sandwich)
library(ggplot2)
graphics.off()
source("R/ar_learning.R") #create a function like uc

rates <- c("1y","10y")
for (i in rates){
  run_ar3_expon_sdr(
    spec = i,
    seed = 616,
    year1 = 1990,
    year2 = 2019,
    alpha = 0.98,
    N = 400,
    plot_figure = TRUE,
    save_results = TRUE,
    results_dir = "results",
    fig_width = 7,
    fig_height = 4.5
  )
}

#PRODUCE TABLE 1
graphics.off()
source("R/table_1.R")
export_sdr_table_tex() #render and export in .tex to the results/tables folder

#PLOT DIFFERENT TERM STRUCTURES OF SDRs ACROSS 
library(ggplot2)
graphics.off()
source("R/figure_2.R") #here we plot - set the previous models to produce only data

#CALCULATE THE ESTIMATES OF THE SOCIAL COST OF CARBON
graphics.off()
library(dplyr)
source("R/scc.R")

export_scc_table_tex() #render and export in .tex to the results/tables folder
