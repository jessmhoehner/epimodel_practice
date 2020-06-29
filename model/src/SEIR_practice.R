#!/usr/bin/env Rscript --vanilla
# set expandtab ft=R ts=4 sw=4 ai fileencoding=utf-7
#
# Author: JR
# Maintainer(s): JR
# License: 2020, GPL v3 or later
#
# -----------------------------------------------------------
# epimodel_practice/model/src/model.R

pacman::p_load("tidyverse", "here", "assertr",
               "readr", "janitor", "EpiModel")

# pull in _state_data and _inc datasets

files <- list(clean_data = here("epimodel_practice/model/input/cendata_clean.csv"),
              
              ppos_facetplot = here("epimodel_practice/graph/input/graph_data.csv")
)

# notes on EpiModel, practice in progress

# read in data 

cen_df <- as.data.frame(read_delim(files$clean_data, delim = "|"))

# total state population
N <- as.integer(cen_df$ga_num[1])

# state pop under 18 
n_sub_18 <- as.integer((cen_df$ga_num[6]/100)*N)

# state pop over 65

n_super_65 = as.integer((cen_df$ga_num[7]/100)*N)

# rest of state pop

n_rest_cp <- as.integer(N-n_sub_18-n_super_65)

# number of initial exposed by 
# first 4 cases identified in GA + 4 est unidentified cases * average household size

n_exp <- 8 * cen_df$ga_num[27]

##################################################

# define SEIR function 

SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + e.num + i.num + r.num
    
    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur
    lambda <- ce * i.num/num
    if (!is.null(parms$inter.eff) && t >= inter.start) {
      lambda <- lambda * (1 - inter.eff)
    }
    
    dS <- -lambda*s.num
    dE <- lambda*s.num - (1/e.dur)*e.num
    dI <- (1/e.dur)*e.num - (1 - cfr)*(1/i.dur)*i.num - cfr*(1/i.dur)*i.num
    dR <- (1 - cfr)*(1/i.dur)*i.num
    
    # Compartments and flows are part of the derivative vector
    # Other calculations to be output are outside the vector, 
    # but within the containing list
    list(c(dS, dE, dI, dR, 
           se.flow = lambda * s.num,
           ei.flow = (1/e.dur) * e.num,
           ir.flow = (1 - cfr)*(1/i.dur) * i.num,
           d.flow = cfr*(1/i.dur)*i.num),
         num = num,
         i.prev = i.num / num,
         ei.prev = (e.num + i.num)/num)
  })
}

# parameterize #
# t0 = Feb 01 2020
# t500 = June 15 20201

# R0 estimate taken from https://rt.live for GA, on 2/26/20 
# recorded on 6/29/2020 at 1:09 pm

# accouting for interventions at 
# t51 = Large gatherings banned
# t79 = shelter in place initiated
# t109 = shelter in place extended for at risk

p <- param.dcm(R0 = 2.73, 
               e.dur = n_exp, 
               i.dur = 8, # number of cases at t0, based on GA DPH + est uncaught cases
               cfr = c(0.2, 0.4, 0.6, 0.8, 1.0))

#set initial state of n susceptible and infected at first time point (t1)

i <- init.dcm(s.num = N,
              e.num = n_exp,
              i.num = 8, 
              r.num = 0, 
              se.flow = 0, 
              ei.flow = 0, 
              ir.flow = 0, 
              d.flow = 0)

# collect model controls, model type and number of time steps (ex: in days) 
#    of the simulation, dt = fraction of time units, default = 1
# can add more controls as necessary for time lagged variables like interventions

c <- control.dcm(new.mod = SEIR, 
                 nsteps = 150, 
                 dt = 1, 
                 dede = TRUE, 
                 nsims = 10)

# produces param, control, and epi
# p: parameters passed to the model through the p data object of class param.dcm
# c: controls passed to the model through the c data object of class control.dcm
# i: a list of data frames, one for each epidemiological output from the model
#     use as.data.frame.dcm to extract results as dataframes, 
#     summarize results with summary.dcm function on the m object of class dcm
#     plot the m object with plot.dcm
#     plot a flow diagram with comp_plot function on dcm object

m <- dcm(p, i, c)

# "Printing the model object provides basic information on model input and output,
# including model parameters and data variables for plotting and analysis." 
#   use as.data.frame to extract results as dataframes, 
#   summarize results with summary.dcm function on the m object of class dcm
#   plot the m object with plot.dcm
#    plot a flow diagram with comp_plot function on dcm object

m_df <- as.data.frame(m)

# custom color palette, green = lower CFR, purple = higher
custpal <- c("#00441b", "#1b7837","#c2a5cf", "#762a83", "#40004b")

#plot disease prevalence
plot(m, 
     y = "i.num",
     col = custpal, 
     main = "CFR Scenarios of Prevalence",
     grid = TRUE, 
     legend="full",
     leg.name = c("CFR = 0.2", "CFR = 0.4", 
                  "CFR = 0.6", "CFR = 0.8", "CFR = 1.0"),
     alpha = 0.7, 
     sim.lines = TRUE, 
     sim.alpha = 0.15, 
     ylim = c(0, 2000)
     )

# done

