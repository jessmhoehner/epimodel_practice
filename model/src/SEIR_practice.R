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
# first cases identified in GA *6

n_0_est <- 4

n_exp <- n_0_est * cen_df$ga_num[27]

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

R0_est = 2.73
# R0 estimate taken from https://rt.live for GA, on 2/26/20 
# recorded on 6/29/2020 at 1:09 pm

# duration of exposed state
# assume a 2 week time of exposure and infectious period
t_dur = 14
t_inf = 14

# accouting for interventions at as documented by the GADPH
# https://dph.georgia.gov/covid-19-daily-status-report

# t51 = Large gatherings banned
# t79 = shelter in place initiated
# t109 = shelter in place extended for at risk
int_t <- c(51, 79)
int_eff_ests <- c(0.01, 0.02)

# potential activity rates 
act_rates <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

#potential infection probabilities
inf_probs <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

# potential cfrs
cfrs <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

# potential recovery rates
rec_rates <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0)

p <- param.dcm(R0 = R0_est, 
               e.dur = t_dur, 
               i.dur = t_inf,
               inter.eff = int_eff_ests, 
               inter.start = int_t,
               cfr = cfrs, 
               act.rate = act_rates, 
               inf.prob = inf_probs,
               rec.rate = rec_rates)

#set initial state of n susceptible and infected at first time point (t1)

i <- init.dcm(s.num = N,
              e.num = n_exp,
              i.num = n_0_est, 
              r.num = 0, 
              se.flow = 0, 
              ei.flow = 0, 
              ir.flow = 0, 
              d.flow = 0)

# collect model controls, model type and number of time steps

c <- control.dcm(new.mod = SEIR, 
                 nsteps = 150, 
                 dt = 1, 
                 dede = TRUE,
                 sens.param = TRUE,
                 nsims = 6)

# produces param, control, and epi
# p: parameters passed to the model through the p data object of class param.dcm
# c: controls passed to the model through the c data object of class control.dcm
# i: a list of data frames, one for each epidemiological output from the model

m <- dcm(p, i, c)

m_df <- as.data.frame(m)

# check if the model is matching observed cases and see how far off it is
m_df$se.flow[51] # close enough at 583.95
m_df$se.flow[79]

#plot disease incidence

plot(m, 
     y = "se.flow",
     grid = TRUE, 
     legend="full",
     alpha = 0.7, 
     sim.lines = TRUE, 
     sim.alpha = 0.15, 
     ylim = c(0, 2000)
     )

# save best looking model as dataframe 

model_df <- as.data.frame(m, run = 4)

# plot best looking model

# normal scale
(out_plot = ggplot(model_df) +
    theme_minimal() + 
    geom_line(aes(x = time, y = ei.flow), color = "#c2a5cf") +
    labs(x ="Time (Days)", y = "Cases", 
         title = "Cases of COVID-19 over First 150 Days", 
         subtitle = "estimated cases") +
    xlim(0, 150) + 
    ylim(0, 2000))
# done

