# notes on EpiModel

pacman::p_load("EpiModel")

# define SEIR function 

SEIR <- function(t, t0, parms) {
  with(as.list(c(t0, parms)), {
    
    # Population size
    num <- s.num + e.num + i.num + r.num
    
    # Effective contact rate and FOI from a rearrangement of Beta * c * D
    ce <- R0 / i.dur
    lambda <- ce * i.num/num
    
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

p <- param.dcm(R0 = 1.0, 
               e.dur = 150, # est 150 people on an international flight into the country 
               i.dur = 1, # should be the number of cases at t0
               cfr = c(0.2, 0.4, 0.6, 0.8, 1.0))

#set initial state of n susceptible and infected at first time point (t1)
# let's say 500 people are susceptible, 1 are infected, and 0
# are recovered/removed from susceptibility
# initial state of 0 flowing from one compartment to another
# country wide model

i <- init.dcm(s.num = 144733050,
              e.num = 144733050,
              i.num = 1, 
              r.num = 0, 
              se.flow = 0, 
              ei.flow = 0, 
              ir.flow = 0, 
              d.flow = 0)

# collect model controls, model type and number of time steps (ex: in days) 
#    of the simulation, dt = fraction of time units, default = 1
# can add more controls as necessary for time lagged variables like interventions

c <- control.dcm(new.mod = SEIR, 
                 nsteps = 500, 
                 dt = 1, 
                 dede = TRUE)

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
     alpha = 0.7)

# network plot
comp_plot(m)
