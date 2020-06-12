# notes on EpiModel

pacman::p_load("EpiModel")

# parameterize #
# let's say our disease has a likelihood of transmission/single act of 0.3
#                           a rate of acts per person per time of 10 sneezes/ 60 minutes or 1 sneeze/6 minutes or 0.167
#                           a recovery rate of 1/14 days of disease duration
#                           an arrival rate/entry/birth rate of 2/100 (let's imagine people enter/exit the population at the same rate)
#                           a death rate of 2/100 among susceptible people /underlying pop death rate
#                           a death rate of 4/100  among infected people
#                           a departure rate of 2/100 among recovered people

p <- param.dcm(inf.prob = 0.3, act.rate = 1, rec.rate = 1/14, 
               a.rate = 2/100, ds.rate = 2/100,
               di.rate = 4/100, dr.rate = 2/100)

#set initial state of n susceptible and infected at first time point (t1)
# let's say 500 people are susceptible, 1 are infected, and 0 are recovered/removed from susceptibility

i <- init.dcm(s.num = 500, i.num = 1, r.num = 0)

# collect model controls, model type and number of time steps (ex: in days) of the simulation, dt = fraction of time units, default = 1
# can add more controls as necessary for time lagged variables like interventions

c <- control.dcm(type = "SIR", nsteps = 500)

# deterministic compartmental model #

# uses the deq solver in deSolve 
# produces param, control, and epi
# param: parameters passed to the model through the p data object of class param.dcm
# control: controls passed to the model through the c data object of class control.dcm
# epi: a list of data frames, one for each epidemiological output from the model
#     use as.data.frame.dcm to extract results as dataframes, 
#     summarize results with summary.dcm function on the m object of class dcm
#     plot the m object with plot.dcm
#     plot a flow diagram with comp_plot function on dcm object

m <- dcm(p, i, c)

# "Printing the model object provides basic information on model input and output,
# including model parameters and data variables for plotting and analysis." -tutorial
#     use as.data.frame to extract results as dataframes, 
#     summarize results with summary.dcm function on the m object of class dcm
#     plot the m object with plot.dcm
#     plot a flow diagram with comp_plot function on dcm object

m_df <- as.data.frame(m)

summary(m, at = 250)

#plot disease prevalence
# plot argument col takes values from RColorBrewer
plot(m, at = 10, 
     col = "RdYlBu", grid = TRUE, legend="lim")

# Plot each compartment alone
# prevalence of susceptibles
plot(m, y = "s.num", popfrac = TRUE, 
     col = "RdYlBu", grid = TRUE, legend="lim", 
     leg.name = "Prevalence of n Susceptible")

# Plot number of susceptibles
plot(m, y = "s.num", popfrac = FALSE, 
     col = "RdYlBu", grid = TRUE, legend="lim", 
     leg.name = "n Susceptible")

# Plot multiple runs of multiple compartments together
plot(m, y = c("s.num", "i.num"),
     run = 1, xlim = c(0, 50), 
     col = "RdYlBu", grid = TRUE, legend="lim")

# add = TRUE adds the graph created below to the graph created above
plot(m, y = c("s.num", "i.num"),
     run = 1, lty = 2, add = TRUE, col = "RdYlBu", grid = TRUE, 
     legend="n")













(m_compplot <- comp_plot(m))
