pacman::p_load("EpiModel")

epiweb("icm")

# Replicate the Poker Chip SI model
# initial parameters
p1 <-param.icm(inf.prob = 1, act.rate = 0.2)
i1<- init.icm(s.num = 9, i.num = 1)
c1 <- control.icm(type = "SI", nsims = 1, nsteps = 50)
m1 <- icm(p1, i1, c1)

# number of susceptibles, number infected, and si flow change with each run 
# while all other parameters stay the same
# the model in general changed each time I ran it because the number of susceptible people
# who became infected did so more quickly or more slowly with each run

# 25 simulations
p2 <-param.icm(inf.prob = 1, act.rate = 0.2)
i2 <- init.icm(s.num = 9, i.num = 1)
c2 <- control.icm(type = "SI", nsims = 1, nsteps = 50)
m2 <- icm(p2, i2, c2)

# This set of simulations could be summarized by stating that the number of infected
# individuls begins to exceed the number of susceptible individuals between days 10 and 20 since
# identification of the index case. By day 40, all susceptible individuals will have been infected
# according to this model

# The relative std deviation of the incidence rise with the mean of the incidence. This indicates that
# the uncertainty in the outcome increases as incidence rises

# Stochastic SI model
# 25 simulations
p2 <-param.icm(inf.prob = 1, act.rate = 0.2)
i2 <- init.icm(s.num = 90, i.num = 10)
c2 <- control.icm(type = "SI", nsims = 25, nsteps = 50)
m2 <- icm(p2, i2, c2)

# The std. deviation of the incidence is 17% of the mean of the incidence
# This tells us that there are between 42-58 people infected at this time step
# and since the standard deviation appears to decrease after this point, 
# this suggests that there is less uncertainty around this measurement than
# in the model with fewer people


# done
