library(deSolve)
# library(dde)   # Not used here, but could try

params = set_parameters()
IC = set_IC_2()
params = add_IC_to_params(params, IC)

times <- seq(0, 200, by = 0.1)

system.time(
  yout <- dede(y = IC, 
               times = times, 
               func = rhs_within_host_deSolve, 
               parms = params)
)
