library(deSolve)

params = set_parameters()
IC = set_IC()
params = add_IC_to_params(params, IC)

times <- seq(0, 40, by = 0.1)

system.time(
  yout <- dede(y = IC, 
               times = times, 
               func = rhs_within_host, 
               parms = params)
)
