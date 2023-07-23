rhs_within_host_deSolve = function(t, x, p) {
  with(as.list(c(x, p)),{
    # Variables are (in order) V,S,I,R,D,F_U,F_B
    # So to get the lagvalue values...
    if (t<tau_I) {
      V_t = V0
      S_t = S0
      I_t = 0
      R_t = 0
    } else {
      V_t = lagvalue(t-tau_I,1)
      S_t = lagvalue(t-tau_I,2)
      I_t = lagvalue(t-tau_I,3)
      R_t = lagvalue(t-tau_I,4)
    }
    dV = p*I-d_V*V
    dS = lambda_S*(1-(S+I+D+R)/S_max)*S-beta*S*V
    dI = beta*S_t*V_t*(1-F_B/(epsilon_FI+F_B))*A_I-d_I*I
    dR = lambda_S*(1-(S+I+D+R)/S_max)*R +
      beta*S_t*V_t*(F_B/(epsilon_FI+F_B))*A_R
    dD = d_I*I-d_D*D
    dF_U = psi_F_prod+p_FI*I/(I+eta_FI) - 
      k_lin_f*F_U-k_B_F*((T_star+I)*A_F-F_B)*F_U+k_U_F*F_B
    dF_B = -k_int_f*F_B+k_B_F*((T_star+I)*A_F-F_B)*F_U-k_U_F*F_B
    dA_I = delta*A_I*(I_t-I)
    dA_R = delta*A_R*(R_t-R)
    return(list(c(dV, dS, dI, dR, dD, dF_U, dF_B,dA_I,dA_R)))
  })
}

set_IC = function() {
  # | Variable | Definition         | Value  | Unit              |
  # |----------|--------------------|--------|-------------------|
  # | V        | Viral load         | 4.5    | log_10(copies/ml) |
  # | S        | Susceptible cells  | 0.16   | 10^9 cells/ml     |
  # | I        | Infected cells     | 0      | 10^9 cells/ml     |
  # | R        | Resistant cells    | 0      | 10^9 cells/ml     |
  # | D        | Dead cells         | 0      | 10^9 cells/ml     |
  # | F_U      | Unbound interferon | 0.015  | pg/ml             |
  # | F_B      | Bound interferon   | 1.1e-8 | pg/ml             |
  IC = c(V = 1, S = 0.16, I = 0, R = 0, D = 0, 
         F_U = 0.015, F_B = 1.1e-8, A_I = 1, A_R = 1)
  return(IC)
}

set_parameters = function() {
  # | Parameter  | Definition | Value | Unit |
  # | lambda_S   |  Proliferation of susceptible cells   | 0.74  |  $day^{-1}$  |
  # | S_max      |  Target cell concentration   | 0.16  |  $10^9$cells/ml  |
  # | d_I        |  Death rate of infected cells   | 0.1  |  $day^{-1}$  |
  # | d_D        |  Degradation rate of dead cells   | 8  |  $day^{-1}$  |
  # | tau_I      |  Eclipse time   | 0.17  |  $day$ |
  # | beta       |  Viral infection rate   | 0.3 (SD:0.1994)  |  $day^{-1}cop/ml$  |
  # | d_V        | Viral decay rate    |  8.4 (SD:0.67)  |  $day^{-1}$ |
  # | p          |  Viral production rate  |   394 (SD:158.65) |  $day^{-1}(cop/10^9cells)$
  # | K_U_F      |  IFN unbinding rate   |  6.072  | $day^{-1}$  |
  # | p_FI       |  IFN production by infected cells   |  2.8235 (SD:1.8741) |  $day^{-1}(pg/ml)$  |
  # | psi_F_prod | IFN production by macrophages and monocytes  |  0.25 |  $day^{-1} (pg/ml)$  |
  # | k_B_F      | IFN binding rate  |  0.0107 |  $day^{-1} (ml/pg)$  |
  # | k_lin_f    |  IFN renal clearance rate  | 16.635 (SD:2.49)  |  $day^{-1}$  |
  # | k_int_f    | IFN internalization rate  |  16.968 (SD:2.54) |  $day^{-1}$  |
  # | epsilon_FI |  Half maximal response   | 2E-4 |  $10^9 cell/ml$  |
  # | eta_FI     | Half-maximal response  | 0.022  |  $10^9$cells/ml  |
  # | T_star     | Initial CD8+ T cells |  1.104E-4 |  $10^9$cells/ml |
  params = c(
    lambda_S = 0.74,
    S_max = 0.16,
    d_I = 0.1,
    d_D = 8,
    tau_I = 0.17,
    beta = 0.3,
    beta_stddev = 0.1994,
    d_V = 8.4,
    d_V_stddev = 0.67,
    p = 394,
    p_stddev = 158.65,
    k_U_F = 6.072,
    p_FI = 2.8235,
    p_FI_stddev = 1.8741,
    psi_F_prod = 0.25,
    k_B_F = 0.0107,
    k_lin_f = 16.635,
    k_lin_f_stddev = 2.49,
    k_int_f = 16.968,
    k_int_f_stddev = 2.54,
    epsilon_FI = 2e-4,
    eta_FI = 0.022328,
    T_star = 1.104e-4,
    delta = 0.1,
    avo=6.02214e23,
    MM_F = 19000,
    R_F_T = 1000,
    R_F_I = 1300
  )
  params = c(params,
             A_F = as.numeric(
               (params["MM_F"]/params["avo"]) *
                 (params["R_F_I"]+params["R_F_T"]) *
                 (1/5000)*(10^9*1e12)
             ))
  return(params)
}

add_IC_to_params = function(params, IC) {
  params = c(params, 
             V0 = as.numeric(IC["V"]), 
             S0 = as.numeric(IC["S"]),
             I0 = as.numeric(IC["I"]), 
             R0 = as.numeric(IC["R"]))
return(params)
}

# Given parameters, find those with a standard deviation
# given and generate a table with sampled values of these parameters, 
# regular values of the others, as well as initial conditions. 
# This way, we have all that's needed for the cohort.
generate_params_patients = function(params, n = 1000) {
  # The whole list of parameters, including values of std dev for some
  names_params = names(params)
  # Which are the parameters that contain std dev information
  idx_stddev = grep("stddev", names_params)
  # Which parameters is that std dev info for
  params_with_stddev = gsub("_stddev", 
                            "", 
                            names_params[idx_stddev])
  idx_params_with_stddev = which(names_params %in% params_with_stddev)
  OUT = 
    data.frame(
      mat.or.vec(nr = n,
                 nc = length(names_params)-length(idx_stddev))
    )
  colnames(OUT) = names_params[setdiff(1:length(names_params),
                                       idx_stddev)]
  for (curr_col in colnames(OUT)) {
    if (!(curr_col %in% params_with_stddev)) {
      # This is a "regular" parameter, we just replicate it n times
      OUT[[curr_col]] = rep(params[curr_col], n)
    } else {
      # This a parameter we must sample. As indicated, we sample from a 
      # normal distribution with mean the value given and std dev 3 times
      # the given std dev
      OUT[[curr_col]] = rnorm(n = n,
                              mean = params[curr_col],
                              sd = 3*params[sprintf("%s_stddev",
                                                    curr_col)])
      # This can give us negative values. If so, resample..
      idx_negative = which(OUT[[curr_col]]<0)
      if (length(idx_negative)>0) {
        repeat{
          # writeLines(paste0("Nb negative ", length(idx_negative)))
          OUT[[curr_col]][idx_negative] = 
            rnorm(n = length(idx_negative),
                  mean = params[curr_col],
                  sd = 3*params[sprintf("%s_stddev",
                                        curr_col)])
          idx_negative = which(OUT[[curr_col]]<0)
          if(length(idx_negative) == 0) {
            break
          }
        }
      }
      OUT[[curr_col]][which(OUT[[curr_col]]<0)] = params[curr_col]
    }
  }
  return(OUT)
}

# Given a patient index idx, a parameters data frame params.df and 
# initial conditions IC, run the simulation of the within host model for
# this patient
run_one_patient = function(idx = 1, 
                           params.df, 
                           IC) {
  writeLines(paste0("patient index = ", idx))
  params_tmp = params.df[idx,]
  params_tmp = add_IC_to_params(params_tmp, IC)
  times <- c(seq(0, ceiling(params_tmp$tau_I), by = 0.01), 
             seq(ceiling(params_tmp$tau_I), 200, by = 0.1))
  yout <- dede(y = IC, 
               times = times, 
               func = rhs_within_host_deSolve, 
               parms = params_tmp)
  return(yout)
}

# Function that sets up the data frame for plotting (changes names, formats, etc.)
# By default, prepares data with mean and 2.5 and 97.5 percentiles. Change as needed.
format_df = function(time,
                     data,
                     line_plotted = "mean", 
                     lower = "2.5%",
                     upper = "97.5%") {
  df = data.frame(time = time,
                  lower = data[,lower],
                  line = data[,line_plotted],
                  upper = data[,upper])
  return(df)
}


# From http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/#Helper%20functions
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}