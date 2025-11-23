# Configuration for the project

config <- list(
  cases = c("correlation"), # pick from c("lag", "correlation", "highrho_simple", "dif", "simple")
  methods = c("fe", "spj"), # pick from c("dgmm", "fe", "spj", "db", "le")
  twoway_cluster = c(TRUE), # pick from c(FALSE, TRUE)
  cluster_method = "driscoll_kraay", # pick from c("driscoll_kraay", "twoway")
  num_iterations = 1000,
  ar_coefficients = list(
    correlation = c(0, 0.2, 0.5, 0.8),
    lag = c(0.2, 0.4, 0.5, 0.6),
    highrho_simple = c(0.9, 0.95), 
    dif = c(get_empirical_arcoef()),
    simple = c(0, 0.2, 0.5, 0.8)
  ),
  beta = list(
    correlation = -0.6, 
    lag = -0.25,
    highrho_simple = -0.6,
    dif = -0.25,
    simple = -0.6
  ),
  parameter_combinations = list(
    correlation = list(c(30, 60), c(30, 120), c(50, 120)),
    lag = list(c(30, 60), c(30, 120), c(50, 120)),
    highrho_simple = list(c(30, 60), c(30, 120), c(50, 120)), 
    dif = list(c(82, 36)), # from empirical data
    simple = list(c(30, 60), c(30, 120), c(50, 120))
  ),
  A = list(
    correlation = 0, 
    lag = 0,
    highrho_simple = 0,
    dif = 0,
    simple = 0
  )
)