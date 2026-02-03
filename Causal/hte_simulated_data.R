simulate = function() {

  # dimensions
  N = 250
  T = 35

  # factor loadings
  λ_1i = rnorm(N) 
  λ_2i = rnorm(N)

  # factors with stochastic drift
  a_t = rnorm(T)
  for (k in 2:length(a_t)) {
    a_t[k] = a_t[k - 1] + a_t[k]
  }
  f_1t = a_t + 0.1 * 1:T + 3 
  f_2t = rnorm(T)

  # time fixed effects with stochastic drift
  a_t = rnorm(T)
  for (k in 2:length(a_t)) {
    a_t[k] = a_t[k - 1] + a_t[k]
  }
  ξ_t = a_t + 0.1 * 1:T + 3 

  # pre-treatment periods
  α_i = rnorm(N) 
  ω_i = rnorm(N)
  tr_i = λ_1i + λ_2i + α_i + ω_i 
  tr_i = rank(tr_i) / length(tr_i)
  T_0i = fcase(tr_i <= .5, 35,
               .5 < tr_i & tr_i <= .6, 32,
               .6 < tr_i & tr_i <= .7, 29,
               .7 < tr_i & tr_i <= .8, 26,
               .8 < tr_i & tr_i <= .9, 23,
               .9 < tr_i, 20)

  # unit-specific data
  unit = data.table(
    unit = 1:N,
    T_0i,
    α_i,
    ω_i,
    λ_1i,
    λ_2i)

  # time-specific data
  time = data.table(
    time = 1:T,
    ξ_t,
    f_1t,
    f_2t)

  # panel data
  dat = unit[CJ(unit, time = 1:T), on = .(unit)]
  dat = time[dat, on = .(time)]

  dat[, ε  := rnorm(N * T)][
      , X1 := rnorm(N * T)][
      , X2 := rnorm(N * T)][
      , D  := fifelse(T_0i < time, 1, 0)][
      , δ  := 0.2 * (time - T_0i) + rnorm(N * T)][
      , Y  := 5 + 
              δ * D + 
              1 * X1 + 
              3 * X2 +
              λ_1i * f_1t +
              λ_2i * f_2t +
              α_i +
              ξ_t +
              ε]

  dat[, `Time after treatment` := time - T_0i]

  return(dat)
}