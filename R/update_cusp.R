update_cusp = function(alpha, omega, nu, a, b, min_var, H) {

  J = nrow(cbind(alpha))

  if(H > 1) {

    l_vec = colSums(dnorm(alpha, mean = 0, sd = sqrt(min_var), log = TRUE))
    u_vec = LaplacesDemon::dmvt(t(alpha), mu = rep(0, J), S = (b/a)*diag(J), df = 2*a, log = T)

    ind_logprob = matrix(0, nrow = H, ncol = H)

    for(h in 1:H) {

      ind_logprob[h, ] = log(omega) + ifelse(1:H <= h, l_vec[h], u_vec[h])
      ind_logprob[h, ] = ind_logprob[h, ] - mclust::logsumexp(ind_logprob[h, ])

    }

    ind = extraDistr::rcatlp(n = H, log_prob = ind_logprob) + 1

    sum_ind = table(factor(ind, levels = 1:H))[-H]
    sum_greater_ind = rowSums(ind > matrix(1:(H-1), nrow = H-1, ncol = H, byrow = FALSE))

    v = rbeta(H-1, 1 + sum_ind, nu + sum_greater_ind)
    v = c(v, 1)
    omega = v * c(1, cumprod(1 - v[-length(v)]))

    prec = ifelse(
      ind <= 1:H,
      1/min_var,
      rgamma(H, shape = a + 0.5*J, rate = b + 0.5*colSums(alpha^2))
    )


  }else{


    ind = 1
    v = rbeta(1, 1 + 1, nu)
    omega = 1
    prec = rgamma(1, shape = a + 0.5*J, rate = b + 0.5*colSums(alpha^2))


  }


  out = list(

    prec = prec,
    v = v,
    omega = omega,
    ind = ind

  )

  return(out)

}
