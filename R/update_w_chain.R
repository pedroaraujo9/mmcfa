update_w_chain = function(z,
                          beta,
                          w,
                          pw,
                          model_data) {

  beta = update_beta(
    beta = beta,
    z = z,
    w = w,
    model_data = model_data
  )

  # w = update_w(
  #   beta = beta,
  #   z = z,
  #   pw = pw,
  #   model_data = model_data
  # )

  pw = update_pw(
    w = w,
    model_data = model_data
  )

  out = list(
    beta = beta,
    w = w,
    pw = pw
  )

  return(out)

}
