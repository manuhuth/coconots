addJuliaFunctions <- function(){

  if (JuliaConnectoR::juliaEval("isdefined(Main, :cocoReg)")) {
  } else {
    #----------------Add Julia functions-----------------------------------
    JuliaConnectoR::juliaEval('
    using Random
    using Distributions
    using Base.Threads
    using ForwardDiff
    using Optim
    using StatsBase
    using LineSearches
    using LinearAlgebra

    function exponential_function(x)
      """
      Wrapper function for the exponential function. Used as link function.
    
      # Arguments
      - `float::x`: function input.
      ...
      """
      return exp.(x)
    end
    
    function logistic_function(x)
      """
      Computes value of the logistic function.
    
      # Arguments
      - `Float::x`: value at which the logistic function should be evaluated.
      ...
      """
      return 1 ./ (1 .+ exp.(-x))
    end
    
    #---------------------------------Distributions---------------------------------
    function poisson_distribution(x, lambda)
        return exp(-lambda) * lambda^x / factorial(Int(x))
    end
    
    function generalized_poisson_distribution(x, lambda, eta)
        if x <= 20
            return exp(-lambda-x*eta) * lambda*(lambda+x*eta)^(x-1) / factorial(Int(x))
        else
            return exp(-lambda-x*eta) * lambda*(lambda+x*eta)^(x-1) / factorial(big(Int(x)))
        end
    end
    
    function bivariate_generalized_poisson(y, z, lambda, alpha1, alpha2, alpha3, eta)
        U = compute_U(alpha1, alpha2, alpha3)
        beta3 = compute_beta_i(lambda, U, alpha3)
        beta2 = compute_beta_i(lambda, U, alpha2)
        beta1 = compute_beta_i(lambda, U, alpha1)
    
        sum = 0.0
        if max(y,z) >= 20
            for j in 0:min(y,z)
                sum = sum + (lambda*U*(1-alpha1-alpha3) + eta*(y-j))^(y-j-1) *
                            (lambda*U*(1-alpha1-alpha3) + eta*(z-j))^(z-j-1) *
                            (lambda*U*(alpha1+alpha3) + eta*(j))^(j-1) /
                            factorial(big(Int(j))) / factorial(big(Int(y-j))) / factorial(big(Int(z-j))) * exp(j*eta)
            end
        else
            for j in 0:min(y,z)
                sum = sum + (lambda*U*(1-alpha1-alpha3) + eta*(y-j))^(y-j-1) *
                            (lambda*U*(1-alpha1-alpha3) + eta*(z-j))^(z-j-1) *
                            (lambda*U*(alpha1+alpha3) + eta*(j))^(j-1) /
                            factorial(Int(j)) / factorial(Int(y-j)) / factorial(Int(z-j)) * exp(j*eta)
            end
        end
    
        return sum * (beta1+beta3) * (lambda*U*(1-alpha1-alpha3))^2 * exp(-(beta1+beta3)-2*(lambda*U*(1-alpha1-alpha3))-y*eta-z*eta)
    end
    
    #--------------------------------Draw random variable---------------------------
    function draw_random_generalized_poisson_variable(u, lambda, eta)
        sum = 0.0
        i = 0
        while sum < u
            sum = sum + generalized_poisson_distribution(i, lambda, eta)
            i = i + 1
        end
        return i - 1
    end
    
    function draw_random_g_r_y_z_variable(u, y, z, lambda, alpha1, alpha2, alpha3, eta, max)
        sum = 0.0
        i = 0
        while sum < u
            sum = sum + compute_g_r_y_z(i, y, z, lambda, alpha1, alpha2, alpha3, eta, max)
            i = i + 1
        end
        return i - 1
    end
    
    function draw_random_g_r_y_variable(u, y, alpha, eta, lambda)
        sum = 0.0
        i = 0
        while sum < u
            sum = sum + compute_g_r_y(y, i, alpha, eta, lambda)
            i = i + 1
        end
        return i - 1
    end
    
    function draw_random_g(order, u, y, z, lambda, alpha1, alpha2, alpha3, alpha, eta, max)
        if order == 1
            return draw_random_g_r_y_variable(u, y, alpha, eta, lambda)
        elseif order == 2
            return draw_random_g_r_y_z_variable(u, y, z, lambda, alpha1, alpha2, alpha3, eta, max)
        end
    end
    
    function cocoSim(type, order, parameter, n, covariates=nothing,
                                link_function=exponential_function, n_burn_in=200,
                                x=zeros(Int(n + n_burn_in)))
        if !isnothing(covariates)
            lambda = link_function.(repeat(covariates * parameter[Int((end-(size(covariates)[2]-1))):Int(end)], Int(ceil(1 + n_burn_in / n)))[Int((end-n_burn_in-n+1)):Int(end)])
        else
            lambda = repeat([last(parameter)], Int(n + n_burn_in))
        end
    
        if order == 2
            alpha1 = parameter[1]
            alpha2 = parameter[2]
            alpha3 = parameter[3]
            alpha = nothing
    
            if type == "GP"
                eta = parameter[4]
            else
                eta = 0
            end
        else
            alpha1 = nothing
            alpha2 = nothing
            alpha3 = nothing
            alpha = parameter[1]
            if type == "GP"
                eta = parameter[2]
            else
                eta = 0
            end
        end
    
        for t in 3:(Int(n + n_burn_in))
            x[t] = draw_random_g(order, rand(Uniform(0,1),1)[1], Int(x[t-1]), Int(x[t-2]), lambda[t], alpha1, alpha2, alpha3, alpha, eta, nothing) +
                    draw_random_generalized_poisson_variable(rand(Uniform(0,1),1)[1], lambda[t], eta)
    
        end
    
        return x[Int((end-n+1)):Int(end)]
    end
    
    function compute_distribution_convolution_x_r_y(x, y, lambda, alpha, eta)
        if x < 0
            return 0
        end
        return sum([compute_convolution_x_r_y(i, y, lambda, alpha, eta) for i in 0:x])
    end
    
    function compute_distribution_convolution_x_r_y_z(x, y, z, alpha1, alpha2, alpha3,
                                                lambda, eta, max_loop=nothing)
        if x < 0
            return 0
        end
        return sum([compute_convolution_x_r_y_z(i, y, z, lambda,
                                alpha1, alpha2, alpha3, eta,
                                max_loop) for i in 0:x])
    end
    
    function cocoPit(cocoReg_fit, n_bins=21)
        cocoReg_fit["data"] = Int.(cocoReg_fit["data"])
        u = collect( range(0, stop = 1, length = Int(n_bins+1)) )
        
        lambda = get_lambda(cocoReg_fit, false)
        if isnothing(cocoReg_fit["covariates"])
          lambda = repeat([lambda], length(cocoReg_fit["data"]) )
        end
    
        if cocoReg_fit["order"] == 1
            if cocoReg_fit["type"] == "Poisson"
                eta = 0
            else
                eta = cocoReg_fit["parameter"][2]
            end
            Px = [compute_distribution_convolution_x_r_y(cocoReg_fit["data"][t],
                        cocoReg_fit["data"][t-1], lambda[t], cocoReg_fit["parameter"][1],
                         eta) for t in 2:length(cocoReg_fit["data"])]
            Pxm1 = [compute_distribution_convolution_x_r_y(cocoReg_fit["data"][t]-1,
                    cocoReg_fit["data"][t-1], lambda[t], cocoReg_fit["parameter"][1],
                    eta) for t in 2:length(cocoReg_fit["data"])]
        elseif cocoReg_fit["order"] == 2
            if cocoReg_fit["type"] == "Poisson"
                eta = 0
            else
                eta = cocoReg_fit["parameter"][4]
            end
            Px = [compute_distribution_convolution_x_r_y_z(cocoReg_fit["data"][t],
                                    cocoReg_fit["data"][t-1], cocoReg_fit["data"][t-2],
                                    cocoReg_fit["parameter"][1], cocoReg_fit["parameter"][2],
                                    cocoReg_fit["parameter"][3], lambda[t],
                                    eta, cocoReg_fit["max_loop"]) for t in 3:length(cocoReg_fit["data"])]
            Pxm1 = [compute_distribution_convolution_x_r_y_z(cocoReg_fit["data"][t]-1,
                        cocoReg_fit["data"][t-1], cocoReg_fit["data"][t-2],
                        cocoReg_fit["parameter"][1], cocoReg_fit["parameter"][2],
                        cocoReg_fit["parameter"][3], lambda[t],
                        eta, cocoReg_fit["max_loop"]) for t in 3:length(cocoReg_fit["data"])]
        end
    
        uniform_distribution = [get_pit_value(Px, Pxm1, u[s]) for s in 1:(Int(n_bins+1))]
    
        return Dict("Pit_values" => [uniform_distribution[s] - uniform_distribution[s-1] for s in 2:(Int(n_bins+1))],
                    "bins" => u[2:end])
    end
    
    function get_pit_value(Px, Pxm1, u)
        value = (u .- Pxm1) ./ (Px .- Pxm1)
        value[value .<  0] .= 0
        value[value .>  1] .= 1
    
        return mean(value)
    end
    
    function compute_scores(cocoReg_fit)
        
        lambda = get_lambda(cocoReg_fit, false)
        
        if isnothing(cocoReg_fit["covariates"])
          lambda = repeat([lambda], length(cocoReg_fit["data"]) )
        end
        
        if Int(cocoReg_fit["order"]) == 1
            if cocoReg_fit["type"] == "Poisson"
                eta = 0
            else
                eta = cocoReg_fit["parameter"][2]
            end
    
            probabilities = [compute_convolution_x_r_y(cocoReg_fit["data"][t],
                             cocoReg_fit["data"][t-1], lambda[t], cocoReg_fit["parameter"][1],
                             eta) for t in 2:length(cocoReg_fit["data"])]
            h_index = [compute_h_index_1(cocoReg_fit["data"][t-1], lambda[t],
                            cocoReg_fit["parameter"][1],
                             eta) for t in 2:length(cocoReg_fit["data"])]
    
            rbs = [compute_ranked_probability_helper_1(cocoReg_fit["data"][t],
                             cocoReg_fit["data"][t-1], lambda[t], cocoReg_fit["parameter"][1],
                             eta) for t in 2:length(cocoReg_fit["data"])]
    
        elseif Int(cocoReg_fit["order"]) == 2
            if cocoReg_fit["type"] == "Poisson"
                eta = 0
            else
                eta = cocoReg_fit["parameter"][4]
            end
    
            probabilities = [compute_convolution_x_r_y_z(cocoReg_fit["data"][t],
                            cocoReg_fit["data"][t-1], cocoReg_fit["data"][t-2], lambda[t],
                            cocoReg_fit["parameter"][1],
                            cocoReg_fit["parameter"][2], cocoReg_fit["parameter"][3],
                             eta, cocoReg_fit["max_loop"]) for t in 3:length(cocoReg_fit["data"])]
    
            h_index = [compute_h_index_2(cocoReg_fit["data"][t-1], cocoReg_fit["data"][t-2],
                                        cocoReg_fit["parameter"][1],
                            cocoReg_fit["parameter"][2], cocoReg_fit["parameter"][3],
                            lambda[t],
                             eta, cocoReg_fit["max_loop"]) for t in 3:length(cocoReg_fit["data"])]
    
            rbs = [compute_ranked_probability_helper_2(cocoReg_fit["data"][t],
                            cocoReg_fit["data"][t-1], cocoReg_fit["data"][t-2],
                            cocoReg_fit["parameter"][1],
                            cocoReg_fit["parameter"][2], cocoReg_fit["parameter"][3],
                            lambda[t],
                            eta, cocoReg_fit["max_loop"]) for t in 3:length(cocoReg_fit["data"])]
        end
    
        return Dict("logarithmic_score" => Float64(- sum(log.(probabilities)) / length(probabilities)),
                    "quadratic_score" => Float64(sum(- 2 .* probabilities .+ h_index) / length(probabilities)),
                    "ranked_probability_score" => Float64(sum(rbs) / length(probabilities))
                    )
    end
    
    function compute_h_index_1(y, lambda, alpha, eta)
        return sum(([compute_convolution_x_r_y(s, y, lambda, alpha, eta) for s in 0:30]).^2)
    end
    
    function compute_ranked_probability_helper_1(x, y, lambda, alpha, eta)
        dist_val = [compute_distribution_convolution_x_r_y(s, y, lambda, alpha, eta) for s in 0:30]
        return sum(((x .<= 0:(length(dist_val)-1)) .* (1 .- dist_val) .+ (x .> 0:(length(dist_val)-1)) .* dist_val).^2)
    end
    
    function compute_h_index_2(y, z, alpha1, alpha2, alpha3, lambda, eta, max_loop=nothing)
        return sum(([compute_convolution_x_r_y_z(s, y, z, lambda, alpha1, alpha2, alpha3,
                                                     eta, max_loop) for s in 0:30]).^2)
    end
    
    function compute_ranked_probability_helper_2(x, y, z, alpha1, alpha2, alpha3,
                                                lambda, eta, max_loop=nothing)
        dist_val = [compute_distribution_convolution_x_r_y_z(s, y, z, alpha1, alpha2, alpha3,
                                                    lambda, eta, max_loop) for s in 0:30]
        return sum(((x .<= 0:(length(dist_val)-1)) .* (dist_val .- 1) .+ (x .> 0:(length(dist_val)-1)) .* dist_val).^2)
    end
    #-------------------------------helper functions--------------------------------
    function compute_beta_i(lambda, U, alpha_i)
        return lambda*U*alpha_i
    end
    
    function compute_U(alpha1, alpha2, alpha3)
        return 1 / (1-alpha1-alpha2-alpha3)
    end
    
    function compute_zeta(lambda, U, alpha1, alpha3)
        return lambda * U * (1-2*alpha1-alpha3)
    end
    
    
    #--------------------------Likelihood relevant----------------------------------
    function compute_g_r_y_z(r, y, z, lambda, alpha1, alpha2, alpha3, eta, max)
        U = compute_U(alpha1, alpha2, alpha3)
        beta3 = compute_beta_i(lambda, U, alpha3)
        beta2 = compute_beta_i(lambda, U, alpha2)
        beta1 = compute_beta_i(lambda, U, alpha1)
        zeta = compute_zeta(lambda, U, alpha1, alpha3)
    
        if isnothing(max)
            max = y
        end
    
        sum = 0.0
        for s in 0:max, v in 0:max, w in 0:max
                    if ((r-s-v) >= 0) & ((z-r+v-w) >= 0) & ((y-s-v-w) >= 0)
                        sum = sum + generalized_poisson_distribution(s, beta3, eta) *
                                    generalized_poisson_distribution(v, beta1, eta) *
                                    generalized_poisson_distribution(w, beta1, eta) *
                                    generalized_poisson_distribution(r-s-v, beta2, eta) *
                                    generalized_poisson_distribution(z-r+v-w, lambda, eta) *
                                    generalized_poisson_distribution(y-s-v-w, zeta, eta)
                    end
        end
    
        return sum / bivariate_generalized_poisson(y, z, lambda, alpha1, alpha2, alpha3, eta)
    end
    
    function compute_g_r_y(y, r, alpha, eta, lambda)
        psi = eta * (1-alpha) / lambda
        if y < 20
            return factorial(Int(y)) / factorial(Int(r)) / factorial(Int(y-r)) * alpha *  (1-alpha) *
                    (alpha + psi*r)^(r-1) * (1-alpha+psi*(y-r))^(y-r-1) /
                    (1+ psi*y)^(y-1)
        else
            return factorial(big(Int(y))) / factorial(big(Int(r))) / factorial(big(Int(y-r))) * alpha *  (1-alpha) *
                    (alpha + psi*r)^(r-1) * (1-alpha+psi*(y-r))^(y-r-1) /
                    (1+ psi*y)^(y-1)
        end
    end
    
    function compute_convolution_x_r_y(x, y, lambda, alpha, eta)
        sum = 0.0
        for r in 0:min(x, y)
          if (y >= r) 
            sum = sum + compute_g_r_y(y, r, alpha, eta, lambda) * generalized_poisson_distribution(x-r, lambda, eta)
          end
        end
        return sum
    end
    
    function compute_negative_log_likelihood_GP1(lambdas, alpha, eta, data)
    
        sum = 0.0
        for t in 2:length(data)
            sum = sum - log(compute_convolution_x_r_y(data[t], data[t-1], lambdas[t], alpha, eta))
        end
    
        return sum
    end
    
    function compute_negative_log_likelihood_GP2(lambdas, alpha1, alpha2, alpha3, eta, data, max=nothing)
    
        sum = 0.0
        for t in 3:length(data)
            sum = sum - log(compute_convolution_x_r_y_z(data[t], data[t-1], data[t-2], lambdas[t], alpha1, alpha2, alpha3, eta, max))
        end
    
        return sum
    end
    
    function compute_convolution_x_r_y_z(x, y, z, lambda, alpha1, alpha2, alpha3, eta, max=nothing)
        sum = 0.0
        for r in 0:min(x, y+z)
            sum = sum + compute_g_r_y_z(r, y, z, lambda, alpha1, alpha2, alpha3, eta, max) * generalized_poisson_distribution(x-r, lambda, eta)
        end
        return sum
    end
    
    function get_lambda(cocoReg_fit, last_val=false)
        cocoReg_fit["link"] = exponential_function
        
        if isnothing(cocoReg_fit["covariates"])
            return last(cocoReg_fit["parameter"])
        else
            if last_val
                
                return cocoReg_fit["link"]( sum(cocoReg_fit["covariates"][end, :] .* cocoReg_fit["parameter"][(end-size(cocoReg_fit["covariates"])[2]+1):end]))
            else
                return cocoReg_fit["link"](cocoReg_fit["covariates"] * cocoReg_fit["parameter"][(end-size(cocoReg_fit["covariates"])[2]+1):end])
            end
        end
    end
    
    function cocoPredict(cocoReg_fit, x=0:10, covariates=nothing,
                          safe_array = Array{Float64}(undef, length(x)))

        lambda = get_lambda(cocoReg_fit, true)
        
        if (!isnothing(covariates))
          lambda = exp(sum(covariates[end,:] .* cocoReg_fit["parameter"][(end-size(cocoReg_fit["covariates"])[2]+1):end]))
        end
    
        if cocoReg_fit["order"] == 2
            if cocoReg_fit["type"] == "Poisson"
                eta = 0
            else
                eta = cocoReg_fit["parameter"][4]
            end
    
            output = Dict("probabilities" => [compute_convolution_x_r_y_z(i, Int(cocoReg_fit["data"][end]),
                                              Int(cocoReg_fit["data"][end-1]), lambda,
                                        cocoReg_fit["parameter"][1], cocoReg_fit["parameter"][2],
                                        cocoReg_fit["parameter"][3], eta,
                                        cocoReg_fit["max_loop"]) for i in x],
                         "prediction_mode" => -3.0)
        elseif cocoReg_fit["order"] == 1
            if cocoReg_fit["type"] == "Poisson"
                eta = 0
            else
                eta = cocoReg_fit["parameter"][2]
            end
    
            output = Dict("probabilities" => [compute_convolution_x_r_y(i, Int(cocoReg_fit["data"][end]),
                                                lambda, cocoReg_fit["parameter"][1], eta) for i in x],
                          "prediction_mode" => -3.0)
        end
    
        output["x"] = x
        output["prediction_mode"] = x[argmax(output["probabilities"])]
        output["prediction_median"] = x[findfirst(cumsum(output["probabilities"]) .>= 0.5)]
    
        return output
    end
    
    
    #------------------------Helper-------------------------------------------------
    function compute_inverse_matrix(M)
    
      return inv(M)
    end
    
    function compute_hessian(f, x)
      return ForwardDiff.hessian(f, x)
    end
    
    function compute_mu_hat_gmm(data)
        sum = 0.0
        x_bar = mean(data)
        for t in 3:length(data)
            sum = sum + (data[t] - x_bar) * (data[t-1] - x_bar) * (data[t-2] - x_bar)#
        end
        return sum / length(data)
    end
    
    function compute_autocorrelation(data, order)
        x_t = data[1:(length(data)-order)]
        x_lag = data[(order+1):length(data)]
    
        return cor(x_t, x_lag) #/ (var(x_t) * var(x_lag))^0.5
    end
    
    function set_to_unit_interval(x)
        return max(min(x, 0.9999), 0.0001)
    end
    
    function reparameterize_alpha(parameter)
        alpha3 = parameter[1] * logistic_function(parameter[2])
        alpha1 = (parameter[1] - alpha3) / 2
        alpha2 = (1-alpha1-alpha3) * logistic_function(parameter[3])
    
        return alpha1, alpha2, alpha3
    end
    
    #-------------------------Optimization relevant----------------------------------
    function minimize_pars_reparameterization_GP2(theta, data, covariates=nothing,
         link_function=exponential_function, max=nothing)
    
        lambdas = repeat([theta[5]], length(data))
        if  !isnothing(covariates)
            lambdas = link_function.(covariates * theta[5:5+(size(covariates, 2)-1)])
        end
    
        alpha1, alpha2, alpha3 = reparameterize_alpha(theta)
    
        return compute_negative_log_likelihood_GP2(lambdas, alpha1, alpha2, alpha3, theta[4], data, max)
    end
    
    function minimize_pars_reparameterization_Poisson2(theta, data, covariates=nothing,
         link_function=exponential_function, max=nothing)
    
        lambdas = repeat([theta[4]], length(data))
        if  !isnothing(covariates)
            lambdas = link_function.(covariates * theta[4:4+(size(covariates, 2)-1)])
        end
    
        alpha1, alpha2, alpha3 = reparameterize_alpha(theta)
    
        return compute_negative_log_likelihood_GP2(lambdas, alpha1, alpha2, alpha3, 0, data, max)
    end
    
    function minimize_pars_GP2(theta, data, covariates=nothing, link_function=exponential_function, max=nothing)
    
        lambdas = repeat([theta[5]], length(data))
        if  !isnothing(covariates)
            lambdas = link_function.(covariates * theta[5:5+(size(covariates, 2)-1)])
        end
    
        return compute_negative_log_likelihood_GP2(lambdas, theta[1], theta[2], theta[3], theta[4], data, max)
    end
    
    function minimize_pars_Poisson2(theta, data, covariates=nothing, link_function=exponential_function, max=nothing)
    
        lambdas = repeat([theta[4]], length(data))
        if  !isnothing(covariates)
            lambdas = link_function.(covariates * theta[4:4+(size(covariates, 2)-1)])
        end
    
        return compute_negative_log_likelihood_GP2(lambdas, theta[1], theta[2], theta[3], 0, data, max)
    end
    
    
    function minimize_pars_GP1(theta, data, covariates=nothing, link_function=exponential_function, max=nothing)
    
        lambdas = repeat([theta[3]], length(data))
        if  !isnothing(covariates)
            lambdas = link_function.(covariates * theta[3:3+(size(covariates, 2)-1)])
        end
    
        return compute_negative_log_likelihood_GP1(lambdas, theta[1], theta[2], data)
    end
    
    function minimize_pars_Poisson1(theta, data, covariates=nothing, link_function=exponential_function, max=nothing)
    
        lambdas = repeat([theta[2]], length(data))
        if  !isnothing(covariates)
            lambdas = link_function.(covariates * theta[2:2+(size(covariates, 2)-1)])
        end
    
        return compute_negative_log_likelihood_GP1(lambdas, theta[1], 0, data)
    end
    
    
    #------------------------Starting values -------------------------------------
    function get_starting_values!(type, order, data, covariates, starting_values, n_parameter_without)
        if isnothing(starting_values)
    
            if (type == "GP") & (order == 2)
                starting_values = compute_starting_values_GP2_reparameterized("GP", data)
                if !isnothing(covariates)
                    starting_values[n_parameter_without+1] = 0
                    starting_values = vcat(starting_values, repeat([0], size(covariates,2)-1))
                end
            end
    
            if (type == "Poisson") & (order == 2)
            starting_values = compute_starting_values_GP2_reparameterized("Poisson", data)
                if !isnothing(covariates)
                    starting_values[n_parameter_without+1] = 0
                    starting_values = vcat(starting_values, repeat([0], size(covariates,2)-1))
                end
            end
    
            if (type == "GP") & (order == 1)
                starting_values = compute_starting_values_GP1("GP", data)
                if !isnothing(covariates)
                    starting_values[n_parameter_without+1] = 0
                    starting_values = vcat(starting_values, repeat([0], size(covariates,2)-1))
                end
            end
    
            if (type == "Poisson") & (order == 1)
            starting_values = compute_starting_values_GP1("Poisson", data)
                if !isnothing(covariates)
                    starting_values[n_parameter_without+1] = 0
                    starting_values = vcat(starting_values, repeat([0], size(covariates,2)-1))
                end
            end
        end
        return starting_values
    end
    
    function compute_starting_values_GP2_reparameterized(type, data)
        if type == "GP"
            eta = compute_eta_starting_value(data)
        elseif type == "Poisson"
            eta = 0
        end
        alpha3 = set_to_unit_interval(compute_mu_hat_gmm(data) / mean(data) / (1+ 2* eta) * (1-eta)^4)
        alpha1 = set_to_unit_interval(compute_autocorrelation(data, 1))
        alpha2 = set_to_unit_interval(compute_autocorrelation(data, 2))
    
        if alpha3 + alpha2 + alpha1 >= 1
            while alpha3 + alpha2 + alpha1 >= 1
                alpha3 = alpha3 * 0.8
                alpha2 = alpha2 * 0.8
                alpha1 = alpha2 * 0.8
            end
        end
    
        if alpha3 + 2*alpha1 >= 1
            while alpha3 + 2*alpha1 >= 1
                alpha3 = alpha3 * 0.8
                alpha1 = alpha2 * 0.8
            end
        end
    
        lambda = mean(data) * (1-alpha1-alpha2-alpha3) * (1-eta)
    
        if type == "GP"
            return [2*alpha1 + alpha3, log(alpha3 / (2*alpha1) ), log(alpha2 / (1-alpha1-alpha2-alpha3)), eta, lambda]
        elseif type == "Poisson"
            return [2*alpha1 + alpha3, log(alpha3 / (2*alpha1) ), log(alpha2 / (1-alpha1-alpha2-alpha3)), lambda]
        end
    end
    
    function compute_starting_values_GP1(type, data)
    
        alpha = compute_autocorrelation(data, 1)
        lambda = mean(data) * (1-alpha)
    
        if type == "GP"
            eta = compute_eta_starting_value(data)
            return [alpha, eta, lambda]
        elseif type == "Poisson"
            eta = 0
            return [alpha, lambda]
        end
    end
    
    function compute_eta_starting_value(data)
        eta = 1 - (mean(data) / var(data))^0.5
        if eta < 0
            eta = 0.0001
        end
        return eta
    end
    
    #---------------------Bounds----------------------------------------------------
    function get_bounds_GP2(type, covariates)
    
        lambda_lower = [0]
        lambda_upper = [Inf]
        if !isnothing(covariates)
            lambda_lower = repeat([-Inf], size(covariates, 2) )
            lambda_upper = repeat([Inf], size(covariates, 2) )
        end
    
        if type == "GP"
            lower = vcat( [0,-10,-10], [0],  lambda_lower)
            upper = vcat( [1, 10,10], [1],  lambda_upper)
        elseif type == "Poisson"
            lower = vcat( [0,-10,-10], lambda_lower)
            upper = vcat( [1, 10,10],  lambda_upper)
        end
        return lower, upper
    end
    
    function get_bounds_GP1(type, covariates)
    
        lambda_lower = [0]
        lambda_upper = [Inf]
        if !isnothing(covariates)
            lambda_lower = repeat([-Inf], size(covariates, 2) )
            lambda_upper = repeat([Inf], size(covariates, 2) )
        end
    
        if type == "GP"
            lower = vcat( [0], [0],  lambda_lower)
            upper = vcat( [1], [1],  lambda_upper)
        elseif type == "Poisson"
            lower = vcat( [0], lambda_lower)
            upper = vcat( [1],  lambda_upper)
        end
        return lower, upper
    end
    
    function get_bounds(order, type, covariates)
        if order == 1
            return get_bounds_GP1(type, covariates)
        elseif order == 2
            return get_bounds_GP2(type, covariates)
        end
    end
    
    #-----------------------cocoreg---------------------------------------------
    function cocoReg(type, order, data, covariates=nothing, starting_values=nothing,
                    link_function=exponential_function, max_loop=nothing,
                    optimizer=Fminbox(LBFGS()))
    
        #-------------------------start dependent on type----------------------------------------------------------
        if order == 2
            if type == "GP"
                starting_values = get_starting_values!(type, order, Int.(data), covariates, starting_values, 4)
                fn = OnceDifferentiable(theta -> minimize_pars_reparameterization_GP2(theta, Int.(data),
                                                                                covariates,
                                                                                link_function,
                                                                                max_loop),
                                        starting_values,
                                        autodiff=:forward)
                f_alphas = theta -> minimize_pars_GP2(theta, Int.(data), covariates,
                                                    link_function, max_loop)
    
            elseif type == "Poisson"
                starting_values = get_starting_values!(type, order, Int.(data), covariates, starting_values, 3)
                fn = OnceDifferentiable(theta -> minimize_pars_reparameterization_Poisson2(theta, Int.(data),
                                                                                covariates,
                                                                                link_function,
                                                                                max_loop),
                                        starting_values,
                                        autodiff=:forward)
                f_alphas = theta -> minimize_pars_Poisson2(theta, Int.(data), covariates,
                                        link_function, max_loop)
            end
        end
        #-----------------------------------------Order 1 models-------------------------
        if order == 1
            if type == "GP"
                starting_values = get_starting_values!(type, order, Int.(data), covariates, starting_values, 2)
                fn = OnceDifferentiable(theta -> minimize_pars_GP1(theta, Int.(data),
                                                                                covariates,
                                                                                link_function),
                                        starting_values,
                                        autodiff=:forward)
                f_alphas = theta -> minimize_pars_GP1(theta,Int.(data), covariates,
                                                    link_function)
    
            elseif type == "Poisson"
                starting_values = get_starting_values!(type, order, Int.(data), covariates, starting_values, 1)
                fn = OnceDifferentiable(theta -> minimize_pars_Poisson1(theta, Int.(data),
                                                                                covariates,
                                                                                link_function),
                                        starting_values,
                                        autodiff=:forward)
                f_alphas = theta -> minimize_pars_Poisson1(theta, Int.(data), covariates,
                                        link_function)
            end
        end
        #--------------------------end dependent on type
    
        #write down constraints
        lower, upper = get_bounds(order, type, covariates)# -> make dependent on type and order
    
        #obtain fit
        fit = optimize(fn, lower, upper, starting_values,
                                            optimizer)#, Optim.Options(show_trace = true, show_every = 1,
                                            #iterations=1000, g_tol = 10^(-5), f_tol = 10^(-20)))
        parameter = Optim.minimizer(fit)
    
        #get alphas from reparameterized results
        if order == 2
            alpha1, alpha2, alpha3 = reparameterize_alpha(parameter)
            parameter[1:3] = [alpha1, alpha2, alpha3]
        end
    
        #construct output
        out = Dict("parameter" => parameter,
                   "covariance_matrix" => compute_inverse_matrix(compute_hessian(f_alphas, parameter)),
                   "log_likelihood" => -f_alphas(parameter),
                   "type" => type,
                   "order" => order,
                   "data" => data,
                   "covariates" => covariates,
                   "link" => link_function,
                   "starting_values" => starting_values,
                   "optimizer" => optimizer,
                   "lower_bounds" => lower,
                   "upper_bounds" => upper,
                   "optimization" => fit,
                   "max_loop" => max_loop
                   )
    
        #compute standard errors
        out["se"] = diag(out["covariance_matrix"]).^0.5
        return out
    end
    
    #----------------------cocoBoot-----------------------------------------------
    function compute_partial_autocorrelation(x, lags)
        return autocov(x, lags, demean=true)
    end
    
    function cocoBoot(cocoReg_fit, lags=[1:1:21;],
                     n_bootstrap=400, alpha=0.05, n_burn_in=200,
                     store_matrix = Array{Float64}(undef, length(lags), 2))
    
        pacfs = compute_random_pacfs(cocoReg_fit, lags, Int(n_bootstrap), n_burn_in,
                                    Array{Float64}(undef, Int(n_bootstrap), length(lags)))
        for i in 1:length(lags)
            store_matrix[ i, :] = quantile!(pacfs[:,i], [alpha/2, (1-alpha)/2])
        end
    
        pacf_data = compute_partial_autocorrelation(Int.(cocoReg_fit["data"]), lags)
    
        return Dict("upper" => store_matrix[:, 1],
                    "lower" => store_matrix[:, 2],
                    "in_interval" => (store_matrix[:, 1] .< pacf_data) .& (store_matrix[:, 2] .> pacf_data),
                    "pacf_data" => pacf_data,
                    "lags" => lags
             )
    end
    
    function compute_random_pacfs(cocoReg_fit, lags,
                     n_bootstrap=400, n_burn_in=200,
                     pacfs=Array{Float64}(undef, n_bootstrap, length(lags)))
        for b in 1:Int(n_bootstrap)
            pacfs[b,:] = compute_partial_autocorrelation(cocoSim(cocoReg_fit["type"], Int(cocoReg_fit["order"]), cocoReg_fit["parameter"],
                    length(cocoReg_fit["data"]), cocoReg_fit["covariates"], cocoReg_fit["link"],
                    n_burn_in, zeros(length(cocoReg_fit["data"]) + n_burn_in)), lags)
        end
        return pacfs
    end
    
    function create_julia_dict(keys, values)
      D = Dict()
      for i in 1:length(keys)
        D[keys[i]] = values[i]
      end
      return D
    end
    ')
    #-------------------------------------
  }
}