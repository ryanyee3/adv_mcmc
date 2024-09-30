# Date: 09/30/2024
# Author: Christian Varner
# Purpose: We need to compare MALA for Bayesian Logistic Regression
# versus HMC implemented in STAN.

using LinearAlgebra
using Plots
using Distributions

######################################
# Generate the data
######################################

# dimension parameters
n = 10000
p = 2

# data and true beta 
Ip = Diagonal(ones(p))
mean = zeros(p)
dist = MvNormal(mean, Ip)
beta_star = rand(dist)

Σ = Matrix(.25 * Diagonal(ones(p)))
mean = zeros(p)
dist = MvNormal(mean, Σ)
X = rand(dist, n)'

# generate observations
logit = X * beta_star
prob = zeros(n)
Y = zeros(n)
for i in 1:n
    prob[i] = 1 / (1 + exp(-logit[i]))
    Y[i] = prob[i] >= rand(1)[1]
end

######################################
# Implementation of MALA 
######################################

"""
Computes log( p(y_1,...,y_n | beta) )
"""
function log_likelihood(Y :: AbstractVector, X :: AbstractMatrix, beta :: AbstractVector)
    return sum(-log.(1 .+ exp.(-X * beta))) + dot( (Y .- 1), X * beta )
end

"""
Compute ∇_β log( p(y_1,...,y_n | beta) )
"""
function grad(Y :: AbstractVector, X :: AbstractMatrix, beta :: AbstractVector)
    g = zeros(size(X)[2])
    Yhat = 1 ./ (1 .+ exp.(-X * beta)) 
    g .= X' * (Y - Yhat)
    return g
end

"""
Computes the proposal β + .5 * λ^2 * ∇_β log( p(beta | Y) ) + ϵ 
"""
function proposal(Y :: AbstractVector, X :: AbstractMatrix, beta :: AbstractVector, λ :: Float64)
    g = grad(Y, X, beta)

    Ip = Matrix(λ^2 * Diagonal(ones(p)))
    mean = zeros(p)
    dist = MvNormal(mean, Ip)
    epsilon = rand(dist)

    return beta + .5 .* (λ^2) .* (g - beta) + epsilon
end

function MALA(Y :: AbstractVector, X :: AbstractMatrix, beta :: AbstractVector, λ :: Float64, NMCMC :: Int64)
    S = zeros(size(beta)[1], NMCMC)
    S[:, 1] .= beta
    beta_prev = zeros(size(beta)[1])
    for i in 2:NMCMC

        # get previous beta
        beta_prev .= S[:, i-1]

        # get proposal
        beta_proposal = proposal(Y, X, beta_prev, λ)

        # get acceptance ratio
        α = min(0, log_likelihood(Y, X, beta_proposal) - log_likelihood(Y, X, beta_prev) ) 
        U = log(rand(1)[1])

        # assign value to S
        if U <= α
            S[:, i] .= beta_proposal
        else
            S[:, i] .= beta_prev
        end
    end
    return S
end

######################################
# Simple Tests
######################################
using Zygote

beta = randn(2)
autog = Zygote.gradient((beta) -> log_likelihood(Y, X, beta), beta)[1]
g0 = grad(Y, X, beta) 
@assert autog ≈ g0

######################################
# Sample using MALA
######################################

P = 1 ./ (1 .+ exp.(-X * beta))
P = Diagonal(P)
Ip = Diagonal(ones(n))
H = -X' * P * (I - P) * X
e, _ = eigen(H)
λ = sqrt(2/maximum(abs.(e)))

beta = zeros(p)
S = MALA(Y, X, beta, λ, 20000)

######################################
# Visualizations -- MALA (2 dim)
######################################

p1 = histogram(S[1, 5000:20000], label = "Posterior for Beta0")
vline!([beta_star[1]])
p2 = histogram(S[2, 5000:20000], label = "Posterior for Beta1")
vline!([beta_star[2]])
plot(p1, p2, layout=(2, 1))

println("Beta estimated from MALA $( sum(S[1, 10000:20000])/10000), $(sum(S[2, 10000:20000])/10000)")

######################################
# Try MLE -- sanity check
######################################

beta0 = zeros(p)
for i in 1:5000
    beta0 .+= .5 * λ^2 * grad(Y, X, beta0)
end
println("Beta estimated from MLE $(beta0)")
println("Ending gradient norm $(norm(grad(Y, X, beta0)))")

######################################
# save data files
######################################

using DelimitedFiles
writedlm("Y.csv", Y, ",")
writedlm("X.csv", X, ",")
writedlm("beta_star.csv", beta_star, ",")
writedlm("S.csv", S', ",")
savefig("posterior_samples.png")