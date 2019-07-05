using Pkg, Random, Distributions

include("ML_functions.jl")
include("train_NN.jl")

# Data prep
const m, n, epoch, α, λ, L, s, s_output =
    1000, 20, UInt16(2500), Float16(0.01), Float16(0), UInt8(4), UInt(3), UInt8(1)
const X, y = dummy_data(m, n, y_family = "binomial");

# train NN
@time nn1 = train_NN(X, y, epoch=epoch, verbose=false);
