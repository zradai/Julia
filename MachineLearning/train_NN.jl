###############################################
##### Neural network: logistic regression #####
###############################################

function train_NN(X, y;
    n::UInt16=UInt16(size(X)[2]),
    m::UInt32=UInt32(size(X)[1]),
    epoch::UInt16=UInt16(1000),
    α::Float16=Float16(0.01),
    λ::Float16=Float16(0),
    L::UInt8=UInt8(4),
    s::UInt8=UInt8(3),
    s_output::UInt8=UInt8(1),
    verbose=true ) # additional argument: activation function

    #include("/media/zozoo/3191D18644760ED0/Tan/MachineLearning/julia_ML/julia_ML_nested_functions.jl")
    include("ML_functions.jl")

    #n, m = UInt8(size(X)[2]), UInt16(size(X)[1])
    tt = build_matrices(L, s, n, s_output)
    costs = Vector{Float16}(undef, epoch)

    for e = 1:epoch

        if verbose
            println("Progress: ", e, " / ", epoch)
        end

        # calculating cost 'J'
        costs[e] = cost_function(X, y, tt)

        # Backpropagation: calculating node (unit within layer) values, then errors
        D_list = backprop(X, y, tt)

        # Gradient: changing theta matrices with calculated gradients, using the errors
        tt = apply_gradient(X, tt, D_list)

    end

    return tt, costs

end
