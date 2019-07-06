
function dummy_data(m, n; y_family="gaussian")

    X = convert(Array{Float16,2}, rand(Normal(0, 1), m, n) )
    B = rand(Float16, size(X)[2])

    if y_family == "gaussian"
        y = X * B
    elseif y_family == "binomial"
        y = sigm_opt(X, B, true)
    end
    return X, y, B
end

#---------------------------------------------------------

function identity_link(X, θ) # activation function
    return X*θ
end

#---------------------------------------------------------

function sigm_opt(X, θ, round_y=false)
    # vectorized
#    if round_y ## round to 0 and 1
#        z0 = Int8.(round.(1 ./(1 .+ exp.(-z0))))
#    else ## not to round
#        z0 = Float16.(1 ./(1 .+ exp.(-z0)))
#    end

    z0 = X * θ
    if round_y ## round to 0 and 1
        for i = 1:size(X)[1]
            z0[i] = round(1 /(1 + exp(-z0[i])))
        end
    else ## not to round
        for i = 1:size(X)[1]
            z0[i] = 1 /(1 + exp(-z0[i]))
        end
    end
    return reshape(z0, size(X)[1], 1)
end

#---------------------------------------------------------

function build_matrices(L, s, s_input, s_output)
    if length(s) == 1 && L>3
        s = repeat([s], (L-2))
    end
    s = vcat(s_input, s, s_output)

    # building L-1 Theta matrices
    TT = Vector{Any}(undef, L-1)
    for t = 1:(L-1)
        dims = [s[t+1], s[t]+1] # row, col
        TT[t] = rand(Float32, dims[1], dims[2])
    end
    return TT
end

#---------------------------------------------------------

function predict_NN(tt, X, round_predy=false)
    L = length(tt)+1
    #s_input = UInt16(size(tt[1])[2] - 1)
    #s_output = UInt16(size(tt[ length(tt) ])[1])

    aL_s = []
    push!(aL_s, X)

    # creating matrices containing unit values
    for l = 2:L
        n_s = size(tt[l-1])[1]
        M_l = reshape(zeros(n_s*size(X)[1]), ( (size(X)[1]), n_s ) )
        for s = 1:n_s
            X_s = hcat(ones(size(aL_s[l-1])[1]), aL_s[l-1])
            tt_s = tt[l-1][s,:]
            M_l[:,s] = sigm_opt(X_s, tt_s, false)
        end
        push!(aL_s, M_l)
    end
    if round_predy
        return round.(aL_s[L])
    else
        return (aL_s[L])
    end

end

#---------------------------------------------------------

function Theta_no_bias(tt)
    tt_noBias = []
    for t = 1:length(tt)
        tt_t = convert(Array{Float32,2}, hcat(zeros( size(tt[t])[1] ), tt[t][:,2:end]) )
        push!(tt_noBias, tt_t)
    end
    return tt_noBias
end

#---------------------------------------------------------

function cost_function(X, y, tt, λ=0)

    m = size(X)[1]
    h0 = Array{Float32}(undef, size(y)[1])
    y_pred = predict_NN(tt, X, false)
    for i = 1:size(y)[1]
        h0[i] = y[i]*log(y_pred[i]) + (1-y[i])*log(1-y_pred[i])
    end
    #h0 = y .* log.(y_pred) .+ (1 .- y) .* log.(1 .- y_pred)

    J = (-1/m)*sum(h0)

    squared_thetas = 0
    tt_noBias = Theta_no_bias(tt)
    for b = 1:length(tt_noBias)
        squared_thetas += sum(tt_noBias[b] .^ 2)
    end

    regpar = (λ/(2m)) * squared_thetas
    cost = J + regpar

    return cost

end



#---------------------------------------------------------

function backprop(X, y, tt)

    L = UInt8(length(tt)+1)

    # array for the 'D' error matrices
    D_list = []
    for k in 1:length(tt)
        #push!(D_list, Array{Float64}(undef, size(tt[length(tt)-(k-1)]) ))
        push!(D_list, convert(Array{Float32, 2}, zeros(size(tt[length(tt)-(k-1)]) ) ) )
    end

    for i = 1:size(X)[1]

        # Forward propagation
        aL_s = Any[]
        push!(aL_s, reshape(vcat(1, X[i,:]), 1, size(X)[2]+1 ) )
        for l = 2:L

            n_s = size(tt[l-1])[1]
            M_l = reshape(zeros(n_s), ( 1, n_s ) )

            for s = 1:n_s
                X_s = aL_s[l-1]
                tt_s = tt[l-1][s,:]
                M_l[:,s] = sigm_opt(X_s, tt_s, false)
            end
            if l != L
                push!(aL_s, hcat(1, M_l))
            else
                push!(aL_s, M_l)
            end
        end

        # Backpropagation
        d_list = []
        push!(d_list, aL_s[L][1] - y[i])
        for d = 2:(L-1)
            ai = L-(d-1)


            if size(d_list[d-1])==()
                d_d = ( (tt[ai][:,2:end] * d_list[d-1]') .* (aL_s[ai][:,2:end].*(1 .- aL_s[ai][:,2:end]) ) )
            else
                d_d = ( (tt[ai][:,2:end] * d_list[d-1]')' .* (aL_s[ai][:,2:end].*(1 .- aL_s[ai][:,2:end]) ) )
            end

            push!(d_list, d_d)

        end
        temporary_dlist_1 = Array{Float64}(undef, 1, length(d_list[1]))
        temporary_dlist_1[1] = d_list[1]
        d_list[1] = temporary_dlist_1

        # Accumulate errors
        for D = 1:length(D_list)
            ai = L - D
            D_list[D] .+= (d_list[D]' * aL_s[ai])
        end

    end

    return D_list

end

#---------------------------------------------------------

function apply_gradient(X, tt, D_list, α=Float16(0.01), λ=Float16(0))

    m = size(X)[1]
    L = UInt8(length(tt)+1)

    # Gradients
    grads = []
    for r = 1:length(D_list)
        D_r = D_list[r]
        ai = L-r

        grad_noReg = D_r/m
        tt_r = tt[ai][:,2:end]
        zrs = zeros(size(tt_r)[1], 1)
        regM = hcat(zrs, (λ/m)*tt_r)

        push!(grads, (grad_noReg + regM) ) # originally was: grads = push!(grads, (grad_noReg + regM) ) ## I think this is redundant with 'grads='
    end

    # Apply gradients
    for g = 1:length(grads)
        ai = L-g
        tt[ai] = tt[ai] - α*grads[g]
    end

    return tt

end
