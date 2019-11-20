function rbinom(n, size, prob)
    # using Distributions

    d = Binomial(size, prob)
    v = rand(d, n)
    return v

end
