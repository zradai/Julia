function simA(mortalities, p, N, Y)

  # using Pkg, Random, Distributions, CSV, DataFrames, Statistics, LinearAlgebra

  # p: probability of good year
  # N: number of females in the simulated population
  # number of years to simulate    

  # mortalities should be a n×1 array (vector), containing
  # mortality values in the following order):
  ## mTot.slow.good, mTot.slow.bad, m1.rapid.good, m1.rapid.bad

  l = 0.5
  f_slow = 44*l
  f_rapid = 22*l
  f0 = f_slow

  mTot_slow_good = mortalities[1]
  mTot_slow_bad = mortalities[2]
  m1_rapid_good = mortalities[3]
  m1_rapid_bad = mortalities[4]

  m2_slow_good = mTot_slow_good
  m2_slow_bad = mTot_slow_bad
  m2_rapid_good = 1 - (1 - m1_rapid_good)*(0.2)
  m2_rapid_bad = 1 - (1 - m1_rapid_bad)*(0.2)

  fecvar = zeros((Y+1), N)
  fecvar[1,:] = rand(Uniform(0,1), N)

  for y in 1:Y

    s = rbinom(1, 1, p)[1]

    if s==1 # good year
      for n in 1:N

        dn = fecvar[1,n]

        t2_slow = rbinom(1, Int64(floor( (1-dn)*f0 )), (1-mTot_slow_good))[1] #
        t1 = rbinom(1, Int64(floor( dn*f0 )), (1-m1_rapid_good))[1] #
        t1_slow = rbinom(1, Int64(floor( t1*(1-dn)*f_rapid )), (1-m2_slow_good))[1] #
        t1_rapid = rbinom(1, Int64(floor( t1*dn*f_rapid )), (1-m2_rapid_good))[1] #
        t2_rapid = t1_slow + t1_rapid

        fecvar[(y+1),n] = t2_slow + t2_rapid

      end
    else # bad year
      for n in 1:N

        dn = fecvar[1,n]

        t2_slow = rbinom(1, Int64(floor( (1-dn)*f0 )), (1-mTot_slow_bad))[1] #
        t1 = rbinom(1, Int64(floor( dn*f0 )), (1-m1_rapid_bad))[1] #
        t1_slow = rbinom(1, Int64(floor( t1*(1-dn)*f_rapid )), (1-m2_slow_bad))[1] #
        t1_rapid = rbinom(1, Int64(floor( t1*dn*f_rapid )), (1-m2_rapid_bad))[1] #
        t2_rapid = t1_slow + t1_rapid

        fecvar[(y+1),n] = t2_slow + t2_rapid

      end
    end

  end

  arithmetic_means = zeros(N)
  standdevs = zeros(N)
  geometric_means = zeros(N)
  for i in 1:N
    arithmetic_means[i] = mean(fecvar[2:end,i])
    standdevs[i] = std(fecvar[2:end,i])
    geometric_means[i] = exp(sum(log.(fecvar[2:end,i] .+ 1))/Y)
  end

  fecs = DataFrame(female = 1:N, d = fecvar[1,:], mean = arithmetic_means,
  sd = standdevs, geom_mean = geometric_means)

  return fecs

end
