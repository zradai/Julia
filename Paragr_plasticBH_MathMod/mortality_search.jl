using Pkg, Compat, Random, Distributions, DataFrames, CSV, DelimitedFiles

#----------------------------------------

morts = zeros(4082925,4)

mvals = 0:0.01:1
global nr = 0
for mTot_slow_good in mvals
    for mTot_slow_bad in mvals
        for m1_rapid_good in mvals
            for m1_rapid_bad in mvals

                if mTot_slow_good > m1_rapid_good &&
                    mTot_slow_bad < m1_rapid_bad &&
                    mTot_slow_good < mTot_slow_bad &&
                    m1_rapid_good < m1_rapid_bad

                    global nr += 1
                    morts[nr,1] = mTot_slow_good
                    morts[nr,2] = mTot_slow_bad
                    morts[nr,3] = m1_rapid_good
                    morts[nr,4] = m1_rapid_bad

                end
            end
        end
    end
end

writedlm("Mortality_search_results.csv", morts, ";")
println("Done!")
