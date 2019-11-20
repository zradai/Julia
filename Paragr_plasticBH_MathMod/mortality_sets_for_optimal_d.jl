using Pkg, Compat, Random, Distributions, DataFrames, CSV, DelimitedFiles

mm = readdlm("Mortality_table_batches/Mortality_batch_V1.csv", ';')

##################################
###----- Run with p = 0.5 -----###
##################################

savepath_A = "p05/All_1000_sims_A/"
savepath_B = "p05/All_1000_sims_B/"
savepath_C = "p05/All_1000_sims_C/"

for i in 1:size(mm)[1]
    println(i)
    resA = simA(mm[i,:], 0.5, 1000, 1000)
    resB = simB(mm[i,:], 0.5, 1000, 1000)
    resC = simC(mm[i,:], 0.5, 1000, 1000)

    CSV.write(string(savepath_A, "FECS_simA_Mort",i,".csv"), resA)
    CSV.write(string(savepath_B, "FECS_simB_Mort",i,".csv"), resB)
    CSV.write(string(savepath_C, "FECS_simC_Mort",i,".csv"), resC)

end

#----------------------------------------------------------------------

######################################
###----- Re-run with p = 0.25 -----###
######################################

savepath_A = "p025/All_1000_sims_A_p025/"
savepath_B = "p025/All_1000_sims_B_p025/"
savepath_C = "p025/All_1000_sims_C_p025/"

for i in 1:size(mm)[1]
    println(i)
    resA = simA(mm[i,:], 0.25, 1000, 1000)
    resB = simB(mm[i,:], 0.25, 1000, 1000)
    resC = simC(mm[i,:], 0.25, 1000, 1000)

    CSV.write(string(savepath_A, "FECS_simA_Mort",i,"_p025.csv"), resA)
    CSV.write(string(savepath_B, "FECS_simB_Mort",i,"_p025.csv"), resB)
    CSV.write(string(savepath_C, "FECS_simC_Mort",i,"_p025.csv"), resC)

end

#----------------------------------------------------------------------

######################################
###----- Re-run with p = 0.75 -----###
######################################

savepath_A = "p075/All_1000_sims_A_p075/"
savepath_B = "p075/All_1000_sims_B_p075/"
savepath_C = "p075/All_1000_sims_C_p075/"

for i in 1:size(mm)[1]
    println(i)
    resA = simA(mm[i,:], 0.75, 1000, 1000)
    resB = simB(mm[i,:], 0.75, 1000, 1000)
    resC = simC(mm[i,:], 0.75, 1000, 1000)

    CSV.write(string(savepath_A, "FECS_simA_Mort",i,"_p075.csv"), resA)
    CSV.write(string(savepath_B, "FECS_simB_Mort",i,"_p075.csv"), resB)
    CSV.write(string(savepath_C, "FECS_simC_Mort",i,"_p075.csv"), resC)

end
