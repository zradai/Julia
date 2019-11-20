

scenarios = "A","B","C"
basic_path = "p05/All_1000_sims_"

max_records = 3000
strg = []
for t in 1:max_records
  append!(strg, randstring(1))
end
scens_intrest = DataFrame(MortPatt = zeros(max_records), GMF_max = zeros(max_records), d_gmfmax = zeros(max_records), scen = strg, strat = strg,
  B1=zeros(max_records), B_d=zeros(max_records), B_d2=zeros(max_records), B_strat=zeros(max_records), B_strat_d=zeros(max_records), B_strat_d2=zeros(max_records) )
global sc = 1

for s in scenarios
  for m in 1:1000

    println(s," -- ", m, " / 1000")

    readpath = string(basic_path, s, "/")
    mtab_sm = CSV.read(string(readpath, "FECS_sim", s, "_Mort",m,".csv"))

    if s=="C"
      XX = zeros(N,6)
      YY = zeros(N)
      for i in 1:N
        XX[i,1] = 1 # intercept
        XX[i,2] = mtab_sm.d[i] # d
        XX[i,3] = mtab_sm.d[i]^2 # d^2

        if mtab_sm.strat_A[i] == 1
          XX[i,4] = 1 # strat_A = 1
          XX[i,5] = mtab_sm.d[i] # d:strat_A=1
          XX[i,6] = mtab_sm.d[i]^2 # d^2:strat_A=1
        end
        YY[i] = mtab_sm.geom_mean[i]
      end
    else
      XX = zeros(N,3)
      YY = zeros(N)
      for i in 1:N
        XX[i,1] = 1 # intercept
        XX[i,2] = mtab_sm.d[i] # d
        XX[i,3] = mtab_sm.d[i]^2 # d^2
        YY[i] = mtab_sm.geom_mean[i]
      end
    end

    θ = pinv(XX' * XX) * XX' * YY # calculating linear regression model parameters

    if s == "C"
      newdata = zeros(202,6)
      for i in 1:size(newdata)[1]
        if i>101
          strat_i = 1
          mns = 101
        else
          strat_i = 0
          mns = 0
        end
        newdata[i,1] = 1
        newdata[i,2] = ((i-mns)-1)/100
        newdata[i,3] = (((i-mns)-1)/100)^2
        newdata[i,4] = strat_i
        newdata[i,5] = strat_i * ((i-mns)-1)/100
        newdata[i,6] = strat_i * (((i-mns)-1)/100)^2
      end
    else
      newdata = zeros(101,3)
      for i in 1:size(newdata)[1]
        newdata[i,1] = 1
        newdata[i,2] = (i-1)/100
        newdata[i,3] = ((i-1)/100)^2
      end
    end

    YY_newdata = newdata * θ

    if s == "C"
      newdata_df = DataFrame(d=newdata[:,2], GM_fitness=YY_newdata, strat_A=newdata[:,4])
    else
      newdata_df = DataFrame(d=newdata[:,2], GM_fitness=YY_newdata)
    end

    d_maxGM = newdata_df.d[newdata_df.GM_fitness.==findmax(newdata_df.GM_fitness)[1]][1]
    if d_maxGM > 0 &&
       d_maxGM < 1

       if s=="A"
         global strat_sm = "static"
       elseif s=="B"
         global strat_sm = "plastic"
       elseif s=="C"
         strat_pre = newdata_df.strat_A[newdata_df.GM_fitness.==findmax(newdata_df.GM_fitness)[1]]
         if strat_pre==1
           global strat_sm = "static"
         elseif strat_pre==0
           global strat_sm = "plastic"
         end
       end

       scens_intrest[sc,1] = m
       scens_intrest[sc,2] = findmax(newdata_df.GM_fitness)[1]
       scens_intrest[sc,3] = d_maxGM
       scens_intrest[sc,4] = s
       scens_intrest[sc,5] = strat_sm

       scens_intrest[sc,6] = θ[1]
       scens_intrest[sc,7] = θ[2]
       scens_intrest[sc,8] = θ[3]

       if s=="C"
         scens_intrest[sc,9] = θ[4]
         scens_intrest[sc,10] = θ[5]
         scens_intrest[sc,11] = θ[6]
       end

       global sc+=1

    end

  end

end

scens_intrest_v2 = scens_intrest[1:(sc-1),:]
CSV.write("p05/All_sims_Analyzed_V1.csv",
  scens_intrest_v2)

#----------------------------------------------------------------------

######################################
###----- Re-run with p = 0.25 -----###
######################################

scenarios = "A","B","C"
basic_path = "p025/All_1000_sims_"

max_records = 3000
strg = []
for t in 1:max_records
  append!(strg, randstring(1))
end
scens_intrest = DataFrame(MortPatt = zeros(max_records), GMF_max = zeros(max_records), d_gmfmax = zeros(max_records), scen = strg, strat = strg,
  B1=zeros(max_records), B_d=zeros(max_records), B_d2=zeros(max_records), B_strat=zeros(max_records), B_strat_d=zeros(max_records), B_strat_d2=zeros(max_records) )
global sc = 1

for s in scenarios
  for m in 1:1000

    println(s," -- ", m, " / 1000")

    readpath = string(basic_path, s,"_p025", "/")
    mtab_sm = CSV.read(string(readpath, "FECS_sim", s, "_Mort",m,"_p025.csv"))

    if s=="C"
      XX = zeros(N,6)
      YY = zeros(N)
      for i in 1:N
        XX[i,1] = 1 # intercept
        XX[i,2] = mtab_sm.d[i] # d
        XX[i,3] = mtab_sm.d[i]^2 # d^2

        if mtab_sm.strat_A[i] == 1
          XX[i,4] = 1 # strat_A = 1
          XX[i,5] = mtab_sm.d[i] # d:strat_A=1
          XX[i,6] = mtab_sm.d[i]^2 # d^2:strat_A=1
        end
        YY[i] = mtab_sm.geom_mean[i]
      end
    else
      XX = zeros(N,3)
      YY = zeros(N)
      for i in 1:N
        XX[i,1] = 1 # intercept
        XX[i,2] = mtab_sm.d[i] # d
        XX[i,3] = mtab_sm.d[i]^2 # d^2
        YY[i] = mtab_sm.geom_mean[i]
      end
    end

    θ = pinv(XX' * XX) * XX' * YY

    if s == "C"
      newdata = zeros(202,6)
      for i in 1:size(newdata)[1]
        if i>101
          strat_i = 1
          mns = 101
        else
          strat_i = 0
          mns = 0
        end
        newdata[i,1] = 1
        newdata[i,2] = ((i-mns)-1)/100
        newdata[i,3] = (((i-mns)-1)/100)^2
        newdata[i,4] = strat_i
        newdata[i,5] = strat_i * ((i-mns)-1)/100
        newdata[i,6] = strat_i * (((i-mns)-1)/100)^2
      end
    else
      newdata = zeros(101,3)
      for i in 1:size(newdata)[1]
        newdata[i,1] = 1
        newdata[i,2] = (i-1)/100
        newdata[i,3] = ((i-1)/100)^2
      end
    end

    YY_newdata = newdata * θ

    if s == "C"
      newdata_df = DataFrame(d=newdata[:,2], GM_fitness=YY_newdata, strat_A=newdata[:,4])
    else
      newdata_df = DataFrame(d=newdata[:,2], GM_fitness=YY_newdata)
    end

    d_maxGM = newdata_df.d[newdata_df.GM_fitness.==findmax(newdata_df.GM_fitness)[1]][1]
    if d_maxGM > 0 &&
       d_maxGM < 1

       if s=="A"
         global strat_sm = "static"
       elseif s=="B"
         global strat_sm = "plastic"
       elseif s=="C"
         strat_pre = newdata_df.strat_A[newdata_df.GM_fitness.==findmax(newdata_df.GM_fitness)[1]]
         if strat_pre==1
           global strat_sm = "static"
         elseif strat_pre==0
           global strat_sm = "plastic"
         end
       end

       scens_intrest[sc,1] = m
       scens_intrest[sc,2] = findmax(newdata_df.GM_fitness)[1]
       scens_intrest[sc,3] = d_maxGM
       scens_intrest[sc,4] = s
       scens_intrest[sc,5] = strat_sm

       scens_intrest[sc,6] = θ[1]
       scens_intrest[sc,7] = θ[2]
       scens_intrest[sc,8] = θ[3]

       if s=="C"
         scens_intrest[sc,9] = θ[4]
         scens_intrest[sc,10] = θ[5]
         scens_intrest[sc,11] = θ[6]
       end

       global sc+=1

    end

  end

end

scens_intrest_v2 = scens_intrest[1:(sc-1),:]
CSV.write("/media/zozoo/10C75B35626FCA66/Tan/PhD/A_P_agr/manuscript/Par_agr_MathMod/pBH_julia/p025/All_sims_Analyzed_p025_V1.csv",
  scens_intrest_v2)

#----------------------------------------------------------------------

######################################
###----- Re-run with p = 0.75 -----###
######################################

scenarios = "A","B","C"
basic_path = "p075/All_1000_sims_"

max_records = 3000
strg = []
for t in 1:max_records
  append!(strg, randstring(1))
end
scens_intrest = DataFrame(MortPatt = zeros(max_records), GMF_max = zeros(max_records), d_gmfmax = zeros(max_records), scen = strg, strat = strg,
  B1=zeros(max_records), B_d=zeros(max_records), B_d2=zeros(max_records), B_strat=zeros(max_records), B_strat_d=zeros(max_records), B_strat_d2=zeros(max_records) )
global sc = 1

for s in scenarios
  for m in 1:1000

    println(s," -- ", m, " / 1000")

    readpath = string(basic_path, s,"_p075", "/")
    mtab_sm = CSV.read(string(readpath, "FECS_sim", s, "_Mort",m,"_p075.csv"))

    if s=="C"
      XX = zeros(N,6)
      YY = zeros(N)
      for i in 1:N
        XX[i,1] = 1 # intercept
        XX[i,2] = mtab_sm.d[i] # d
        XX[i,3] = mtab_sm.d[i]^2 # d^2

        if mtab_sm.strat_A[i] == 1
          XX[i,4] = 1 # strat_A = 1
          XX[i,5] = mtab_sm.d[i] # d:strat_A=1
          XX[i,6] = mtab_sm.d[i]^2 # d^2:strat_A=1
        end
        YY[i] = mtab_sm.geom_mean[i]
      end
    else
      XX = zeros(N,3)
      YY = zeros(N)
      for i in 1:N
        XX[i,1] = 1 # intercept
        XX[i,2] = mtab_sm.d[i] # d
        XX[i,3] = mtab_sm.d[i]^2 # d^2
        YY[i] = mtab_sm.geom_mean[i]
      end
    end

    θ = pinv(XX' * XX) * XX' * YY

    if s == "C"
      newdata = zeros(202,6)
      for i in 1:size(newdata)[1]
        if i>101
          strat_i = 1
          mns = 101
        else
          strat_i = 0
          mns = 0
        end
        newdata[i,1] = 1
        newdata[i,2] = ((i-mns)-1)/100
        newdata[i,3] = (((i-mns)-1)/100)^2
        newdata[i,4] = strat_i
        newdata[i,5] = strat_i * ((i-mns)-1)/100
        newdata[i,6] = strat_i * (((i-mns)-1)/100)^2
      end
    else
      newdata = zeros(101,3)
      for i in 1:size(newdata)[1]
        newdata[i,1] = 1
        newdata[i,2] = (i-1)/100
        newdata[i,3] = ((i-1)/100)^2
      end
    end

    YY_newdata = newdata * θ

    if s == "C"
      newdata_df = DataFrame(d=newdata[:,2], GM_fitness=YY_newdata, strat_A=newdata[:,4])
    else
      newdata_df = DataFrame(d=newdata[:,2], GM_fitness=YY_newdata)
    end

    d_maxGM = newdata_df.d[newdata_df.GM_fitness.==findmax(newdata_df.GM_fitness)[1]][1]
    if d_maxGM > 0 &&
       d_maxGM < 1

       if s=="A"
         global strat_sm = "static"
       elseif s=="B"
         global strat_sm = "plastic"
       elseif s=="C"
         strat_pre = newdata_df.strat_A[newdata_df.GM_fitness.==findmax(newdata_df.GM_fitness)[1]]
         if strat_pre==1
           global strat_sm = "static"
         elseif strat_pre==0
           global strat_sm = "plastic"
         end
       end

       scens_intrest[sc,1] = m
       scens_intrest[sc,2] = findmax(newdata_df.GM_fitness)[1]
       scens_intrest[sc,3] = d_maxGM
       scens_intrest[sc,4] = s
       scens_intrest[sc,5] = strat_sm

       scens_intrest[sc,6] = θ[1]
       scens_intrest[sc,7] = θ[2]
       scens_intrest[sc,8] = θ[3]

       if s=="C"
         scens_intrest[sc,9] = θ[4]
         scens_intrest[sc,10] = θ[5]
         scens_intrest[sc,11] = θ[6]
       end

       global sc+=1

    end

  end

end

scens_intrest_v2 = scens_intrest[1:(sc-1),:]
CSV.write("/media/zozoo/10C75B35626FCA66/Tan/PhD/A_P_agr/manuscript/Par_agr_MathMod/pBH_julia/p075/All_sims_Analyzed_p075_V1.csv",
  scens_intrest_v2)
