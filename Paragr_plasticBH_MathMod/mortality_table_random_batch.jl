morts = readdlm("Mortality_search_results.csv", ';')

batchsize = 1000 # number of rows to sample
random_rows_to_extract = rand(1:size(morts)[1], batchsize)

extr_mort = morts[random_rows_to_extract,:]
writedlm("Mortality_table_batches/Mortality_batch_V1.csv", extr_mort, ";")
