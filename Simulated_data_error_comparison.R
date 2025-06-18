

col_names = c("L1","L2","L_inf","MSE","RMSE","L1_rel","L2_rel","L_inf_rel")


####-- OUR MODEL SIM 1 --####
error_matrix_ourModel_sim1_mis = matrix(nrow = 0, ncol = 8)
error_matrix_ourModel_sim1_mis_smoot = matrix(nrow = 0, ncol = 8)
error_matrix_ourModel_sim1_al = matrix(nrow = 0, ncol = 8)
error_matrix_ourModel_sim1_al_smoot = matrix(nrow = 0, ncol = 8)

colnames(error_matrix_ourModel_sim1_mis) = col_names
colnames(error_matrix_ourModel_sim1_mis_smoot) = col_names
colnames(error_matrix_ourModel_sim1_al) = col_names
colnames(error_matrix_ourModel_sim1_al_smoot) = col_names

for(seed in 1:100){
  
  file_name <- paste0("errore_sim1_OURmodel_seed", seed, ".RData")
  load(file_name)
  
  new_row_mis = c(errors$curves$y_mis$err$norm_1, errors$curves$y_mis$err$norm_2, errors$curves$y_mis$err$norm_inf,
                  errors$curves$y_mis$err$MSE, errors$curves$y_mis$err$RMSE, errors$curves$y_mis$err_rel$norm_1_rel,
                  errors$curves$y_mis$err_rel$norm_2_rel, errors$curves$y_mis$err_rel$norm_inf_rel)
  
  new_row_al = c(errors$curves$y_al$err$norm_1, errors$curves$y_al$err$norm_2, errors$curves$y_al$err$norm_inf,
                 errors$curves$y_al$err$MSE, errors$curves$y_al$err$RMSE, errors$curves$y_al$err_rel$norm_1_rel,
                 errors$curves$y_al$err_rel$norm_2_rel, errors$curves$y_al$err_rel$norm_inf_rel)
  
  #print(errors$curves$y_al$err$MSE)
  
  new_row_mis_smoot = c(errors$curves$y_mis_smoot$err$norm_1, errors$curves$y_mis_smoot$err$norm_2, errors$curves$y_mis_smoot$err$norm_inf,
                        errors$curves$y_mis_smoot$err$MSE, errors$curves$y_mis_smoot$err$RMSE, errors$curves$y_mis_smoot$err_rel$norm_1_rel,
                        errors$curves$y_mis_smoot$err_rel$norm_2_rel, errors$curves$y_mis_smoot$err_rel$norm_inf_rel)
  
  new_row_al_smoot = c(errors$curves$y_al_smoot$err$norm_1, errors$curves$y_al_smoot$err$norm_2, errors$curves$y_al_smoot$err$norm_inf,
                       errors$curves$y_al_smoot$err$MSE, errors$curves$y_al_smoot$err$RMSE, errors$curves$y_al_smoot$err_rel$norm_1_rel,
                       errors$curves$y_al_smoot$err_rel$norm_2_rel, errors$curves$y_al_smoot$err_rel$norm_inf_rel)
  
  error_matrix_ourModel_sim1_mis = rbind(error_matrix_ourModel_sim1_mis, new_row_mis)
  error_matrix_ourModel_sim1_mis_smoot = rbind(error_matrix_ourModel_sim1_mis_smoot, new_row_mis_smoot)
  error_matrix_ourModel_sim1_al = rbind(error_matrix_ourModel_sim1_al, new_row_al)
  error_matrix_ourModel_sim1_al_smoot = rbind(error_matrix_ourModel_sim1_al_smoot, new_row_al_smoot)
  
  rownames(error_matrix_ourModel_sim1_mis)[nrow(error_matrix_ourModel_sim1_mis)] = seed
  rownames(error_matrix_ourModel_sim1_mis_smoot)[nrow(error_matrix_ourModel_sim1_mis_smoot)] = seed
  rownames(error_matrix_ourModel_sim1_al)[nrow(error_matrix_ourModel_sim1_al)] = seed
  rownames(error_matrix_ourModel_sim1_al_smoot)[nrow(error_matrix_ourModel_sim1_al_smoot)] = seed
  
}

####--- OUR MODEL SIM 2 --####
error_matrix_ourModel_sim2_mis = matrix(nrow = 0, ncol = 8)
error_matrix_ourModel_sim2_mis_smoot = matrix(nrow = 0, ncol = 8)
error_matrix_ourModel_sim2_al = matrix(nrow = 0, ncol = 8)
error_matrix_ourModel_sim2_al_smoot = matrix(nrow = 0, ncol = 8)

colnames(error_matrix_ourModel_sim2_mis) = col_names
colnames(error_matrix_ourModel_sim2_mis_smoot) = col_names
colnames(error_matrix_ourModel_sim2_al) = col_names
colnames(error_matrix_ourModel_sim2_al_smoot) = col_names

for(seed in 1:100){
  
  file_name <- paste0("errore_sim2_OURmodel_seed", seed, ".RData")
  load(file_name)
  
  new_row_mis = c(errors$curves$y_mis$err$norm_1, errors$curves$y_mis$err$norm_2, errors$curves$y_mis$err$norm_inf,
                  errors$curves$y_mis$err$MSE, errors$curves$y_mis$err$RMSE, errors$curves$y_mis$err_rel$norm_1_rel,
                  errors$curves$y_mis$err_rel$norm_2_rel, errors$curves$y_mis$err_rel$norm_inf_rel)
  
  new_row_mis_smoot = c(errors$curves$y_mis_smoot$err$norm_1, errors$curves$y_mis_smoot$err$norm_2, errors$curves$y_mis_smoot$err$norm_inf,
                        errors$curves$y_mis_smoot$err$MSE, errors$curves$y_mis_smoot$err$RMSE, errors$curves$y_mis_smoot$err_rel$norm_1_rel,
                        errors$curves$y_mis_smoot$err_rel$norm_2_rel, errors$curves$y_mis_smoot$err_rel$norm_inf_rel)
  
  new_row_al = c(errors$curves$y_al$err$norm_1, errors$curves$y_al$err$norm_2, errors$curves$y_al$err$norm_inf,
                 errors$curves$y_al$err$MSE, errors$curves$y_al$err$RMSE, errors$curves$y_al$err_rel$norm_1_rel,
                 errors$curves$y_al$err_rel$norm_2_rel, errors$curves$y_al$err_rel$norm_inf_rel)
  
  new_row_al_smoot = c(errors$curves$y_al_smoot$err$norm_1, errors$curves$y_al_smoot$err$norm_2, errors$curves$y_al_smoot$err$norm_inf,
                       errors$curves$y_al_smoot$err$MSE, errors$curves$y_al_smoot$err$RMSE, errors$curves$y_al_smoot$err_rel$norm_1_rel,
                       errors$curves$y_al_smoot$err_rel$norm_2_rel, errors$curves$y_al_smoot$err_rel$norm_inf_rel)
  
  error_matrix_ourModel_sim2_mis = rbind(error_matrix_ourModel_sim2_mis, new_row_mis)
  error_matrix_ourModel_sim2_mis_smoot = rbind(error_matrix_ourModel_sim2_mis_smoot, new_row_mis_smoot)
  error_matrix_ourModel_sim2_al = rbind(error_matrix_ourModel_sim2_al, new_row_al)
  error_matrix_ourModel_sim2_al_smoot = rbind(error_matrix_ourModel_sim2_al_smoot, new_row_al_smoot)
  
  rownames(error_matrix_ourModel_sim2_mis)[nrow(error_matrix_ourModel_sim2_mis)] = seed
  rownames(error_matrix_ourModel_sim2_mis_smoot)[nrow(error_matrix_ourModel_sim2_mis_smoot)] = seed
  rownames(error_matrix_ourModel_sim2_al)[nrow(error_matrix_ourModel_sim2_al)] = seed
  rownames(error_matrix_ourModel_sim2_al_smoot)[nrow(error_matrix_ourModel_sim2_al_smoot)] = seed
  
}


####--- TEL MODEL SIM 1 --####

error_matrix_telModel_sim1_mis = matrix(nrow = 0, ncol = 8)
error_matrix_telModel_sim1_mis_smoot = matrix(nrow = 0, ncol = 8)
error_matrix_telModel_sim1_al = matrix(nrow = 0, ncol = 8)
error_matrix_telModel_sim1_al_smoot = matrix(nrow = 0, ncol = 8)

colnames(error_matrix_telModel_sim1_mis) = col_names
colnames(error_matrix_telModel_sim1_mis_smoot) = col_names
colnames(error_matrix_telModel_sim1_al) = col_names
colnames(error_matrix_telModel_sim1_al_smoot) = col_names

for(seed in 1:100){
  
  file_name_group1 <- paste0("errore_sim1_gruppo1_TELmodel_seed", seed, ".RData")
  load(file_name_group1)
  
  errors_mis_group1 = c(errors$curves$y_mis$err$norm_1, errors$curves$y_mis$err$norm_2, errors$curves$y_mis$err$norm_inf,
                        errors$curves$y_mis$err$MSE, errors$curves$y_mis$err$RMSE, errors$curves$y_mis$err_rel$norm_1_rel,
                        errors$curves$y_mis$err_rel$norm_2_rel, errors$curves$y_mis$err_rel$norm_inf_rel)
  
  
  errors_mis_smoot_group1 = c(errors$curves$y_mis_smoot$err$norm_1, errors$curves$y_mis_smoot$err$norm_2, errors$curves$y_mis_smoot$err$norm_inf,
                              errors$curves$y_mis_smoot$err$MSE, errors$curves$y_mis_smoot$err$RMSE, errors$curves$y_mis_smoot$err_rel$norm_1_rel,
                              errors$curves$y_mis_smoot$err_rel$norm_2_rel, errors$curves$y_mis_smoot$err_rel$norm_inf_rel)
  
  errors_al_group1 = c(errors$curves$y_al$err$norm_1, errors$curves$y_al$err$norm_2, errors$curves$y_al$err$norm_inf,
                       errors$curves$y_al$err$MSE, errors$curves$y_al$err$RMSE, errors$curves$y_al$err_rel$norm_1_rel,
                       errors$curves$y_al$err_rel$norm_2_rel, errors$curves$y_al$err_rel$norm_inf_rel)
  
  errors_al_smoot_group1 = c(errors$curves$y_al_smoot$err$norm_1, errors$curves$y_al_smoot$err$norm_2, errors$curves$y_al_smoot$err$norm_inf,
                             errors$curves$y_al_smoot$err$MSE, errors$curves$y_al_smoot$err$RMSE, errors$curves$y_al_smoot$err_rel$norm_1_rel,
                             errors$curves$y_al_smoot$err_rel$norm_2_rel, errors$curves$y_al_smoot$err_rel$norm_inf_rel)
  
  file_name_group2 <- paste0("errore_sim1_gruppo2_TELmodel_seed", seed, ".RData")
  load(file_name_group2)
  
  errors_mis_group2 = c(errors$curves$y_mis$err$norm_1, errors$curves$y_mis$err$norm_2, errors$curves$y_mis$err$norm_inf,
                        errors$curves$y_mis$err$MSE, errors$curves$y_mis$err$RMSE, errors$curves$y_mis$err_rel$norm_1_rel,
                        errors$curves$y_mis$err_rel$norm_2_rel, errors$curves$y_mis$err_rel$norm_inf_rel)
  
  errors_mis_smoot_group2 = c(errors$curves$y_mis_smoot$err$norm_1, errors$curves$y_mis_smoot$err$norm_2, errors$curves$y_mis_smoot$err$norm_inf,
                              errors$curves$y_mis_smoot$err$MSE, errors$curves$y_mis_smoot$err$RMSE, errors$curves$y_mis_smoot$err_rel$norm_1_rel,
                              errors$curves$y_mis_smoot$err_rel$norm_2_rel, errors$curves$y_mis_smoot$err_rel$norm_inf_rel)
  
  errors_al_group2 = c(errors$curves$y_al$err$norm_1, errors$curves$y_al$err$norm_2, errors$curves$y_al$err$norm_inf,
                       errors$curves$y_al$err$MSE, errors$curves$y_al$err$RMSE, errors$curves$y_al$err_rel$norm_1_rel,
                       errors$curves$y_al$err_rel$norm_2_rel, errors$curves$y_al$err_rel$norm_inf_rel)
  
  errors_al_smoot_group2 = c(errors$curves$y_al_smoot$err$norm_1, errors$curves$y_al_smoot$err$norm_2, errors$curves$y_al_smoot$err$norm_inf,
                             errors$curves$y_al_smoot$err$MSE, errors$curves$y_al_smoot$err$RMSE, errors$curves$y_al_smoot$err_rel$norm_1_rel,
                             errors$curves$y_al_smoot$err_rel$norm_2_rel, errors$curves$y_al_smoot$err_rel$norm_inf_rel)
  
  w1 = 0.5
  w2 = 0.5
  
  new_row_mis =  (errors_mis_group1 * w1 + errors_mis_group2 * w2) / (w1 + w2)
  new_row_mis_smoot =  (errors_mis_smoot_group1 * w1 + errors_mis_smoot_group2 * w2) / (w1 + w2)
  new_row_al =  (errors_al_group1 * w1 + errors_al_group2 * w2) / (w1 + w2)
  new_row_al_smoot =  (errors_al_smoot_group1 * w1 + errors_al_smoot_group2 * w2) / (w1 + w2)
  
  error_matrix_telModel_sim1_mis = rbind(error_matrix_telModel_sim1_mis, new_row_mis)
  error_matrix_telModel_sim1_mis_smoot = rbind(error_matrix_telModel_sim1_mis_smoot, new_row_mis_smoot)
  error_matrix_telModel_sim1_al = rbind(error_matrix_telModel_sim1_al, new_row_al)
  error_matrix_telModel_sim1_al_smoot = rbind(error_matrix_telModel_sim1_al_smoot, new_row_al_smoot)
  
  rownames(error_matrix_telModel_sim1_mis)[nrow(error_matrix_telModel_sim1_mis)] = seed
  rownames(error_matrix_telModel_sim1_mis_smoot)[nrow(error_matrix_telModel_sim1_mis_smoot)] = seed
  rownames(error_matrix_telModel_sim1_al)[nrow(error_matrix_telModel_sim1_al)] = seed
  rownames(error_matrix_telModel_sim1_al_smoot)[nrow(error_matrix_telModel_sim1_al_smoot)] = seed
  
}


####--- TEL MODEL SIM 2 --####

error_matrix_telModel_sim2_mis = matrix(nrow = 0, ncol = 8)
error_matrix_telModel_sim2_mis_smoot = matrix(nrow = 0, ncol = 8)
error_matrix_telModel_sim2_al = matrix(nrow = 0, ncol = 8)
error_matrix_telModel_sim2_al_smoot = matrix(nrow = 0, ncol = 8)

colnames(error_matrix_telModel_sim2_mis) = col_names
colnames(error_matrix_telModel_sim2_mis_smoot) = col_names
colnames(error_matrix_telModel_sim2_al) = col_names
colnames(error_matrix_telModel_sim2_al_smoot) = col_names

for(seed in 1:100){
  
  file_name <- paste0("errors_sim2_telMODEL_seed", seed, ".RData")
  load(file_name)
  
  new_row_mis = c(errors$curves$y_mis$err$norm_1, errors$curves$y_mis$err$norm_2, errors$curves$y_mis$err$norm_inf,
                  errors$curves$y_mis$err$MSE, errors$curves$y_mis$err$RMSE, errors$curves$y_mis$err_rel$norm_1_rel,
                  errors$curves$y_mis$err_rel$norm_2_rel, errors$curves$y_mis$err_rel$norm_inf_rel)
  
  new_row_mis_smoot = c(errors$curves$y_mis_smoot$err$norm_1, errors$curves$y_mis_smoot$err$norm_2, errors$curves$y_mis_smoot$err$norm_inf,
                        errors$curves$y_mis_smoot$err$MSE, errors$curves$y_mis_smoot$err$RMSE, errors$curves$y_mis_smoot$err_rel$norm_1_rel,
                        errors$curves$y_mis_smoot$err_rel$norm_2_rel, errors$curves$y_mis_smoot$err_rel$norm_inf_rel)
  
  new_row_al = c(errors$curves$y_al$err$norm_1, errors$curves$y_al$err$norm_2, errors$curves$y_al$err$norm_inf,
                 errors$curves$y_al$err$MSE, errors$curves$y_al$err$RMSE, errors$curves$y_al$err_rel$norm_1_rel,
                 errors$curves$y_al$err_rel$norm_2_rel, errors$curves$y_al$err_rel$norm_inf_rel)
  
  new_row_al_smoot = c(errors$curves$y_al_smoot$err$norm_1, errors$curves$y_al_smoot$err$norm_2, errors$curves$y_al_smoot$err$norm_inf,
                       errors$curves$y_al_smoot$err$MSE, errors$curves$y_al_smoot$err$RMSE, errors$curves$y_al_smoot$err_rel$norm_1_rel,
                       errors$curves$y_al_smoot$err_rel$norm_2_rel, errors$curves$y_al_smoot$err_rel$norm_inf_rel)
  
  error_matrix_telModel_sim2_mis = rbind(error_matrix_telModel_sim2_mis, new_row_mis)
  error_matrix_telModel_sim2_mis_smoot = rbind(error_matrix_telModel_sim2_mis_smoot, new_row_mis_smoot)
  error_matrix_telModel_sim2_al = rbind(error_matrix_telModel_sim2_al, new_row_al)
  error_matrix_telModel_sim2_al_smoot = rbind(error_matrix_telModel_sim2_al_smoot, new_row_al_smoot)
  
  rownames(error_matrix_telModel_sim2_mis)[nrow(error_matrix_telModel_sim2_mis)] = seed
  rownames(error_matrix_telModel_sim2_mis_smoot)[nrow(error_matrix_telModel_sim2_mis_smoot)] = seed
  rownames(error_matrix_telModel_sim2_al)[nrow(error_matrix_telModel_sim2_al)] = seed
  rownames(error_matrix_telModel_sim2_al_smoot)[nrow(error_matrix_telModel_sim2_al_smoot)] = seed
  
}


##-- Aligned curves errors --##

col = "MSE"  # actually the MSE is the L2 
data_names = c("OUR model sim 1", "TEL model sim 1 ", "OUR model sim 2", "TEL model sim 2")
ind = which(col_names == col)

data1 = error_matrix_ourModel_sim1_al[,ind]
data2 = error_matrix_telModel_sim1_al[,ind]
data3 = error_matrix_ourModel_sim2_al[,ind]
data4 = error_matrix_telModel_sim2_al[,ind]

df1 = data.frame(data1)
df2 = data.frame(data2)
df3 = data.frame(data3)
df4 = data.frame(data4)

colnames(df1) = col
colnames(df2) = col
colnames(df3) = col
colnames(df4) = col

combined_data_sim1 <- cbind(df1, df2)
combined_data_sim2 <- cbind(df3, df4)

colnames(combined_data_sim1) = c("","")
colnames(combined_data_sim2) = c("","")

x11()
boxplot(combined_data_sim2,
        col = c("blue", "red"),
        ylim = c(0, 2.8),
        names = c("A", "B"))




