library(stringr) 
library(randomForest)
library(ggplot2)
library("rstudioapi")

get_abundance_table <- function(reads_table) {
   reads_table_abundance = reads_table
   reads_table_abundance = sapply(1: ncol(reads_table), function(j) (reads_table_abundance[,j] <- reads_table[,j] / colSums(reads_table)[j] ))
   colnames(reads_table_abundance) = colnames(reads_table)
   row.names(reads_table_abundance) = row.names(reads_table)
   
   return(reads_table_abundance)
}

##### get and prepare data #########
path_input = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path_input)
load("./MV_BV_model.RData")
RF2_input_colnames_list = read.csv('RF2_input_colnames_list.csv', header = T, row.names = 1)
RF2_input_colnames_list = RF2_input_colnames_list$x

reads_table_all = read.csv('40168_2017_305_MOESM2_ESM_reads_table_genus.csv')
metadata_all = read.csv('40168_2017_305_MOESM2_ESM_metadata.csv')
metadata_all$ID = colnames(reads_table_all)[-1]

# convert to genus level  
genus_list = reads_table_all[,1]
genus_list[genus_list == "Lactobacillus"] = "Other_Lactobacillus"
genus_list = str_replace(genus_list,'xxx','_')

genus_list_2 = unique(genus_list)

reads_table_all_genus = as.data.frame(matrix(data = 0, ncol = ncol(reads_table_all)-1, nrow = length(genus_list_2)))
colnames(reads_table_all_genus) = colnames(reads_table_all)[-1]
row.names(reads_table_all_genus) = genus_list_2
reads_table_all = reads_table_all[,-1]
for (a in 1: nrow(reads_table_all_genus)) {
   n = which(genus_list == genus_list_2[a])
   reads_table_all_genus[a,] = colSums(reads_table_all[which(genus_list == genus_list_2[a]),])
}

out_liers = read.csv('outlier_boxplot.stats_genus.csv')
keep = row.names(reads_table_all_genus) %in% out_liers$out_genus
sum(keep)
reads_table_all_genus = reads_table_all_genus[!keep,]
write.csv(reads_table_all_genus,'reads_table_all_genus_rm_outlier.csv')

##### predict using urine directly #####
library(pROC)
library(doSNOW)
library(parallel)
library(matrixStats)
set.seed(2022)
reads_table_all_genus_abundance = reads_table_all_genus
reads_table_all_genus_abundance = sapply(1: ncol(reads_table_all_genus), function(j) (reads_table_all_genus_abundance[,j] <- reads_table_all_genus[,j] / colSums(reads_table_all_genus)[j] ))
colnames(reads_table_all_genus_abundance) = colnames(reads_table_all_genus)
row.names(reads_table_all_genus_abundance) = row.names(reads_table_all_genus)

# remove samples with total predicted abundance less than x% ?
reads_table_all_genus_abundance_2 = reads_table_all_genus_abundance[row.names(reads_table_all_genus_abundance) %in% RF2_input_colnames_list,]
y = as.numeric(colSums(reads_table_all_genus_abundance_2))

abundance_impact_2_ML = data.frame(Total_abundance_in_model = c(seq(0,0.9,by=0.33),'all'), auROC = NA, pvalue = NA, case_number = NA, difference_in_median_values = NA)

for (z in 1:4) {
   reads_table = reads_table_all_genus[row.names(reads_table_all_genus) %in% RF2_input_colnames_list,]
   metadata = metadata_all
   
   if (z %in% c(1:3)) {
      x = 0.33*z - 0.33
      keep = y> x & y <= x+0.33
      
      reads_table = reads_table_all_genus[row.names(reads_table_all_genus) %in% RF2_input_colnames_list,keep]
      metadata = metadata_all[keep,]
   }
   
   
   keep = colSums(reads_table) >= 500
   sum(keep)
   reads_table = reads_table[,keep]
   metadata = metadata[keep,]
   
   # prediction
   reads_table = get_abundance_table(reads_table) 
   reads_table = as.matrix(t(reads_table))
   
   keep = !is.na(metadata$Nugent.Score)
   sum(keep)
   reads_table = reads_table[keep,]
   metadata = metadata[keep,]
   real_value = metadata$Nugent.Score
   keep = real_value >= 7
   real_value[keep] = 'Yes'
   real_value[!keep] = 'No'
   
   if (length(unique(real_value)) <2) {
      next
   }
   
   abundance_impact_2_ML$case_number[z] = paste0(length(real_value)," (",sum(real_value == 'Yes'),')')
   
   model_input = as.data.frame(matrix(data = 0, ncol = length(RF2_input_colnames_list), nrow = nrow(reads_table)))
   colnames(model_input) = RF2_input_colnames_list
   row.names(model_input) = row.names(reads_table)
   
   for (a in 1: ncol(model_input)) {
      n = which(colnames(reads_table) == colnames(model_input)[a])
      if (length(n) == 0) {
         next
      }
      model_input[,a] = reads_table[,n]
   }
   
   ppp = sapply(1:nrow(model_input), function(j) (predict(RF2, model_input[j,], type='prob')[1] ))
   
   roc = roc(real_value, ppp)
   abundance_impact_2_ML$auROC[z] = roc$auc
   
   abundance_impact_2_ML$pvalue[z] = kruskal.test(real_value, ppp, paired = F)$p.value
   
   abundance_impact_2_ML$difference_in_median_values[z] = median(ppp[real_value == 'No']) - median(ppp[real_value == 'Yes'])
   
   if (z ==0) {
      data = cbind(roc$sensitivities,roc$specificities)
      data = as.data.frame(data)
      colnames(data) = c('Sensitivities','Specificities')
      data = data[order(data$Specificities, decreasing = F),]
      
      ggplot(data, aes(x = Sensitivities, y = Specificities)) +
         geom_point(size = 1)+geom_step()+scale_x_reverse()+
         theme(axis.title = element_text(size = 6), 
               axis.text = element_text(size = 6), 
               legend.text = element_text(size = 6), 
               legend.title = element_text(size = 6)) + theme_bw() 
      ggsave('Sensitivities_Specificities_final_models_dx_bv_input_urine.pdf',width=3, height=3)
      
      data = data.frame(Real_value=real_value, Predicted_value= ppp)
      
      ggplot(data, aes(x=Real_value, y=Predicted_value,fill=Real_value)) +
         geom_boxplot()+ theme_bw() 
      ggsave('model_final_dx_bv_input_urine.pdf',width=2.5, height=3)
   }
   
}

write.csv(abundance_impact_2_ML,'model_final_dx_bv_input_urine.csv')


abundance_impact_2_ML = data.frame(Total_abundance_in_model = c(seq(0,0.9,by=0.33),'all'), auROC = NA, pvalue = NA, case_number = NA, difference_in_median_values = NA)

for (z in 1:4) {
   reads_table = reads_table_all_genus[row.names(reads_table_all_genus) %in% RF2_input_colnames_list,]
   metadata = metadata_all
   
   if (z %in% c(1:3)) {
      x = 0.33*z - 0.33
      keep = y> x & y <= x+0.33
      
      reads_table = reads_table_all_genus[row.names(reads_table_all_genus) %in% RF2_input_colnames_list,keep]
      metadata = metadata_all[keep,]
   }
   
   keep = colSums(reads_table) >= 500
   sum(keep)
   reads_table = reads_table[,keep]
   metadata = metadata[keep,]
   
   # prediction
   reads_table = get_abundance_table(reads_table) 
   reads_table = as.matrix(t(reads_table))
   
   keep = !is.na(metadata$Nugent.Score)
   sum(keep)
   reads_table = reads_table[keep,]
   metadata = metadata[keep,]
   real_value = metadata$Biofilm.on.Epithelial.cells
   real_value[real_value == 'pos'] = 'Yes'
   real_value[real_value == 'neg'] = 'No'
   
   if (length(unique(real_value)) <2) {
      next
   }
   
   abundance_impact_2_ML$case_number[z] = paste0(length(real_value)," (",sum(real_value == 'Yes'),')')
   
   model_input = as.data.frame(matrix(data = 0, ncol = length(RF2_input_colnames_list), nrow = nrow(reads_table)))
   colnames(model_input) = RF2_input_colnames_list
   row.names(model_input) = row.names(reads_table)
   
   for (a in 1: ncol(model_input)) {
      n = which(colnames(reads_table) == colnames(model_input)[a])
      if (length(n) == 0) {
         next
      }
      model_input[,a] = reads_table[,n]
   }
   
   ppp = sapply(1:nrow(model_input), function(j) (predict(RF2, model_input[j,], type='prob')[1] ))
   
   roc = roc(real_value, ppp)
   abundance_impact_2_ML$auROC[z] = roc$auc
   
   abundance_impact_2_ML$pvalue[z] = kruskal.test(real_value, ppp, paired = F)$p.value
   
   abundance_impact_2_ML$difference_in_median_values[z] = median(ppp[real_value == 'No']) - median(ppp[real_value == 'Yes'])
   
   if (z ==0) {
      data = cbind(roc$sensitivities,roc$specificities)
      data = as.data.frame(data)
      colnames(data) = c('Sensitivities','Specificities')
      data = data[order(data$Specificities, decreasing = F),]
      
      ggplot(data, aes(x = Sensitivities, y = Specificities)) +
         geom_point(size = 1)+geom_step()+scale_x_reverse()+
         theme(axis.title = element_text(size = 6), 
               axis.text = element_text(size = 6), 
               legend.text = element_text(size = 6), 
               legend.title = element_text(size = 6)) + theme_bw() 
      ggsave('Sensitivities_Specificities_final_models_dx_bv_input_urine_output_biofilm.pdf',width=3, height=3)
      
      data = data.frame(Real_value=real_value, Predicted_value= ppp)
      
      ggplot(data, aes(x=Real_value, y=Predicted_value,fill=Real_value)) +
         geom_boxplot() + theme_bw() 
      ggsave('model_final_dx_bv_input_urine_output_biofilm.pdf',width=2.5, height=3)
   }
   
}

write.csv(abundance_impact_2_ML,'model_final_dx_bv_input_urine_output_biofilm.csv')