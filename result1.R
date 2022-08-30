#### result1: summary if study Chinese cohort #####
library(vioplot)
setwd('/Users/liangbingbing/Documents/turing/3.bioinformatics/1.project/1.wanren_xiehe/10.文章/10家医院数据文章/文章1-中美比对/rscript')
set.seed(123)

qbm_mNGS_rename_clean_batch_filter_factor = read.delim('qbm_mNGS_English_filter_questionaire_917s.txt')

# age, BMI, Symptoms - 917s
pdf(file ='plot/01_age_vioplot_917.pdf' , width = 4,height = 6)
vioplot(qbm_mNGS_rename_clean_batch_filter_factor$age, main = "Age",col.axis = NA) 
dev.off() 
summary(qbm_mNGS_rename_clean_batch_filter_factor$age)
sd(qbm_mNGS_rename_clean_batch_filter_factor$age,na.rm=T)

pdf(file ='plot/01_BMI_vioplot_917.pdf' , width = 4,height = 6)
vioplot(qbm_mNGS_rename_clean_batch_filter_factor$BMI, main = "BMI",col.axis = NA) 
dev.off() 
summary(qbm_mNGS_rename_clean_batch_filter_factor$BMI, na.rm=T)
sd(qbm_mNGS_rename_clean_batch_filter_factor$BMI, na.rm = T)


# Other question
summary_factors = c('age_cut','race','city','longterm_smoking','longterm_drinking','suffers_insomnia','BMI_cut',
                 'meat','vegan','heavy_salt_oil','spicy_deit','sweet_deit','intake_of_yogurt','disease', 
                 'menstrual_cycle','menstrual_status','has_children','abnormal_pregnancy','frequency_of_sex','sex48h',
                 'contraception','education_level','income_level')

count_df_917 = data.frame()
count_df_917_2 = data.frame()
for(i in 1:length(summary_factors)){
  factor = summary_factors[i]
  table_df =as.data.frame(table(qbm_mNGS_rename_clean_batch_filter_factor[factor]))
  table_df$count = paste(table_df$Var1, table_df$Freq,sep=': ')
  table_df$per = round(table_df$Freq/length(unique(qbm_mNGS_rename_clean_batch_filter_factor$sample_id))*100,2)
  table_df$factor = factor
  count_df_sub = data.frame(effect = factor, n = paste(table_df$count, collapse = ', '))
  count_df_917 = rbind(count_df_917, count_df_sub)
  count_df_917_2 = rbind(count_df_917_2, table_df)
}

write.table(count_df_917_2, 'other_results/factor_count_917s_percetage.txt',sep='\t',row.names = F, quote = F)
