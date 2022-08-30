#### results 4	Factors associated with CSTs and microbiota ######
###############   result 4.1 multinom logsistic regresion for CSTs  ------
library(reshape2)

### ---- Part 9: Multinomial Logistic Regression Analysis -------
# reshape feature
ingroup_info_CST = read.delim('ingroup_info_CST.txt')
ingroup_morph = read.delim('ingroup_morph.txt')


feature_df1 = ingroup_info_CST[,-(which(names(ingroup_info_CST) %in% c('batch_id','age_cut','BMI_cut','disease','sex48h','contraception','symptom_q','self_report_sym','Nugent.score','Lactobacillus','Community.type','Self.reported.symptoms',''))),]
for(i in 1:(ncol(feature_df1)-1)){
  print(names(feature_df1)[i+1])
  feature_df1[,names(feature_df1)[i+1]] = as.character(feature_df1[,names(feature_df1)[i+1]])
  feature_df1[which(feature_df1[,names(feature_df1)[i+1]] == 'missing'),names(feature_df1)[i+1]] = NA
  feature_df1[which(feature_df1[,names(feature_df1)[i+1]] == 'yes'),names(feature_df1)[i+1]] = 1
  feature_df1[which(feature_df1[,names(feature_df1)[i+1]] == 'no'),names(feature_df1)[i+1]] = 0
}

feature_df1$age = as.numeric(feature_df1$age)
feature_df1$BMI = as.numeric(feature_df1$BMI)

feature_df2 = feature_df1
feature_df2[which(feature_df2$symptoms_group == 'No_symptoms'),'symptoms_group']=0
feature_df2[which(feature_df2$symptoms_group == 'Any_symptoms'),'symptoms_group']=1
feature_df2[which(feature_df2$menstrual_status == 'regularity'),'menstrual_status']=0
feature_df2[which(feature_df2$menstrual_status == 'irregularity'),'menstrual_status']=1
feature_df2$race = as.factor(feature_df2$race)
#"Han"=1   "other" =4  "Tibetan"=2 "Uighur" =3
levels(feature_df2$race) = c(1,4,2,3)
feature_df2$city = as.factor(feature_df2$city)
#"Beijing"      "Chongqing"    "Guangdong"    "Guangxi"      "Heilongjiang" "Hunan"        "Shanghai"     "Shanxi"       "Tibet"      "Xinjiang"    
levels(feature_df2$city)=c(1,2,3,4,5,6,7,8,9,10)
#follicular phase=1     luteal phase =3  menstrual period =4 ovulation phase = 2
feature_df2$menstrual_cycle = as.factor(feature_df2$menstrual_cycle)
levels(feature_df2$menstrual_cycle) = c(1,3,4,2)
# "once_a_month"  =3     "once_a_week"  =2        "once_few_month"=4       "other"=5  "several_times_a_week" =1
feature_df2$frequency_of_sex = as.factor(feature_df2$frequency_of_sex)
levels(feature_df2$frequency_of_sex) = c(3,2,4,5,1)
#"bachelor" =3  ,   "bachelor_below" =2  "docter"  =5     "hs_below" =1    "master" =4
feature_df2$education_level = as.factor(feature_df2$education_level)
levels(feature_df2$education_level) = c(3,2,5,1,4)
#  ">50"=5   "15-50"=4   "5"=1     "5-8"=2   "8-15"=3
feature_df2$income_level = as.factor(feature_df2$income_level)
levels(feature_df2$income_level) = c(5,4,1,2,3)
for(i in 1:(ncol(feature_df2)-1)){
  print(names(feature_df2)[i+1])
  feature_df2[,names(feature_df2)[i+1]] = as.character(feature_df2[,names(feature_df2)[i+1]])
  feature_df2[,names(feature_df2)[i+1]] = as.numeric(feature_df2[,names(feature_df2)[i+1]])
}

## imputation
# missing value : Multiple Imputation using mice

library(mice)
row.names(feature_df2) =feature_df2$sample_id

imp<-mice(feature_df2[,c(2:22)],seed=1234)
feature_df3 = complete(imp)

feature_df4=feature_df3
feature4 = names(feature_df4)[-c(2,8)]

for(i in 1:(length(feature4))){
  print(feature4[i])
  feature_df4[,feature4[i]] = as.factor(feature_df4[,feature4[i]])
}


levels(feature_df4$race)  = c('Han','Tibetan','Uighur','Other')
levels(feature_df4$city) = c('Beijing','Chongqing','Guangdong','Guangxi','Heilongjiang','Hunan','Shanghai','Shanxi','Tibet','Xinjiang')
levels(feature_df4$menstrual_cycle)  = c('follicular phase','ovulation phase','luteal phase','menstrual period')
levels(feature_df4$frequency_of_sex)  = c('several_times_a_week','once_a_week','once_a_month','once_few_month','other')
levels(feature_df4$education_level) =c('hs_below','bachelor_below','bachelor','master','docter')
levels(feature_df4$income_level) =c('5','5-8','8-15','15-50','>50')


# run model
factors_inmodel = c('age','BMI','suffers_insomnia','symptoms_group','city','menstrual_cycle','frequency_of_sex',
                    'has_children','abnormal_pregnancy',"meat","vegan", "heavy_salt_oil","spicy_deit","sweet_deit","intake_of_yogurt",
                    'education_level','income_level')
ingroup_info_CST_multinom = feature_df4[,factors_inmodel]

f = names(ingroup_info_CST_multinom)[-c(1,2)]
ref_f = c('0','0','Beijing','follicular phase','once_a_week','1','0',
          '0','0','0','0','0','0','bachelor','8-15')
for (i in 1:length(f)){
  print(paste(c("============> ", f[i], " : ", ref_f[i], " <==========="),collapse = ''))
  ingroup_info_CST_multinom[,f[i]] = as.factor(ingroup_info_CST_multinom[,f[i]])
  ingroup_info_CST_multinom[,f[i]] = relevel(ingroup_info_CST_multinom[,f[i]],ref=ref_f[i])
}

ingroup_info_CST_multinom$sample_id = row.names(ingroup_info_CST_multinom)
write.table(ingroup_info_CST_multinom, 'ingroup_info_CST_multinom_imputation.txt',sep='\t',row.names=F, quote = F)

library(nnet)

## ---- 9.1 CSTs vs metadata ------
CSTs_df = merge(ingroup_info_CST_multinom,ingroup_info_CST[,c('sample_id','Community.type')], by='sample_id')
CSTs_df$Community.type=as.factor(CSTs_df$Community.type)
CSTs_df$Community.type =  relevel(CSTs_df$Community.type, ref = "CST I")
CST_multinom =  multinom(Community.type ~ . , data = CSTs_df[-1])
summary(CST_multinom)
CST_OR=round(exp(coef(CST_multinom)),2) #计算OR 
CST_OR95 = as.data.frame(round(exp(confint(CST_multinom)),2)) #计算OR的95%CI
CST_OR95$factors = row.names(CST_OR95)
# 计算p value
CST_z <- summary(CST_multinom)$coefficients/summary(CST_multinom)$standard.errors
CST_p <-round(((1 - pnorm(abs(CST_z), 0, 1)) * 2),5)

# 森林图
library(forestplot)
CST_OR_melt = melt(CST_OR)
names(CST_OR_melt) = c('CST','factors','OR')
CST_p_melt = melt(CST_p)
names(CST_p_melt) = c('CST','factors','p value')
### CST  III
CST3_OR_df = CST_OR_melt[which(CST_OR_melt$CST == 'CST III'),]
CST3_OR_melt_merge = merge(CST3_OR_df, CST_OR95[,c('factors','2.5 %.CST III','97.5 %.CST III')], by='factors')
CST3_OR_melt_merge = merge(CST3_OR_melt_merge, CST_p_melt, by=c('CST','factors'))
for(i in 1:nrow(CST3_OR_melt_merge)){
  CST3_OR_melt_merge[i,'OR_95OR'] = paste( c(round(CST3_OR_melt_merge[i,'OR'],2),' (', round(CST3_OR_melt_merge[i,'2.5 %.CST III'],2), '-' , round(CST3_OR_melt_merge[i,'97.5 %.CST III'],2),')'), collapse = '')
}
#write.table(CST3_OR_melt_merge, 'other_results/09_CST3_OR_melt_merge.txt',sep='\t',row.names = F, quote = F)
CST3_OR_melt_merge = read.delim('other_results/09_CST3_OR_melt_merge.txt',na.strings = "")

forestplot(labeltext = as.matrix(CST3_OR_melt_merge[,c(3,9)]), 
           mean = CST3_OR_melt_merge$OR,
           lower = CST3_OR_melt_merge$X2.5...CST.III, 
           upper = CST3_OR_melt_merge$X97.5...CST.III, 
           #graph.pos=3,
           zero = 1, 
           clip = c(0.1,20), 
           graph.pos = 2, 
           xticks = c(0.1,1,10),
           boxsize = 0.3, 
           xlog = TRUE,
           graphwidth = unit(0.3,"npc"),
           col=fpColors(box="royalblue",line="darkblue"),
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           #箱线图中基准线的位置
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.2)


### CST II
CST2_OR_df = CST_OR_melt[which(CST_OR_melt$CST == 'CST II'),]
CST2_OR_melt_merge = merge(CST2_OR_df, CST_OR95[,c('factors','2.5 %.CST II','97.5 %.CST II')], by='factors')
CST2_OR_melt_merge = merge(CST2_OR_melt_merge, CST_p_melt, by=c('CST','factors'))
for(i in 1:nrow(CST2_OR_melt_merge)){
  CST2_OR_melt_merge[i,'OR_95OR'] = paste( c(round(CST2_OR_melt_merge[i,'OR'],2),' (', round(CST2_OR_melt_merge[i,'2.5 %.CST II'],2), '-' , round(CST2_OR_melt_merge[i,'97.5 %.CST II'],2),')'), collapse = '')
}
#write.table(CST2_OR_melt_merge, 'other_results/09_CST2_OR_melt_merge.txt',sep='\t',row.names = F, quote = F)
CST2_OR_melt_merge = read.delim('other_results/09_CST2_OR_melt_merge.txt',na.strings = "")

forestplot(labeltext = as.matrix(CST2_OR_melt_merge[,c(3,9)]), 
           mean = CST2_OR_melt_merge$OR,
           lower = CST2_OR_melt_merge$X2.5...CST.II, 
           upper = CST2_OR_melt_merge$X97.5...CST.II, 
           #graph.pos=3,
           zero = 1, 
           clip = c(0.1,20), 
           graph.pos = 2, 
           xticks = c(0.1,1,10),
           boxsize = 0.3, 
           xlog = TRUE,
           graphwidth = unit(0.3,"npc"),
           col=fpColors(box="royalblue",line="darkblue"),
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           #箱线图中基准线的位置
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.2)

### CST IV
CST4_OR_df = CST_OR_melt[which(CST_OR_melt$CST == 'CST IV'),]
CST4_OR_melt_merge = merge(CST4_OR_df, CST_OR95[,c('factors','2.5 %.CST IV','97.5 %.CST IV')], by='factors')
CST4_OR_melt_merge = merge(CST4_OR_melt_merge, CST_p_melt, by=c('CST','factors'))
for(i in 1:nrow(CST4_OR_melt_merge)){
  CST4_OR_melt_merge[i,'OR_95OR'] = paste( c(round(CST4_OR_melt_merge[i,'OR'],2),' (', round(CST4_OR_melt_merge[i,'2.5 %.CST IV'],2), '-' , round(CST4_OR_melt_merge[i,'97.5 %.CST IV'],2),')'), collapse = '')
}
#write.table(CST4_OR_melt_merge, 'other_results/09_CST4_OR_melt_merge.txt',sep='\t',row.names = F, quote = F)
CST4_OR_melt_merge = read.delim('other_results/09_CST4_OR_melt_merge.txt',na.strings = "")


forestplot(labeltext = as.matrix(CST4_OR_melt_merge[,c(3,9)]), 
           mean = CST4_OR_melt_merge$OR,
           lower = CST4_OR_melt_merge$X2.5...CST.IV, 
           upper = CST4_OR_melt_merge$X97.5...CST.IV, 
           #graph.pos=3,
           zero = 1, 
           clip = c(0.1,20), 
           graph.pos = 2, 
           xticks = c(0.1,1,10),
           boxsize = 0.3, 
           xlog = TRUE,
           graphwidth = unit(0.3,"npc"),
           col=fpColors(box="royalblue",line="darkblue"),
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           #箱线图中基准线的位置
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.2)



### CST V
CST5_OR_df = CST_OR_melt[which(CST_OR_melt$CST == 'CST V'),]
CST5_OR_melt_merge = merge(CST5_OR_df, CST_OR95[,c('factors','2.5 %.CST V','97.5 %.CST V')], by='factors')
CST5_OR_melt_merge = merge(CST5_OR_melt_merge, CST_p_melt, by=c('CST','factors'))
for(i in 1:nrow(CST5_OR_melt_merge)){
  CST5_OR_melt_merge[i,'OR_95OR'] = paste( c(round(CST5_OR_melt_merge[i,'OR'],2),' (', round(CST5_OR_melt_merge[i,'2.5 %.CST V'],2), '-' , round(CST5_OR_melt_merge[i,'97.5 %.CST V'],2),')'), collapse = '')
}
#write.table(CST5_OR_melt_merge, 'other_results/09_CST5_OR_melt_merge.txt',sep='\t',row.names = F, quote = F)
CST5_OR_melt_merge = read.delim('other_results/09_CST5_OR_melt_merge.txt',na.strings = "")

forestplot(labeltext = as.matrix(CST5_OR_melt_merge[,c(3,9)]), 
           mean = CST5_OR_melt_merge$OR,
           lower = CST5_OR_melt_merge$X2.5...CST.V, 
           upper = CST5_OR_melt_merge$X97.5...CST.V, 
           #graph.pos=3,
           zero = 1, 
           clip = c(0.1,20), 
           graph.pos = 2, 
           xticks = c(0.1,1,10),
           boxsize = 0.3, 
           xlog = TRUE,
           graphwidth = unit(0.3,"npc"),
           col=fpColors(box="royalblue",line="darkblue"),
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           #箱线图中基准线的位置
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱线图两端添加小竖线，高度
           ci.vertices=TRUE, ci.vertices.height = 0.2)


############### result 4.2 RDA analysis  ------
##--------------------------- CCA and RDA ------------------------------------
# process env data
CCA_env1 = feature_df4[,c(1,2,8:14,16:18)]
CCA_env1$sample_id = row.names(CCA_env1)
CCA_env2 = feature_df4[,c(4,15,19:21)]
CCA_env2$sample_id = row.names(CCA_env2)
for(i in 1:5){
  CCA_env2_sub=CCA_env2[,c('sample_id',names(CCA_env2)[i])]
  names(CCA_env2_sub)[2] = 'variable'
  CCA_env2_sub$value = 1
  CCA_env2_sub_d = dcast(data =CCA_env2_sub, sample_id ~ variable, value.var = 'value')
  CCA_env2_sub_d[is.na(CCA_env2_sub_d)] = 0
  CCA_env1 = merge(CCA_env1, CCA_env2_sub_d, by='sample_id')
}

# process microbiota data
filter_virgo_abdper_melt_select = read.delim('filter_virgo_abdper_melt_select.txt')
filter0.001_abd_per_mean =  aggregate(data = filter_virgo_abdper_melt_select,abd_per ~ species, mean)
filter_cutoff1 = filter0.001_abd_per_mean[which(filter0.001_abd_per_mean$abd_per >=0.001),'species']
filter0.001_abd_per = filter_virgo_abdper_melt_select[which(filter_virgo_abdper_melt_select$species %in% filter_cutoff1 & filter_virgo_abdper_melt_select$abd_per !=0),]

# group
CST_group = ingroup_info_CST[,c('sample_id','Community.type')]
names(CST_group)[2] = 'group'
row.names(CST_group) = CST_group$sample_id
group  = CST_group
#env = ingroup_info_CST[,c(1,2,4,6:11,13:28)]
env = CCA_env1
write.table(env, 'ingroup_info_CST_multinom_imputation_01_binary_variable.txt',sep='\t',row.names=F, quote = F)


for(i in 1:(ncol(env)-1)){
  print(names(env)[i+1])
  env[,names(env)[i+1]] = as.character(env[,names(env)[i+1]])
  env[,names(env)[i+1]] = as.numeric(env[,names(env)[i+1]])
}

row.names(env) = env$sample_id
library(ggrepel)
sabd_sub = merge(filter0.001_abd_per, env, by='sample_id')
sabd_sub_d = dcast(sabd_sub, sample_id ~ species, value.var = 'abd_per')
sabd_sub_d[is.na(sabd_sub_d)]=0
row.names(sabd_sub_d) = sabd_sub_d$sample_id
# group
#group2 <- as.list(group[,c('sample_id','group')])
###hellinger transform
sabd_sub_d.hell <- decostand(sabd_sub_d[-1], "hellinger")

decorana(veg = sabd_sub_d.hell)

####------ RDA  -------
rda_tb<-rda(sabd_sub_d.hell~.,env[,-which(names(env) == 'sample_id')])
result<-summary(rda_tb)
rda_coef <- coef(rda_tb)

rda_tb.scaling1 <- summary(rda_tb, scaling = 1)
rda_tb.site <- data.frame(rda_tb.scaling1$sites)[1:2]
rda_tb.env <- data.frame(rda_tb.scaling1$biplot)[1:2]
#ggplot2 
rda_tb.site$sample <- rownames(rda_tb.site)
rda_tb.site <- merge(rda_tb.site, group, by = 'row.names')
rda_tb.env$sample <- rownames(rda_tb.env)

show_factors = c('symptoms_group','BMI','abnormal_pregnancy','Beijing','Chongqing','Guangdong','Guangxi','Heilongjiang','Hunan','Shanghai',
                 'Shanxi','Tibet','`ovulation phase`','`luteal phase`','hs_below','`5`')
rda_tb.env_show = rda_tb.env[which(rda_tb.env$sample %in% show_factors),]



library(ggplot2)
color=c( "#3C5488B2","#00A087B2", 
         "#F39B7FB2","#91D1C2B2", 
         "#8491B4B2", "#DC0000B2", 
         "#7E6148B2","yellow", 
         "darkolivegreen1", "lightskyblue", 
         "darkgreen", "deeppink", "khaki2", 
         "firebrick", "brown1", "darkorange1", 
         "cyan1", "royalblue4", "darksalmon", 
         "darkgoldenrod1", "darkseagreen", "darkorchid")

p_rda <- ggplot(rda_tb.site, aes(RDA1, RDA2)) +
  geom_point(aes(color = group,shape = group)) +
  stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE, linetype = 2) +
  scale_color_manual(values = color[1:length(unique(group$group))]) +
  #theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5), legend.key = element_rect(fill = 'transparent')) + 
  
  theme(panel.grid.major = element_line(color = 'gray', size = 0.1), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'), plot.title = element_text(hjust = 0.5)) + 
  #labs(x = 'RDA1 (42.91%)', y = 'RDA2 (9.80%)') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb.env_show, aes(x = 0, y = 0, xend = RDA1,yend = RDA2), arrow = arrow(length = unit(0.1, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = rda_tb.env_show, aes(RDA1 * 1.1, RDA2 * 1.1, label = sample), color = 'blue', size = 3)
#geom_label_repel(aes(label =sample, color = group), size = 3, box.padding = unit(0, 'lines'), show.legend = FALSE)

p_rda
ggsave(p_rda, filename ='plot/RDA_plot.png',width = 10, height = 7)

rda.perm=permutest(rda_tb,permu=999)
rda.perm
## 检查每个环境因子与群落变化的相关性
rda.env=envfit(rda_tb,env,permu=999)
rda_pvalue =  rda.env$vectors            

