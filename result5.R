#### results 5	Subgroup of L.iners, L.cripatus, L.gasseri, L.jensenii, G,vaginalis ######

set.seed(123)
library(reshape2)
source('CA_funciton_20220827.R')
ingroup_info_CST = read.delim('ingroup_info_CST.txt')
names(ingroup_info_CST)[c(31,33,34)] = c("Nugent score", "Community type" , "Self-reported symptoms")

#---------------------- Part 5.1  : Gene heatmap in different factors ---------------------
##----- Part 8.1 : Subgroup in five species ----
## Lactobacillus_iners: VOG = 5216, total gene = 1,200, 1236
## Lactobacillus_jensenii: VOG =5529， total gene =  1,686, 1566
## Lactobacillus_gasseri: VOG = 4563, total gene = 1,747, 1830
## Lactobacillus_crispatus: VOG = 6607, total gene = 2,056, 2064
## Gardnerella_vaginalis: VOG = 32234, total gene = 1,454, 1329
# L.ienrs
L.iners_gene = read.delim('other_results/gene_heatmap/L.iners_gene_presence.txt',header=F)
names(L.iners_gene) = c('sample_id','gene_id')
L.iners_gene$value =1
L.iners_count = aggregate(data = L.iners_gene, value ~ sample_id, sum)
L.iners_count$per = L.iners_count$value/length(unique(L.iners_gene$gene_id))
L.iners_gene_select1 = L.iners_gene[which(L.iners_gene$sample_id %in% L.iners_count[which(L.iners_count$value >= 1236*0.8),'sample_id']),]
L.iners_gene_select2 = merge(ingroup_info_CST['sample_id'],L.iners_gene_select1,by='sample_id')
L.iners_gene_select_d = dcast(L.iners_gene_select2, sample_id ~ gene_id,value.var = 'value')
L.iners_gene_select_d[is.na(L.iners_gene_select_d)] = 0
row.names(L.iners_gene_select_d) = L.iners_gene_select_d$sample_id

CST_group_Li = ingroup_info_CST[which(ingroup_info_CST$sample_id %in% unique(L.iners_gene_select2$sample_id)),c('sample_id','Community type')]
row.names(CST_group_Li) = CST_group_Li$sample_id
L.iners_subgroup = gene_heatmap_func(L.iners_gene_select_d,CST_group_Li, 'plot/09_gene_heatmap_L.iners.png',2)
L.iners_subgroup_meta = merge(L.iners_subgroup, ingroup_info_CST, by='sample_id')
write.table(L.iners_subgroup_meta, 'other_results/gene_heatmap/L.iners_subgroup_meta.txt',sep='\t',row.names = F, quote = F)

#L.crispatus
L.crispatus_gene = read.delim('other_results/gene_heatmap/L.crispatus_gene_presence.txt',header=F)
names(L.crispatus_gene) = c('sample_id','gene_id')
L.crispatus_gene$value =1
L.crispatus_count = aggregate(data = L.crispatus_gene, value ~ sample_id, sum)
L.crispatus_count$per = L.crispatus_count$value/length(unique(L.crispatus_gene$gene_id))
L.crispatus_gene_select1 = L.crispatus_gene[which(L.crispatus_gene$sample_id %in% L.crispatus_count[which(L.crispatus_count$value >= 2064*0.8),'sample_id']),]
L.crispatus_gene_select2 = merge(ingroup_info_CST['sample_id'],L.crispatus_gene_select1,by='sample_id')
L.crispatus_gene_select_d = dcast(L.crispatus_gene_select2, sample_id ~ gene_id,value.var = 'value')
L.crispatus_gene_select_d[is.na(L.crispatus_gene_select_d)] = 0
row.names(L.crispatus_gene_select_d) = L.crispatus_gene_select_d$sample_id

CST_group_Lc = ingroup_info_CST[which(ingroup_info_CST$sample_id %in% unique(L.crispatus_gene_select2$sample_id)),c('sample_id','Community type')]
row.names(CST_group_Lc) = CST_group_Lc$sample_id
L.crispatus_subgroup = gene_heatmap_func(L.crispatus_gene_select_d,CST_group_Lc, 'plot/09_gene_heatmap_L.crispatus.png',2)
L.crispatus_subgroup_meta = merge(L.crispatus_subgroup, ingroup_info_CST, by='sample_id')
write.table(L.crispatus_subgroup_meta, 'other_results/gene_heatmap/L.crispatus_subgroup_meta.txt',sep='\t',row.names = F, quote = F)


#L.jensenii
L.jensenii_gene = read.delim('other_results/gene_heatmap/L.jensenii_gene_presence.txt',header=F)
names(L.jensenii_gene) = c('sample_id','gene_id')
L.jensenii_gene$value =1
L.jensenii_count = aggregate(data = L.jensenii_gene, value ~ sample_id, sum)
L.jensenii_count$per = L.jensenii_count$value/length(unique(L.jensenii_gene$gene_id))
L.jensenii_gene_select1 = L.jensenii_gene[which(L.jensenii_gene$sample_id %in% L.jensenii_count[which(L.jensenii_count$value >= 1566*0.8),'sample_id']),]
L.jensenii_gene_select2 = merge(ingroup_info_CST['sample_id'],L.jensenii_gene_select1,by='sample_id')
L.jensenii_gene_select_d = dcast(L.jensenii_gene_select2, sample_id ~ gene_id,value.var = 'value')
L.jensenii_gene_select_d[is.na(L.jensenii_gene_select_d)] = 0
row.names(L.jensenii_gene_select_d) = L.jensenii_gene_select_d$sample_id

CST_group_Lj = ingroup_info_CST[which(ingroup_info_CST$sample_id %in% unique(L.jensenii_gene_select2$sample_id)),c('sample_id','Community type')]
row.names(CST_group_Lj) = CST_group_Lj$sample_id
L.jensenii_subgroup=gene_heatmap_func(L.jensenii_gene_select_d,CST_group_Lj, 'plot/09_gene_heatmap_L.jensenii.png',2)
L.jensenii_subgroup_meta = merge(L.jensenii_subgroup, ingroup_info_CST, by='sample_id')
write.table(L.jensenii_subgroup_meta, 'other_results/gene_heatmap/L.jensenii_subgroup_meta.txt',sep='\t',row.names = F, quote = F)

#L.gasseri
L.gasseri_gene = read.delim('other_results/gene_heatmap/L.gasseri_gene_presence.txt',header=F)
names(L.gasseri_gene) = c('sample_id','gene_id')
L.gasseri_gene$value =1
L.gasseri_count = aggregate(data = L.gasseri_gene, value ~ sample_id, sum)
L.gasseri_count$per = L.gasseri_count$value/length(unique(L.gasseri_gene$gene_id))
L.gasseri_gene_select1 = L.gasseri_gene[which(L.gasseri_gene$sample_id %in% L.gasseri_count[which(L.gasseri_count$value >= 1830*0.8),'sample_id']),]
L.gasseri_gene_select2 = merge(ingroup_info_CST['sample_id'],L.gasseri_gene_select1,by='sample_id')
L.gasseri_gene_select_d = dcast(L.gasseri_gene_select2, sample_id ~ gene_id,value.var = 'value')
L.gasseri_gene_select_d[is.na(L.gasseri_gene_select_d)] = 0
row.names(L.gasseri_gene_select_d) = L.gasseri_gene_select_d$sample_id

CST_group_Lg = ingroup_info_CST[which(ingroup_info_CST$sample_id %in% unique(L.gasseri_gene_select2$sample_id)),c('sample_id','Community type')]
row.names(CST_group_Lg) = CST_group_Lg$sample_id
L.gasseri_subgroup = gene_heatmap_func(L.gasseri_gene_select_d,CST_group_Lg, 'plot/09_gene_heatmap_L.gasseri.png',2)
L.gasseri_subgroup_meta = merge(L.gasseri_subgroup, ingroup_info_CST, by='sample_id')
write.table(L.gasseri_subgroup_meta, 'other_results/gene_heatmap/L.gasseri_subgroup_meta.txt',sep='\t',row.names = F, quote = F)



#G.vaginal
G.vaginal_gene = read.delim('other_results/gene_heatmap/G.vaginalis_gene_presence.txt',header=F)
names(G.vaginal_gene) = c('sample_id','gene_id')
G.vaginal_gene$value =1
G.vaginal_count = aggregate(data = G.vaginal_gene, value ~ sample_id, sum)
G.vaginal_count$per = G.vaginal_count$value/length(unique(G.vaginal_gene$gene_id))
G.vaginal_gene_select1 = G.vaginal_gene[which(G.vaginal_gene$sample_id %in% G.vaginal_count[which(G.vaginal_count$value >= 1329*0.8),'sample_id']),]
G.vaginal_gene_select2 = merge(ingroup_info_CST['sample_id'],G.vaginal_gene_select1,by='sample_id')
G.vaginal_gene_select_d = dcast(G.vaginal_gene_select2, sample_id ~ gene_id,value.var = 'value')
G.vaginal_gene_select_d[is.na(G.vaginal_gene_select_d)] = 0
row.names(G.vaginal_gene_select_d) = G.vaginal_gene_select_d$sample_id
write.table(G.vaginal_gene_select_d, 'other_results/gene_heatmap/G.vaginalis_gene_select_d_filter.txt',sep='\t',row.names=F,quote = F)

CST_group_Gv = ingroup_info_CST[which(ingroup_info_CST$sample_id %in% unique(G.vaginal_gene_select2$sample_id)),c('sample_id','Community type')]
row.names(CST_group_Gv) = CST_group_Gv$sample_id
G.vaginal_subgroup = gene_heatmap_func(G.vaginal_gene_select_d,CST_group_Gv, 'plot/09_gene_heatmap_G.vaginal.png',4)
#G.vaginal_subgroup = read.delim('other_results/gene_heatmap/G.vaginal_subgroup.txt')
G.vaginal_subgroup_meta = merge(G.vaginal_subgroup, ingroup_info_CST, by='sample_id')
write.table(G.vaginal_subgroup_meta, 'other_results/gene_heatmap/G.vaginal_subgroup_meta.txt',sep='\t',row.names = F, quote = F)

##----- Part 5.2 : association with metadata and subgroup ----
# 
#ingroup_info_CST = read.delim('ingroup_info_CST.txt')
factors= c('race','suffers_insomnia','Self.reported.symptoms','city','menstrual_cycle','frequency_of_sex',
                    'has_children','abnormal_pregnancy',"meat","vegan", "heavy_salt_oil","spicy_deit","sweet_deit","intake_of_yogurt",
                    'education_level','income_level','Community.type')


L.iners_subgroup_meta = read.delim('other_results/gene_heatmap/L.iners_subgroup_meta.txt')
L.iners_asso = subgroup_chisq(L.iners_subgroup_meta,factors)
summary(aov(age~group, L.iners_subgroup_meta))
summary(aov(BMI~group, L.iners_subgroup_meta))

L.crispatus_subgroup_meta = read.delim('other_results/gene_heatmap/L.crispatus_subgroup_meta.txt')
dim(L.crispatus_subgroup_meta)
table(L.crispatus_subgroup_meta$group)
L.crispatu_asso = subgroup_chisq(L.crispatus_subgroup_meta,factors)
summary(aov(age~group, L.crispatus_subgroup_meta))
summary(aov(BMI~group, L.crispatus_subgroup_meta))


L.jensenii_subgroup_meta = read.delim('other_results/gene_heatmap/L.jensenii_subgroup_meta.txt')
L.jensenii_asso = subgroup_chisq(L.jensenii_subgroup_meta,factors)
summary(aov(age~group, L.jensenii_subgroup_meta))
summary(aov(BMI~group, L.jensenii_subgroup_meta))
table(L.jensenii_subgroup_meta$Community.type,L.jensenii_subgroup_meta$group)[-c(2,4),]
chisq.test(table(L.jensenii_subgroup_meta$Community.type,L.jensenii_subgroup_meta$group)[-c(2,4),])


L.gasseri_subgroup_meta = read.delim('other_results/gene_heatmap/L.gasseri_subgroup_meta.txt')
L.gasseri_asso = subgroup_chisq(L.gasseri_subgroup_meta,factors)
summary(aov(age~group, L.gasseri_subgroup_meta))
summary(aov(BMI~group, L.gasseri_subgroup_meta))


G.vaginal_subgroup_meta = read.delim('other_results/gene_heatmap/G.vaginal_subgroup_meta.txt')
G.vaginal_asso = subgroup_chisq(G.vaginal_subgroup_meta,factors)
summary(aov(age~group, G.vaginal_subgroup_meta))
summary(aov(BMI~group, G.vaginal_subgroup_meta))



write.table(L.iners_asso, 'other_results/gene_heatmap/L.iners_asso.txt',sep='\t',row.names=F, quote = F)
write.table(L.crispatu_asso, 'other_results/gene_heatmap/L.crispatu_asso.txt',sep='\t',row.names=F, quote = F)
write.table(L.jensenii_asso, 'other_results/gene_heatmap/L.jensenii_asso.txt',sep='\t',row.names=F, quote = F)
write.table(G.vaginal_asso, 'other_results/gene_heatmap/G.vaginal_asso.txt',sep='\t',row.names=F, quote = F)


#---------------------- Part 5.2  : assoiciation between Gene heatmap group and different factors ---------------------
library(nnet)
ingroup_info_CST_multinom = read.delim('ingroup_info_CST_multinom_imputation.txt')
f = names(ingroup_info_CST_multinom)[-c(which(names(ingroup_info_CST_multinom) %in% c('sample_id','age','BMI')))]
ref_f = c('0','0','Beijing','follicular phase','once_a_week','1','0',
          '0','0','0','0','0','0','bachelor','8-15')
for (i in 1:length(f)){
  print(paste(c("============> ", f[i], " : ", ref_f[i], " <==========="),collapse = ''))
  ingroup_info_CST_multinom[,f[i]] = as.factor(ingroup_info_CST_multinom[,f[i]])
  ingroup_info_CST_multinom[,f[i]] = relevel(ingroup_info_CST_multinom[,f[i]],ref=ref_f[i])
}
## logistic regression

L.iners_OR_total = taxa_subgroup_logistic_func(L.iners_subgroup_meta, ingroup_info_CST_multinom, 'group1')

L.crispatus_OR_total = taxa_subgroup_logistic_func(L.crispatus_subgroup_meta, ingroup_info_CST_multinom, 'group1')

L.jensenii_OR_total = taxa_subgroup_logistic_func(L.jensenii_subgroup_meta, ingroup_info_CST_multinom, 'group1')

# Pvalue = 1
# L.gasseri_OR_total = taxa_subgroup_logistic_func(L.gasseri_subgroup_meta, ingroup_info_CST_multinom)

#multinom logistic regression in the subgroup in gene heatmap
G.vaginal_OR_total = taxa_subgroup_multinom_func(G.vaginal_subgroup_meta, ingroup_info_CST_multinom, 'group2')

write.table(L.iners_OR_total, 'other_results/gene_heatmap/L.iners_OR_total.txt',sep='\t',row.names=F, quote = F)
write.table(L.crispatus_OR_total, 'other_results/gene_heatmap/L.crispatus_OR_total.txt',sep='\t',row.names=F, quote = F)
#write.table(L.gasseri_OR_total, 'other_results/gene_heatmap/L.gasseri_OR_total.txt',sep='\t',row.names=F, quote = F)
write.table(L.jensenii_OR_total, 'other_results/gene_heatmap/L.jensenii_OR_total.txt',sep='\t',row.names=F, quote = F)
write.table(G.vaginal_OR_total, 'other_results/gene_heatmap/G.vaginal_OR_total.txt',sep='\t',row.names=F, quote = F)

### Forest plot ----------------
Lc_OR_melt_merge = read.delim('other_results/gene_heatmap/L.crispatus_OR_total.txt',na.strings = "")

library(forestplot)
forestplot(labeltext = as.matrix(Lc_OR_melt_merge[,c(3,9)]), 
           mean = Lc_OR_melt_merge$OR,
           lower = Lc_OR_melt_merge$lower2.5. ,
           upper = Lc_OR_melt_merge$upper97.5., 
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



Li_OR_melt_merge = read.delim('other_results/gene_heatmap/L.iners_OR_total.txt',na.strings = "")

forestplot(labeltext = as.matrix(Li_OR_melt_merge[,c(3,9)]), 
           mean = Li_OR_melt_merge$OR,
           lower = Li_OR_melt_merge$lower2.5. ,
           upper = Li_OR_melt_merge$upper97.5., 
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

Lj_OR_melt_merge = read.delim('other_results/gene_heatmap/L.jensenii_OR_total.txt',na.strings = "")

forestplot(labeltext = as.matrix(Lj_OR_melt_merge[,c(3,9)]), 
           mean = Lj_OR_melt_merge$OR,
           lower = Lj_OR_melt_merge$lower2.5. ,
           upper = Lj_OR_melt_merge$upper97.5., 
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


Gv1_OR_melt_merge = read.delim('other_results/gene_heatmap/G.vaginal_OR_total_group1.txt',na.strings = "")

forestplot(labeltext = as.matrix(Gv1_OR_melt_merge[,c(3,9)]), 
           mean = Gv1_OR_melt_merge$OR,
           lower = Gv1_OR_melt_merge$lower2.5. ,
           upper = Gv1_OR_melt_merge$upper97.5., 
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

