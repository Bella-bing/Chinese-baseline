#### results 3: Microbiome composition and CST classification in Chinese women  ###############
library(reshape2)
library(phyloseq)
library(stringr)
source('CA_funciton_20220827.R')

set.seed(123)

#----------------------------------------------- 3.1 taxa filter : read abd count df and filter taxa ----------------------------------
ingroup_info = read.delim('ingroup_info_CST.txt')
ingroup_morph = read.delim('ingroup_morph.txt')
ingroup_samples_list = ingroup_info$sample_id
# virgo
summ.Abundance = read.delim('virgo_results/summary.Abundance.txt')
names(summ.Abundance)[1] = 'sample_id'
summ.Abundance_filter1 = filter_taxa_read_func(summ.Abundance)
summ.Abundance_filter2 = summ.Abundance_filter1[,!names(summ.Abundance_filter1) %in% c('Burkholderia_mallei','Coprobacillus_sp.',
                                                                                       'Thermobifida_cellulosilytica','Paracoccus_denitrificans',
                                                                                       'Neisseria_gonorrhoeae','Actinobacillus_pleuropneumoniae')]

# for NG
ng = read.delim('virgo_results/Neisseria_gonorrhoeae_reads_count.txt')
summ.Abundance_filter3 = merge(ng, summ.Abundance_filter2, by='sample_id',all=T)
summ.Abundance_filter3[is.na(summ.Abundance_filter3)]=0
virgo_abd_all = summ.Abundance_filter3
row.names(virgo_abd_all)= virgo_abd_all$sample_id

#------------------------------------ input data to phyloseq tool 
# otu_table
virgo_abd_all_matrix = t(virgo_abd_all[,-1])
otu = otu_table(virgo_abd_all_matrix,taxa_are_rows = TRUE)
# tax
taxdf= data.frame(Species = names(virgo_abd_all)[-1])
taxdf$Genus = NA
for (i in 1:nrow(taxdf)){
  taxdf[i,2] = str_split(taxdf[i,1],'_')[[1]][1]
}
row.names(taxdf) = taxdf$Species
tax = tax_table(as.matrix(taxdf))
physeq_virgo = phyloseq(otu, tax)

physeq_virgo_abdper = transform_sample_counts(physeq_virgo, function(x) x/sum(x))
filter_virgo_abdper = as.data.frame(otu_table(physeq_virgo_abdper))
# ingroup mcirobiome 

filter_virgo_abdper$species = row.names(filter_virgo_abdper)
filter_virgo_abdper_melt = melt(filter_virgo_abdper, id= 'species')
names(filter_virgo_abdper_melt) = c('species','sample_id','abd_per')
filter_virgo_abdper_melt_select = filter_virgo_abdper_melt[which(filter_virgo_abdper_melt$sample_id %in% ingroup_samples_list),]
col_check0  = aggregate(data = filter_virgo_abdper_melt_select, abd_per ~ species,  sum)
del_0species = col_check0[which(col_check0==0),'species']
filter_virgo_abdper_melt_select_d = dcast(filter_virgo_abdper_melt_select, sample_id ~ species, value.var = 'abd_per')
row.names(filter_virgo_abdper_melt_select_d) = filter_virgo_abdper_melt_select_d$sample_id

print(paste('Number of species pass filter : ', length(names(filter_virgo_abdper_melt_select_d))-1, sep=' '))
write.table(filter_virgo_abdper_melt_select,'filter_virgo_abdper_melt_select.txt',sep='\t',row.names=F, quote = F)
# genus level
species_genus_mapping = data.frame(species = unique(filter_virgo_abdper_melt_select$species),genus =NA)
for(i in 1:nrow(species_genus_mapping)){
  species_genus_mapping[i,'genus'] = str_split(species_genus_mapping[i,'species'], '_')[[1]][1]
}
filter_virgo_abdper_melt_select_genus_tmp = merge(filter_virgo_abdper_melt_select,species_genus_mapping, by='species')
filter_virgo_abdper_melt_select_genus = aggregate(abd_per ~ sample_id +genus, data=filter_virgo_abdper_melt_select_genus_tmp, sum)


# microbiome composition bar plot: abd_per >0.3 in 15% samples, three plot, normal WBC dysbiosis
pbarplot2_group = ingroup_morph[,c('sample_id','microecologial_group')]
pbarplot2_group[which(pbarplot2_group$microecologial_group %in% c('WBC≥10','Normal microecology')),'group'] = 'Normal microecology'
pbarplot2_group[which(pbarplot2_group$microecologial_group == 'Dysbiosis'),'group'] = 'Dysbiosis'
pbarplot2_group$group = factor(pbarplot2_group$group, levels=c('Normal microecology','Dysbiosis'))
pbarplot2 = abd_filter_plot_func_CHN(filter_virgo_abdper_melt_select, pbarplot2_group) 
ggsave(pbarplot2, filename = 'plot/02_p2_all_morph_2barplot_917s_test.png',width = 12, height = 5)


#--- Part 3.2 shannon diversity ----
Shannon_index_917s=as.data.frame(vegan::diversity(filter_virgo_abdper_melt_select_d[-1], index = "shannon"))
names(Shannon_index_917s) = 'shannon_index'
Shannon_index_917s$sample_id = row.names(Shannon_index_917s)
shannon_nugent = merge(Shannon_index_917s, ingroup_morph, by='sample_id')
shannon_nugent$Nugent评分 = as.factor(shannon_nugent$Nugent评分)
shannon_nugent_plot =  ggplot(data=shannon_nugent, aes(x = Nugent评分, y=shannon_index, group=Nugent评分)) + 
  geom_boxplot()+geom_jitter(width = 0.1,alpha = 0.2,size=1.5)+
  labs(x='Nugent score', y='Shannon diveristy index')
ggsave(shannon_nugent_plot, filename = 'plot/03_nugent_shannon_boxplot_917s.png',width = 6, height = 4)

# shannon index
ingroup_info_CST_shannon = merge(ingroup_info, Shannon_index_917s, by='sample_id')

shanno_mean = aggregate(data=ingroup_info_CST_shannon, shannon_index ~ Community.type, FUN = mean)
shanno_sd = aggregate(data=ingroup_info_CST_shannon, shannon_index ~ Community.type, FUN = sd)
shanno_mean_sd = merge(shanno_mean, shanno_sd, by='Community.type')
shanno_mean_sd$shannon_index.x = round(shanno_mean_sd$shannon_index.x,2)
shanno_mean_sd$shannon_index.y = round(shanno_mean_sd$shannon_index.y,2)
shanno_mean_sd$mean_sd = paste(shanno_mean_sd$shannon_index.x, shanno_mean_sd$shannon_index.y,sep = '±')
summary(aov(shannon_index~Community.type, ingroup_info_CST_shannon) )

#---- Part 3.3 CST count -----
# input_VALENCIA_CST
V_cst = read.csv('VALENCIA_Chinese_917s_output.csv')
V_cst_sub = V_cst[,c('sampleID','subCST','CST')]
names(V_cst_sub)[1] = 'sample_id'


table(V_cst_sub$CST)
round(table(V_cst_sub$CST)/917*100,2)

#---- Part 3.3 pheatmap -----
filter_virgo_abdper_melt_select_genus_group=merge(filter_virgo_abdper_melt_select_genus, ingroup_morph[,c('sample_id','Nugent_3c')], by='sample_id')

filter_virgo_abdper_melt_select_genus_group[which(filter_virgo_abdper_melt_select_genus_group$genus == 'Lactobacillus' & filter_virgo_abdper_melt_select_genus_group$abd_per > 0.8),'Lactobacillus'] = 'Above 80%'
filter_virgo_abdper_melt_select_genus_group[which(filter_virgo_abdper_melt_select_genus_group$genus == 'Lactobacillus' & filter_virgo_abdper_melt_select_genus_group$abd_per < 0.2),'Lactobacillus'] = 'Below 20%'
Lactobacillus_group = unique.data.frame(filter_virgo_abdper_melt_select_genus_group[which(filter_virgo_abdper_melt_select_genus_group$genus == 'Lactobacillus'),c('sample_id','Lactobacillus')])
Lactobacillus_group[is.na(Lactobacillus_group$Lactobacillus),'Lactobacillus'] = '20%-80%'

heatmap_group = merge( ingroup_morph[,c('sample_id','Nugent_3c')],Lactobacillus_group, by='sample_id' )
heatmap_group = merge(heatmap_group, V_cst_sub, by='sample_id')
heatmap_group[which(heatmap_group$CST %in% c('IV-A','IV-B','IV-C')),'CST'] = 'IV'
heatmap_group = heatmap_group[,c("sample_id", "Nugent_3c", "Lactobacillus","CST")]
names(heatmap_group) = c('sample_id','Nugent score','Lactobacillus','Community state types')
heatmap_group$`Community state types` = as.factor(heatmap_group$`Community state types`)
levels(heatmap_group$`Community state types`) = c('CST I','CST II','CST III','CST IV','CST V')
row.names(heatmap_group) = heatmap_group$sample_id
heatmap_group$Lactobacillus = factor(heatmap_group$Lactobacillus,levels =c('Below 20%','20%-80%','Above 80%') )
heatmap_group = merge(heatmap_group, ingroup_info[,c('sample_id','symptoms_group')], by='sample_id')
names(heatmap_group)[5] = 'Self-reported symptoms'

library(RColorBrewer)
ann_colors = list(
  `Community type`=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"),
  `Nugent score`= c("#E5F5E0", "#A1D99B", "#31A354"),
  Lactobacillus = c("#FEE0D2", "#FC9272", "#DE2D26"),
  `Self-reported symptoms` = c("#20B2AA","#FFA500")
)


pheatmap_Chinese_917s = heatmap_filter_plot_func(filter_virgo_abdper_melt_select,heatmap_group, heatmap_group)
ggsave(pheatmap_Chinese_917s, filename = 'plot/02_p3_healthy_all_917s_pheatmap.png',width = 10, height=7)



    
