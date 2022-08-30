###################### result6: Chinese VS North American ################################

library(vegan)
library(ggplot2)
library(reshape2)
source('CA_funciton_20220827.R')

set.seed(123)
#----------------------- 6.1 CSTs between Chinese vs North American: study cohort1  -------------------------------
# Chinese CSTs
V_cst = read.csv('VALENCIA_Chinese_917s_output.csv')
V_cst_sub = V_cst[,c('sampleID','subCST','CST')]
names(V_cst_sub)[1] = 'sample_id'
V_cst_sub$region = 'Chinese'



# NA CTTs cohort1 : Proportion of CSTs
usa_virgo = read.csv('other_results/USA_from_Mabing/all_samples_taxonomic_composition_data[31].csv')
names(usa_virgo)[1] = 'sample_id'
row.names(usa_virgo) = usa_virgo$sample_id

# meta data
usa_meta = read.csv('other_results/USA_from_Mabing/all_samples_metadata[80]_uniq.csv')
names(usa_meta)[1] = 'sample_id'
usa_meta = merge(usa_meta, usa_virgo[,c(1:7)], by=c('sample_id','Subject_number'))
usa_meta_select = usa_meta[which(usa_meta$age >=21 & usa_meta$age <= 45),]
usa_sample_select = unique(usa_meta_select$sample_id)
usa_meta_select = usa_meta[which(usa_meta$sample_id %in% usa_sample_select),]

# age
#--------------------------- sample  summary ------------------------------------------

usa_meta_select[which(usa_meta_select$city == 'Burmingham'),'city'] = 'Birmingham'
table(usa_meta_select$Race)
# city , race
usa_meta_select_subject = unique.data.frame(usa_meta_select[,c('Subject_number','city','Race','age')])
usa_meta_select_subject = usa_meta_select_subject[-which(duplicated(usa_meta_select_subject$Subject_number)),]
table(usa_meta_select_subject$Race)
table(usa_meta_select_subject$city,usa_meta_select_subject$Race)

# age
summary(usa_meta_select_subject$age)
sd(usa_meta_select_subject$age)
# age, race
age_mean = aggregate(data=usa_meta_select_subject, age ~ Race, mean)
age_sd = aggregate(data=usa_meta_select_subject, age ~ Race, sd)
age_mean_sd = merge(age_mean, age_sd, by='Race')
age_mean_sd$age.x = round(age_mean_sd$age.x,2)
age_mean_sd$age.y = round(age_mean_sd$age.y,2)
age_mean_sd$mean_sd = paste(age_mean_sd$age.x,age_mean_sd$age.y,sep='±')
age_mean_sd

# Nugent score
table(usa_meta_select$nugent_cat, usa_meta_select$Race)

# CSTs 
V_cst_sub
usa_cst_cohort1 = usa_meta_select[,c('sample_id','Race','Val_CST','Val_subCST')]
names(usa_cst_cohort1) = c('sample_id','region','CST','subCST')

CA_CSTs_cohort1 = rbind(V_cst_sub, usa_cst_cohort1)
table( CA_CSTs_cohort1$CST, CA_CSTs_cohort1$region)/5600
table(CA_CSTs_cohort1$subCST,CA_CSTs_cohort1$region)

#-------- 6.2 MG between Chinese vs North American  -------------------------------
#-------- 6.2.1 : Data processing of North American  -------------
# north American women taxonomy
usa_virgo_tmp = read.delim('other_results/USA_from_Mabing/MG_NA_082322_summary_abd.txt')
usa_virgo = dcast(usa_virgo_tmp, sample_id~Taxa, value.var = 'reads')
usa_virgo = usa_virgo[,-2]
usa_virgo[is.na(usa_virgo)] = 0
usa_virgo_filter1 = filter_taxa_read_func(usa_virgo)
row.names(usa_virgo_filter1) = usa_virgo_filter1$sample_id
usa_virgo_filter2 = usa_virgo_filter1[,!names(usa_virgo_filter1) %in% c('Burkholderia_mallei','Coprobacillus_sp.',
                                                                        'Thermobifida_cellulosilytica','Paracoccus_denitrificans',
                                                                        'Neisseria_gonorrhoeae','Actinobacillus_pleuropneumoniae')]
usa_taxdf = data.frame(species = names(usa_virgo_filter2)[-1])
usa_taxdf$genus = NA
for (i in 1:nrow(usa_taxdf)){
  usa_taxdf[i,2] = str_split(usa_taxdf[i,1],'_')[[1]][1]
}
row.names(usa_taxdf) = usa_taxdf$species

usa_physeq_virgo = phyloseq(otu_table(t(usa_virgo_filter2[,-1]), taxa_are_rows = TRUE), tax_table(as.matrix(usa_taxdf)))
usa_physeq_virgo_abdper = transform_sample_counts(usa_physeq_virgo, function(x) x/sum(x))
filter_usa_virgo_abdper = as.data.frame(otu_table(usa_physeq_virgo_abdper))
filter_usa_virgo_abdper$species = row.names(filter_usa_virgo_abdper)
usa_abdper_melt = melt(filter_usa_virgo_abdper,id='species')

# North American: 1508 samples ,254 species
usa_sample_select = unique(usa_abdper_melt$sample_id)
names(usa_abdper_melt)=c('species','sample_id','abd_per')
usa_abdper_melt_s = usa_abdper_melt[which(usa_abdper_melt$sample_id %in% usa_sample_select),]

#--------------- Part 6.2.2 USA  CSTs ------------------------------

# input_VALENCIA_CST
usa_V_cst = read.csv('VALENCIA_USA_1508s_output.csv')
usa_V_cst_sub = usa_V_cst[,c('sampleID','subCST','CST')]
names(usa_V_cst_sub)[1] = 'sample_id'
usa_V_cst_sub$region = 'North_American'


#-------- Part 6.2.3 : Data processing of Chinese and North American  -------------
# Chinese: 917 samples 
## merge two cohort and filter
filter_virgo_abdper_melt_select = read.delim('filter_virgo_abdper_melt_select.txt')
CA_abdper_merge_raw = rbind(usa_abdper_melt_s,filter_virgo_abdper_melt_select)
CA_abdper_merge_raw_d = dcast(CA_abdper_merge_raw, sample_id ~ species,value.var = 'abd_per')
CA_abdper_merge_raw_d[is.na(CA_abdper_merge_raw_d)] = 0
row.names(CA_abdper_merge_raw_d) = CA_abdper_merge_raw_d$sample_id
apply(CA_abdper_merge_raw_d[,-1],2,sum)
CA_abdper_merge_raw_melt = melt(CA_abdper_merge_raw_d,id='sample_id')
names(CA_abdper_merge_raw_melt) = c('sample_id','species','abd_per')
# filter  mean_abdper > 0.001 propertion > 0.1, 44 species pass filter
CA_mean_species =  aggregate(data = CA_abdper_merge_raw_melt,abd_per ~ species, mean)
CA_abdper_melt_cutoff1 = CA_mean_species[which(CA_mean_species$abd_per >=0.001),'species']
CA_abdper_melt_non0 = CA_abdper_merge_raw[which(CA_abdper_merge_raw$abd_per !=0),]
CA_taxa_table = as.data.frame(table(CA_abdper_melt_non0$species)/length(unique(CA_abdper_melt_non0$sample_id)))
CA_taxa_table$Freq = round(CA_taxa_table$Freq,1)
CA_species_cutoff2 = CA_taxa_table[which(CA_taxa_table$Freq >= 0.05),'Var1']
CA_species_cutoff_final = intersect(CA_abdper_melt_cutoff1, CA_species_cutoff2)
CA_abdper_melt_filter = CA_abdper_merge_raw[which(CA_abdper_merge_raw$species %in% CA_species_cutoff_final ),]

usa_final_species = unique(CA_abdper_melt_filter[which(CA_abdper_melt_filter$sample_id %in% usa_sample_select),'species'])
chn_final_species = unique(CA_abdper_melt_filter[which(CA_abdper_melt_filter$sample_id %in% unique(filter_virgo_abdper_melt_select$sample_id)),'species'])
CA_abdper_filter_d = dcast(CA_abdper_melt_filter, sample_id ~ species, value.var = 'abd_per')
row.names(CA_abdper_filter_d) = CA_abdper_filter_d$sample_id

setdiff(usa_final_species,chn_final_species)
setdiff(chn_final_species,usa_final_species)
intersect(chn_final_species,chn_final_species)
# USA and Chinese before filter taxa
library(VennDiagram)
venn.diagram(
  x = list(unique(usa_abdper_melt_s$species),
           unique(filter_virgo_abdper_melt_select$species)),
  category.names = c("North American","Chinese"),
  fill = c("#d37a20","#dbcb09"),
  filename = 'plot/08_USA_CHN_taxa_before_filter_venn_MGonly.png',
  cex=.5,cat.cex=.5,width=1600,height=1600
)
all_taxa_USA_CHN = merge(data.frame(taxa = unique(usa_abdper_melt_s$species), from='North American'), data.frame(taxa =  unique(filter_virgo_abdper_melt_select$species), from = 'Chinese'), by='taxa',all=T)
write.table(all_taxa_USA_CHN, 'other_results/08_all_taxa_USA_CHN_MGonly.txt',sep='\t',row.names=F)

# USA and Chinese after filter taxa
venn.diagram(
  x = list(usa_final_species,
           chn_final_species),
  category.names = c("North American","Chinese"),
  fill = c("#3a9cbc","#a30019"),
  filename = 'plot/08_USA_CHN_taxa_after_filter_venn_MGonly.png',
  cex=.5,cat.cex=.5,width=1600,height=1600
)
all_taxa_USA_CHN_after_filter = merge(data.frame(taxa = unique(usa_abdper_melt_s$species), from='North American'), data.frame(taxa =  unique(filter_virgo_abdper_melt_select$species), from = 'Chinese'), by='taxa',all=T)
write.table(all_taxa_USA_CHN_after_filter, 'other_results/08_all_taxa_USA_CHN_after_filter_MGonly.txt',sep='\t',row.names=F,quote = F)


#-------- Part 6.2.4: CA : CST and UMAP : NA cohort2 -------
library(umap)
library(uwot)
CA_CST = rbind(V_cst_sub,usa_V_cst_sub)
names(CA_CST)[4] = 'group'
table(CA_CST$group,CA_CST$CST)
table(CA_CST$group,CA_CST$subCST)

CA_umap_res = umap(CA_abdper_merge_raw_d[-c(1)])

#umap_res = umap(jsd_564s)
CA_umap_res = as.data.frame(CA_umap_res)
names(CA_umap_res)[c(1,2)]=c('UMAP1','UMAP2')
CA_umap_res$sample_id = row.names(CA_umap_res)
CA_umap_taxa_df_vcst  = merge(CA_umap_res, CA_CST, by='sample_id')
CA_umap_taxa_df_vcst$Val_subCST = paste('CST ', CA_umap_taxa_df_vcst$subCST, sep='')

CA_CST_count = data.frame(table(CA_umap_taxa_df_vcst$CST,CA_umap_taxa_df_vcst$group))
names(CA_CST_count) = c('CST','group','count')
CA_subCST_count = data.frame(table(CA_umap_taxa_df_vcst$subCST,CA_umap_taxa_df_vcst$group))
names(CA_subCST_count) = c('subCST','group','count')

CA_subCST_plot = ggplot(data = CA_subCST_count, aes(x=subCST, y= count, fill = group)) +
  geom_bar(stat = 'identity',width = .7, position = 'dodge') +
  labs(x='Community state types', y='Number of sample')+
  scale_fill_brewer(palette="Dark2")+
  ylim(0, 600)
CA_subCST_plot
ggsave(CA_subCST_plot, filename = 'plot/08_CA_subCST_count_barplot_MGonly.png',width=8, height = 6)

CA_umap_CST = ggplot(CA_umap_taxa_df_vcst,aes(x=UMAP1,y=UMAP2))+
  geom_point(aes(colour = Val_subCST),size=1)+
  scale_colour_brewer(palette = "Set1")+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black")) +
  guides(color=guide_legend(override.aes = list(size=2)))

CA_umap_CST
ggsave(CA_umap_CST,filename = 'plot/08_CA_umap_euclidean_subCST.png', width=6, height=4)



CA_umap_region = ggplot(CA_umap_taxa_df_vcst,aes(UMAP1,UMAP2))+
  geom_point(aes(colour = group ),size=0.3)+
  scale_colour_brewer(palette = "Set1")+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw()+
  facet_wrap(~ Val_subCST , scales = "free", ncol = 3)+
  theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=1,colour="black")) +
  guides(color=guide_legend(override.aes = list(size=2)))
CA_umap_region

ggsave(CA_umap_region,filename = 'plot/08_CA_2022_umap_euclidean_region_subCST.png', width=8, height=6)

# same CST, Chinese vs North American adonis
CA_taxadf = data.frame(species=names(CA_abdper_merge_raw_d)[-1])
row.names(CA_taxadf) = CA_taxadf$species

# adonis analysis in the same subCST
adonis_IA = CA_adonis_func('I-A')
adonis_IB = CA_adonis_func('I-B')
adonis_II = CA_adonis_func('II')
adonis_IIIA = CA_adonis_func('III-A')
adonis_IIIB = CA_adonis_func('III-B')
adonis_IVA = CA_adonis_func('IV-A')
adonis_IVB = CA_adonis_func('IV-B')
adonis_IVC0 = CA_adonis_func('IV-C0')
adonis_V = CA_adonis_func('V')
adonis_all = rbind(adonis_IA, adonis_IB, adonis_II, adonis_IIIA, adonis_IIIB, adonis_IVA, adonis_IVB, adonis_IVC0, adonis_V)
write.table(adonis_all,'other_results/09_CA_CST_adonis_MGonly.txt',sep='\t',row.names = F, quote = F)

### ------------ Part 6.2.5 CST I, II, III, IV ，V lefse analysis ----------------------
# lefse input
CSTI_sample_df = CA_CST[which(CA_CST$CST == 'I'),]
CSTI_lefse_input = lefse_input_func(CA_abdper_melt_filter, CSTI_sample_df[,c('sample_id','group')])
write.table(CSTI_lefse_input, 'other_results/lefse_res/CA_lefse_input/CSTI_lefse_input.txt',sep='\t',row.names = T, col.names=F, quote = F)

CSTII_sample_df = CA_CST[which(CA_CST$CST == 'II'),]
CSTII_lefse_input = lefse_input_func(CA_abdper_melt_filter, CSTII_sample_df[,c('sample_id','group')])
write.table(CSTII_lefse_input, 'other_results/lefse_res/CA_lefse_input/CSTII_lefse_input.txt',sep='\t',row.names = T, col.names=F, quote = F)

CSTIII_sample_df = CA_CST[which(CA_CST$CST == 'III'),]
CSTIII_lefse_input = lefse_input_func(CA_abdper_melt_filter, CSTIII_sample_df[,c('sample_id','group')])
write.table(CSTIII_lefse_input, 'other_results/lefse_res/CA_lefse_input/CSTIII_lefse_input.txt',sep='\t',row.names = T, col.names=F, quote = F)

CSTV_sample_df = CA_CST[which(CA_CST$CST == 'V'),]
CSTV_lefse_input = lefse_input_func(CA_abdper_melt_filter, CSTV_sample_df[,c('sample_id','group')])
write.table(CSTV_lefse_input, 'other_results/lefse_res/CA_lefse_input/CSTV_lefse_input.txt',sep='\t',row.names = T, col.names=F, quote = F)

CSTIV_sample_df = CA_CST[which(CA_CST$CST %in% c('IV-A','IV-B','IV-C')),]
CSTIV_lefse_input = lefse_input_func(CA_abdper_melt_filter, CSTIV_sample_df[,c('sample_id','group')])
write.table(CSTIV_lefse_input, 'other_results/lefse_res/CA_lefse_input/CSTIV_lefse_input.txt',sep='\t',row.names = T, col.names=F, quote = F)


# barplot

boxplot_taxa_func <- function(CST_species_list, cst){
  cst_df = CA_umap_taxa_df_vcst[which(CA_umap_taxa_df_vcst$CST %in% cst),c('sample_id','group')]
  abd_per_vstsub_df = merge(CA_abdper_merge_raw[which(CA_abdper_merge_raw$species %in% CST_species_list),], cst_df, by='sample_id')
  abd_per_vstsub_df$log_abd_per = log10(abd_per_vstsub_df$abd_per)
  
  pcst <- ggboxplot(abd_per_vstsub_df, x = "species", y = "log_abd_per",
                    color = "group", palette = "jama"
                    #add = "jitter"
  )
  # palette可以按照期刊选择相应的配色，如"npg"等
  pcst_boxp = pcst + stat_compare_means(aes(group = group),method = "wilcox.test", label = "p.signif") +
    labs(x='Species', y='Relative Abundance (%)') +
    theme(axis.text.x = element_text(angle = 30, hjust = 1,size = 10))
  return(pcst_boxp)
}

CSTI_species_list = c('Enterococcus_faecalis','Lactobacillus_jensenii','Lactobacillus_gasseri','Bifidobacterium_breve','Bifidobacterium_longum','Propionibacterium_sp','Enterococcus_faecium')
CSTI_taxa_plot = boxplot_taxa_func(CSTI_species_list, 'I')
ggsave(CSTI_taxa_plot, filename = 'plot/08_CA_CSTI_lefse_taxa_boxplot.png', width = 8, height= 6)

CSTII_species_list = c('Bifidobacterium_longum','Gardnerella_vaginalis','Lactobacillus_jensenii','Sneathia_amnii','Prevotella_melaninogenica','Bifidobacterium_breve','Enterococcus_faecium','Lactobacillus_johnsonii')
CSTII_taxa_plot = boxplot_taxa_func(CSTII_species_list, 'II')
ggsave(CSTII_taxa_plot, filename = 'plot/08_CA_CSTII_lefse_taxa_boxplot.png', width = 8, height= 6)

CSTIII_species_list = c('Gardnerella_vaginalis','Lactobacillus_jensenii','Streptococcus_agalactiae','BVAB1','Lactobacillus_gasseri','Propionibacterium_sp_','Streptococcus_anginosus','Prevotella_bivia','Enterococcus_faecium','Escherichia_coli','Enterococcus_faecalis','Bifidobacterium_breve','Streptococcus_bovis','Lactobacillus_iners')
CSTIII_taxa_plot = boxplot_taxa_func(CSTIII_species_list, 'III')
ggsave(CSTIII_taxa_plot, filename = 'plot/08_CA_CSTIII_lefse_taxa_boxplot.png', width = 8, height= 6)

CSTIV_species_list = c('BVAB1','Lactobacillus_iners','Ruminococcus_lactaris','Megasphaera_genomosp_','Prevotella_amnii','Prevotella_sp_','Mobiluncus_mulieris','Sneathia_amnii','Prevotella_buccalis','Peptostreptococcus_anaerobius','Prevotella_timonensis','Mageeibacillus_indolicus','Anaerococcus_tetradius','Lactobacillus_jensenii','Sneathia_sanguinegens','Prevotella_disiens','Porphyromonas_uenonis','Prevotella_melaninogenica','Propionibacterium_sp_','Lactobacillus_gasseri','Bifidobacterium_longum','Enterococcus_faecium','Bifidobacterium_breve','Escherichia_coli','Prevotella_bivia','Streptococcus_anginosus')
CSTIV_taxa_plot = boxplot_taxa_func(CSTIV_species_list, c('IV-A','IV-B','IV-C0'))
ggsave(CSTIV_taxa_plot, filename = 'plot/08_CA_CSTIV_lefse_taxa_boxplot.png', width = 8, height= 6)

CSTV_species_list = c('Lactobacillus_crispatus', 'Enterococcus_faecium')
CSTV_taxa_plot = boxplot_taxa_func(CSTV_species_list, c('V'))
ggsave(CSTV_taxa_plot, filename = 'plot/08_CA_CSTV_lefse_taxa_boxplot.png', width = 8, height= 6)
