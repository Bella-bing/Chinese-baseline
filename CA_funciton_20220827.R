# function


## filter taxa based on reads  -----------------------------------------------------------
#If a species is present in <5 reads and only present in <5% samples, then remove the species.

filter_taxa_read_func <- function(summ.abd=summary.Abundance){
  summ.abd_melt = melt(summ.abd, id='sample_id')
  names(summ.abd_melt) = c('sample_id','taxa','reads')
  summ.abd_melt$reads = as.numeric(summ.abd_melt$reads)
  sample_num = length(summ.abd$sample_id)
  max_df = aggregate(reads ~ taxa,data=summ.abd_melt, max)
  names(max_df)[2] = 'max_reads'
  mean_df = aggregate(reads ~ taxa, data=summ.abd_melt, mean)
  names(mean_df)[2] = 'mean_reads'
  stats_df = merge(mean_df,max_df, by='taxa')
  summ.abd_melt2 = summ.abd_melt[which(summ.abd_melt$reads>0),]
  taxa_count_df = as.data.frame(table(summ.abd_melt2$taxa))
  names(taxa_count_df) = c('taxa','number_of_samples')
  taxa_count_df$per = round(taxa_count_df$number_of_samples/sample_num*100,0)
  m = merge(stats_df,taxa_count_df, by='taxa' )
  m_filter = m[which(m$per < 5 & m$mean_reads <5),'taxa']
  m_filter2 = m_filter[-which(m_filter %in% c('Mycoplasma_genitalium','Neisseria_gonorrhoeae','Chlamydia_trachomatis','Ureaplasma_parvum','Ureaplasma_urealyticum'))]
  summ.abd_after_filter = summ.abd[,!names(summ.abd) %in% m_filter2]
  return(summ.abd_after_filter)
}


# # filter data and plot ： percentage of relative abundance %  -----------------------------------------------------------------
### barplot of microbiome

abd_filter_plot_func_CHN <- function(abd_per_df, group){
  abd_df_filter_sample = abd_per_df[which(abd_per_df$sample_id %in% unique(group$sample_id)),]
  species_plot_df = abd_df_filter_sample[which(abd_df_filter_sample$abd_per >= 0.3),] 
  group_samples = unique(group$sample_id)
  species_table = as.data.frame(table(species_plot_df$species)/length(group_samples))
  #species_table = species_table[which(species_table$Freq > 0.01),]
  if(nrow(species_table) > 8){
    order_species = species_table[order(-species_table$Freq),'Var1'][1:8]
  } else {
    order_species = species_table[order(-species_table$Freq),'Var1']
  }
  #order_species = c(order_species[2],order_species[1],order_species[3:8])
  abd_df_group_species_filter = abd_df_filter_sample[which(abd_df_filter_sample$species %in% order_species),]
  n_d =  dcast(abd_df_group_species_filter[,c('sample_id','species','abd_per')], sample_id ~ species, value.var = 'abd_per')
  n_d[is.na(n_d)] = 0
  n_d$other =  1-apply(n_d[,-1],1,sum)
  n_d = n_d[,c('sample_id',as.character(order_species),'other')]
  
  order_sample = unique(n_d[order(-n_d[,as.character(order_species[1])],-n_d[,as.character(order_species[2])],-n_d[,as.character(order_species[3])],-n_d[,as.character(order_species[4])],-n_d[,as.character(order_species[5])]),'sample_id'])
  abd_df_group_species_filter2 = melt(n_d, id = "sample_id" )
  names(abd_df_group_species_filter2) = c('sample_id','species','abd_per')
  abd_df_group2 = merge(abd_df_group_species_filter2, group, by='sample_id')
  abd_df_group2$species = factor(abd_df_group2$species, levels = c(as.character(order_species),'other'))
  abd_df_group2$sample_id = factor(abd_df_group2$sample_id, levels =order_sample )

   p= ggplot(data=abd_df_group2, mapping=aes(x =sample_id, y = abd_per,fill=species))+
    geom_bar(stat="identity",position="stack",width = 1) + 
    labs(x='Samples', y='Relative Abundance')+
    scale_y_continuous(expand=c(0, 0))+
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    facet_grid(~group, scales = 'free')+
    scale_fill_brewer(palette = "Set1")
  return(p)
}

# heatmap and cluster : percentage of relative abundance % ------------------------------------------------------------------
library(pheatmap)
func_pheatmap <- function(tax_df, group_df, taxa_anno_df){
  # tax_df and group_df must by with row.names. tax_df row.names is speices, and group_df row.names is sample_id 
  color = colorRampPalette(c("blue", "white", "red"))(10)
  #ann_colors = list(
  #  `Community type`=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"),
  #  `Nugent score`= c("#E5F5E0", "#A1D99B", "#31A354"),
  #  Lactobacillus = c("#FEE0D2", "#FC9272", "#DE2D26")
  #)
  #ann_colors=list(Group=c(XJ="#EE7AE9",YN="#00447E",Cr="#F34800"))
  pheatmap(as.matrix(tax_df) , 
           scale = "none", #对行归一化 
           show_rownames=T, #显示行名-基因名
           cluster_rows=F, #是否对列进行聚类
           cluster_cols=T, #是否对行进行聚类
           #treeheight_row = 15, #调整纵向数的高度
           #annotation_row=taxa_anno_df, ##将行分组信息加入到列注释信息,如taxa的注释
           annotation_col=group_df, #将列分组信息加入到列注释信息
           clustering_distance_rows = "correlation", #聚类线长度优化
           clustering_distance_cols = "maximum",
           #color = colorRampPalette(c("yellow", "blue", "red"))(11),
           #annotation_colors=ann_colors, #给注释条块的颜色赋值
           #cutree_rows=5, #按聚类树将所有行分为5组
           #cellheight=2,cellwidth=8)#每个各自的大小
           #fontsize_col=6, #列名的字号
           border_color=NA,
           #annotation_names_row=T, #因为taxa anno的名字太长,设置给隐藏掉
           show_colnames = F,
           fontsize_row = 10, #行名的字号
           fontsize_annotation=8 )
  #main = "HSPs")#命名主标题
  #filename="heatmap_HSPs_new.png") #可保存为File type should be: pdf, png, bmp, jpg, tiff
}

heatmap_filter_plot_func <- function(abd_per_df, group, taxa_anno){
  abd_df_filter_sample = abd_per_df[which(abd_per_df$sample_id %in% unique(group$sample_id)),]
  species_abd_sum = aggregate(abd_per ~ species, data =abd_df_filter_sample , sum)
  order_species = species_abd_sum[order(-species_abd_sum$abd),'species'][1:15]
  abd_df_species_filter = abd_df_filter_sample[which(abd_df_filter_sample$species %in% order_species),]
  abd_df_species_filter$species = factor(abd_df_species_filter$species, levels = as.character(order_species))
  abd_df_species_filter_dcast = dcast(abd_df_species_filter, species ~ sample_id, value.var = 'abd_per')
  abd_df_species_filter_dcast[is.na(abd_df_species_filter_dcast)] = 0
  row.names(abd_df_species_filter_dcast) = abd_df_species_filter_dcast$species
  row.names(group) = group$sample_id
  row.names(taxa_anno) = taxa_anno$species
  taxa_anno_filter = taxa_anno[which(taxa_anno$species %in% order_species),]
  p=func_pheatmap(abd_df_species_filter_dcast[-1], group[-1],taxa_anno_filter[-1])
  return(p)
}

#---- gene_pheatmap ------------
gene_heatmap_func <- function(gene_d, group_df,file_name, tree_k){
  
  # hclust using jaccard distance
  row.names(gene_d) = gene_d$sample_id
  sum_col= data.frame(gene_id = names(gene_d)[-1],sum= apply(gene_d[-1],2,sum))
  a = sum_col[which(sum_col$sum == 0),'gene_id']
  if (length(a) == 0){
    gene_d_filter = gene_d
  } else {
    gene_d_filter = gene_d[,!(names(gene_d) %in% a)]
  }
  sample_dist = vegdist(as.matrix(gene_d_filter[-1]), method = "jaccard")
  sample_hclust =  hclust(sample_dist, method = "ward.D")
  
  # regroup
  sample_regroup = as.data.frame(cutree(sample_hclust, k=tree_k))
  names(sample_regroup) = 'group'
  sample_regroup$group = factor(sample_regroup$group)
  levels(sample_regroup$group) = paste('group', 1:tree_k,sep='')
  sample_regroup$sample_id = row.names(sample_regroup)
  group_df_2 = merge(group_df, sample_regroup, by='sample_id')
  row.names(group_df_2) = group_df_2$sample_id
  #plot(sample_hclust)
  gene_dist = vegdist(t(as.matrix(gene_d_filter[-1])), method = "jaccard")
  gene_hclust = hclust(gene_dist, method = "ward.D")
  #plot(gene_hclust)
  
  # sort samples
  gene_d_filter$sample_id = factor(gene_d_filter$sample_id, levels = row.names(group_df))
  matrix_data = t(as.matrix(gene_d_filter[-1]))
  matrix_data = matrix_data[,row.names(group_df)]
  # pheatmap
  ann_colors = list(
    `Community type`=c(`CST I` ="#E41A1C",`CST II` = "#377EB8",`CST III` ="#4DAF4A",`CST IV`="#984EA3",`CST V` ="#FF7F00"),
    group = c(group1 = "#66C2A5", group2 = "#FC8D62", group3 = "#8DA0CB", group4 = "#E78AC3", group5 ="#FEE0D2")[1:tree_k],
    race= c(Han="#66C2A5", Tibetan="#FC8D62",Uighur= "#8DA0CB",other="#E78AC3"),
    `Self-reproted symptoms` = c(Any_symptoms = "#66C2A5",No_symptoms = "#FC8D62"),
    city = c(Beijing="#E41A1C", Chongqing= "#377EB8",Heilongjiang= "#4DAF4A", Tibet="#984EA3", Guangxi="#FF7F00", Shanghai="#FFFF33", Guangdong="#A65628", Shanxi="#F781BF",Xinjiang= "#B3DE69", Hunan="#FCCDE5")
    #Lactobacillus = c("#FEE0D2", "#FC9272", "#DE2D26")
  )
  gene_heatmap = pheatmap(matrix_data , 
                          scale = "none", #对行归一化 
                          show_rownames=F, #显示行名-基因名
                          cluster_rows=gene_hclust, #是否对行进行聚类
                          cluster_cols=sample_hclust, #是否对列进行聚类
                          #cluster_cols=F,
                          color = colorRampPalette(colors = c("white","darkblue"))(100),
                          #treeheight_row = 15, #调整纵向数的高度
                          #annotation_row=group_df['group_row'], ##将行分组信息加入到列注释信息,如taxa的注释
                          annotation_col=group_df_2[-1], #将列分组信息加入到列注释信息
                          #clustering_distance_rows = "correlation", #聚类线长度优化
                          #clustering_distance_cols = "maximum",
                          annotation_colors=ann_colors, #给注释条块的颜色赋值
                          #cutree_rows=5, #按聚类树将所有行分为5组
                          #cellheight=2,cellwidth=8)#每个各自的大小
                          #fontsize_col=6, #列名的字号
                          border_color=NA,
                          #annotation_names_row=T, #因为taxa anno的名字太长,设置给隐藏掉
                          show_colnames = F,
                          fontsize_row = 10, #行名的字号
                          fontsize_annotation=8 ,
                          filename=file_name, width=10, height = 4)
  return(sample_regroup)
}


## chisq test or fisher's extact test in the subgroup in gene heatmap 
subgroup_chisq <- function(meta, factors){
  asso_subgroup_meta_res = data.frame()
  for (i in 1:length(factors)){
    print(paste(c("============> Factors : ", factors[i], " <============"), collapse = ''))
    group_factor_df = meta[,c('sample_id',factors[i],'group')]
    group_factor_df2 = group_factor_df[which(group_factor_df[,factors[i]] != 'missing'),]
    group_factor_df2[,factors[i]] = as.character(group_factor_df2[,factors[i]])
    print(table(group_factor_df2[,factors[i]],group_factor_df2$group))
    fisher_res = fisher.test(table(group_factor_df2[,factors[i]],group_factor_df2$group),simulate.p.value=TRUE)
    print(fisher_res)
    chisq_res = chisq.test(table(group_factor_df2[,factors[i]],group_factor_df2$group))
    print(chisq_res)
    df = as.data.frame(table(group_factor_df2[,factors[i]],group_factor_df2$group))
    names(df) = c('level','group','count')
    df$factor = factors[i]
    df$fisher_p_value = fisher_res$p.value
    df$chis_p_value = chisq_res$p.value
    asso_subgroup_meta_res = rbind(asso_subgroup_meta_res,df)
  }
  return(asso_subgroup_meta_res)
}

### logsistic regresion in the subgroup in gene heatmap 
taxa_subgroup_logistic_func <- function(taxa_subgroup_meta, ingroup_info_CST_multinom, ref_group){
  taxa_meta= merge(taxa_subgroup_meta[,c('sample_id','group')], ingroup_info_CST_multinom, by='sample_id')
  taxa_meta$group = as.factor(taxa_meta$group)
  taxa_meta$group = relevel(taxa_meta$group ,ref=ref_group)
  table(taxa_meta$group)
  row.names(taxa_meta) = taxa_meta$sample_id
  fit = glm(group ~ age+BMI+suffers_insomnia+symptoms_group+city+menstrual_cycle+
              frequency_of_sex+has_children+abnormal_pregnancy+meat+vegan+
              heavy_salt_oil+spicy_deit+sweet_deit+intake_of_yogurt+education_level+income_level , 
            data = taxa_meta,family = binomial(link='logit'))
  
  #P value
  p<-summary(fit)$coefficients[,4]
  #waldvalue
  wald<-summary(fit)$coefficients[,3]^2
  #B value
  valueB<-coef(fit)
  #OR 
  valueOR<-exp(coef(fit))
  #OR 95%CI
  confitOR<-exp(confint(fit))
  OR_total =  data.frame(
    factors = names(p),
    OR=round(valueOR,2),
    pvalue=round(p,5),
    `lower2.5%` = round(confitOR[,1],2),`upper97.5%` = round(confitOR[,2],2),
    OR95CI=paste(round(valueOR,2)," (",round(confitOR[,1],2),"-",round(confitOR[,2],2),")",sep=""))
  return(OR_total)
}



### multinom logsistic regresion in the subgroup in gene heatmap 
taxa_subgroup_multinom_func <- function(taxa_subgroup_meta, ingroup_info_CST_multinom, ref_group){
  taxa_subgroup_meta_new = merge(taxa_subgroup_meta[,c('sample_id','group')], ingroup_info_CST_multinom, by='sample_id')
  taxa_subgroup_meta_new$group  = as.factor(taxa_subgroup_meta_new$group)
  taxa_subgroup_meta_new$group =  relevel(taxa_subgroup_meta_new$group, ref = ref_group)
  row.names(taxa_subgroup_meta) = taxa_subgroup_meta$sample_id
  taxa_multinom =  multinom(group ~ . , data = taxa_subgroup_meta_new[-1])
  summary(taxa_multinom)
  taxa_OR=round(exp(coef(taxa_multinom)),2) #计算OR 
  taxa_OR95 = as.data.frame(round(exp(confint(taxa_multinom)),2)) #计算OR的95%CI
  taxa_OR95$factors = row.names(taxa_OR95)
  # 计算p value
  taxa_z <- summary(taxa_multinom)$coefficients/summary(taxa_multinom)$standard.errors
  taxa_p <-round(((1 - pnorm(abs(taxa_z), 0, 1)) * 2),5)
  taxa_OR_melt = melt(taxa_OR)
  names(taxa_OR_melt) = c('taxa','factors','OR')
  taxa_OR95_melt = melt(taxa_OR95)
  for(i in 1:nrow(taxa_OR95_melt)){
    taxa_OR95_melt[i,'taxa'] = str_split(taxa_OR95_melt[i,'variable'],'%.')[[1]][2]
    taxa_OR95_melt[i,'95CI'] = str_split(taxa_OR95_melt[i,'variable'],'%.')[[1]][1]
  }
  taxa_OR95_melt_2 =dcast(taxa_OR95_melt,taxa+factors ~ `95CI`, value.var = 'value' )
  names(taxa_OR95_melt_2)[c(3,4)] = c('lower2.5%','upper97.5%')
  taxa_p_melt = melt(taxa_p)
  names(taxa_p_melt) = c('taxa','factors','p value')
  
  
  taxa_OR_total = merge(taxa_OR_melt, taxa_p_melt, by=c('taxa','factors'))
  taxa_OR_total = merge(taxa_OR_total, taxa_OR95_melt_2, by=c('taxa','factors'))
  for(i in 1:nrow(taxa_OR_total)){
    taxa_OR_total[i,'OR95CI'] = paste(c(taxa_OR_total[i,'OR'], ' (', taxa_OR_total[i,'lower2.5%'], '-', taxa_OR_total[i,'upper97.5%'], ')'),collapse = '')
  }
  return(taxa_OR_total)
}

## CA adonis analysis in the same CST
CA_adonis_func <- function(subcst){
  cst_sub = CA_umap_taxa_df_vcst[which(CA_umap_taxa_df_vcst$subCST == subcst),c('sample_id','group')]
  cst_sub_abdper_d = CA_abdper_merge_raw_d[which(CA_abdper_merge_raw_d$sample_id %in% cst_sub$sample_id),]
  species_sum = apply(cst_sub_abdper_d[-1], 2, sum)
  cst_sub_abdper_d_2 = cst_sub_abdper_d[,-(which(names(cst_sub_abdper_d) %in% names(which(species_sum == 0)))), ]
  cst_sub_taxdf = data.frame(species=names(cst_sub_abdper_d)[-1])
  row.names(cst_sub_taxdf) = cst_sub_taxdf$species
  CA_subCST_pyloseq = phyloseq(otu_table(t(cst_sub_abdper_d_2[,-1]), taxa_are_rows = TRUE), tax_table(as.matrix(cst_sub_taxdf)))
  CA_subcst_jsd = distance(CA_subCST_pyloseq, method='jsd')
  adonis_subcst = as.data.frame(adonis2(CA_subcst_jsd ~ group , data = cst_sub,permutations = 1000,na.action=na.omit))
  adonis_subcst$subCST = subcst
  return(adonis_subcst)
}

