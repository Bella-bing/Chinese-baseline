#### result2.	Self-reported symptoms and morphological results #########################


ques_split_df = read.delim('ques_selfreported_symptoms_split_df.txt')
ingroup_morph = read.delim('ingroup_morph.txt')
ingroup_info = read.delim('ingroup_info_CST.txt')

# morhology vs with or without sym
ques_symptoms_group = unique.data.frame(ques_split_df[,c('sample_id','symptoms_group')])
table(ques_symptoms_group$symptoms_group)
table(ingroup_morph$microecologial_group,ingroup_morph$symptoms_group)
table(ingroup_morph$rclass,ingroup_morph$symptoms_group)
table(ingroup_morph$Nugent_3c,ingroup_morph$symptoms_group)
table(ingroup_morph$rclass)

# no symptoms
length(unique(ques_split_df[-which(is.na(ques_split_df$No_symptoms)),'sample_id']))

# Any unusual vaginal discharge
length(unique(ques_split_df[-which(is.na(ques_split_df$unusual_vaginal_discharge)),'sample_id']))
unusual_vaginal_discharge_df = merge(unique.data.frame(ques_split_df[-which(is.na(ques_split_df$unusual_vaginal_discharge)),c('sample_id','unusual_vaginal_discharge')]),ingroup_morph, by='sample_id')
table(unusual_vaginal_discharge_df$unusual_vaginal_discharge,unusual_vaginal_discharge_df$microecologial_group)

# Genital itching or burning or pain
length(unique(ques_split_df[-which(is.na(ques_split_df$unusual_genital)),'sample_id']))
unusual_genital_df = merge(unique.data.frame(ques_split_df[-which(is.na(ques_split_df$unusual_genital)),c('sample_id','unusual_genital')]),ingroup_morph, by='sample_id')
table(unusual_genital_df$unusual_genital,unusual_genital_df$microecologial_group)

# Pain during urination or pain during sex
length(unique(ques_split_df[-which(is.na(ques_split_df$panin_urination_sex)),'sample_id']))
panin_urination_sex_df = merge(unique.data.frame(ques_split_df[-which(is.na(ques_split_df$panin_urination_sex)),c('sample_id','panin_urination_sex')]),ingroup_morph, by='sample_id')
table(panin_urination_sex_df$panin_urination_sex,panin_urination_sex_df$microecologial_group)

all_sys = merge(ques_split_df[,c('sample_id','symptom_q_split')],ingroup_morph, by='sample_id')
table(all_sys$symptom_q_split, all_sys$microecologial_group)
chisq.test(matrix(c(497,101,227,74),2,2))

# qPCR
qPCR=read.delim('qPCR_0918.txt')
qPCR$species = NA
names(qPCR)[2] = 'species'
names(qPCR)[1] = 'sample_id'
qPCR_d= dcast(qPCR, sample_id ~ species, value.var = 'pclass')
qPCR_d_ingroup_morph = merge(qPCR_d,ingroup_morph, by='sample_id')
table(qPCR_d_ingroup_morph$microecologial_group,qPCR_d_ingroup_morph$symptoms_group)
# CT
table(qPCR_d_ingroup_morph$Chlamydia_trachomatis,qPCR_d_ingroup_morph$microecologial_group,qPCR_d_ingroup_morph$symptoms_group)
table(qPCR_d_ingroup_morph$Neisseria_gonorrhoeae,qPCR_d_ingroup_morph$microecologial_group,qPCR_d_ingroup_morph$symptoms_group)
#

