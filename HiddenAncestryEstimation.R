# load libraries
library(tidyverse)
library(Summix)
# Fst library
library(BEDASSLE)
library(RColorBrewer)

# colorblind color values
brewerval = 7
colorval = c(brewer.pal(9, 'Purples')[brewerval],
             brewer.pal(9, 'Blues')[brewerval],
             brewer.pal(9, 'Oranges')[brewerval - 1],
             brewer.pal(9, 'Greens')[brewerval],
             brewer.pal(9, 'Reds')[brewerval])

# load in Summix paper genome data by chromosome, and merge
gframe = data.frame()

for (i in 1:22){
  g_l_dat = read.delim(paste0("~/Summix_genome_data_byCHR/Summix_genome_data_CHR", i, ".txt.gz"))
  g_l_dat2 = g_l_dat[,c(1:11, 92:98)]
  gframe = rbind(gframe, g_l_dat2)
}

# Sample 100K SNPS
gsample = gframe %>% 
  sample_n(100000)

# Simulate admixed population of 10,000
# Manually change numbers to hit 10k (CHANGE)

### EUR POP
#CHANGE
eur_name = "ref_AF_eur_1000G"
eur_pop = 4000

eur_geno = t(sapply(gsample[[eur_name]], function(x){x2<-as.numeric(x); rmultinom(1, eur_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
eur_dat = data.frame(AC = 2*eur_geno[,1] + eur_geno[,2],
                     AN = 2*eur_pop) %>% 
  mutate(AF = AC / AN)

### AFR POP
#CHANGE
afr_name = "ref_AF_afr_1000G"
afr_pop = 200

afr_geno = t(sapply(gsample[[afr_name]], function(x){x2<-as.numeric(x); rmultinom(1, afr_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
afr_dat = data.frame(AC = 2*afr_geno[,1] + afr_geno[,2],
                     AN = 2*afr_pop) %>% 
  mutate(AF = AC / AN)

### EAS POP
#CHANGE
eas_name = "ref_AF_eas_1000G"
eas_pop = 3300

eas_geno = t(sapply(gsample[[eas_name]], function(x){x2<-as.numeric(x); rmultinom(1, eas_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
eas_dat = data.frame(AC = 2*eas_geno[,1] + eas_geno[,2],
                     AN = 2*eas_pop) %>% 
  mutate(AF = AC / AN)

### SAS POP
#CHANGE
sas_name = "ref_AF_sas_1000G"
sas_pop = 1500

sas_geno = t(sapply(gsample[[sas_name]], function(x){x2<-as.numeric(x); rmultinom(1, sas_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
sas_dat = data.frame(AC = 2*sas_geno[,1] + sas_geno[,2],
                     AN = 2*sas_pop) %>% 
  mutate(AF = AC / AN)

### IAM POP
#CHANGE
iam_name = "ref_AF_iam_1000G"
iam_pop = 1000

iam_geno = t(sapply(gsample[[iam_name]], function(x){x2<-as.numeric(x); rmultinom(1, iam_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
iam_dat = data.frame(AC = 2*iam_geno[,1] + iam_geno[,2],
                     AN = 2*iam_pop) %>% 
  mutate(AF = AC / AN)

# combine data into final admixed genotype
combdat = data.frame(AC = afr_dat$AC + eur_dat$AC + eas_dat$AC + sas_dat$AC + iam_dat$AC,
                     AN = afr_dat$AN + eur_dat$AN + eas_dat$AN + sas_dat$AN + iam_dat$AN) %>% 
  mutate(AF = AC / AN)

# set reference and observed AFs
#CHANGE
sumframe = data.frame(ref_eur = gsample$ref_AF_eur_1000G,
                      ref_afr = gsample$ref_AF_afr_1000G,
                      ref_eas = gsample$ref_AF_eas_1000G,
                      ref_sas = gsample$ref_AF_sas_1000G,
                      ref_iam = gsample$ref_AF_iam_1000G,
                      obs = combdat$AF)

# run summmix
sumres = summix(sumframe,
                c("ref_eur", "ref_afr", "ref_eas", "ref_sas", "ref_iam"),
                "obs")

# extract allele count for fst
ac_frame = t(as.matrix(data.frame(eur_dat$AC,
                                  afr_dat$AC,
                                  eas_dat$AC,
                                  sas_dat$AC,
                                  iam_dat$AC,
                                  combdat$AC)))

# extract allele number for fst
an_frame = t(as.matrix(data.frame(eur_dat$AN,
                                  afr_dat$AN,
                                  eas_dat$AN,
                                  sas_dat$AN,
                                  iam_dat$AN,
                                  combdat$AN)))
# test order to name output fst matrix
testorder = c("ref_eur", "ref_afr", "ref_eas", "ref_sas", "ref_iam", "obs")

# calculate fst
pairfst = calculate.all.pairwise.Fst(
  ac_frame,
  an_frame
)

# display results
fstout = data.frame(pairfst)
names(fstout) = testorder; rownames(fstout) = testorder
fstout
sumres

# remove one reference ancestry for algorithm ("hidden" ancestry)
#CHANGE
br_test = gsample %>%
  select(ref_AF_eur_1000G, ref_AF_afr_1000G, ref_AF_eas_1000G, ref_AF_iam_1000G) %>%
  mutate(obs_AF = combdat$AF,
         ref_U = runif(100000))

# initialize variables
emiter = 0
finalframe_small = data.frame()
# proportion of reference ancestry removed
trueprop = 0.15

# block relaxation algorithm
while (TRUE){
  # count iterations
  emiter = emiter + 1
  # summix
  sumres = ancestr(br_test[,c(1:4,6)], br_test[,5])
  
  # Update allele frequency after summix has run
  #CHANGE
  update_AF = (br_test$obs_AF - sumres[1]*br_test$ref_AF_eur_1000G - sumres[2]*br_test$ref_AF_afr_1000G - sumres[3]*br_test$ref_AF_eas_1000G - sumres[4]*br_test$ref_AF_iam_1000G)/sumres[5]
  update_AF[update_AF > 1] = 1
  update_AF[update_AF < 0] = 0
  br_test$ref_U = update_AF
  
  # save iteration output
  outframe = data.frame(AFRp = sumres[2],
                        EURp = sumres[1],
                        EASp = sumres[3],
                        IAMp = sumres[4],
                        Up = sumres[5],
                        objective = sumres[6],
                        iter = emiter)
  # bind to output frame
  finalframe_small = rbind(finalframe_small, outframe)
  
  # terminate loop when hidden ancestry gets within 0.5% of true proportion
  if (abs(sumres[5] - trueprop) < 0.005) return()
}

# plot difference betwen true and unknown AFs
#CHANGE
afdif = br_test$ref_U - gsample$ref_AF_sas_1000G
hist(afdif)

plot(br_test$ref_U,
     gsample$ref_AF_sas_1000G)

# plot BR algo output
ggplot(finalframe_small) +
  geom_line(aes(iter, EURp), color = colorval[3], size = 1)+
  geom_line(aes(iter, AFRp), color = colorval[1], size = 1)+
  geom_line(aes(iter, EASp), color = colorval[2], size = 1)+
  geom_line(aes(iter, IAMp), color = colorval[4], size = 1)+
  geom_line(aes(iter, Up), color = "black", size = 1) +
  geom_hline(aes(yintercept = 0.4), color = colorval[3], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.02), color = colorval[1], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.33), color = colorval[2], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.1), color = colorval[4], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.15), color = "black", linetype = "dashed", size = 1) +
  labs(x = "Iterations",
       y = "Ancestry Proportion") +
  theme_bw()

# BR algo output
finalframe_small

##### IAM GNOMAD TEST - REMOVE IAM
# Same algorithm from before, using real gnomad data.
# Remove IAM from reference, and see if BR algo returns close proportions to Summix

br_test = gframe %>%
  select(ref_AF_eur_1000G, ref_AF_afr_1000G, ref_AF_eas_1000G, ref_AF_sas_1000G, gnomad_AF_amr)%>% 
  drop_na() %>%
  mutate(ref_U = runif(582155)) 

emiter = 0
finalframe_small = data.frame()

while (TRUE){
  emiter = emiter + 1
  
  sumres = ancestr(br_test[,c(1:4,6)], br_test[,5])
  
  update_AF = (br_test$gnomad_AF_amr - sumres[1]*br_test$ref_AF_eur_1000G - sumres[2]*br_test$ref_AF_afr_1000G - sumres[3]*br_test$ref_AF_eas_1000G - sumres[4]*br_test$ref_AF_sas_1000G)/sumres[5]
  update_AF[update_AF > 1] = 1
  update_AF[update_AF < 0] = 0
  br_test$ref_U = update_AF
  
  outframe = data.frame(AFRp = sumres[2],
                        EURp = sumres[1],
                        EASp = sumres[3],
                        SASp = sumres[4],
                        Up = sumres[5],
                        objective = sumres[6],
                        iter = emiter)
  
  finalframe_small = rbind(finalframe_small, outframe)
  # terminate when summix objective value is < 5.8
  if (sumres[6] < 5.8) return()
}

which(is.na(gframe$gnomad_AF_amr))
iam_refdat = gframe$ref_AF_iam_1000G[-which(is.na(gframe$gnomad_AF_amr))]

#CHANGE
afdif = br_test$ref_U - iam_refdat
hist(afdif)

plot(br_test$ref_U,
     iam_refdat)

ggplot(finalframe_small) +
  geom_line(aes(iter, EURp), color = colorval[3], size = 1)+
  geom_line(aes(iter, AFRp), color = colorval[1], size = 1)+
  geom_line(aes(iter, EASp), color = colorval[2], size = 1)+
  geom_line(aes(iter, SASp), color = colorval[5], size = 1)+
  geom_line(aes(iter, Up), color = "black", size = 1) +
  geom_hline(aes(yintercept = 0.505), color = colorval[3], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.058), color = colorval[1], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.038), color = colorval[2], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.019), color = colorval[5], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.380), color = "black", linetype = "dashed", size = 1)+
  labs(x = "Iterations",
       y = "Ancestry Proportion") +
  theme_bw()

# GNOMAD AFR FULL ANC TEST
# Use all 5 reference ancestries, and 6th unknown with gnomad AFR data.

br_test = gframe %>%
  select(ref_AF_eur_1000G, ref_AF_afr_1000G, ref_AF_eas_1000G, ref_AF_sas_1000G, ref_AF_iam_1000G, gnomad_AF_afr) %>%
  mutate(ref_U = runif(582156))

emiter = 0
finalframe_small = data.frame()

while (TRUE){
  emiter = emiter + 1
  
  sumres = ancestr(br_test[,c(1:5,7)], br_test[,6])
  
  update_AF = (br_test$gnomad_AF_afr - sumres[1]*br_test$ref_AF_eur_1000G - sumres[2]*br_test$ref_AF_afr_1000G - sumres[3]*br_test$ref_AF_eas_1000G - sumres[4]*br_test$ref_AF_sas_1000G - sumres[5]*br_test$ref_AF_iam_1000G)/sumres[6]
  update_AF[update_AF > 1] = 1
  update_AF[update_AF < 0] = 0
  br_test$ref_U = update_AF
  
  outframe = data.frame(AFRp = sumres[2],
                        EURp = sumres[1],
                        EASp = sumres[3],
                        SASp = sumres[4],
                        IAMp = sumres[5],
                        Up = sumres[6],
                        objective = sumres[7],
                        iter = emiter)
  
  finalframe_small = rbind(finalframe_small, outframe)
  
  if (sumres[7] < 5.8) return()
}

ggplot(finalframe_small) +
  geom_line(aes(iter, EURp), color = colorval[3], size = 1)+
  geom_line(aes(iter, AFRp), color = colorval[1], size = 1)+
  geom_line(aes(iter, EASp), color = colorval[2], size = 1)+
  geom_line(aes(iter, SASp), color = colorval[5], size = 1)+
  geom_line(aes(iter, IAMp), color = colorval[4], size = 1)+
  geom_line(aes(iter, Up), color = "black", size = 1) +
  geom_hline(aes(yintercept = 0.157), color = colorval[3], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.825), color = colorval[1], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.005), color = colorval[2], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.005), color = colorval[5], linetype = "dashed", size = 1)+
  geom_hline(aes(yintercept = 0.008), color = colorval[4], linetype = "dashed", size = 1)+
  labs(x = "Iterations",
       y = "Ancestry Proportion") +
  theme_bw()

# Pseudo FST estimates using an arbitrary population numbers
# true Allele Numbers can be extracted and used
fstpopnum = 1000

ac_frame = t(as.matrix(data.frame(floor(br_test$ref_AF_eur_1000G*fstpopnum),
                                  floor(br_test$ref_AF_afr_1000G*fstpopnum),
                                  floor(br_test$ref_AF_eas_1000G*fstpopnum),
                                  floor(br_test$ref_AF_sas_1000G*fstpopnum),
                                  floor(br_test$ref_AF_iam_1000G*fstpopnum),
                                  floor(br_test$ref_U*fstpopnum))))

an_frame = t(as.matrix(data.frame(eur_dat=rep(fstpopnum, 582156),
                                  afr_dat=rep(fstpopnum, 582156),
                                  eas_dat=rep(fstpopnum, 582156),
                                  sas_dat=rep(fstpopnum, 582156),
                                  iam_dat=rep(fstpopnum, 582156),
                                  combdat=rep(fstpopnum, 582156))))

testorder = c("ref_eur", "ref_afr", "ref_eas", "ref_sas", "ref_iam", "UNKNOWN")

pairfst = calculate.all.pairwise.Fst(
  ac_frame,
  an_frame
)

fstout = data.frame(pairfst)
names(fstout) = testorder; rownames(fstout) = testorder