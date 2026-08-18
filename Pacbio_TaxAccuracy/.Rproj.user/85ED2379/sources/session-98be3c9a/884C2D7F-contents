# Determine Pacbio Fidelity and Accuracy
# Manually changed in excel:
# 1. replace all space with _
# 2. unassigned -> NA
# 3. ITS XXXIncertaesedis across each taxonomic level or XXXsp as species -> NA
# 4. ITS species that are genusspecies, keep only species

library(tidyverse)

#### LOAD DATA #####

# read matching of curated DB LSU vs SSU and ITS:
match <- read.csv("data/SSU_ITS_LSU_AMFASV_long-read_database_regioncomparison_July12026.csv")

# read final curated DB & join match :1122
dat <- read.csv("data/SSU_ITS_LSU_AMFASV_long-read_database_July12026.csv") %>%
  left_join(match, by = "ASV_updated")
length(unique(dat$ASV_culture)) # 869

##########################################
########## Get True LSU-based AMF ########
##########################################

# non-contaminated AMF asvs
length(unique(dat$ASV_updated)) #1122
length(unique(dat$LSU_family)) # 11
length(unique(dat$LSU_genus)) # 20
length(unique(dat$ASV_genus)) # 20
length(unique(dat$ASV_species)) # 57
length(unique(dat$ASV_culture)) # 869

##########################################
######### percent taxonomy match #########
######## What is primer accuracy? ########
##########################################

# For SSU:
dat.SSU.kin <- dat %>% filter(SSU_kingdom == "Fungi")
dat.SSU.phy <- dat %>% filter(SSU_phylum == "Mucoromycota")
dat.SSU.fam <- dat %>% filter(SSU_family == LSU_family)
dat.SSU.gen <- dat %>% filter(SSU_genus == LSU_genus)
# check against culture species
dat.SSU.sp <- dat %>% filter(SSU_species == ASV_species)

length(unique(dat.SSU.kin$ASV_updated))/length(unique(dat$ASV_updated)) # 0.9955556
length(unique(dat.SSU.phy$ASV_updated))/length(unique(dat$ASV_updated)) # 0.9955556
length(unique(dat.SSU.fam$ASV_updated))/length(unique(dat$ASV_updated)) # 0.536
length(unique(dat.SSU.gen$ASV_updated))/length(unique(dat$ASV_updated)) # 0.4053333
length(unique(dat.SSU.sp$ASV_updated))/length(unique(dat$ASV_updated)) # 0 # VT so not really comparable.

# per LSU family SSUfamily matches
SSUfam_perfam <- dat %>%
  mutate(all_match = SSU_family == LSU_family) %>%
  group_by(LSU_family) %>%
  summarise(n_total = n(),n_match = sum(all_match, na.rm = TRUE),prop_match = n_match / n_total,.groups = "drop")

# per LSU family SSU genus matches
SSUgen_perfam <- dat %>%
  mutate(all_match = SSU_genus == LSU_genus) %>%
  group_by(LSU_family) %>%
  summarise(n_total = n(),n_match = sum(all_match, na.rm = TRUE),prop_match = n_match / n_total,.groups = "drop")

# For ITS:
dat.ITS.kin <- dat %>% filter(ITS_kingdom == "Fungi")
dat.ITS.phy <- dat %>% filter(ITS_phylum == LSU_phylum)
dat.ITS.fam <- dat %>% filter(ITS_family == LSU_family)
dat.ITS.gen <- dat %>% filter(ITS_genus == LSU_genus)
# check against culture species
dat.ITS.sp <- dat %>% filter(ITS_species == ASV_species)

length(unique(dat.ITS.kin$ASV_updated))/length(unique(dat$ASV_updated)) # 0.9937778
length(unique(dat.ITS.phy$ASV_updated))/length(unique(dat$ASV_updated)) # 0.9937778
length(unique(dat.ITS.fam$ASV_updated))/length(unique(dat$ASV_updated)) # 0.3946667
length(unique(dat.ITS.gen$ASV_updated))/length(unique(dat$ASV_updated)) # 0.7342222
length(unique(dat.ITS.sp$ASV_updated))/length(unique(dat$ASV_updated)) # 0.05422222

# per LSU familyITS family matches
ITSfam_perfam <- dat %>%
  mutate(all_match = ITS_family == LSU_family) %>%
  group_by(LSU_family) %>%
  summarise(n_total = n(),n_match = sum(all_match, na.rm = TRUE),prop_match = n_match / n_total,.groups = "drop")

# per LSU family SSU genus matches
ITSgen_perfam <- dat %>%
  mutate(all_match = ITS_genus == LSU_genus) %>%
  group_by(LSU_family) %>%
  summarise(n_total = n(),n_match = sum(all_match, na.rm = TRUE),prop_match = n_match / n_total,.groups = "drop")

# For Funguild
dat.ITS.fg <- dat %>% filter(ft_primary_lifestyle == "arbuscular_mycorrhizal")
length(unique(dat.ITS.fg$ASV_updated))/length(unique(dat$ASV_updated)) # 0.7662222

##########################################
############# plot results ###############
##########################################

# proportion match SSU ITS

match.for.prop <- match %>% 
  rename(SSU_family = SSU_FAMILY_x_LSU_family, SSU_genus = SSU_genus_x_LSU_genus, ITS_family = ITS_family_x_LSU_family, ITS_genus = ITS_genus_x_LSU_genus, ITS_species = ITS_species_x_LSU_species) %>%
  mutate(across(everything(), ~replace_na(.x, "None")))

prop_match <- match.for.prop %>%
  pivot_longer(
    cols = c(SSU_family, SSU_genus, ITS_family, ITS_genus, ITS_species),
    names_to = "Variable",
    values_to = "Category"
  ) %>%
  filter(Category != "LSU_NA") %>%
  count(Variable, Category) %>%
  group_by(Variable) %>%
  mutate(Proportion = n / sum(n))

plot <- ggplot(prop_match,
       aes(x = Variable,
           y = Proportion,
           fill = Category)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c(
    "Correct" = "#66C2A5",  # soft green
    "Wrong"   = "#FC8D62",  # soft coral/red
    "Outdated"    = "#FFD92F",  # soft yellow
    "None"    = "#D9D9D9"   # light grey
  )) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = "Taxonomic assignment"
  ) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "bottom")

# save fig
png("figures/Taxonomic_assignment_by_region.jpg", width = 16, height = 10, units = 'in', res = 300)
plot
dev.off()

# proportion of each family LSU that is outdated

dat.for.prop <- match.for.prop %>%
  select(ASV_updated, SSU_family, ITS_family) %>%
  pivot_longer(
    cols = c(SSU_family, ITS_family),
    names_to = "Variable",
    values_to = "Category"
  ) %>%
  filter(Category %in% c("Outdated"))

dat.for.prop.LSU <- dat %>% select("ASV_updated", "LSU_family") 

dat.for.prop.all <- dat.for.prop %>% left_join(dat.for.prop.LSU, by = "ASV_updated") %>% filter(!is.na(LSU_family))

prop_match_LSUfam <- dat.for.prop.all %>%
  count(Variable, LSU_family) %>%
  group_by(Variable) %>%
  mutate(Proportion = n / sum(n))

plot <- ggplot(prop_match_LSUfam,
       aes(x = Variable,
           y = Proportion,
           fill = LSU_family)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "Dominikiaceae" = "#5B8E7D",
      "Entrophosporaceae" = "#A7C957",
      "Gigasporaceae" = "#F2CC8F",
      "Sclerocystaceae_Kamienskiaceae" = "#81B29A",
      "Septoglomeraceae" = "#E07A5F"
    ),
    labels = function(x) gsub("_", " ", x)
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = "LSU family"
  ) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none")

# save fig
png("figures/Outdated_LSU_family.jpg", width = 9.5, height = 10, units = 'in', res = 300)
plot
dev.off()

# proportion of each genus LSU that is outdated/wrong

dat.for.prop <- match.for.prop %>%
  select(ASV_updated, SSU_genus, ITS_genus) %>%
  pivot_longer(
    cols = c(SSU_genus, ITS_genus),
    names_to = "Variable",
    values_to = "Category"
  ) %>%
  filter(Category %in% c("Outdated"))

dat.for.prop.LSU <- dat %>% select("ASV_updated", "LSU_genus") 

dat.for.prop.all <- dat.for.prop %>% left_join(dat.for.prop.LSU, by = "ASV_updated") %>% filter(!is.na(LSU_genus))

prop_match_LSUgen <- dat.for.prop.all %>%
  count(Variable, LSU_genus) %>%
  group_by(Variable) %>%
  mutate(Proportion = n / sum(n))

plot <- ggplot(prop_match_LSUgen,
               aes(x = Variable,
                   y = Proportion,
                   fill = LSU_genus)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
  "Entrophospora" = "#6D8F71",
 "Rhizophagus" = "#9CBF88",
 "Septoglomus" = "#C2B280",
 "Cetraspora" = "#DDB892",
 "Dentiscutata" = "#C97C5D",
 "Funneliformis" = "#7FA99B",
 "Oehlia" = "#B7C9A8",
 "Racocetra" = "#A68A64",
 "Sclerocystis" = "#C88C7A"
    ),
    labels = function(x) gsub("_", " ", x)
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = "LSU genus"
  ) +
  theme_minimal(base_size = 30) +
  theme(legend.position = "none")

# save fig
png("figures/Outdated_LSU_genus.jpg", width = 9.5, height = 10, units = 'in', res = 300)
plot
dev.off()
