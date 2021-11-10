# chamando as libs necessárias

pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, RAM, permute, dplyr, stringr)

df <- read_csv('analise.csv') %>% 
  select(-Gênero, - Tribo) %>%
  filter(abundância > 0 )

top_species <- df %>%
  group_by(Espécies) %>%
  summarise(abundância = sum(abundância)) %>%
  arrange(-abundância) %>% filter(abundância >= 12)

# Por Coleta
top_specie <- df %>%
  group_by(Espécies, Coleta) %>%
  summarise(abundância = sum(abundância)) %>%
  arrange(-abundância) %>% filter(Espécies %in% top_species$Espécies)

ggplot(top_specie, 
       aes(x = Coleta, y = abundância, color = Espécies, group = Espécies)) + 
  geom_point() + geom_line() + theme_light() + xlab("Distância (m)") +
  ylab("Abundância")

# Trilha
top_specie <- df %>%
  group_by(Espécies, Trilha) %>%
  summarise(abundância = sum(abundância)) %>%
  arrange(-abundância) %>% filter(Espécies %in% top_species$Espécies)

ggplot(top_specie, 
       aes(x = Trilha, y = abundância, color = Espécies, group = Espécies)) + 
  geom_point() + geom_line() + theme_light() + ylab("Abundância")

# NMDS

rich <- df %>%
  group_by(Trilha, Coleta) %>%
  summarise(rich = n_distinct(Espécies))

data_nmds <- df %>%
  group_by(Trilha, Coleta, Espécies) %>%
  summarise(abundância = sum(abundância)) 


data_nmds %<>%
  pivot_wider(names_from = Espécies, values_from = abundância) %>%
  mutate(
    across(everything(), ~ replace_na(.x, 0))
  ) %>%
  left_join(rich, by = c("Trilha", "Coleta"))


run_nmds <- data_nmds

run_nmds$Trilha <- NULL
run_nmds$Coleta <- NULL
run_nmds$rich <- NULL

dist_bray <- vegdist(run_nmds, method = "bray", binary = TRUE)

nmds <- metaMDS(dist_bray)

scores(nmds) %>%
  as_tibble() %>%
  cbind(data_nmds) %>%
  as_tibble() %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(size = Trilha, color = Coleta)) +
  stat_ellipse(geom = "polygon", aes(group = rich), alpha = 0.3) +
  annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw()+
  scale_color_discrete(breaks=c("0m","100m","300m", "700m", "1700m"))+
  scale_color_discrete(name = "Distância (m)")

# Permanova

adonis(dist_bray~data_nmds$Trilha, permutations = 1000)


# Riqueza por Trilha

df %>%
  group_by(Trilha) %>%
  summarise(rich = n_distinct(Espécies)) %>%
  ggplot(aes(x=Trilha, y=rich)) +
  geom_bar(stat="identity", fill="lightgray") +
  geom_text(aes(label=rich), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_light() +
  ylab("Riqueza")

# Riqueza por Coleta

df %>%
  group_by(Coleta) %>%
  summarise(rich = n_distinct(Espécies)) %>%
  arrange(-rich) %>%
  ggplot(aes(x=reorder(Coleta,+rich), y=rich)) +
  geom_bar(stat="identity", fill="lightgray") +
  geom_text(aes(label=rich), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_light() +
  ylab("Riqueza") +
  xlab("Distância (m)")

# Diversidade 

abund <- df %>% group_by(Coleta, Espécies) %>%
  summarise(abundância = sum(abundância)) %>%
  pivot_wider(names_from = Coleta, values_from = abundância) %>%
  mutate(
    across(everything(), ~ replace_na(.x, 0))
  ) %>%
  column_to_rownames(var = "Espécies") 

abund %<>% select('0m', "100m", "300m", '700m', "1700m")

resultados <- iNEXT(abund,
  q = 2,
  datatype = "abundance",
  endpoint = 57
)

ggiNEXT(resultados, type = 1) + theme_light() + 
  scale_color_discrete(labels = paste(c('0m', "100m", "300m", '700m', "1700m")))+
  scale_fill_discrete(labels = paste(c('0m', "100m", "300m", '700m', "1700m"))) +
  ylab("Diversidade de espécies")+ 
  xlab("Número de indivíduos")+
  guides(fill = guide_legend(title="Guia"),
         color = guide_legend(title="Guia"),
         linetype = guide_legend(title="Método"),
         shape = guide_legend(title="Guia"))
