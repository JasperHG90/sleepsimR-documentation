## Legacy code

# Unique combinations 
e<- d %>%
  select(emiss_var_short,
         zeta) %>%
  mutate(expandvar = paste0(emiss_var_short, "_", zeta)) %>%
  select(expandvar) %>%
  distinct() %>%
  pull() %>% 
  expand.grid(., .) %>%
  filter(Var1 != Var2) %>%
  apply(1, function(x) sort(x)) %>%
  t() %>%
  data.frame(., stringsAsFactors = FALSE) %>%
  distinct() %>%
  separate(X1, into=c("type", "stat", "channel", "zeta_1"), sep = "_") %>%
  mutate(emiss_var_short_1 = paste0(type, "_", stat, "_", channel)) %>%
  select(-stat, -channel, -type) %>%
  separate(X2, into=c("type", "stat", "channel", "zeta_2"), sep = "_") %>%
  mutate(emiss_var_short_2 = paste0(type, "_", stat, "_", channel)) %>%
  select(-stat, -channel, -type)
