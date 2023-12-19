## code to prepare `basicRelationships` dataset goes here

library(tidyverse)

data = list(
  unrelated                = list(label = "UN", kappa = c(1, 0, 0), pos = 1),
  `parent-offspring`       = list(label = "PO", kappa = c(0, 1, 0), pos = 1),
  `monozygotic twins`      = list(label = "MZ", kappa = c(0, 0, 1), pos = 4),
  `full siblings`          = list(label = "S",  kappa = c(.25, 0.5, .25), pos = 4),
  `half siblings`          = list(label = "H",  kappa = c(.5, .5, 0), pos = 1),
  `avuncular`              = list(label = "A",  kappa = c(.5, .5, 0), pos = 1),
  `grandparent-grandchild` = list(label = "G",  kappa = c(.5, .5, 0), pos = 1),
  `halfsib/uncle/grandp`   = list(label = "H,U,G", kappa = c(.5, .5, 0), pos = 1),
  `first cousins`          = list(label = "FC", kappa = c(.75, .25, 0), pos = 1),
  `second cousins`         = list(label = "SC", kappa = c(15/16, 1/16, 0), pos = 1),
  `double first cousins`   = list(label = "DFC", kappa = c(9/16, 6/16, 1/16), pos = 3),
  `quad half first cousins`= list(label = "Q",  kappa = c(17/32, 14/32, 1/32), pos = 4))

relsTibble = enframe(data) |>
  unnest_wider(value) |>
  unnest_wider(kappa, names_sep = "") |>
  rename(relationship = name, kappa0 = kappa1, kappa1 = kappa2, kappa2 = kappa3) |>
  mutate(phi = kappa2/2 + kappa1/4, .after = label) |>
  print()

basicRelationships = as.data.frame(relsTibble)

usethis::use_data(basicRelationships, overwrite = TRUE)
