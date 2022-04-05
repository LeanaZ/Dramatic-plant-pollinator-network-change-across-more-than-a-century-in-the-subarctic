# Overview

This repository contains the code and data needed to reproduce the results and develop the figures from:

**Zoller, L., Bennett, J.M., Knight, T.M.**; *Dramatic plant-pollinator network change across more than a century in the subarctic*

[Download][1] or clone the repository then run the scripts using the `A_century_of_plant_pollinator_interactions.Rproj` file ([R][2] and [R Studio][3] are needed)

[1]: https://github.com/LeanaZ/Dramatic-plant-pollinator-network-change-across-more-than-a-century-in-the-subarctic/archive/master.zip
[2]: https://www.r-project.org/
[3]: https://www.rstudio.com/products/rstudio/download/


# Data analysis

Code for data analysis and producing figures is in `analysis/main_analysis_Zoller_et_al.R`. Reads data from `data/InteractionData_Zoller_et_al.csv`.

# R Session information
Information about the current R Session at the time of running the analysis:

> sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 16299)

Matrix products: default

locale:
[1] LC_COLLATE=English_Germany.1252  LC_CTYPE=English_Germany.1252    LC_MONETARY=English_Germany.1252 LC_NUMERIC=C                    
[5] LC_TIME=English_Germany.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
 [1] sna_2.6              pillar_1.6.0         compiler_4.0.3       iterators_1.0.13     tools_4.0.3          dotCall64_1.0-1     
 [7] rle_0.9.2            lifecycle_1.0.0      tibble_3.1.1         gtable_0.3.0         nlme_3.1-152         lattice_0.20-44     
[13] mgcv_1.8-35          pkgconfig_2.0.3      rlang_0.4.11         foreach_1.5.1        Matrix_1.2-18        igraph_1.2.6        
[19] rstudioapi_0.13      cli_2.5.0            DBI_1.1.1            yaml_2.2.1           parallel_4.0.3       spam_2.6-0          
[25] coda_0.19-4          dplyr_1.0.6          cluster_2.1.2        maps_3.3.0           fields_11.6          generics_0.1.0      
[31] vctrs_0.3.8          grid_4.0.3           cowplot_1.1.1        tidyselect_1.1.1     glue_1.4.2           data.table_1.14.0   
[37] bipartite_2.16       R6_2.5.0             fansi_0.4.2          tidyr_1.1.3          ggplot2_3.3.5        purrr_0.3.4         
[43] magrittr_2.0.1       matrixStats_0.58.0   codetools_0.2-18     scales_1.1.1         ellipsis_0.3.2       MASS_7.3-54         
[49] splines_4.0.3        assertthat_0.2.1     bootstrapnet_1.0.0   permute_0.9-5        colorspace_2.0-1     utf8_1.2.1          
[55] network_1.16.1       doParallel_1.0.16    munsell_0.5.0        statnet.common_4.4.1 crayon_1.4.1         vegan_2.5-7         
