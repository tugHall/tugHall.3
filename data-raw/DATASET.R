## code to prepare `tugHall_dataset` dataset goes here

tugHall_dataset  =  simulation_example( verbose = FALSE, to_plot = FALSE )

usethis::use_data( tugHall_dataset, overwrite = TRUE )

devtools::install( build_vignettes = TRUE )


