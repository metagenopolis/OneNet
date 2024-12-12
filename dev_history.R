#initial setup
usethis::use_mit_license("Camille Champion")
usethis::use_build_ignore("dev_history.R")
attachment::att_amend_desc()
usethis::use_package_doc()
usethis::use_readme_rmd()
devtools::build_manual()
#vignettes
usethis::use_vignette("Demo_offline", title = "Demo_offline")
#Rcpp
usethis::use_rcpp()
Rcpp::compileAttributes()
roxygen2::roxygenize(roclets="rd")
#imports
usethis::use_package("EMtree")
usethis::use_package("SPRING")
usethis::use_package("rMAGMA")
usethis::use_package("SpiecEasi")
usethis::use_package("PLNmodels")
usethis::use_package("RColorBrewer")
usethis::use_package("Rfit")
usethis::use_package("ggm")
usethis::use_package("dplyr")
usethis::use_package("ggplot2")
usethis::use_package("grDevices")
usethis::use_package("huge")
usethis::use_package("magrittr")
usethis::use_package("mixedCCA")
usethis::use_package("parallel")
usethis::use_package("pulsar")
usethis::use_package("RColorBrewer")
usethis::use_package("reshape2")
usethis::use_package("tibble")
usethis::use_package("blockmodels")
usethis::use_package("igraph")
# data set
usethis::use_data_raw("liver")
#worflow
devtools::document()
devtools::load_all()
devtools::run_examples()


# VIGNETTES
# build vignetes with RStudio:
# dans Build>configure build tools > configure > check vignette
# warning : check deletes inst/doc !
# to move vignette to ins/doc
dir.create("inst/doc")
file.copy(dir("vignettes", full.names=TRUE), "inst/doc", overwrite=TRUE)
