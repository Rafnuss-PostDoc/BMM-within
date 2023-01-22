library(bioRad)
library(rjson)


flr <- "../data/error_high_elev_target/"

dd_name <- fromJSON(file = "../test_error_high_elev/dd_name.json")

for (i_r in seq_len(length(dd_name))){
  par(mfrow=c(4,2))
  for (i_d in seq_len(length(dd_name[[i_r]]))){

    vpfile <- calculate_vp(file = paste0(flr,substr(dd_name[[i_r]][[i_d]],17,10000)),
                           range_max=150000,
                           h_layer = 100,
                           n_layer = 3000/100)
    plot(vpfile)

    pvolfile <- read_pvolfile(file = paste0(flr,substr(dd_name[[i_r]][[i_d]],17,10000)))
    #ppi <- project_as_ppi(pvolfile$scans[[1]], range_max = 150000)
    #plot(ppi)

    corrected_ppi <- integrate_to_ppi(pvolfile, vpfile, res = 2000, xlim = c(-150000, 150000), ylim = c(-150000, 150000))
    plot(corrected_ppi)
  }
}


my_scan <- example_pvol$scans[[1]]
# project it as a PPI on the ground:
my_ppi <- project_as_ppi(my_scan, range_max = 100000)

plot(my_ppi)
