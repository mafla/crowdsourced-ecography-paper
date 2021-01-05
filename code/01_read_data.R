# shape file of the MTBs
print("reading shape file with ordnance survey maps")
MTB = readOGR("data/mtbschnitt.shp")

# shape file of Germany
print("reading shape file of Germany")
germany_0 = readOGR("data/gadm36_DEU_0.shp")

# read German census data 
print("reading census data")
zensus = read_csv("data/mtbschnitt_zensus2011.csv")

# read florkart; in the following fk and  flora inkognita; in the following fi
print("reading Florkart and Flora Incognita and filtering to common species and grid cells")
tibble_fk  = read_csv("data/florkart_table.csv")
tibble_fi  = read_csv("data/fia_table.csv")
tibble_fk  = rename(tibble_fk, NAME = Name)

# find common MTBs (MTB = Messtischblatt, traditional German geogrid)
idx_mtb_fk = unique(tibble_fk$NAME)
idx_mtb_fi = unique(tibble_fi$NAME)
idx_mtb    = sort(intersect(idx_mtb_fk, idx_mtb_fi))

# retain common locations only
tibble_fk  = filter(tibble_fk, idx_mtb_fk %in% idx_mtb) 
tibble_fi  = filter(tibble_fi, idx_mtb_fi %in% idx_mtb)

# sort locations to be able to compare
tibble_fk  = arrange(tibble_fk, NAME)
tibble_fi  = arrange(tibble_fi, NAME)

# no we can drop the column with MTBs - we need clean matrices
tibble_fk  = dplyr::select(tibble_fk, -NAME)
tibble_fi  = dplyr::select(tibble_fi, -NAME)

# find common species
idx_spc_fk = unique(colnames(tibble_fk))
idx_spc_fi = unique(colnames(tibble_fi))
idx_spc    = intersect(idx_spc_fk, idx_spc_fi)

# retain common species only
tibble_fk  = dplyr::select(tibble_fk, idx_spc)
tibble_fi  = dplyr::select(tibble_fi, idx_spc)

# convert to presence/absence
tibble_fk[tibble_fk > 0] = 1
tibble_fi[tibble_fi > 0] = 1