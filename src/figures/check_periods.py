from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

table = NasaExoplanetArchive.query_object("Kepler-8 b")
print(table['pl_orbper'])

table2 = NasaExoplanetArchive.query_object("Kepler-7 b")
print(table2['pl_orbper'])