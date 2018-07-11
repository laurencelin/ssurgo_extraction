arg=commandArgs(T)

library(rgrass7)

rast = readRAST('soil_ssurgo')


soil_cat_mukey = arg[1] # product from "v.db.select map=ssurgo separator=comma file=$PROJDIR/rhessys/ssurgo_cat_mukey.csv"
soil_mukey_texture = arg[2] # product from "Rscript ssurgo_extraction.R $downloadedSSURGO_directory"

cat_mukey = read.csv(soil_cat_mukey, stringsAsFactors=F)
mukey_texture = read.csv(soil_mukey_texture, stringsAsFactors=F)
	
	cat_texture = mukey_texture[match(cat_mukey[,'MUKEY'],mukey_texture[,'mukey']),'texture']
	rast$soil_texture = cat_texture[match(rast@data[[1]], cat_mukey[,'cat'])]
	
	writeRAST(rast,'soil_texture',zcol='soil_texture')
	
