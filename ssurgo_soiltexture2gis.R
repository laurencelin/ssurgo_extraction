arg=commandArgs(T)


library(rgrass7)

rast = readRAST('soil_ssurgo')
soil_cat_mukey = arg[1] 
soil_mukey_texture = arg[2]

cat_mukey = read.csv(soil_cat_mukey, stringsAsFactors=F)
mukey_texture = read.csv(soil_mukey_texture, stringsAsFactors=F)
	
	# matching SSURGO dataset product to the mapping mukey
	cond = match(cat_mukey[,'MUKEY'],mukey_texture[,'mukey'])
	cat_texture = mukey_texture[cond,'rhessys_soilid']
	
	# sorting matched mukey to the map organization
	cond2 = match(rast@data[[1]], cat_mukey[,'cat'])
	rast$soil_texture = as.integer(cat_texture[cond2]); writeRAST(rast,'soil_texture',zcol='soil_texture',overwrite=T)
	
	## more variable
	cat_texture = mukey_texture[cond,'ksat_0']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_ksat_0',zcol='soil_texture',overwrite=T)
	
	cat_texture = mukey_texture[cond,'ksat_decay']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_ksat_decay',zcol='soil_texture',overwrite=T)
	
	cat_texture = mukey_texture[cond,'por_0']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_por_0',zcol='soil_texture',overwrite=T)
	
	cat_texture = mukey_texture[cond,'por_decay']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_por_decay',zcol='soil_texture',overwrite=T)
	
	cat_texture = mukey_texture[cond,'sand']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_sand',zcol='soil_texture',overwrite=T)

	cat_texture = mukey_texture[cond,'silt']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_silt',zcol='soil_texture',overwrite=T)
	
	cat_texture = mukey_texture[cond,'clay']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_clay',zcol='soil_texture',overwrite=T)
	
	cat_texture = mukey_texture[cond,'BD']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_BD',zcol='soil_texture',overwrite=T)

    cat_texture = mukey_texture[cond,'partdensity']
    rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_horizonDepth',zcol='soil_texture',overwrite=T)

	cat_texture = mukey_texture[cond,'soilDepth']
	rast$soil_texture = cat_texture[cond2]; writeRAST(rast,'soil_horizonDepth',zcol='soil_texture',overwrite=T)


	
