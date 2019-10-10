arg=commandArgs(T)

library(rgrass7)

rast = readRAST('soil_ssurgo')
soil_rhessys = read.csv(arg[1], stringsAsFactors=F)

	
	# sorting matched mukey to the map organization
	cond2 = match(rast@data[[1]], soil_rhessys[,'cat'])
	
    # texture class
    rast$soil_texture = as.integer(soil_rhessys[cond2,'texture']); writeRAST(rast,'soil_texture',zcol='soil_texture',overwrite=T)
	
	# other variables
	rast$soil_texture = soil_rhessys[cond2,'ksat0']; writeRAST(rast,'soil_ksat_0',zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'ksatdecay']; writeRAST(rast,'soil_ksat_decay',zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'por0']; writeRAST(rast,'soil_por_0',zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'pordecay']; writeRAST(rast,'soil_por_decay',zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'sand']; writeRAST(rast,'soil_sand',zcol='soil_texture',overwrite=T)

	rast$soil_texture = soil_rhessys[cond2,'silt']; writeRAST(rast,'soil_silt',zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'clay']; writeRAST(rast,'soil_clay',zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'bulkdensity']; writeRAST(rast,'soil_bulkdensity',zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'particledensity']; writeRAST(rast,'soil_particledensity',zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'om']; writeRAST(rast,'soil_soilc',zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'omdecay']; writeRAST(rast,'soil_omdecay',zcol='soil_texture',overwrite=T)

	rast$soil_texture = soil_rhessys[cond2,'soildepth']; writeRAST(rast,'soil_soildepth',zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'activedepth']; writeRAST(rast,'soil_activedepth',zcol='soil_texture',overwrite=T)
    
    rast$soil_texture = soil_rhessys[cond2,'maxrootdepth']; writeRAST(rast,'soil_maxrootdepth',zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'albedo']; writeRAST(rast,'soil_albedo',zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'por_size_index']; writeRAST(rast,'soil_por_size_index',zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'psi_air_entry']; writeRAST(rast,'soil_psi_air_entry',zcol='soil_texture',overwrite=T)

	
