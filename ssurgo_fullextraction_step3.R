arg=commandArgs(T)

library(sp)
library(XML)
library(rgrass7)
library(rgdal)
tryCatch({ use_sp() },error=function(cond){message(cond)},warning=function(cond){message(cond)},finally={message("Please update the rgrass7 package on R")})

rast = readRAST(arg[1])
soil_rhessys = read.csv(arg[2], stringsAsFactors=F)

	
	# sorting matched mukey to the map organization
	cond2 = match(rast@data[[1]], soil_rhessys[,'cat'])
	if(length(arg)>2){ suffix = arg[3]; } else { suffix=''; }
    # texture class
    rast$soil_texture = as.integer(soil_rhessys[cond2,'texture']); writeRAST(rast,paste('soil_texture',suffix,sep=''),zcol='soil_texture',overwrite=T)
	
	# other variables
	rast$soil_texture = soil_rhessys[cond2,'ksat0']; writeRAST(rast,paste('soil_ksat_0',suffix,sep=''),zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'ksatdecay']; writeRAST(rast,paste('soil_ksat_decay',suffix,sep=''),zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'por0']; writeRAST(rast,paste('soil_por_0',suffix,sep=''),zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'pordecay']; writeRAST(rast,paste('soil_por_decay',suffix,sep=''),zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'sand']; writeRAST(rast,paste('soil_sand',suffix,sep=''),zcol='soil_texture',overwrite=T)

	rast$soil_texture = soil_rhessys[cond2,'silt']; writeRAST(rast,paste('soil_silt',suffix,sep=''),zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'clay']; writeRAST(rast,paste('soil_clay',suffix,sep=''),zcol='soil_texture',overwrite=T)
	
	rast$soil_texture = soil_rhessys[cond2,'bulkdensity']; writeRAST(rast,paste('soil_bulkdensity',suffix,sep=''),zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'particledensity']; writeRAST(rast,paste('soil_particledensity',suffix,sep=''),zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'om']; writeRAST(rast,paste('soil_soilc',suffix,sep=''),zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'omdecay']; writeRAST(rast,paste('soil_omdecay',suffix,sep=''),zcol='soil_texture',overwrite=T)

	rast$soil_texture = soil_rhessys[cond2,'soildepth']; writeRAST(rast,paste('soil_soildepth',suffix,sep=''),zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'activedepth']; writeRAST(rast,paste('soil_activedepth',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    rast$soil_texture = soil_rhessys[cond2,'maxrootdepth']; writeRAST(rast,paste('soil_maxrootdepth',suffix,sep=''),zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'albedo']; writeRAST(rast,paste('soil_albedo',suffix,sep=''),zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'por_size_index']; writeRAST(rast,paste('soil_por_size_index',suffix,sep=''),zcol='soil_texture',overwrite=T)

    rast$soil_texture = soil_rhessys[cond2,'psi_air_entry']; writeRAST(rast,paste('soil_psi_air_entry',suffix,sep=''),zcol='soil_texture',overwrite=T)

	
