arg=commandArgs(T)

library(sp)
library(XML)
library(rgrass7)
library(rgdal)
tryCatch({ use_sp() },error=function(cond){message(cond)},warning=function(cond){message(cond)},finally={message("Please update the rgrass7 package on R")})

soil_rhessys = read.csv(arg[2], stringsAsFactors=F)
if(length(arg)>2){ suffix = arg[3]; } else { suffix=''; }

if(is.null(soil_rhessys$map)){
    rast = readRAST(arg[1])
    
    # sorting matched mukey to the map organization
    cond2 = match(rast@data[[1]], soil_rhessys[,'cat'])
    
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
    
}else{
    ssurgoSuffix = unique(soil_rhessys$map)
    rast = readRAST(paste(arg[1],ssurgoSuffix,sep=''))
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = as.integer(soil_rhessys[cond2[cond3],'texture']);
    };  writeRAST(rast,paste('soil_texture',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'ksat0'];
    };  writeRAST(rast,paste('soil_ksat_0',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'por0'];
    };  writeRAST(rast,paste('soil_por_0',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'pordecay'];
    };  writeRAST(rast,paste('soil_por_decay',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'sand'];
    };  writeRAST(rast,paste('soil_sand',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'silt'];
    };  writeRAST(rast,paste('soil_silt',suffix,sep=''),zcol='soil_texture',overwrite=T)
 
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'clay'];
    };  writeRAST(rast,paste('soil_clay',suffix,sep=''),zcol='soil_texture',overwrite=T)
 
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'bulkdensity'];
    };  writeRAST(rast,paste('soil_bulkdensity',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'particledensity'];
    };  writeRAST(rast,paste('soil_particledensity',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'om'];
    };  writeRAST(rast,paste('soil_soilc',suffix,sep=''),zcol='soil_texture',overwrite=T)
 
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'omdecay'];
    };  writeRAST(rast,paste('soil_omdecay',suffix,sep=''),zcol='soil_texture',overwrite=T)
 
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'soildepth'];
    };  writeRAST(rast,paste('soil_soildepth',suffix,sep=''),zcol='soil_texture',overwrite=T)
 
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'activedepth'];
    };  writeRAST(rast,paste('soil_activedepth',suffix,sep=''),zcol='soil_texture',overwrite=T)
 
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'maxrootdepth'];
    };  writeRAST(rast,paste('soil_maxrootdepth',suffix,sep=''),zcol='soil_texture',overwrite=T)
 
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'albedo'];
    };  writeRAST(rast,paste('soil_albedo',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'por_size_index'];
    };  writeRAST(rast,paste('soil_por_size_index',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
    for(j in seq_along(ssurgoSuffix)){
        cond2 = match(rast@data[[j]], soil_rhessys[soil_rhessys$map==ssurgoSuffix[j],'cat'])
        cond3 = !is.na(cond2); rast$soil_texture[cond3] = soil_rhessys[cond2[cond3],'psi_air_entry'];
    };  writeRAST(rast,paste('soil_psi_air_entry',suffix,sep=''),zcol='soil_texture',overwrite=T)
    
}# end of if else


	
	

	
