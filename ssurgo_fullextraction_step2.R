arg=commandArgs(T)

replaceNA0 = function(x,value){
	x[is.na(x) | x<=0] = value;
	return <- x;	
}#function
weightedAverage = function(xx,weight){
    new_weight = weight;
    new_weight[is.na(xx)] = 0;
    if( sum(new_weight)>0 ){
        return <- sum(xx * new_weight/sum(new_weight), na.rm=T)
    }else{
        return <- NA
    }
}#function
horizonOrder = c('O','A','E','B','C','R')
htitle = c('horizonThinkness','horizonMeanDepth','horizonTopDepth','horizonBottomDepth','sand','silt','clay','ksat','POR','BD','om','density','por_size_index','psi_air_entry','albedo')
horizonThicknessDefault = c(2, 5, 5, 10, 10, 10)*0.0254 #inches to meter; https://www.sheffield.ac.uk/ssa/soil-facts/horizons
soildepth_rtz_ratio = sum(horizonThicknessDefault)/sum(horizonThicknessDefault[1:5])
soildepth_act_ratio = sum(horizonThicknessDefault)/sum(horizonThicknessDefault[1:4])



horizons = read.csv(paste(arg[1],'/soil_mukey_texture.csv',sep=''),stringsAsFactors=F)
cat_mukey = read.csv(arg[2], stringsAsFactors=F, header=F)# use command "v.out.ascii"
colnames(cat_mukey) = c('x','y','cat','mukey')
LOCATION_cat_mukey = match(cat_mukey[,'mukey'], horizons[,'mukey']) ## cat(LOCATION) <-- mukey


what = data.frame(
	catid = cat_mukey[,'cat'],
	mukey = horizons[LOCATION_cat_mukey,'mukey'],
	# -- full	
	soil_texture = horizons[LOCATION_cat_mukey,'texture'], 
	soildepth = horizons[LOCATION_cat_mukey,'R_horizonBottomDepth'],
	soil_a_meanz = horizons[LOCATION_cat_mukey,'A_horizonMeanDepth'],
	soil_a_BD = horizons[LOCATION_cat_mukey,'A_BD'],
	soil_a_PD = horizons[LOCATION_cat_mukey,'A_density'],
	soil_a_POR = horizons[LOCATION_cat_mukey,'A_POR'],
	soil_a_ksat = horizons[LOCATION_cat_mukey,'A_ksat'],
	soil_a_por_size_index = horizons[LOCATION_cat_mukey,'A_por_size_index'],
	soil_a_psi_air_entry = horizons[LOCATION_cat_mukey,'A_psi_air_entry'],
	soil_b_meanz = horizons[LOCATION_cat_mukey,'B_horizonMeanDepth'],
	soil_b_BD = horizons[LOCATION_cat_mukey,'B_BD'],
	soil_b_PD = horizons[LOCATION_cat_mukey,'B_density'],
	soil_b_POR = horizons[LOCATION_cat_mukey,'B_POR'],
	soil_b_ksat = horizons[LOCATION_cat_mukey,'B_ksat'],
	soil_b_por_size_index = horizons[LOCATION_cat_mukey,'B_por_size_index'],
	soil_b_psi_air_entry = horizons[LOCATION_cat_mukey,'B_psi_air_entry'],
	soil_c_meanz = horizons[LOCATION_cat_mukey,'C_horizonMeanDepth'],
	soil_c_BD = horizons[LOCATION_cat_mukey,'C_BD'],
	soil_c_PD = horizons[LOCATION_cat_mukey,'C_density'],
	soil_c_POR = horizons[LOCATION_cat_mukey,'C_POR'],
	soil_c_ksat = horizons[LOCATION_cat_mukey,'C_ksat'],
	soil_c_por_size_index = horizons[LOCATION_cat_mukey,'C_por_size_index'],
	soil_c_psi_air_entry = horizons[LOCATION_cat_mukey,'C_psi_air_entry'],
	soil_c_thickness = horizons[LOCATION_cat_mukey,'C_horizonThinkness'],
	#soil_BD
	#soil_PD
	# soil_ksat0
	# soil_ksatdecay
	# soil_POR0
	# soil_PORdecay
	# -- upper layers
	soil_maxrootdepth = horizons[LOCATION_cat_mukey,'C_horizonBottomDepth'],
	soil_activedepth = horizons[LOCATION_cat_mukey,'B_horizonBottomDepth'],
	soil_a_albedo = horizons[LOCATION_cat_mukey,'A_albedo'],
	soil_a_thickness = horizons[LOCATION_cat_mukey,'A_horizonThinkness'],
	soil_a_om = horizons[LOCATION_cat_mukey,'A_om'],
	soil_a_sand = horizons[LOCATION_cat_mukey,'A_sand'],
	soil_a_clay = horizons[LOCATION_cat_mukey,'A_clay'],
	soil_a_silt = horizons[LOCATION_cat_mukey,'A_silt'],
	soil_b_thickness = horizons[LOCATION_cat_mukey,'B_horizonThinkness'],
	soil_b_om = horizons[LOCATION_cat_mukey,'B_om'],
	soil_b_sand = horizons[LOCATION_cat_mukey,'B_sand'],
	soil_b_clay = horizons[LOCATION_cat_mukey,'B_clay'],
	soil_b_silt = horizons[LOCATION_cat_mukey,'B_silt'],
    soil_c_om = horizons[LOCATION_cat_mukey,'C_om'],
	soil_r_thickness = horizons[LOCATION_cat_mukey,'R_horizonThinkness'],
	soil_r_active_ratio = (horizons[LOCATION_cat_mukey,'A_horizonThinkness']+horizons[LOCATION_cat_mukey,'B_horizonThinkness']+horizons[LOCATION_cat_mukey,'C_horizonThinkness']+horizons[LOCATION_cat_mukey,'R_horizonThinkness'])/(horizons[LOCATION_cat_mukey,'A_horizonThinkness']+horizons[LOCATION_cat_mukey,'B_horizonThinkness']),
    soil_r_maxrtz_ratio = (horizons[LOCATION_cat_mukey,'A_horizonThinkness']+horizons[LOCATION_cat_mukey,'B_horizonThinkness']+horizons[LOCATION_cat_mukey,'C_horizonThinkness']+horizons[LOCATION_cat_mukey,'R_horizonThinkness'])/(horizons[LOCATION_cat_mukey,'A_horizonThinkness']+horizons[LOCATION_cat_mukey,'B_horizonThinkness']+horizons[LOCATION_cat_mukey,'C_horizonThinkness'])
)# data.frame
cond = !is.na(what$soildepth) & !is.na(what$soil_maxrootdepth) & (what$soildepth < what$soil_maxrootdepth);sum(cond)
cond = !is.na(what$soildepth) & !is.na(what$soil_activedepth) & (what$soildepth < what$soil_activedepth);sum(cond)
cond = !is.na(what$soil_maxrootdepth) & !is.na(what$soil_activedepth) & (what$soil_maxrootdepth < what$soil_activedepth);sum(cond)
cond = !is.na(what$soildepth); sum(cond); dim(what)[1]
# what[cond,c('mukey','soildepth','soil_maxrootdepth','soil_activedepth')][1:10,]

what$check = sapply(seq_len(dim(what)[1]),function(ii){ sum(is.na(what[ii,]))>0 })
what_LEN = seq_len(dim(what)[1])
counting = 1
constrainedDist=1000
while(sum(what$check)>0){	
	finalTable = what
	print(paste('counting', counting, sum(what$check) ))
	for(ii in what_LEN[what$check]){
		dist = sqrt((cat_mukey[ii,'x']-cat_mukey[,'x'])^2 + (cat_mukey[ii,'y']-cat_mukey[,'y'])^2)
		weights = exp(-0.04*dist) * ifelse(dist>0 & dist<constrainedDist,1,0)
		cond = weights>0
        
        if(is.na(what$soil_texture[ii])){
            texturecount = tapply(weights[cond],what$soil_texture[cond],sum)
            finalTable $soil_texture[ii] = as.numeric(names(texturecount))[which.max(texturecount)]
        }#if
            
		if(is.na(what$soildepth[ii])) finalTable $soildepth[ii] = weightedAverage(what$soildepth[cond], weights[cond])
		if(is.na(what$soil_r_thickness[ii])) finalTable $soil_r_thickness[ii] = weightedAverage(what$soil_r_thickness[cond], weights[cond])
		if(is.na(what$soil_r_active_ratio[ii])) finalTable $soil_r_active_ratio[ii] = weightedAverage(what$soil_r_active_ratio[cond], weights[cond])
        if(is.na(what$soil_r_maxrtz_ratio[ii])) finalTable $soil_r_maxrtz_ratio[ii] = weightedAverage(what$soil_r_maxrtz_ratio[cond], weights[cond])
		
		if(is.na(what$soil_a_meanz[ii])) finalTable $soil_a_meanz[ii] = weightedAverage(what$soil_a_meanz[cond], weights[cond])
		if(is.na(what$soil_a_BD[ii])) finalTable $soil_a_BD[ii] = weightedAverage(what$soil_a_BD[cond], weights[cond])
		if(is.na(what$soil_a_PD[ii])) finalTable $soil_a_PD[ii] = weightedAverage(what$soil_a_PD[cond], weights[cond])
		if(is.na(what$soil_a_POR[ii])) finalTable $soil_a_POR[ii] = weightedAverage(what$soil_a_POR[cond], weights[cond])
		if(is.na(what$soil_a_ksat[ii])) finalTable $soil_a_ksat[ii] = weightedAverage(what$soil_a_ksat[cond], weights[cond])
		if(is.na(what$soil_a_por_size_index[ii])) finalTable $soil_a_por_size_index[ii] = weightedAverage(what$soil_a_por_size_index[cond], weights[cond])
		if(is.na(what$soil_a_psi_air_entry[ii])) finalTable $soil_a_psi_air_entry[ii] = weightedAverage(what$soil_a_psi_air_entry[cond], weights[cond])
		
		if(is.na(what$soil_b_meanz[ii])) finalTable $soil_b_meanz[ii] = weightedAverage(what$soil_b_meanz[cond], weights[cond])
		if(is.na(what$soil_b_BD[ii])) finalTable $soil_b_BD[ii] = weightedAverage(what$soil_b_BD[cond], weights[cond])
		if(is.na(what$soil_b_PD[ii])) finalTable $soil_b_PD[ii] = weightedAverage(what$soil_b_PD[cond], weights[cond])
		if(is.na(what$soil_b_POR[ii])) finalTable $soil_b_POR[ii] = weightedAverage(what$soil_b_POR[cond], weights[cond])
		if(is.na(what$soil_b_ksat[ii])) finalTable $soil_b_ksat[ii] = weightedAverage(what$soil_b_ksat[cond], weights[cond])
		if(is.na(what$soil_b_por_size_index[ii])) finalTable $soil_b_por_size_index[ii] = weightedAverage(what$soil_b_por_size_index[cond], weights[cond])
		if(is.na(what$soil_b_psi_air_entry[ii])) finalTable $soil_b_psi_air_entry[ii] = weightedAverage(what$soil_b_psi_air_entry[cond], weights[cond])
		
		if(is.na(what$soil_c_meanz[ii])) finalTable $soil_c_meanz[ii] = weightedAverage(what$soil_c_meanz[cond], weights[cond])
		if(is.na(what$soil_c_BD[ii])) finalTable $soil_c_BD[ii] = weightedAverage(what$soil_c_BD[cond], weights[cond])
		if(is.na(what$soil_c_PD[ii])) finalTable $soil_c_PD[ii] = weightedAverage(what$soil_c_PD[cond], weights[cond])
		if(is.na(what$soil_c_POR[ii])) finalTable $soil_c_POR[ii] = weightedAverage(what$soil_c_POR[cond], weights[cond])
		if(is.na(what$soil_c_ksat[ii])) finalTable $soil_c_ksat[ii] = weightedAverage(what$soil_c_ksat[cond], weights[cond])
		if(is.na(what$soil_c_por_size_index[ii])) finalTable $soil_c_por_size_index[ii] = weightedAverage(what$soil_c_por_size_index[cond], weights[cond])
		if(is.na(what$soil_c_psi_air_entry[ii])) finalTable $soil_c_psi_air_entry[ii] = weightedAverage(what$soil_c_psi_air_entry[cond], weights[cond])
		if(is.na(what$soil_c_thickness[ii])) finalTable $soil_c_thickness[ii] = weightedAverage(what$soil_c_thickness[cond], weights[cond])
		
		if(is.na(what$soil_maxrootdepth[ii])) finalTable $soil_maxrootdepth[ii] = weightedAverage(what$soil_maxrootdepth[cond], weights[cond])
		if(is.na(what$soil_activedepth[ii])) finalTable $soil_activedepth[ii] = weightedAverage(what$soil_activedepth[cond], weights[cond])
		
		if(is.na(what$soil_a_albedo[ii])) finalTable $soil_a_albedo[ii] = weightedAverage(what$soil_a_albedo[cond], weights[cond])
		if(is.na(what$soil_a_thickness[ii])) finalTable $soil_a_thickness[ii] = weightedAverage(what$soil_a_thickness[cond], weights[cond])
		if(is.na(what$soil_a_om[ii])) finalTable $soil_a_om[ii] = weightedAverage(what$soil_a_om[cond], weights[cond])
		if(is.na(what$soil_a_sand[ii])) finalTable $soil_a_sand[ii] = weightedAverage(what$soil_a_sand[cond], weights[cond])
		if(is.na(what$soil_a_clay[ii])) finalTable $soil_a_clay[ii] = weightedAverage(what$soil_a_clay[cond], weights[cond])
		if(is.na(what$soil_a_silt[ii])) finalTable $soil_a_silt[ii] = weightedAverage(what$soil_a_silt[cond], weights[cond])
		
		if(is.na(what$soil_b_thickness[ii])) finalTable $soil_b_thickness[ii] = weightedAverage(what$soil_b_thickness[cond], weights[cond])
		if(is.na(what$soil_b_om[ii])) finalTable $soil_b_om[ii] = weightedAverage(what$soil_b_om[cond], weights[cond])
		if(is.na(what$soil_b_sand[ii])) finalTable $soil_b_sand[ii] = weightedAverage(what$soil_b_sand[cond], weights[cond])
		if(is.na(what$soil_b_clay[ii])) finalTable $soil_b_clay[ii] = weightedAverage(what$soil_b_clay[cond], weights[cond])
		if(is.na(what$soil_b_silt[ii])) finalTable $soil_b_silt[ii] = weightedAverage(what$soil_b_silt[cond], weights[cond])
        
        if(is.na(what$soil_c_om[ii])) finalTable $soil_c_om[ii] = weightedAverage(what$soil_c_om[cond], weights[cond])
        ##----------- fixing rtz and soildepth, 4011
        # finalTable[ii,]
        finaldepth = max(
            finalTable$soil_maxrootdepth[ii] * finalTable$soil_r_maxrtz_ratio[ii],
            finalTable$soil_activedepth[ii] * finalTable$soil_r_active_ratio[ii],
            finalTable$soil_maxrootdepth[ii] * soildepth_rtz_ratio,
            finalTable$soil_activedepth[ii] * soildepth_act_ratio,
            finalTable$soildepth[ii], -1, na.rm=T)
        if(finaldepth>0){
            finalTable$soildepth[ii] = finaldepth
            finalTable$soil_r_thickness[ii] = finaldepth - finalTable$soil_maxrootdepth[ii]
            finalTable$soil_r_maxrtz_ratio[ii] = finalTable$soildepth[ii]/finalTable$soil_maxrootdepth[ii]
            finalTable$soil_r_active_ratio[ii] = finalTable$soildepth[ii]/finalTable$soil_activedepth[ii]
        }
      
    }#ii for
	finalTable $check = sapply(seq_len(dim(finalTable)[1]),function(ii){ sum(is.na(finalTable[ii,]))>0 })
	if(sum(finalTable $check)==sum(what $check)) constrainedDist=constrainedDist*1.5
	what = finalTable
	counting = counting+1
	
}#while
# write.csv(finalTable,'~/Downloads/test.csv',row.names=F) # debug
# sum(finalTable$soildepth < finalTable$soil_maxrootdepth)
# sum(finalTable$soildepth < finalTable$soil_activedepth)





active_wa = finalTable$soil_a_thickness/(finalTable$soil_a_thickness+finalTable$soil_b_thickness)
active_wb = 1 - active_wa

full_wa = finalTable$soil_a_thickness/(finalTable$soil_a_thickness+finalTable$soil_b_thickness+finalTable$soil_c_thickness)
full_wb = finalTable$soil_b_thickness/(finalTable$soil_a_thickness+finalTable$soil_b_thickness+finalTable$soil_c_thickness)
full_wc = 1 - full_wa - full_wb

	# debug:
	# range(finalTable$soil_a_meanz)
	# range(finalTable$soil_a_ksat)
	# range(finalTable$soil_a_POR)
	# range(finalTable$soil_a_om)
	
	# range(finalTable$soil_b_meanz)
	# range(finalTable$soil_b_ksat)
	# range(finalTable$soil_b_POR)
	# range(finalTable$soil_b_om) ##<---- 0 om
	
	# range(finalTable$soil_c_meanz)
	# range(finalTable$soil_c_ksat)
	# range(finalTable$soil_c_POR)
	# range(finalTable$soil_c_om) ##<---- 0 om
	
profileRates = do.call(rbind, lapply(what_LEN, function(ii){
	profile = data.frame(
		z = c(finalTable$soil_a_meanz[ii],finalTable$soil_b_meanz[ii],finalTable$soil_c_meanz[ii]),
		ksat = c(finalTable$soil_a_ksat[ii], finalTable$soil_b_ksat[ii], finalTable$soil_c_ksat[ii]),
		por = c(finalTable$soil_a_POR[ii], finalTable$soil_b_POR[ii], finalTable$soil_c_POR[ii]),
        om = c(finalTable$soil_a_om[ii], finalTable$soil_b_om[ii], finalTable$soil_c_om[ii])
	)#
	if(sum(profile$ksat<=0)>0){ print(paste('soil profile',ii,'has 0 ksat')); profile$ksat[profile$ksat<=0]=0.0001; }
	if(sum(profile$por <=0)>0){ print(paste('soil profile',ii,'has 0 por' )); profile$por[profile$por <=0]=0.0001; }
	if(sum(profile$om<=0)>0){ profile$om[profile$om<=0]=0.00001; }
	ksat_log = log(profile$ksat)
	por_log = log(profile$por)
	om_log = log(profile$om)
	
	ksat_cof = lm(ksat_log~z,profile)$coefficient
	por_cof = lm(por_log~z,profile)$coefficient
    om_cof = lm(om_log~z,profile)$coefficient
	return <- unlist(list(
		ksat0 = exp(ksat_cof[1]),
		ksatdecay = ifelse(abs(1/ksat_cof[2])>4000,4000,-1/ksat_cof[2]),
		por0 = exp(por_cof[1]),
		pordecay = ifelse(abs(1/por_cof[2])>4000,4000,-1/por_cof[2]),
        om0 = exp(om_cof[1]),
        omdecay = ifelse(abs(1/om_cof[2])>4000,4000,-1/om_cof[2])
	))# 
}))#
mean_ksat = finalTable$soil_a_ksat*full_wa + finalTable$soil_b_ksat*full_wb + finalTable$soil_c_ksat*full_wc
mean_por = finalTable$soil_a_POR*full_wa + finalTable$soil_b_POR*full_wb + finalTable$soil_c_POR*full_wc

# adjust for rhessys
#cond = profileRates[,2]<0; profileRates[cond,1]= mean_ksat[cond]; # cbind(profileRates[cond,1], mean_ksat[cond])
#cond = profileRates[,4]<0;
#profileRates[cond,3]= mean_por[cond]; # cbind(profileRates[cond,3], mean_por[cond])
#profileRates[cond,4]= 4000


rhessys_soil_table = data.frame(
	cat = finalTable$catid,
	mukey = finalTable$mukey,
	texture = finalTable$soil_texture,
	soildepth = finalTable$soildepth,
	albedo = finalTable$soil_a_albedo,
	activedepth = finalTable$soil_activedepth,
	maxrootdepth = finalTable$soil_maxrootdepth,
	
	sand = finalTable$soil_a_sand*active_wa + finalTable$soil_b_sand*active_wb,
	clay = finalTable$soil_a_clay*active_wa + finalTable$soil_b_clay*active_wb,
	silt = finalTable$soil_a_silt*active_wa + finalTable$soil_b_silt*active_wb,
	om = finalTable$soil_a_om + finalTable$soil_b_om, # kgC/m2
    omdecay = profileRates[,6],
	
	bulkdensity = finalTable$soil_a_BD*full_wa + finalTable$soil_b_BD*full_wb + finalTable$soil_c_BD*full_wc,
	por_size_index = finalTable$soil_a_por_size_index*full_wa + finalTable$soil_b_por_size_index*full_wb + finalTable$soil_c_por_size_index*full_wc,
	psi_air_entry = finalTable$soil_a_psi_air_entry*full_wa + finalTable$soil_b_psi_air_entry*full_wb + finalTable$soil_c_psi_air_entry*full_wc,
	particledensity = finalTable$soil_a_PD*full_wa + finalTable$soil_b_PD*full_wb + finalTable$soil_c_PD*full_wc,
	
	ksat0 = profileRates[,1],
	ksatdecay = profileRates[,2],
	por0 = profileRates[,3],
	pordecay = profileRates[,4]
	)# data.frame
	


write.csv(rhessys_soil_table,arg[3], row.names=F)


























