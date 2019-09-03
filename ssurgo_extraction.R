mapunitTitle=c('musym','muname','mukind','mustatus','muacres','mapunitlfw_l','mapunitlfw_r','mapunitlfw_h','mapunitpfa_l','mapunitpfa_r','mapunitpfa_h','farmlndcl','muhelcl','muwathelcl','muwndhelcl','interpfocus','invesintens','iacornsr','nhiforsoigrp','nhspiagr','vtsepticsyscl','mucertstat','lkey','mukey')
compTitle=c('comppct_l','comppct_r','comppct_h','compname','compkind','majcompflag','otherph','localphase','slope_l','slope_r','slope_h','slopelenusle_l','slopelenusle_r','slopelenusle_h','runoff','tfact','wei','weg','erocl','earthcovkind1','earthcovkind2','hydricon','hydricrating','drainagecl','elev_l','elev_r','elev_h','aspectccwise','aspectrep','aspectcwise','geomdesc','albedodry_l','albedodry_r','albedodry_h','airtempa_l','airtempa_r','airtempa_h','map_l','map_r','map_h','reannualprecip_l','reannualprecip_r','reannualprecip_h','ffd_l','ffd_r','ffd_h','nirrcapcl','nirrcapscl','nirrcapunit','irrcapcl','irrcapscl','irrcapunit','cropprodindex','constreeshrubgrp','wndbrksuitgrp','rsprod_l','rsprod_r','rsprod_h','foragesuitgrpid','wlgrain','wlgrass','wlherbaceous','wlshrub','wlconiferous','wlhardwood','wlwetplant','wlshallowwat','wlrangeland','wlopenland','wlwoodland','wlwetland','soilslippot','frostact','initsub_l','initsub_r','initsub_h','totalsub_l','totalsub_r','totalsub_h','hydgrp','corcon','corsteel','taxclname','taxorder','taxsuborder','taxgrtgroup','taxsubgrp','taxpartsize','taxpartsizemod','taxceactcl','taxreaction','taxtempcl','taxmoistscl','taxtempregime','soiltaxedition','castorieindex','flecolcomnum','flhe','flphe','flsoilleachpot','flsoirunoffpot','fltemik2use','fltriumph2use','indraingrp','innitrateleachi','misoimgmtgrp','vasoimgtgrp','mukey','cokey')
chorizonTitle = c('hzname','desgndisc','desgnmaster','desgnmasterprime','desgnvert','hzdept_l','hzdept_r','hzdept_h','hzdepb_l','hzdepb_r','hzdepb_h','hzthk_l','hzthk_r','hzthk_h','fraggt10_l','fraggt10_r','fraggt10_h','frag3to10_l','frag3to10_r','frag3to10_h','sieveno4_l','sieveno4_r','sieveno4_h','sieveno10_l','sieveno10_r','sieveno10_h','sieveno40_l','sieveno40_r','sieveno40_h','sieveno200_l','sieveno200_r','sieveno200_h','sandtotal_l','sandtotal_r','sandtotal_h','sandvc_l','sandvc_r','sandvc_h','sandco_l','sandco_r','sandco_h','sandmed_l','sandmed_r','sandmed_h','sandfine_l','sandfine_r','sandfine_h','sandvf_l','sandvf_r','sandvf_h','silttotal_l','silttotal_r','silttotal_h','siltco_l','siltco_r','siltco_h','siltfine_l','siltfine_r','siltfine_h','claytotal_l','claytotal_r','claytotal_h','claysizedcarb_l','claysizedcarb_r','claysizedcarb_h','om_l','om_r','om_h','dbtenthbar_l','dbtenthbar_r','dbtenthbar_h','dbthirdbar_l','dbthirdbar_r','dbthirdbar_h','dbfifteenbar_l','dbfifteenbar_r','dbfifteenbar_h','dbovendry_l','dbovendry_r','dbovendry_h','partdensity','ksat_l','ksat_r','ksat_h','awc_l','awc_r','awc_h','wtenthbar_l','wtenthbar_r','wtenthbar_h','wthirdbar_l','wthirdbar_r','wthirdbar_h','wfifteenbar_l','wfifteenbar_r','wfifteenbar_h','wsatiated_l','wsatiated_r','wsatiated_h','lep_l','lep_r','lep_h','ll_l','ll_r','ll_h','pi_l','pi_r','pi_h','aashind_l','aashind_r','aashind_h','kwfact','kffact','caco3_l','caco3_r','caco3_h','gypsum_l','gypsum_r','gypsum_h','sar_l','sar_r','sar_h','ec_l','ec_r','ec_h','cec7_l','cec7_r','cec7_h','ecec_l','ecec_r','ecec_h','sumbases_l','sumbases_r','sumbases_h','ph1to1h2o_l','ph1to1h2o_r','ph1to1h2o_h','ph01mcacl2_l','ph01mcacl2_r','ph01mcacl2_h','freeiron_l','freeiron_r','freeiron_h','feoxalate_l','feoxalate_r','feoxalate_h','extracid_l','extracid_r','extracid_h','extral_l','extral_r','extral_h','aloxalate_l','aloxalate_r','aloxalate_h','pbray1_l','pbray1_r','pbray1_h','poxalate_l','poxalate_r','poxalate_h','ph2osoluble_l','ph2osoluble_r','ph2osoluble_h','ptotal_l','ptotal_r','ptotal_h','excavdifcl','excavdifms','cokey','chkey')
chtexgrpTitle=c('texture','stratextsflag','rvindicator','texdesc','chkey','chtgkey')
chtexturTitle=c('texcl','lieutex','chtgkey','chtkey')
chporesTitle=c('poreqty','poresize','porecont','poreshp','rvindicator','chkey','chporeskey')


arg=commandArgs(T);
arg=c('/Users/laurencelin/Library/Mobile Documents/com~apple~CloudDocs/Workspace/current_projects/baisman_new_setup/raw_data/MD005')
target = paste(arg[1],'/tabular',sep='') # e.g., 'VA113/tabular'
print(target)


mapunit = read.table(paste(target,'/mapunit.txt',sep=''),sep='|',stringsAsFactors=F) #mukey
colnames(mapunit) = mapunitTitle
	## 1 mukey to many componments [cokey]
comp = read.table(paste(target,'/comp.txt',sep=''),sep='|',stringsAsFactors=F) #cokey
colnames(comp) = compTitle
	## 1 cokey to many horizons [chkey]
chorizon = read.table(paste(target,'/chorizon.txt',sep=''),sep='|',stringsAsFactors=F) #chkey
colnames(chorizon) = chorizonTitle

chtexgrp = read.table(paste(target,'/chtexgrp.txt',sep=''),sep='|',stringsAsFactors=F) #chkey
colnames(chtexgrp) = chtexgrpTitle

chtextur = read.table(paste(target,'/chtextur.txt',sep=''),sep='|',stringsAsFactors=F) #chkey
colnames(chtextur) = chtexturTitle

# chpores = read.table(paste(target,'/chpores.txt',sep=''),sep='|',stringsAsFactors=F) #chkey
# colnames(chpores) = chporesTitle
#chpores # new table
	# chkey -> chorizon
	# chporeskey ->	

## mapunit:component:horizon
## extract (weighted average[thickness] of soil horizon in a soil component; then weighted [%composition] by componments):
	# 1 permeability
	# 2 watercapacity
	# 3 bulkdensity
	# 4 saturatedhydraulicconductivity
	# 5 erodibility
	# 6 field capacity
	# 7 porosity
	# 8 soilthickness
	# 9 organ matter 
    # AWC?

	# relationshup
	# 1 mukey: N cokey (componments) : N chkey (horizons)
 	#        : fraction "comppct_r"  : fraction "hzdepb_r" -> thickness

replaceNA0 = function(x,value){
	x[is.na(x) | x<=0] = value;
	return <- x;	
}#function

horizonOrder = c('O','A','E','B','C','R')	# there is horizon 'H' H horizons or layers: Layers dominated by organic material, formed from accumulations of undecomposed or partially decomposed organic material at the soil surface which may be underwater. All H horizons are saturated with water for prolonged periods or were once saturated but are now artificially drained. An H horizon may be on top of mineral soils or at any depth beneath the surface if it is buried. (http://www.fao.org/3/w8594e/w8594e0g.htm)
    # H - B microbial activities, but no evidence of such activities in C or deeper
    # plant roots can penetrate C horizon
    # R is hard bedrock

horizonThicknessDefault = c(2, 5, 5, 10, 10, 10)*0.0254 #inches to meter; https://www.sheffield.ac.uk/ssa/soil-facts/horizons
horizonUpperZDefault = c(0, cumsum(horizonThicknessDefault)[1:(length(horizonThicknessDefault)-1)])
horizonlowerZDefault = cumsum(horizonThicknessDefault)[1:length(horizonThicknessDefault)]

	# this part need to be read in from a file!
	soilscoreNames = c('Clay','Silty clay','Silty clay loam','Sandy clay','Sandy clay loam','Clay loam','Silt','Silt loam','Loam','Sand','Loamy sand','Sandy loam')
	soilscore = c(1,2,3,4,5,6,7,8,9,10,11,12); names(soilscore)=soilscoreNames
	
	soilscore_db = data.frame(
		sand = c(0.2,0.05,0.1,0.5,0.6,0.34,0.05,0.17,0.43,0.92,0.85,0.7),
		silt = c(0.2,0.475,0.55,0.05,0.1,0.33,0.9,0.7,0.39,0.05,0.1,0.2),
		clay = c(0.6,0.475,0.35,0.45,0.3,0.33,0.05,0.13,0.18,0.03,0.05,0.1),
		ksat = c(0.111,0.089,0.147,0.188,0.544,0.212,0.622,0.622,0.6,15.2064,13.5,3),
		por = c(0.482,0.492,0.477,0.426,0.42,0.476,0.485,0.485,0.451,0.395,0.41,0.435),
		por_size_index = c(0.088,0.96,0.129,0.096,0.14,0.117,0.189,0.189,0.204,0.204,0.228,0.204),
		psi_air_entry = c(0.405,0.49,0.356,0.153,0.299,0.63,0.786,0.786,0.478,0.121,0.09,0.218),
		albedo = c(0.23,0.229,0.238,0.253,0.265,0,0.253,0.253,0.28,0.28,0.29,0.28)
	)#data.frame

mukey = mapunit[,'mukey']	
mukeyHold = matrix(NA,length(mukey), 19);
profileDEPTH = seq(0,10,0.01)
profile_ksat_table = data.frame(z=profileDEPTH)
profile_POR_table = data.frame(z=profileDEPTH)
for(i in seq_along(mukey)){
	## for each mukey, we have N cokey
    ## mukey soil texture class
    
    
	cond1 = comp[,'mukey']== mukey[i]; sum(cond1)
	
	## find overall soil type and develop default parameter values; not component/horizon specific
	overall_NODATA = F
    cond25 = chorizon[,'cokey'] %in% comp[cond1,'cokey']; sum(cond25)
    if(sum(cond25)>0){
    	cond3 = chtexgrp[,'chkey'] %in% chorizon[cond25,'chkey']; sum(cond3)
	    cond4 = chtextur[,'chtgkey'] %in% chtexgrp[cond3,'chtgkey']; sum(cond4)
        default_overallsoil = table(soilscore[match(chtextur[cond4,1], soilscoreNames)]) # could have chtextur[cond4,1]='' problem!
        
	    	default_overallsoil_index = as.numeric(names(default_overallsoil))
	        default_overallsoil_weight = default_overallsoil/sum(default_overallsoil)
	        default_overallsoil_char = unlist(list(
	            sand = sum(soilscore_db$sand[default_overallsoil_index]* default_overallsoil_weight),
	            silt = sum(soilscore_db$silt[default_overallsoil_index]* default_overallsoil_weight),
	            clay = sum(soilscore_db$clay[default_overallsoil_index]* default_overallsoil_weight),
				por_size_index = sum(soilscore_db$por_size_index[default_overallsoil_index]* default_overallsoil_weight),
				psi_air_entry = sum(soilscore_db$psi_air_entry[default_overallsoil_index]* default_overallsoil_weight),
				albedo = sum(soilscore_db$albedo[default_overallsoil_index]* default_overallsoil_weight),
	            ksat = sum(soilscore_db$ksat[default_overallsoil_index]* default_overallsoil_weight),
	            POR = sum(soilscore_db$por[default_overallsoil_index]* default_overallsoil_weight)
	            ))#
	    soilID_ = ifelse(length(default_overallsoil)>0, as.numeric(names(which.max(default_overallsoil))), NA)
	    soilID_name = ifelse(length(default_overallsoil)>0, soilscoreNames[soilID_], NA)
	    # default_overallsoil_char
	    # soilID_name
	}else{
		overall_NODATA = T
		default_overallsoil_char = unlist(list(
	            sand = NA, silt = NA, clay = NA,
				por_size_index = NA, psi_air_entry = NA, albedo = NA,
	            ksat = NA, POR = NA
	            ))#
		soilID_name = NA
    }# cond25
	
	
	## find all the components and horizons and extract specific information
	if(! overall_NODATA){
		
		## developing component (vertical partition) > horizon (horizonal partition)
		componentSummary_zdown = list()
		componentSummary_zup = list()
		componentSummary_ksat = list()
		componentSummary_por = list()
        componentSummary_ksat_profile = list()
        componentSummary_por_profile = list()
		componentSummaryTable = as.data.frame( do.call(rbind, lapply( comp[cond1,'cokey'], function(cokey){
			# linking horizons to component
			
	        # cokey = comp[cond1,'cokey'][1] # for debugging
            cond2 = chorizon[,'cokey']== cokey; sum(cond2) # finding component-specific horizons
			if(sum(cond2)>0){
                # found matched component-specific horizons
		        component_NODATA = F
                
		        # sort the component specific soil types --> develop component-specific default parameters
		        cond3 = chtexgrp[,'chkey'] %in% chorizon[cond2,'chkey']; sum(cond3)
		        cond4 = chtextur[,'chtgkey'] %in% chtexgrp[cond3,'chtgkey']; sum(cond4)
		        default_soil = table(soilscore[match(chtextur[cond4,1], soilscoreNames)]); length(default_soil)
		        	if(length(default_soil)==0){
		        		# cannot find the component soil name
		        		component_NODATA = T
		        		default_soil = default_overallsoil;
		        	}#if
		        default_soil_index = as.numeric(names(default_soil))
		        default_soil_weight = default_soil/sum(default_soil)
		        default_soil = unlist(list(
		            sand = sum(soilscore_db$sand[default_soil_index]*default_soil_weight),
		            silt = sum(soilscore_db$silt[default_soil_index]*default_soil_weight),
		            clay = sum(soilscore_db$clay[default_soil_index]*default_soil_weight),
                    por_size_index = sum(soilscore_db$por_size_index[default_soil_index]*default_soil_weight),  # <<-- not from surrgo
                    psi_air_entry = sum(soilscore_db$psi_air_entry[default_soil_index]*default_soil_weight),    # <<-- not from surrgo
                    albedo = sum(soilscore_db$albedo[default_soil_index]*default_soil_weight),                  # <<-- not from surrgo
		            ksat = sum(soilscore_db$ksat[default_soil_index]*default_soil_weight),
		            POR = sum(soilscore_db$por[default_soil_index]*default_soil_weight)
		            ))#
	
    
		        horizonSummaryTable = data.frame(id = chorizon[cond2,'chkey'])
                horizonSummaryTable$horizonThinkness = NA#(chorizon[cond2,'hzdepb_r']-chorizon[cond2,'hzdept_r'])*0.01 #cm -> m
                horizonSummaryTable$horizonMeanDepth = NA#0.5*(chorizon[cond2,'hzdepb_r']+chorizon[cond2,'hzdept_r'])*0.01 # m (real depth)
                horizonSummaryTable$horizonTopDepth = chorizon[cond2,'hzdept_r']*0.01 #m
                horizonSummaryTable$horizonBottomDepth = chorizon[cond2,'hzdepb_r']*0.01 #m <<-- for microbial activities depth & max root depth
                # microbial activities stop at the bottom of B
                # max plant root could be at the bottom of C
		        horizonSummaryTable$horizonName = sapply(chorizon[cond2,'hzname'],function(nn){
                    tmp = sapply(horizonOrder, function(x){grepl(x,nn)})
                    if(sum(tmp)>0) return <- horizonOrder[tmp]
                    else return <- nn
                })#
                horizonSummaryTable$horizonIndex = match(horizonSummaryTable$horizonName,horizonOrder)
                horizonSummaryTable$sand = replaceNA0(chorizon[cond2,'sandtotal_r']*0.01, default_soil['sand']); #% -> [0,1]
		        horizonSummaryTable$silt = replaceNA0(chorizon[cond2,'silttotal_r']*0.01, default_soil['silt']); #% -> [0,1]
		        horizonSummaryTable$clay = replaceNA0(chorizon[cond2,'claytotal_r']*0.01, default_soil['clay']); #% -> [0,1]
                
                horizonSummaryTable$BD = chorizon[cond2,'dbovendry_r'] #g/cm3 bulk density [0.02 - 2.6]; the oven dry weight of the less than 2 mm soil material per unit volume of soil exclusive of the desication cracks, measured on a costed clod;
                horizonSummaryTable$fc = chorizon[cond2,'wthirdbar_r']*0.01 #% -> [0,1]; field capacity; the volumetric content of soil water retained at the tension of 1/3 bar (33 kPa), expressed as a percentage of the whole soil.
                horizonSummaryTable$awc = chorizon[cond2,'awc_r']*0.01 # cm/cm; Amount of water that an increment of soil depth, inclusive of fragements, can store that is available to plants.
                horizonSummaryTable$kffact = chorizon[cond2,'kffact'] # soil erodibility factor [0-1]; an erodibility factor which quantifies the susceptibility of doil particles to detechment by water
                horizonSummaryTable$om = chorizon[cond2,'om_r']*0.01 #% -> [0,1]; the amount by weight of decomposed plant and animal residue expressed as weight percentage of the less than 2 mm soil material (not the top 2 mm soil!).
                horizonSummaryTable$density = chorizon[cond2,'partdensity'] # g/cm3; [0.01 - 5] mass per unit volume (not including pore space) of the solid particle either mineral or organic.
                
                horizonData_zdown = chorizon[cond2,'hzdepb_r']*0.01
                horizonData_zup = chorizon[cond2,'hzdept_r']*0.01
                horizonData_ksat = chorizon[cond2,'ksat_r']*1e-6*3600*24 # um/s [0-705] --> m/day; the amount of water that would move verticaly through a unit area of saturated soil un inut time under unit hydraulic gradient.
                horizonData_por = chorizon[cond2,'wsatiated_r']*0.01; # %; estimated volumetric soil water content at or near zero bar tension, expressed as a percentage of the whole soil.
                horizonSummaryTable$ksat = replaceNA0(chorizon[cond2,'ksat_r']*1e-6*3600*24, default_soil['ksat']); # m/day;
                    horizonSummaryTable$ksat[horizonSummaryTable$ksat<=0] = 0.001 # m/day
		        horizonSummaryTable$POR = replaceNA0(chorizon[cond2,'wsatiated_r']*0.01, default_soil['POR']); #% -> [0,1]
		            horizonSummaryTable$POR[horizonSummaryTable$POR<=0] = 0.01
		            
                    
                # sorting horizonSummaryTable by horizonName
                horizonNameOrder = order(horizonSummaryTable$horizonIndex)
                horizonSummaryTable = horizonSummaryTable[horizonNameOrder,]
                horizonSummaryTable$horizonThinkness = sapply(seq_along(horizonNameOrder),function(jj){
                    if( is.na(horizonSummaryTable$horizonTopDepth[jj]) | horizonSummaryTable$horizonTopDepth[jj]<0 |
                    is.na(horizonSummaryTable$horizonBottomDepth[jj]) | horizonSummaryTable$horizonBottomDepth[jj]<0){
                        ## no data boundaries
                        return <- horizonThicknessDefault[horizonSummaryTable$horizonIndex[jj]]
                    }else{
                        ## yes data define boundaries
                        tmp = horizonSummaryTable$horizonBottomDepth[jj] - horizonSummaryTable$horizonTopDepth[jj]
                        return <- ifelse(tmp<=0, horizonThicknessDefault[horizonSummaryTable$horizonIndex[jj]], tmp)
                    }# if
                })#sapply
                endpoints = c(0, cumsum(horizonSummaryTable$horizonThinkness))
                horizonSummaryTable$horizonTopDepth = endpoints[1:(length(endpoints)-1)]
                horizonSummaryTable$horizonBottomDepth = endpoints[-1]
                horizonSummaryTable$horizonMeanDepth = 0.5*(horizonSummaryTable$horizonTopDepth+horizonSummaryTable$horizonBottomDepth)
                horizonSummaryTable$microbialActiveDepth = max(sapply(seq_along(horizonSummaryTable$horizonBottomDepth), function(ii){
                        ifelse(horizonSummaryTable$horizonName[ii] %in% horizonOrder[1:4],horizonSummaryTable$horizonBottomDepth[ii],NA)
                    }),na.rm=T)
                horizonSummaryTable$maxRootDepth = max(sapply(seq_along(horizonSummaryTable$horizonBottomDepth), function(ii){
                        ifelse(horizonSummaryTable$horizonName[ii] %in% horizonOrder[5:6],horizonSummaryTable$horizonBottomDepth[ii],1.0 )
                    }),na.rm=T)
                horizonSummaryTable$horizonWeight = horizonSummaryTable$horizonThinkness/sum(horizonSummaryTable$horizonThinkness)
                
                
                ## extracting profile every cm from 0 to 10 m
                numHorizon = dim(horizonSummaryTable)[1]
		        if(numHorizon>1){
                    # more than one horizon layers
                    # developing profile by each components
		            profileFunction = approxfun(horizonSummaryTable$horizonMeanDepth, horizonSummaryTable$ksat); #splinefun
                    profileFunction_boundaryZ = range(horizonSummaryTable$horizonMeanDepth)
                    profileFunction_boundaryY = c(
                        horizonSummaryTable$ksat[1],
                        horizonSummaryTable$ksat[numHorizon]
                    )
		            componentSummary_ksat_profile[[toString(cokey)]]<<- sapply(profileDEPTH, function(z){
		                if(z<profileFunction_boundaryZ[1]) return <- profileFunction_boundaryY[1]
		                else if(z>profileFunction_boundaryZ[2]) return <- profileFunction_boundaryY[2]
		                else profileFunction(z)
		            })
                    # plot(componentSummary_ksat_profile[[toString(cokey)]], profileDEPTH, type='l', ylim=c(2,0));
                    # points(horizonSummaryTable$ksat, horizonSummaryTable$horizonMeanDepth, col='red')
                    # cbind(horizonSummaryTable$horizonMeanDepth, horizonSummaryTable$ksat)
                    
                    
		            profileFunction = approxfun(horizonSummaryTable$horizonMeanDepth, horizonSummaryTable$POR); #splinefun
		            profileFunction_boundaryZ = range(horizonSummaryTable$horizonMeanDepth)
		            profileFunction_boundaryY = c(
                        horizonSummaryTable$POR[1],
                        horizonSummaryTable$POR[numHorizon]
		            )
		            componentSummary_por_profile[[toString(cokey)]]<<- sapply(profileDEPTH, function(z){
		                if(z<profileFunction_boundaryZ[1]) return <- profileFunction_boundaryY[1]
		                else if(z>profileFunction_boundaryZ[2]) return <- profileFunction_boundaryY[2]
		                else profileFunction(z)
		            })
		            # plot(profile_POR, profileDEPTH, type='l', ylim=c(10,0))
                    # points(horizonSummaryTable$POR, horizonSummaryTable$horizonMeanDepth, col='red')
                    
		        }else{
                    # only one horizon layer from the data
		            componentSummary_ksat_profile[[toString(cokey)]]<<- horizonSummaryTable$ksat[1]*exp(-profileDEPTH* 0.12) # 1/8.33333
		            componentSummary_por_profile[[toString(cokey)]]<<- horizonSummaryTable$POR[1]*exp(-profileDEPTH* 0.00025) # decay by 1/4000
		            
		        }
	    
	        	# debugging each component
	        	# horizonSummaryTable
	        	# plot(profile_ksat, profileDEPTH,type='l', ylim=c(10,0), xlim=c(0,1))
	        	# lines(profile_POR, profileDEPTH, lty=2)
	        
                ## aggregating horizons to component (i.e., cokey)
	        	componentSummary_zdown[[toString(cokey)]]<<- horizonData_zdown
				componentSummary_zup[[toString(cokey)]]<<- horizonData_zup
				componentSummary_ksat[[toString(cokey)]]<<- horizonData_ksat
				componentSummary_por[[toString(cokey)]]<<- horizonData_por
                ############ <<---- working here ----->> #################
		        return <- unlist(list(
		        	id = cokey,
		            Thinkness = sum(horizonSummaryTable$horizonThinkness,na.rm=T), # meter
		        	down2Horizon = horizonSummaryTable$horizonIndex[numHorizon],
                    # ... general lumped proporties (not being used in RHESSys)
                    kffact = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$kffact, na.rm=T),
                    fc = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$fc, na.rm=T),
                    # ... lumped proporties solute sorption (in deep soil)
                    partdensity = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$density, na.rm=T),
                    om = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$om, na.rm=T),
                    BD = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$BD, na.rm=T),
                        # PARTICLE_DENSITY 2.65  - soil particle density g/cm3 (Dingman) in RHESSys is
                        # a) used to calculate kg_soil -> soil NO3 density (Parton et al. 1996; µgN/g soil)
                        # b) used to calculate kg_soil -> max_nit_rate = kg_soil * MAX_RATE * 1e-06; // (mgN/kg soil/day) * (kg soil/m2) ==> mgN/m2/day ==> 1e-06 kgN/m2/day (Parton et al. 1996)
                        # c) use to calculate "absorptionConst" = PARTICLE_DENSITY * 1000.0 * N_absorption_rate;
                        # ammonium binding to soil (https://crops.extension.iastate.edu/cropnews/2019/03/understanding-anhydrous-ammonia-application-soil):
                        # dissolved ammonium reacts with soil organic matter and clay, and ammonium ions attach on the soil cation exchange complex; These reactions all tend to limit the movement of ammonia
                        # Cation-exchange capacity is a measure of how many cations can be retained on soil particle surfaces. Negative charges on the surfaces of soil particles bind positively-charged atoms or molecules, but allow these to exchange with other positively charged particles in the surrounding soil water
                        # so what binds ammonium in the soil is not the soil mineral!
                        # However, DOC is usually negative charge (soil anion sorption), but DOC has some functional groupswhich form positive sites.
                        # DOC sorption is mostly via anion exchange in the soil depth (mineral density)
                        # most soil particles are nagetively charged, but some anions (e.g., metal oxides) are also bound
                        # soil anion exchange capacity is much lower than the cation exchange capacity (differen by 20X less)
                        # in general both exchanges are highly related to clay
                        # additionally, DON and DOP behave differently than DOC
                    # ... lumped porporties for biochemistry in soil column (0 - horizon B)...
                    # microbial activities stop at the bottom of B
                    # max plant root could be at the bottom of C
                    clay = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$clay, na.rm=T),
                    sand = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$sand, na.rm=T),
                    silt = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$silt, na.rm=T), # not used for any calculations
                    
    
                    
                    awc = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$awc, na.rm=T),
                    #awc: cm/cm; Amount of water that an increment of soil depth, inclusive of fragements, can store that is available to plants. volume fraction and is commonly estimated as the difference between the water contents at 1/3 and 15 bars (permanent wilting point).
                    #POR: estimated volumetric soil water content at or near zero bar tension, expressed as a percentage of the whole soil.
                    
                    # not reporting here; profile should be studied after merging components
                    # ... if ksat0 and ksat_decay are calculated correctly from the ssurgo dataset
                    # A) last Ksat matching; (last or certain meter depth?)
                    # B) ksat area equal
                    #ksat = sum(horizonSummaryTable$horizonWeight * horizonSummaryTable$ksat, na.rm=T),
                    #vksat = sum(horizonSummaryTable$horizonThinkness,na.rm=T) / sum(horizonSummaryTable$horizonThinkness/horizonSummaryTable$ksat, na.rm=T)
		        ))
                # within cond2 if
	        }else{
	        	# component_NODATA = T
	        	componentSummary_zdown[[toString(cokey)]]<<- NA
				componentSummary_zup[[toString(cokey)]]<<- NA
				componentSummary_ksat[[toString(cokey)]]<<- NA
				componentSummary_por[[toString(cokey)]]<<- NA
                componentSummary_ksat_profile[[toString(cokey)]]<<- NA
                componentSummary_por_profile[[toString(cokey)]]<<- NA
	        	return <- unlist(list(
		        	id = cokey,
		            horizonThinkness = NA, # meter
		        	deepHorizon = NA, # meter
		            ksat = NA,
		            vksat = NA, 
		        	sand = NA,
		        	silt = NA,
		        	clay = NA,
		        	BD = NA,
		        	partdensity = NA,
		        	fc = NA,
		        	awc = NA,
		        	kffact = NA,
		        	om = NA #14
		        ))
	        }#if else 
		})) )# end of construction: componentSummaryTable
		componentSummaryTable$Percent = comp[cond1,'comppct_r']
		componentSummaryTable$Weight = componentSummaryTable$Percent/sum(componentSummaryTable$Percent)
		componentSummaryTable$profileWeight = rep(NA, dim(componentSummaryTable)[1])
		componentSummaryTable$profileWeight[ !is.na(componentSummaryTable[,15]) ] = componentSummaryTable$Weight[ !is.na(componentSummaryTable[,15]) ]
		componentSummaryTable$profileWeight = componentSummaryTable$profileWeight/sum(componentSummaryTable$profileWeight,na.rm=T)
		componentSummaryTable$profileWeight[is.na(componentSummaryTable$profileWeight)] = 0


	    ## summarizing profiles from different components
	    profile_POR_index = seq(15,length.out=1001, by=1) # 15 - 1015
	    profile_POR = apply(componentSummaryTable[,profile_POR_index],2,function(xx){ sum(xx*componentSummaryTable$profileWeight,na.rm=T) } )
	    #plot(profile_POR, profileDEPTH, type='l', ylim=c(10,0))
	    
	    profile_ksat_index = seq(1016,length.out=1001, by=1) # 1016 - 2016
	    profile_ksat = apply(componentSummaryTable[,profile_ksat_index],2,function(xx){ sum(xx*componentSummaryTable$profileWeight,na.rm=T) } )
	    #plot(profile_ksat, profileDEPTH, type='l', ylim=c(10,0))
	   
	   
	    # cbind(componentSummaryTable[,c(1:14)], componentSummaryTable[,c('Percent','Weight','profileWeight')])
	    # componentSummary_zdown
		# componentSummary_zup 
		# componentSummary_ksat
		# componentSummary_por 
        upper_depth = min(sapply(seq_along(componentSummary_zdown),function(jj){
            num_of_sample = sum(!is.na(componentSummary_zdown[[jj]]))
            if(num_of_sample>0){
                # yes samples
                return <- min( componentSummary_zup[[jj]], na.rm=T )
            }else{
                # no samples
                return <- NA
            }#
        }),na.rm=T); upper_depth_index = which.min(abs(profileDEPTH-upper_depth))
        
        
        lower_depth = max(sapply(seq_along(componentSummary_zdown),function(jj){
            num_of_sample = sum(!is.na(componentSummary_zdown[[jj]]))
            if(num_of_sample>0){
                # yes samples
                return <- max( componentSummary_zdown[[jj]],na.rm=T ) 
            }else{
                # no samples
                return <- NA
            }#
        }),na.rm=T); lower_depth_index = which.min(abs(profileDEPTH-lower_depth))
        
        
        
        
        sample_index = upper_depth_index:lower_depth_index
        len = length(sample_index)
       
        sample_ksat_profile = profile_ksat[sample_index];
        tot_ksat = 0.01*0.5*sum(sample_ksat_profile[1:(len-1)]+sample_ksat_profile[2:len]);
        if( sum(is.na(sample_ksat_profile))==0 ){
        	# ... search lamda
            mukey_ksat0 = mean(sample_ksat_profile[1:10])
            search_lamda = seq(1/4000,1,0.0001)
            search_lamda_result = sapply(search_lamda, function(x){mukey_ksat0/x*(exp(-x*profileDEPTH[upper_depth_index])-exp(-x*profileDEPTH[lower_depth_index]))-tot_ksat})
            mukey_lamda = search_lamda[which.min(abs(search_lamda_result))]
            if( search_lamda_result[which.min(abs(search_lamda_result))]<0 & which.min(abs(search_lamda_result))==1 ){
                mukey_lamda=1/4000;
                mukey_ksat0 = tot_ksat*mukey_lamda/(exp(-mukey_lamda*profileDEPTH[upper_depth_index])-exp(-mukey_lamda*profileDEPTH[lower_depth_index]))
            }
            # c(mukey_lamda,mukey_ksat0)
            # plot(sample_ksat_profile,ylim=c(0,3)); lines(sample_index, mukey_ksat0*exp(-mukey_lamda*profileDEPTH[sample_index]), col='red')
        } else {  
        	mukey_ksat0 = NA
        	mukey_lamda = NA
        }
            
        sample_por_profile = profile_POR[sample_index]
        tot_por = 0.01*0.5*sum(sample_por_profile[1:(len-1)]+sample_por_profile[2:len]);
        if( sum(is.na(sample_por_profile))==0 ){
            # ... search lamda
            mukey_por0 = mean(sample_por_profile[1:10])
            search_lamdaP = seq(1/4000,1,0.0001)
            search_lamdaP_result = sapply(search_lamdaP, function(x){mukey_por0/x*(exp(-x*profileDEPTH[upper_depth_index])-exp(-x*profileDEPTH[lower_depth_index]))-tot_por})
            mukey_lamdaP = search_lamdaP[which.min(abs(search_lamdaP_result))]
            if( search_lamdaP_result[which.min(abs(search_lamdaP_result))]<0 & which.min(abs(search_lamdaP_result))==1 ){
                mukey_lamdaP = 1/4000;
                mukey_por0 = tot_por*mukey_lamda/(exp(-mukey_lamda*profileDEPTH[upper_depth_index])-exp(-mukey_lamda*profileDEPTH[lower_depth_index]))
            }
            # c(mukey_lamdaP,mukey_por0)
        }else{
        	mukey_por0 = NA
        	mukey_lamdaP = NA
        }
        
        
	    ## debugging all components
		    # componentSummaryTable[,c(
		    	# 'id','horizonThinkness','deepHorizon','ksat','vksat','sand','silt','clay','BD','partdensity',
		    	# 'fc','awc','kffact','om','Percent','Weight')]
		    	
	    	# plot(as.numeric(componentSummaryTable[1,profile_ksat_index]), profileDEPTH,type='l', ylim=c(10,0), xlim=c(0,10), col='gray')
			# lines(as.numeric(componentSummaryTable[2,profile_ksat_index]), profileDEPTH,col='gray')
	    	# lines(profile_ksat, profileDEPTH)
	    	
	    	# plot(as.numeric(componentSummaryTable[1, profile_POR_index]), profileDEPTH,type='l', ylim=c(10,0), xlim=c(0,1), col='gray')
			# lines(as.numeric(componentSummaryTable[2, profile_POR_index]), profileDEPTH,col='gray')
	    	# lines(profile_POR, profileDEPTH)
	
		## aggregating 
		mukey_sand = sum(componentSummaryTable$Weight*componentSummaryTable$sand,na.rm=T)
		mykey_silt = sum(componentSummaryTable$Weight*componentSummaryTable$silt,na.rm=T)
		mykey_clay = sum(componentSummaryTable$Weight*componentSummaryTable$clay,na.rm=T)
		mykey_total = mukey_sand + mykey_silt + mykey_clay
			mukey_sand = mukey_sand/mykey_total
			mykey_silt = mykey_silt/mykey_total
			mykey_clay = mykey_clay/mykey_total
		
			## fixing componentSummaryTable with NA (horizonThinkness, deepHorizon, ksat)
			mean_component_thickness = mean(componentSummaryTable$horizonThinkness[componentSummaryTable$horizonThinkness>0])
			componentSummaryTable$horizonThinkness[componentSummaryTable$horizonThinkness<=0] = mean_component_thickness
			
			soilhorizonLayers = componentSummaryTable$deepHorizon[order(componentSummaryTable$Weight, decreasing=T)]
			
			mean_component_vksat = mean(componentSummaryTable$vksat[!is.na(componentSummaryTable$vksat)])
			componentSummaryTable$vksat[is.na(componentSummaryTable$vksat)] = mean_component_vksat
			
		tmp = unlist(list(
			mukey = mukey[i],
			rhessys_soilid = soilID_, #soil_name = soilID_name,
			
			soilDepth = sum(componentSummaryTable$Weight*componentSummaryTable$horizonThinkness,na.rm=T), #<<- problem of NA
			soilhorizon = soilhorizonLayers[!is.na(soilhorizonLayers)][1], #<<- problem of NA
			
			ksat = sum(componentSummaryTable$Weight*componentSummaryTable$ksat,na.rm=T), 
			vksat = sum(componentSummaryTable$Weight*componentSummaryTable$vksat,na.rm=T), #<<- problem of NA
            ksat0 = mukey_ksat0,
            ksat0_decay = mukey_lamda,
            por0 = mukey_por0,
            por0_decay = 1/mukey_lamdaP,
       
			sand = mukey_sand, 
			silt = mykey_silt,
			clay = mykey_clay,
			
			BD = sum(componentSummaryTable$Weight*componentSummaryTable$BD,na.rm=T), #<<- problem of NA
			partdensity = sum(componentSummaryTable$Weight*componentSummaryTable$partdensity,na.rm=T), #<<- problem of NA
			fc = sum(componentSummaryTable$Weight*componentSummaryTable$fc,na.rm=T), #<<- problem of NA
			awc = sum(componentSummaryTable$Weight*componentSummaryTable$awc,na.rm=T), #<<- problem of NA
			kffact = sum(componentSummaryTable$Weight*componentSummaryTable$kffact,na.rm=T), #<<- problem of NA
			om = sum(componentSummaryTable$Weight*componentSummaryTable$om,na.rm=T) #<<- problem of NA
		))
        
		if(i==1) colnames(mukeyHold)=names(tmp)
        mukeyHold[i,] = tmp
		profile_ksat_table[,paste(mukey[i])] = profile_ksat
	    profile_POR_table[,paste(mukey[i])] = profile_POR
        # ... within overall_NODATA = F
        # ... components ... #
        
	}else{
		# overall_NODATA=T case (working on...)
		tmp = unlist(list(
			mukey = mukey[i],
			rhessys_soilid = NA, #soil_name = soilID_name,
			
			soilDepth = NA, #<<- problem of NA
			soilhorizon = NA, #<<- problem of NA
			
			ksat = NA, 
			vksat = NA, #<<- problem of NA
            ksat0 = NA,
            ksat0_decay = NA,
            por0 = NA,
            por0_decay = NA,
            
			sand = NA, 
			silt = NA,
			clay = NA,
			
			BD = NA, #<<- problem of NA
			partdensity = NA, #<<- problem of NA
			fc = NA, #<<- problem of NA
			awc = NA, #<<- problem of NA
			kffact = NA, #<<- problem of NA
			om = NA #<<- problem of NA
		))
		mukeyHold[i,] = tmp
		if(i==1) colnames(mukeyHold)=names(tmp)
		profile_ksat_table[,paste(mukey[i])] = rep(NA,1001)
	    profile_POR_table[,paste(mukey[i])] = rep(NA,1001)
	}# end of if
}#i
# write.csv(mukeyHold, paste('~/Downloads/soil_mukey_texture.csv',sep=''),row.names=F)
write.csv(mukeyHold, paste(arg[1],'/soil_mukey_texture.csv',sep=''),row.names=F)
write.csv(profile_ksat_table, paste(arg[1],'/soil_mukey_ksat.csv',sep=''),row.names=F)
write.csv(profile_POR_table, paste(arg[1],'/soil_mukey_por.csv',sep=''),row.names=F)

	## ----- checking, debugging
	#which(apply(mukeyHold,1,function(x){sum(is.na(x))})>0)
	# 46  47 144 166 168 169 170 171 172 182 197 202
	#mukeyHold[46,]
	
	# mukeyHold = read.csv(paste(arg[1],'/soil_mukey_texture.csv',sep=''))
	# which(mukeyHold$mukey == 2403674) # 22
	# which(mukeyHold$mukey == 2403693) # 41
	# which(mukeyHold$mukey == 2403885) # 42
	# which(mukeyHold$mukey == 2403681) # 29
	
	# profile_ksat_table = read.csv(paste(arg[1],'/soil_mukey_ksat.csv',sep=''))
	# plot(profile_ksat_table[,29+1], profileDEPTH,type='l', ylim=c(3,0), xlim=c(0,10) )
	# for(i in 2:204) lines(profile_ksat_table[,1+i], profileDEPTH)
	
	# profile_POR_table = read.csv(paste(arg[1],'/soil_mukey_por.csv',sep=''))
	# plot(profile_POR_table[,1+1], profileDEPTH,type='l', ylim=c(3,0), xlim=c(0,1) )
	# for(i in 2:204) lines(profile_POR_table[,1+i], profileDEPTH)
	
	



