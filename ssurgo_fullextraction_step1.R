mapunitTitle=c('musym','muname','mukind','mustatus','muacres','mapunitlfw_l','mapunitlfw_r','mapunitlfw_h','mapunitpfa_l','mapunitpfa_r','mapunitpfa_h','farmlndcl','muhelcl','muwathelcl','muwndhelcl','interpfocus','invesintens','iacornsr','nhiforsoigrp','nhspiagr','vtsepticsyscl','mucertstat','lkey','mukey')
compTitle=c('comppct_l','comppct_r','comppct_h','compname','compkind','majcompflag','otherph','localphase','slope_l','slope_r','slope_h','slopelenusle_l','slopelenusle_r','slopelenusle_h','runoff','tfact','wei','weg','erocl','earthcovkind1','earthcovkind2','hydricon','hydricrating','drainagecl','elev_l','elev_r','elev_h','aspectccwise','aspectrep','aspectcwise','geomdesc','albedodry_l','albedodry_r','albedodry_h','airtempa_l','airtempa_r','airtempa_h','map_l','map_r','map_h','reannualprecip_l','reannualprecip_r','reannualprecip_h','ffd_l','ffd_r','ffd_h','nirrcapcl','nirrcapscl','nirrcapunit','irrcapcl','irrcapscl','irrcapunit','cropprodindex','constreeshrubgrp','wndbrksuitgrp','rsprod_l','rsprod_r','rsprod_h','foragesuitgrpid','wlgrain','wlgrass','wlherbaceous','wlshrub','wlconiferous','wlhardwood','wlwetplant','wlshallowwat','wlrangeland','wlopenland','wlwoodland','wlwetland','soilslippot','frostact','initsub_l','initsub_r','initsub_h','totalsub_l','totalsub_r','totalsub_h','hydgrp','corcon','corsteel','taxclname','taxorder','taxsuborder','taxgrtgroup','taxsubgrp','taxpartsize','taxpartsizemod','taxceactcl','taxreaction','taxtempcl','taxmoistscl','taxtempregime','soiltaxedition','castorieindex','flecolcomnum','flhe','flphe','flsoilleachpot','flsoirunoffpot','fltemik2use','fltriumph2use','indraingrp','innitrateleachi','misoimgmtgrp','vasoimgtgrp','mukey','cokey')
chorizonTitle = c('hzname','desgndisc','desgnmaster','desgnmasterprime','desgnvert','hzdept_l','hzdept_r','hzdept_h','hzdepb_l','hzdepb_r','hzdepb_h','hzthk_l','hzthk_r','hzthk_h','fraggt10_l','fraggt10_r','fraggt10_h','frag3to10_l','frag3to10_r','frag3to10_h','sieveno4_l','sieveno4_r','sieveno4_h','sieveno10_l','sieveno10_r','sieveno10_h','sieveno40_l','sieveno40_r','sieveno40_h','sieveno200_l','sieveno200_r','sieveno200_h','sandtotal_l','sandtotal_r','sandtotal_h','sandvc_l','sandvc_r','sandvc_h','sandco_l','sandco_r','sandco_h','sandmed_l','sandmed_r','sandmed_h','sandfine_l','sandfine_r','sandfine_h','sandvf_l','sandvf_r','sandvf_h','silttotal_l','silttotal_r','silttotal_h','siltco_l','siltco_r','siltco_h','siltfine_l','siltfine_r','siltfine_h','claytotal_l','claytotal_r','claytotal_h','claysizedcarb_l','claysizedcarb_r','claysizedcarb_h','om_l','om_r','om_h','dbtenthbar_l','dbtenthbar_r','dbtenthbar_h','dbthirdbar_l','dbthirdbar_r','dbthirdbar_h','dbfifteenbar_l','dbfifteenbar_r','dbfifteenbar_h','dbovendry_l','dbovendry_r','dbovendry_h','partdensity','ksat_l','ksat_r','ksat_h','awc_l','awc_r','awc_h','wtenthbar_l','wtenthbar_r','wtenthbar_h','wthirdbar_l','wthirdbar_r','wthirdbar_h','wfifteenbar_l','wfifteenbar_r','wfifteenbar_h','wsatiated_l','wsatiated_r','wsatiated_h','lep_l','lep_r','lep_h','ll_l','ll_r','ll_h','pi_l','pi_r','pi_h','aashind_l','aashind_r','aashind_h','kwfact','kffact','caco3_l','caco3_r','caco3_h','gypsum_l','gypsum_r','gypsum_h','sar_l','sar_r','sar_h','ec_l','ec_r','ec_h','cec7_l','cec7_r','cec7_h','ecec_l','ecec_r','ecec_h','sumbases_l','sumbases_r','sumbases_h','ph1to1h2o_l','ph1to1h2o_r','ph1to1h2o_h','ph01mcacl2_l','ph01mcacl2_r','ph01mcacl2_h','freeiron_l','freeiron_r','freeiron_h','feoxalate_l','feoxalate_r','feoxalate_h','extracid_l','extracid_r','extracid_h','extral_l','extral_r','extral_h','aloxalate_l','aloxalate_r','aloxalate_h','pbray1_l','pbray1_r','pbray1_h','poxalate_l','poxalate_r','poxalate_h','ph2osoluble_l','ph2osoluble_r','ph2osoluble_h','ptotal_l','ptotal_r','ptotal_h','excavdifcl','excavdifms','cokey','chkey')
chtexgrpTitle=c('texture','stratextsflag','rvindicator','texdesc','chkey','chtgkey')
chtexturTitle=c('texcl','lieutex','chtgkey','chtkey')
chporesTitle=c('poreqty','poresize','porecont','poreshp','rvindicator','chkey','chporeskey')


arg=commandArgs(T);
target = paste(arg[1],'/tabular',sep='') 
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
weightedAverage = function(xx,weight){
    new_weight = weight;
    new_weight[is.na(xx)] = 0;
    if( sum(new_weight)>0 ){
        return <- sum(xx * new_weight/sum(new_weight), na.rm=T)
    }else{
        return <- NA
    }
}#function
getDistMatrix = function(x){
    len = dim(x)[1]
    len_index = seq_len(len)
    len_index2 = seq_len(dim(x)[2])
    dm = matrix(NA,len,len)
    for(i in len_index){
        for(j in len_index){
            dm[i,j] = sqrt(sum(sapply(len_index2,function(ii){
                return <- (x[i,ii]-x[j,ii])^2
            }),na.rm=T ))
        }#j
    }#i
    return <- dm
}#function

horizonOrder = c('H','O','A','E','B','C','R')	# there is horizon 'H' H horizons or layers: Layers dominated by organic material, formed from accumulations of undecomposed or partially decomposed organic material at the soil surface which may be underwater. All H horizons are saturated with water for prolonged periods or were once saturated but are now artificially drained. An H horizon may be on top of mineral soils or at any depth beneath the surface if it is buried. (http://www.fao.org/3/w8594e/w8594e0g.htm)
    # H - B microbial activities, but no evidence of such activities in C or deeper
    # plant roots can penetrate C horizon
    # R is hard bedrock
    # ---------------------------------------------
	# make horizon H as O
horizonOrderScore = c(1,1,1,1,2,3,4)
horizonThicknessDefault = c(0, 2, 5, 5, 10, 10, 10)*0.0254 #inches to meter; https://www.sheffield.ac.uk/ssa/soil-facts/horizons
horizonThicknessDefault_acc = cumsum(horizonThicknessDefault)
horizonThicknessDefault_scale2soildepth = max(horizonThicknessDefault_acc)/horizonThicknessDefault_acc
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
#mukeyHold = matrix(NA,length(mukey), 17);
#profileDEPTH = seq(0,10,0.01)
#profile_ksat_table = data.frame(z=profileDEPTH)
#profile_POR_table = data.frame(z=profileDEPTH)

###---------------------------------------------------------- extracting
textureFILE = paste(arg[1],'/soil_mukey_texture.csv',sep='')
for(i in seq_along(mukey)){
	## for each mukey, we have N cokey
    ## mukey soil texture class
    
    
	cond1 = comp[,'mukey']== mukey[i]; sum(cond1)
	
	## 1 find overall soil type and develop default parameter values; not component/horizon specific
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
	
	
    
    
	## 2 find all the components and horizons and extract specific information
	if(! overall_NODATA){
		
        ## 2.1 reoganizing information by horizon
        mukeyHorizonLIST = list()
        for(cokeyii in seq_along(comp[cond1,'cokey'])){
        
            cokey = comp[cond1,'cokey'][cokeyii]
            cond2 = chorizon[,'cokey']== cokey; sum(cond2) # finding component-specific horizons
            if(sum(cond2)>0){
                # found matched component-specific horizons
                component_NODATA = F
                
                ## 2.1.1 ..... default
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
                
                
                ## 2.1.2 ... SSURGO specific
                horizonSummaryTable = data.frame(id = chorizon[cond2,'chkey'])
                horizonSummaryTable$horizonThinkness = NA#(chorizon[cond2,'hzdepb_r']-chorizon[cond2,'hzdept_r'])*0.01 #cm -> m
                horizonSummaryTable$horizonMeanDepth = NA#0.5*(chorizon[cond2,'hzdepb_r']+chorizon[cond2,'hzdept_r'])*0.01 # m (real depth)
                horizonSummaryTable$horizonTopDepth = chorizon[cond2,'hzdept_r']*0.01 #m
                horizonSummaryTable$horizonBottomDepth = chorizon[cond2,'hzdepb_r']*0.01 #m <<-- for microbial activities depth & max root depth
                # microbial activities stop at the bottom of B
                # max plant root could be at the bottom of C
                #horizonSummaryTable$horizonOriginal = chorizon[cond2,'hzname']
                horizonSummaryTable$horizonClass = sapply(chorizon[cond2,'hzname'],function(nn){
                    tmp = sapply(horizonOrder, function(x){grepl(x,nn)})
                    if(sum(tmp)>0) return <- horizonOrder[tmp]
                    else return <- NA
                })#
                horizonSummaryTable$horizonName = chorizon[cond2,'hzname']
                
                horizonSummaryTable$sand = replaceNA0(chorizon[cond2,'sandtotal_r']*0.01, default_soil['sand']); #% -> [0,1]
                horizonSummaryTable$silt = replaceNA0(chorizon[cond2,'silttotal_r']*0.01, default_soil['silt']); #% -> [0,1]
                horizonSummaryTable$clay = replaceNA0(chorizon[cond2,'claytotal_r']*0.01, default_soil['clay']); #% -> [0,1]
                horizonSummaryTable$ksat = replaceNA0(chorizon[cond2,'ksat_r']*1e-6*3600*24, default_soil['ksat']); # m/day;
                horizonSummaryTable$ksat[horizonSummaryTable$ksat<=0] = 0.001 # m/day
                horizonSummaryTable$POR = replaceNA0(chorizon[cond2,'wsatiated_r']*0.01, default_soil['POR']); #% -> [0,1]
                horizonSummaryTable$POR[horizonSummaryTable$POR<=0] = 0.01
                
                horizonSummaryTable$BD = chorizon[cond2,'dbovendry_r'] #g/cm3 bulk density [0.02 - 2.6]; the oven dry weight of the less than 2 mm soil material per unit volume of soil exclusive of the desication cracks, measured on a costed clod;
                horizonSummaryTable$fc = chorizon[cond2,'wthirdbar_r']*0.01 #% -> [0,1]; field capacity; the volumetric content of soil water retained at the tension of 1/3 bar (33 kPa), expressed as a percentage of the whole soil.
                horizonSummaryTable$awc = chorizon[cond2,'awc_r']*0.01 # cm/cm; Amount of water that an increment of soil depth, inclusive of fragements, can store that is available to plants.
                horizonSummaryTable$kffact = chorizon[cond2,'kffact'] # soil erodibility factor [0-1]; an erodibility factor which quantifies the susceptibility of doil particles to detechment by water
                horizonSummaryTable$density = horizonSummaryTable$BD / (1 - horizonSummaryTable$POR)
                #chorizon[cond2,'partdensity']
                # g/cm3; [0.01 - 5] mass per unit volume (not including pore space) of the solid particle either mineral or organic.
                # from text book: % pore space = porosity = (1 - BD / density) * 100 %
                # from my note: density = BD / (1 - POR)
                horizonSummaryTable$om = chorizon[cond2,'om_r']*0.01
                #% -> [0,1]; the amount by weight of decomposed plant and animal residue expressed as weight percentage of the less than 2 mm soil material (not the top 2 mm soil!).
                
                horizonSummaryTable$compWeight = comp[cond1,'comppct_r'][cokeyii]
                
                # try to find horizon specific soil class, but it return multiple!!
                # why there is multiple soil class in one horizon?
                horizonSoilCl_varx = sapply(chorizon[cond2,'chkey'],function(hh){
                    horizonSoilCl = chtextur[ (chtextur[,'chtgkey'] %in% chtexgrp[ (chtexgrp[,'chkey'] == hh),'chtgkey']),'texcl']
                    horizonSoilCl_index = match(horizonSoilCl, soilscoreNames)
                    horizonSoilCl_count = table(horizonSoilCl_index)
                    horizonSoilCl_weight = horizonSoilCl_count/sum(horizonSoilCl_count)
                    horizonSoilCl_index = as.numeric(names(horizonSoilCl_count))
                    return <- c(
                        weightedAverage(soilscore_db$por_size_index[horizonSoilCl_index],horizonSoilCl_weight),
                        weightedAverage(soilscore_db$psi_air_entry[horizonSoilCl_index],horizonSoilCl_weight),
                        weightedAverage(soilscore_db$albedo[horizonSoilCl_index],horizonSoilCl_weight))
                })#
                horizonSummaryTable$por_size_index = horizonSoilCl_varx[1,]
                horizonSummaryTable$psi_air_entry = horizonSoilCl_varx[2,]
                horizonSummaryTable$albedo = horizonSoilCl_varx[3,]
       
       
                ## ... ## ... 2.1.2.1 resolve horizon names :: analyzing
                # a profile containing a buried sequence could be structured O, A1, A2, B2, 2A2, 2B21, 2B22, 2C with the buried profile commencing at 2A2.
                # take another approach for this: keep the seq and find break
                workingPattern_order = order( horizonSummaryTable$horizonTopDepth + horizonSummaryTable$horizonBottomDepth )
                horizonSummaryTable = horizonSummaryTable[workingPattern_order,]
                horizonSummaryTable$horizonIndex = sapply(horizonSummaryTable$horizonClass,function(x){ horizonOrderScore[match(x,horizonOrder)] })
                
                excludedR = sapply(horizonSummaryTable$horizonClass,function(x){ sum(grepl('R',x))==0 })
                if(sum(excludedR)>2){
                    ## 1) sort the horizon by depths
                    ## 2) determine to horizon class --> output A, B, C, R
                    ## 3) how to determine? ksat? POR? om? sand, slit, clay? letter classes
                    ## horizonSummaryTable[,c('horizonClass','horizonName','sand','silt','clay','ksat','POR','om')]
                    ##
                    #  horizonClass horizonName sand silt clay   ksat  POR    om
                    #             A           A 0.44 0.41 0.15 0.7776 0.50 0.100
                    #             B          Bw 0.47 0.37 0.16 0.7776 0.33 0.005
                    #             C          BC 0.67 0.20 0.13 2.4192 0.30 0.002
                    ## how to form dist matrix?
                    horizonPTs = cbind(
                        numclass=sapply(horizonSummaryTable$horizonClass[excludedR],length),
                        indexclass=sapply(horizonSummaryTable$horizonIndex[excludedR],mean),
                        depth=0.5*(horizonSummaryTable$horizonTopDepth[excludedR]+horizonSummaryTable$horizonBottomDepth[excludedR])/max(horizonSummaryTable$horizonBottomDepth),
                        horizonSummaryTable[excludedR,c('sand','silt','clay','POR','om')])#cbind 'ksat',
                    distmatrix = getDistMatrix( horizonPTs )
                    
                    horizonLumpList = list()
                    currentJJ = 1;
                    for(jj in 1:2){
                        JJseq = currentJJ:dim(distmatrix)[1]
                        distinguishingDist = min(1.5, max(distmatrix[JJseq,currentJJ])/(4-jj) )
                         selectedIndex = JJseq[ distmatrix[JJseq,currentJJ]<distinguishingDist ]
                         selectedIndex_cond = sapply(seq_along(selectedIndex),function(yy){ length(selectedIndex[1]:selectedIndex[yy]) == yy })
                         
                         horizonLumpList[[jj]] = selectedIndex[selectedIndex_cond]
                         currentJJ = 1+max(selectedIndex[selectedIndex_cond])
                         if(currentJJ>=dim(distmatrix)[1]){ break; }
                         ## problem: two horizons; jj=1 pass; then what is the second horizon? B or C?
                    }#jj for
                    if(currentJJ<=dim(distmatrix)[1]){ horizonLumpList[[3]] = currentJJ:dim(distmatrix)[1] }
                    
                    outputHorizon = c('A','B','C')
                    horizonSummaryTable$horizonClass = sapply(horizonSummaryTable$horizonClass,function(x){x[1]})
                    horizonSummaryTable$horizonIndex = sapply(horizonSummaryTable$horizonIndex,function(x){x[1]})
                    for(jj in seq_along(horizonLumpList)){
                        horizonSummaryTable$horizonClass[ horizonLumpList[[jj]] ] = outputHorizon[jj]
                        horizonSummaryTable$horizonIndex[ horizonLumpList[[jj]] ] = horizonOrderScore[match(outputHorizon[jj],horizonOrder)]
                    }#jj
                }else{
                    horizonSummaryTable$horizonClass = sapply(horizonSummaryTable$horizonName,function(nn){
                        tmp = sapply(horizonOrder, function(x){grepl(x,substr(nn,1,1) )})
                        if(sum(tmp)>0) return <- horizonOrder[tmp]
                        else return <- NA
                    })#
                    # lumping
                    horizonSummaryTable$horizonClass[horizonSummaryTable$horizonClass=='H'] = 'A'
                    horizonSummaryTable$horizonClass[horizonSummaryTable$horizonClass=='O'] = 'A'
                    horizonSummaryTable$horizonClass[horizonSummaryTable$horizonClass=='E'] = 'A'
                    
                    horizonSummaryTable$horizonIndex = sapply(horizonSummaryTable$horizonClass,function(x){ horizonOrderScore[match(x,horizonOrder)] })
                    
                }# if(dim(horizonSummaryTable)[1]>1)
                
                
                ## ... ## ... 2.1.2.2 sorting horizonSummaryTable by horizonName
                horizonSummaryTable$horizonThinkness = sapply(seq_len(dim(horizonSummaryTable)[1]),function(jj){
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
                
    
                ## ... ## ... 2.1.2.3 convert OM from % to kgC/m2 (total decaying organic)
                # the amount by weight of decomposed plant and animal residue expressed as weight percentage of the less than 2 mm soil material (not the top 2 mm soil!).
                horizonSummaryTable$om = 500*horizonSummaryTable$BD*chorizon[cond2,'om_r']*0.01*horizonSummaryTable$horizonThinkness # kgC/m2
                # g/cm3 = 100*100*100 g/m3 = 100*100*100/1000 kg/m3 => 500 kgC/m3
                # percent * PD * thickness (kg/m3 * m) -> kgC/m2
                # 0.5 kg/m3 = kgC/m3
                # we need to count for the POR in the actual soil core, so we use "BD" instead of "PD"
                
                
                
                ## ... ## ... 2.1.2.4 how to unpack "horizonSummaryTable" to "mukeyHorizonLIST" ?
                ### <<----- lumping O and E to A for simplicity
                #horizonSummaryTable$horizonName[horizonSummaryTable$horizonClass=='H']= 'A'
                #horizonSummaryTable$horizonName[horizonSummaryTable$horizonClass=='O']= 'A'
                #horizonSummaryTable$horizonName[horizonSummaryTable$horizonClass=='E']= 'A'
                horizonSummaryTable_class = match(horizonSummaryTable$horizonClass,horizonOrder)
                horizonSummaryTable_class_index = tapply(seq_along(horizonSummaryTable_class),horizonSummaryTable_class,function(ii){return <- ii})
                
                for(iij in seq_along(horizonSummaryTable_class_index) ){
                    
                    tmp = horizonSummaryTable[horizonSummaryTable_class_index[[iij]],]
                    tmp_name = tmp$horizonClass[1]
                    if(dim(tmp)[1]>1){
                        ## merging repeated horizon
                        temp_weight = tmp$horizonThinkness/sum(tmp$horizonThinkness)
                        
                        finalHorizon = tmp[1,]
                        ## somehow it creates new column below!!!
                        finalHorizon$horizonThinkness = sum(tmp$horizonThinkness, na.rm=T)
                        finalHorizon$horizonTopDepth = min(tmp$horizonTopDepth, na.rm=T)
                        finalHorizon$horizonBottomDepth = max(tmp$horizonBottomDepth, na.rm=T)
                        finalHorizon$horizonMeanDepth = 0.5*(finalHorizon$horizonTopDepth + finalHorizon$horizonBottomDepth)
                        
                        finalHorizon$sand = weightedAverage(tmp$sand, temp_weight)
                        finalHorizon$silt = weightedAverage(tmp$silt, temp_weight)
                        finalHorizon$clay = weightedAverage(tmp$clay, temp_weight)
                        finalHorizon$ksat = weightedAverage(tmp$ksat, temp_weight)
                        finalHorizon$POR = weightedAverage(tmp$POR, temp_weight); #dim(finalHorizon)
                        finalHorizon$BD = weightedAverage(tmp$BD, temp_weight); #dim(finalHorizon)
                        finalHorizon$fc = weightedAverage(tmp$fc, temp_weight); #dim(finalHorizon)
                        finalHorizon$awc = weightedAverage(tmp$awc, temp_weight); #dim(finalHorizon)
                        finalHorizon$kffact = weightedAverage(tmp$kffact, temp_weight); #dim(finalHorizon)
                        finalHorizon$om = weightedAverage(tmp$om, temp_weight); #dim(finalHorizon)
                        
                        finalHorizon$density = finalHorizon$BD / (1 - finalHorizon$POR); #dim(finalHorizon)
                    }else{
                        finalHorizon = tmp
                    }
                    
                    if( is.null(mukeyHorizonLIST[[ tmp_name ]]) ){
                        # need to make data.frame
                        # list$A = table [soil variable]
                        # colname [cokey]
                        mukeyHorizonLIST[[ tmp_name ]] = finalHorizon
                    }else{
                        # add to data.frame
                        mukeyHorizonLIST[[ tmp_name ]] = rbind(mukeyHorizonLIST[[ tmp_name ]], finalHorizon)
                    }
                }# for loop iij
                
            
            }else{
                ## 2.1.2
                # component_NODATA = T :: no information about this component at all :: no info on default_overallsoil
                # do something?
            }# else of sum(cond2)>0
        }# for loop of cokey
        # mukeyHorizonLIST
        
        ## 3 mukeyHorizon Table (row = mukey, column = variation
        # horizonOrder = c('O','A','E','B','C','R')
       	mukeyHorizonLIST_order = order(sapply(seq_along(mukeyHorizonLIST),function(ii){
            return <- mukeyHorizonLIST[[ii]]$horizonIndex[1]
        }))
        aggregatedHorian_name = sapply(mukeyHorizonLIST_order,function(ii){
            return <- mukeyHorizonLIST[[ii]]$horizonClass[1]
        })
        aggregatedHorian_num = sapply(mukeyHorizonLIST_order,function(ii){
        		length(mukeyHorizonLIST[[ii]]$horizonClass)
        	})
        aggregatedHorian_table = as.data.frame(do.call(rbind,lapply(mukeyHorizonLIST_order,function(ii){
            if(sum(mukeyHorizonLIST[[ii]]$compWeight)==0){
            		weight = rep(1,length(mukeyHorizonLIST[[ii]]$compWeight))
            		weight = weight/sum(weight)
            }else{
            		weight = mukeyHorizonLIST[[ii]]$compWeight/sum(mukeyHorizonLIST[[ii]]$compWeight);
            }
            return <- c(
            horizonThinkness = weightedAverage(mukeyHorizonLIST[[ii]]$horizonThinkness,weight),
            horizonMeanDepth = weightedAverage(mukeyHorizonLIST[[ii]]$horizonMeanDepth,weight),
            horizonTopDepth = weightedAverage(mukeyHorizonLIST[[ii]]$horizonTopDepth,weight),
            horizonBottomDepth = weightedAverage(mukeyHorizonLIST[[ii]]$horizonBottomDepth,weight),
            sand = weightedAverage(mukeyHorizonLIST[[ii]]$sand,weight),
            silt = weightedAverage(mukeyHorizonLIST[[ii]]$silt,weight),
            clay = weightedAverage(mukeyHorizonLIST[[ii]]$clay,weight),
            ksat = weightedAverage(mukeyHorizonLIST[[ii]]$ksat,weight),
            POR = weightedAverage(mukeyHorizonLIST[[ii]]$POR,weight),
            
            BD = weightedAverage(mukeyHorizonLIST[[ii]]$BD,weight), # one of it is NA
            om = weightedAverage(mukeyHorizonLIST[[ii]]$om,weight),
            density = weightedAverage(mukeyHorizonLIST[[ii]]$BD,weight) / (1 - weightedAverage(mukeyHorizonLIST[[ii]]$POR,weight)),
            
            por_size_index = weightedAverage(mukeyHorizonLIST[[ii]]$por_size_index,weight),
            psi_air_entry = weightedAverage(mukeyHorizonLIST[[ii]]$psi_air_entry,weight),
            albedo = weightedAverage(mukeyHorizonLIST[[ii]]$albedo,weight)
            )
        })))
        ## ... adjust horizon depths after averaging
        len = dim(aggregatedHorian_table)[1]
        if(len>1){
        
            issuedepth = aggregatedHorian_table$horizonTopDepth != c(0,aggregatedHorian_table$horizonBottomDepth[1:(len-1)]) |
                c(F,diff(aggregatedHorian_table$horizonTopDepth)<0) |
                c(F,diff(c(0,aggregatedHorian_table$horizonBottomDepth[1:(len-1)]))<0)
            #up_adjust = c(F,aggregatedHorian_num[1:(len-1)] < aggregatedHorian_num[2:len]) # make adjust to the up horizon
            #down_adjust = c(aggregatedHorian_num[1:(len-1)] >= aggregatedHorian_num[2:len],F)
	        
            # cbind(aggregatedHorian_table$horizonTopDepth, aggregatedHorian_table$horizonBottomDepth, issuedepth, aggregatedHorian_table$horizonThinkness)
            
            for(jj in seq_along(issuedepth)){
                
                up_adjust = ifelse(jj==1, F, ifelse(jj==len,F, aggregatedHorian_num[jj] < aggregatedHorian_num[jj+1]) )
                
                if(issuedepth[jj] & !up_adjust){
                    # because cannot adjust up, we adjust current
                    newtop = ifelse(jj-1>=1, aggregatedHorian_table$horizonBottomDepth[jj-1], aggregatedHorian_table$horizonTopDepth[1])
                }else{ newtop = aggregatedHorian_table$horizonTopDepth[jj] }
                    
                if(issuedepth[jj] & up_adjust){
                    # because cannot adjust down, we adjust current
                    newbottom = ifelse(jj+1<=len, aggregatedHorian_table$horizonTopDepth[jj+1], aggregatedHorian_table$horizonBottomDepth[len])
                }else{ newbottom = aggregatedHorian_table$horizonBottomDepth[jj] }
                    
                    
                if(newtop > newbottom) newbottom = newtop + aggregatedHorian_table$horizonThinkness[jj]
                
                aggregatedHorian_table$horizonTopDepth[jj] = newtop
                aggregatedHorian_table$horizonBottomDepth[jj] = newbottom
            }# jj for loop
            
        }# if(len>1)
        # cbind(aggregatedHorian_name, aggregatedHorian_num, aggregatedHorian_table)
        # cbind(aggregatedHorian_name, aggregatedHorian_num, aggregatedHorian_table[,c('horizonTopDepth','horizonBottomDepth')])
        # cbind(aggregatedHorian_name, aggregatedHorian_table[,c('horizonTopDepth','horizonBottomDepth')])
        
        ## final check on thickness and top-bottom depths
        #tmp_thickness = aggregatedHorian_table$horizonBottomDepth - aggregatedHorian_table$horizonTopDepth
        #cbind(aggregatedHorian_table$horizonThinkness,tmp_thickness)
        if(len>1) aggregatedHorian_table$horizonThinkness[1:(len-1)] = aggregatedHorian_table$horizonBottomDepth[1:(len-1)] - aggregatedHorian_table$horizonTopDepth[1:(len-1)]
        
        aggregatedHorian_table$horizonBottomDepth[len] = aggregatedHorian_table$horizonTopDepth[len] +aggregatedHorian_table$horizonThinkness[len]
        
        ## 4a output aggregated horizon information, rhessys_soilid
        if(i==1){
        		write( paste(c('mukey','texture',do.call(c,lapply(horizonOrder[c(2,4,5,6)+1],function(x){ paste(x, colnames(aggregatedHorian_table), sep="_") }))),collapse=','), textureFILE)
        	}# if header
        write( paste(c(mukey[i],soilID_,unlist(lapply(horizonOrder[c(2,4,5,6)+1],function(x){
        		index_ = match(x,aggregatedHorian_name)
        		if(is.na(index_)){
        			return <- rep(NA,15)
        		}else{
        			return <- aggregatedHorian_table[index_,]
        		}
        	}))),collapse=','), textureFILE,append=T)
        
               
       
        
	}else{
        # overall_NODATA=T case (absolutely no data at all; do not know the default soil neither )
        
        ## 4a output aggregated horizon information
        if(i==1){
       		write( paste(c('mukey','texture',do.call(c,lapply(horizonOrder[c(2,4,5,6)+1],function(x){ paste(x, colnames(aggregatedHorian_table), sep="_") }))),collapse=','), textureFILE)
        }# if header
        write( paste(c(mukey[i],soilID_,rep(NA,60)), collapse=','), textureFILE, append=T) #90
        


	}# end of if
}# for loop i

