import os,sys,scipy,numpy,pandas,geopandas,pygeohydro,ee,geemap,us,random,ACPTrajectoryNamelist,folium

class ACPDistributionType(scipy.stats.rv_continuous):
    """Generic customizable distribution class - to be used to define nutrient yield and/or effectiveness coefficient distributions"""
    def __init__(self, percentiles, values):
        super().__init__()
        self.percentiles = percentiles
        self.values = values
        self.interp_func = scipy.interpolate.interp1d(percentiles, values, kind='linear', fill_value="extrapolate")
    def _cdf(self, x):
        return numpy.interp(x, self.values, self.percentiles)
    def rvs(self, size=1, random_state=None):
        return self.interp_func(numpy.random.rand(size))
        
class ACPData:
    """Class to hold the data"""
    
    class SpatialData:
        """Spatial data"""
        domain                         = numpy.nan
        domain_huc12s                  = numpy.nan # empty unless overwrite = True
        domain_fields                  = numpy.nan
        domain_counties                = numpy.nan
        
    class NutrientData:
        """Baseline nutrient data"""
        baseline_yields_kgha_ACRE      = numpy.nan
        baseline_distributions_SWAT    = numpy.nan
        swatsims                       = numpy.nan

    class EfficiencyCoefficientData:
        """Effectiveness coefficient (EC) data"""
        ecs_ACRE                       = numpy.nan
        ecs_ACRE_simple                = numpy.nan
        ec_distributions_ACRE          = numpy.nan

    class SWATOutputVariables:
        """SWAT output variable names"""
        n_vars                         = ('NSURQ','ORGN','NLAT','GWQN')
        p_vars                         = ('ORGP','SOLP')
        
    class MiscStats:
        """Miscellaneous"""
        ## need to do - better places /classifications for these things? Do we even need?
        domain_area_ha                 = numpy.nan
        dt_fips_key_fieldid            = numpy.nan
        dt_areaha_key_fieldid          = numpy.nan
        acp_types                      = numpy.nan
        hrugiss_agonly                 = numpy.nan # used?
        field_ids                      = numpy.nan

    def __init__(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Initialize ACPData"""
        if acpnamelist.vars.verbose: print('Initializing ACP data class')
        self._init_vars(acpnamelist)
        self._init_dirs(acpnamelist)
        self._read_spatial_data(acpnamelist)
        self._read_nutrient_and_ec_data(acpnamelist)

    def _init_vars(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Initialize class variables"""
        if acpnamelist.vars.verbose: print('Initializing variables')
        self.miscstats    = ACPData.MiscStats()
        self.spatialdata  = ACPData.SpatialData()
        self.nutrientdata = ACPData.NutrientData()
        self.ecs          = ACPData.EfficiencyCoefficientData()
        self.swatoutvars  = ACPData.SWATOutputVariables()
        
    def _init_dirs(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Make project directories"""
        if acpnamelist.vars.verbose: print('Making project directories')
        if not os.path.isdir(acpnamelist.dirnames.project): 
            os.mkdir(acpnamelist.dirnames.project)
        if not os.path.isdir(acpnamelist.dirnames.county_boundary): 
            os.mkdir(acpnamelist.dirnames.county_boundary)
        if not os.path.isdir(acpnamelist.dirnames.field_boundary): 
            os.mkdir(acpnamelist.dirnames.field_boundary)
        if not os.path.isdir(acpnamelist.dirnames.huc_boundary): 
            os.mkdir(acpnamelist.dirnames.huc_boundary)
        if not os.path.isdir(acpnamelist.dirnames.swat_sims): 
            os.mkdir(acpnamelist.dirnames.swat_sims)
        if not os.path.isdir(acpnamelist.dirnames.output): 
            os.mkdir(acpnamelist.dirnames.output)
        
    def _read_spatial_data(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Make or read spatial data"""
        if acpnamelist.vars.verbose: print('Making or reading project spatial data')
        self._setDomain(acpnamelist)
        self._setCounties(acpnamelist)
        self._setFieldBoundaries(acpnamelist)
        self._setHRUBoundaries(acpnamelist)
        self._evalHRUFieldOverlap(acpnamelist)

    def _read_nutrient_and_ec_data(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Read baseline nutrient data and ACP effectiveness coefficients (ECs)"""
        if acpnamelist.vars.verbose: print('Reading baseline nutrient data and ACP effectiveness coefficients')
        self._readBaselineNutrientYields_ACRE(acpnamelist)
        self._readSWATSimulations(acpnamelist)
        self._setBaselineNutrientYieldDistributions_SWAT(acpnamelist)
        self._readECs_ACRE(acpnamelist)
        self._setECDistributions_ACRE(acpnamelist)
        
    def _setDomain(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Make or read domain spatial boundary"""
        if acpnamelist.vars.verbose: print('Making or reading domain spatial boundary')
        if os.path.isfile(acpnamelist.fnames.domain) and not acpnamelist.vars.overwrite_flag:
            if acpnamelist.vars.verbose: print('    Reading: '+acpnamelist.fnames.domain)
            self.spatialdata.domain = geopandas.read_file(acpnamelist.fnames.domain)
            self.spatialdata.domain[acpnamelist.varnames.area_ha] = self.spatialdata.domain.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.miscstats.domain_area_ha = float(self.spatialdata.domain[acpnamelist.varnames.area_ha].item())
        else:
            if acpnamelist.vars.verbose: print('    Making: '+acpnamelist.fnames.domain)
            self.spatialdata.domain_huc12s = pygeohydro.WBD("huc12").byids("huc12", acpnamelist.vars.input_huc12_ids)
            self.spatialdata.domain_huc12s.to_file(acpnamelist.fnames.hucs,driver='GPKG')
            self.spatialdata.domain_huc12s = geopandas.read_file(acpnamelist.fnames.hucs)
            self.spatialdata.domain = self.spatialdata.domain_huc12s.dissolve()
            self.spatialdata.domain.to_file(acpnamelist.fnames.domain, driver="GPKG")
            self.spatialdata.domain[acpnamelist.varnames.area_ha] = self.spatialdata.domain.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.miscstats.domain_area_ha = float(self.spatialdata.domain[acpnamelist.varnames.area_ha].item())
            
    def _setCounties(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Set/read US county boundaries into gdf"""
        if acpnamelist.vars.verbose: print('Making or reading US county spatial boundary')
        if os.path.isfile(acpnamelist.fnames.counties) and not acpnamelist.vars.overwrite_flag:
            if acpnamelist.vars.verbose: print('    Reading: '+acpnamelist.fnames.counties)
            self.spatialdata.domain_counties = geopandas.read_file(acpnamelist.fnames.counties)
        else:
            if acpnamelist.vars.verbose: print('    Making: '+acpnamelist.fnames.counties)
            ee.Authenticate()
            ee.Initialize(project='simmgmttrajectory')
            us_counties_fc  = ee.FeatureCollection("TIGER/2016/Counties")
            us_counties = geemap.ee_to_gdf(us_counties_fc)
            us_counties = us_counties[us_counties.STATEFP != '11'] # remove wash dc
            us_counties['STATE_NAME'] = [str(us.states.lookup(r['STATEFP']).name).title() for _,r in us_counties.iterrows()]
            domain_county_intersect = geopandas.overlay(self.spatialdata.domain, us_counties, how='intersection')
            self.spatialdata.domain_counties = us_counties[us_counties["GEOID"].isin(list(set(domain_county_intersect['GEOID'])))]
            self.spatialdata.domain_counties.to_file(acpnamelist.fnames.counties)
            
    def _setFieldBoundaries(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Set/read field boundaries into gdf"""
        if acpnamelist.vars.verbose: print('Processing or reading field spatial boundary')
        if os.path.isfile(acpnamelist.fnames.fields) and not acpnamelist.vars.overwrite_flag:
            if acpnamelist.vars.verbose: print('    Processing: '+acpnamelist.fnames.fields)
            self.spatialdata.domain_fields = geopandas.read_file(acpnamelist.fnames.fields)
            self.miscstats.dt_fips_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields['FIPS']))
            self.miscstats.dt_areaha_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields[acpnamelist.varnames.area_ha]))
        else:
            if acpnamelist.vars.verbose: print('    Making: '+acpnamelist.fnames.fields)
            fields_dfs = list()
            for _,r in self.spatialdata.domain_counties.iterrows():
                dirstate = os.path.join(acpnamelist.vars.input_fields,str(us.states.lookup(str(r['STATE_NAME']).title()).abbr).lower())
                fname = os.path.join(dirstate,str(r['NAME']).title()+'.shp')
                fields_dfs.append(geopandas.read_file(fname))
            fields = pandas.concat(fields_dfs)
            self.spatialdata.domain_fields = geopandas.overlay(self.spatialdata.domain, fields, how='intersection')
            self.spatialdata.domain_fields[acpnamelist.varnames.area_ha] = self.spatialdata.domain_fields.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.spatialdata.domain_fields.rename(columns={"name": "county_name"},inplace=True)
            self.spatialdata.domain_fields = geopandas.sjoin(left_df=self.spatialdata.domain_fields,right_df=self.spatialdata.domain_counties,how='inner')
            self.spatialdata.domain_fields['field_id'] = self.spatialdata.domain_fields.index
            self.miscstats.dt_fips_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields['FIPS']))
            self.miscstats.dt_areaha_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields[acpnamelist.varnames.area_ha]))
            self.spatialdata.domain_fields.to_file(acpnamelist.fnames.fields,driver='GPKG')
        self.miscstats.field_ids = tuple(self.miscstats.dt_areaha_key_fieldid.keys())
    
    def _setHRUBoundaries(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Set/read HRU boundaries into gdf"""
        if acpnamelist.vars.verbose: print('Making or reading HRU spatial boundary')
        if os.path.isfile(acpnamelist.fnames.hrus) and not acpnamelist.vars.overwrite_flag:
            if acpnamelist.vars.verbose: print('    Reading: '+acpnamelist.fnames.hrus)
            self.spatialdata.domain_hrus = geopandas.read_file(acpnamelist.fnames.hrus)
        else:
            if acpnamelist.vars.verbose: print('    Making: '+acpnamelist.fnames.hrus)
            hru_dfs = list()
            for fname_swathru in acpnamelist.vars.input_hru_shps:
                df = geopandas.read_file(fname_swathru)
                modelid = "ERROR"
                if fname_swathru.find('LEFW') != -1: modelid = 'LEFW'
                if fname_swathru.find('UEFW') != -1: modelid = 'UEFW'
                df['modelid'] = modelid
                hru_dfs.append(df)
            hrus = pandas.concat(hru_dfs)
            self.spatialdata.domain_hrus = geopandas.overlay(self.spatialdata.domain,hrus.to_crs(self.spatialdata.domain.crs.to_wkt()),how='intersection')
            self.spatialdata.domain_hrus[acpnamelist.varnames.area_ha] = self.spatialdata.domain_hrus.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.spatialdata.domain_hrus.drop(['objectid','shape_Length','Shape_Area'],axis=1,inplace=True)
            self.spatialdata.domain_hrus['HRUGIS'] = self.spatialdata.domain_hrus['modelid'] + self.spatialdata.domain_hrus['HRUGIS']
            self.spatialdata.domain_hrus.to_file(acpnamelist.fnames.hrus,driver='GPKG')
        
    def _evalHRUFieldOverlap(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Determine the areal-dominant HRU for each field"""
        if acpnamelist.vars.verbose: print('Reading or assessing areal dominant HRU for each field')
        if os.path.isfile(acpnamelist.fnames.field_hrugis) and not acpnamelist.vars.overwrite_flag:
            if acpnamelist.vars.verbose: print('    Reading: '+acpnamelist.fnames.field_hrugis)
            ls = list(open(acpnamelist.fnames.field_hrugis,'r'))
            ls = [l.replace('\n','').split(',') for l in ls] 
            self.miscstats.dtHRUGISKeyFieldID = {int(ls[i][0]):str(ls[i][1]) for i in range(1,len(ls))}
        else:
            if acpnamelist.vars.verbose: print('    Making: '+acpnamelist.fnames.field_hrugis)
            domain_fields_overlay = geopandas.overlay(self.spatialdata.domain_fields, self.spatialdata.domain_hrus.to_crs(self.spatialdata.domain_fields.crs), how='intersection')
            max_area_rows = domain_fields_overlay.groupby("field_id")["Area"].idxmax()
            domain_fields_overlay = domain_fields_overlay.loc[max_area_rows]
            self.miscstats.dtHRUGISKeyFieldID = dict(zip(domain_fields_overlay['field_id'],domain_fields_overlay['HRUGIS']))
            with open(acpnamelist.fnames.field_hrugis,'w') as openfile:
                openfile.write('field_id,hrugis\n')
                for field_id in self.miscstats.dtHRUGISKeyFieldID:
                    openfile.write(str(field_id)+','+str(self.miscstats.dtHRUGISKeyFieldID[field_id])+'\n')
    
    def _readBaselineNutrientYields_ACRE(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Read mean annual simulations from ACRE inputs"""
        if acpnamelist.vars.verbose: print('Reading baseline nutrient yields from ACRE data')
        csv_fips = os.path.join(acpnamelist.dirnames.acre_data,'FIPS')+'.csv'
        _, file_ext = os.path.splitext(csv_fips)
        if file_ext.find('csv') != -1: df_fips = pandas.read_csv(filepath_or_buffer=csv_fips)
        dfs = list()
        for i,r in self.spatialdata.domain_counties.iterrows():
            dir_state = os.path.join(acpnamelist.dirnames.acre_data,str(r['STATE_NAME']).title())
            csv_file_county = os.path.join(dir_state,str(r['NAME']).title())+'.csv'
            if acpnamelist.vars.verbose: print('    Reading: '+csv_file_county)
            _, file_ext = os.path.splitext(os.path.join(acpnamelist.dirnames.acre_data,dir_state,csv_file_county))
            if file_ext.find('csv') != -1:
                df = pandas.read_csv(filepath_or_buffer=csv_file_county)
                df['state'] = str(r['STATE_NAME']).title()
                df['county'] = str(r['NAME']).title()
                dfs.append(df)
        self.nutrientdata.baseline_yields_kgha_ACRE = pandas.concat(dfs)
        self.nutrientdata.baseline_yields_kgha_ACRE = pandas.merge(left=self.nutrientdata.baseline_yields_kgha_ACRE,right=df_fips,left_on=['state','county'],right_on=['STATE_NAME','NAME'])
        self.nutrientdata.baseline_yields_kgha_ACRE = self.nutrientdata.baseline_yields_kgha_ACRE.drop(['STATE_NAME', 'NAME'], axis=1)

    def _readSWATSimulations(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist): 
        """Read mean annual simulations from input swat output.hru files"""
        ## TO DO - refactor to exclude study area specific code - i.e., LEFW, UEFW, etc.
        if acpnamelist.vars.verbose: print('Reading baseline nutrient yields from SWAT data')
        if os.path.isfile(acpnamelist.fnames.swatsims):
            if acpnamelist.vars.verbose: print('    Reading: '+acpnamelist.fnames.swatsims)
            self.nutrientdata.swatsims = pandas.read_csv(acpnamelist.fnames.swatsims)
            self.nutrientdata.swatsims.set_index('HRUGIS',inplace=True)
        else:
            dfs = list()
            for fname in acpnamelist.vars.input_swat_outputs:
                if fname.find('notill') == -1:
                    if acpnamelist.vars.verbose: print('    Reading: '+fname)
                    ls = list(open(fname,'r'))
                    dt = {'HRUGIS':[10,20],'MON':[30,34],'AREAkm2':[34,45],'ORGNkg/ha':[114,125],'ORGPkg/ha':[124,135],'NSURQkg/ha':[144,155],'NLATQkg/ha':[154,165],'NO3Lkg/ha':[164,175],'NO3GWkg/ha':[174,185],'SOLPkg/ha':[184,195]}
                    vs = {name:list() for name in dt}
                    for i in range(9,len(ls)):
                        if float(ls[i][dt['MON'][0]:dt['MON'][1]]) < 1900:
                            for name in dt:
                                if name.find('HRUGIS') != -1: v = str(ls[i][dt[name][0]:dt[name][1]]).replace(' ','')
                                else: v = float(ls[i][dt[name][0]:dt[name][1]])
                                vs[name].append(v)
                    df = pandas.DataFrame(vs)
                    modelid = -1
                    if fname.find('LEFW') != -1: modelid = 'LEFW'
                    elif fname.find('UEFW') != -1: modelid = 'UEFW'
                    df['HRUGIS'] = modelid + df['HRUGIS']
                    dfs.append(df)
            self.nutrientdata.swatsims = pandas.concat(dfs)
            self.nutrientdata.swatsims.set_index('HRUGIS',inplace=True)
            self.nutrientdata.swatsims.rename(columns={'ORGNkg/ha':'ORGN','ORGPkg/ha':'ORGP','NSURQkg/ha':'NSURQ','NLATQkg/ha':'NLAT','NO3GWkg/ha':'GWQN','SOLPkg/ha':'SOLP'},inplace=True)
            if acpnamelist.vars.verbose: print('    Making: '+acpnamelist.fnames.swatsims)
            self.nutrientdata.swatsims.to_csv(acpnamelist.fnames.swatsims,index=True)
    
    def _setBaselineNutrientYieldDistributions_SWAT(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Set/create statistical distribution objects using SWAT simulated nutrient yields form agricultural HRUs"""
        if acpnamelist.vars.verbose: print('Making statistical distribution objects for baseline nutrient yields from SWAT simulations')
        self.nutrientdata.baseline_distributions_SWAT = dict()
        if 'LU_CODE' not in list(self.nutrientdata.swatsims.columns):
            if self.nutrientdata.swatsims.index.name != 'HRUGIS': self.nutrientdata.swatsims.set_index('HRUGIS',inplace=True)
            if self.spatialdata.domain_hrus.index.name != 'HRUGIS': self.spatialdata.domain_hrus.set_index('HRUGIS',inplace=True)
            self.nutrientdata.swatsims = self.nutrientdata.swatsims.join(other=self.spatialdata.domain_hrus[['LU_CODE']],on='HRUGIS')
        for name in self.swatoutvars.n_vars + self.swatoutvars.p_vars:
            iqr = self.nutrientdata.swatsims[name].quantile(0.75)-self.nutrientdata.swatsims[name].quantile(0.25)
            upb = self.nutrientdata.swatsims[name].quantile(0.75) + 3 * iqr
            lwb = self.nutrientdata.swatsims[name].quantile(0.25) - 3 * iqr
            val = self.nutrientdata.swatsims[(self.nutrientdata.swatsims[name] >= lwb)&(self.nutrientdata.swatsims[name] <= upb)&(~self.nutrientdata.swatsims['LU_CODE'].isin(['SEPT','WETN','WATR','URHD','URLD','FRSD']))][name] #['SEPT','WETN','WATR','URHD','URLD','FRSD'] #exclude these lu codes
            kde = scipy.stats.gaussian_kde(val)
            xk = numpy.arange(min(val),max(val),(max(val)-min(val))*0.01)
            pk = [kde.pdf(v).item() for v in xk]
            pk_norm = tuple(p/sum(pk) for p in pk)
            self.nutrientdata.baseline_distributions_SWAT[name] =  scipy.stats.rv_discrete(a=0, b=max(val), name=name, values=(xk, pk_norm))
            #random_number = dtdistkeyparm[name].rvs() #how-to draw random number from distribution

    def _readECs_ACRE(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Read effectiveness coefficient distributions from ACRE directory"""
        if acpnamelist.vars.verbose: print('Reading effectiveness coefficients from ACRE data')
        dfs = list()
        for acp in ['CONSERVATION','Contour Farming','Fert 75%','Fert 90%','Filterstrip','MIN_TILL','Ponds 50%','Ponds 75%','Terraces And Waterway','Terraces Only','Waterway Only']:
            csv_acp = os.path.join(acpnamelist.dirnames.acre_data,acp)+'.csv'
            _, file_ext = os.path.splitext(csv_acp)
            if file_ext.find('csv') != -1:
                if acpnamelist.vars.verbose: print('    Reading: '+csv_acp)
                df = pandas.read_csv(filepath_or_buffer=csv_acp)
                df['ACP'] = acp
                dfs.append(df)
        self.ecs.ecs_ACRE = pandas.concat(dfs)
        self.ecs.ecs_ACRE.loc[self.ecs.ecs_ACRE["Parm"]=="DP", "Parm"] = "SOLP" # assume 'DP' = 'SOLP'
        csv_acp = os.path.join(acpnamelist.dirnames.acre_data,'Simple_con')+'.csv'
        _, file_ext = os.path.splitext(csv_acp)
        if file_ext.find('csv') != -1:
            if acpnamelist.vars.verbose: print('    Reading: '+csv_acp)
            self.ecs.ecs_ACRE_simple = pandas.read_csv(filepath_or_buffer=csv_acp)
            self.ecs.ecs_ACRE_simple = self.ecs.ecs_ACRE_simple.rename(columns={"Practice": "ACP","Constituent": "Parm","min": "pct_1","median": "pct_50","max": "pct_100"})
            self.ecs.ecs_ACRE_simple['pct_1'] = self.ecs.ecs_ACRE_simple['pct_1']/100
            self.ecs.ecs_ACRE_simple['pct_50'] = self.ecs.ecs_ACRE_simple['pct_50']/100
            self.ecs.ecs_ACRE_simple['pct_100'] = self.ecs.ecs_ACRE_simple['pct_100']/100
        self.miscstats.acp_types = list(set(list(self.ecs.ecs_ACRE['ACP'].unique())+list(self.ecs.ecs_ACRE_simple['ACP'].unique())))

    def _setECDistributions_ACRE(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Make effectiveness coefficient distribution objects per ACP-type, county, and nutrient form and loss pathway"""
        if acpnamelist.vars.verbose: print('Making statistical distribution objects for effectiveness coefficients from ACRE data')
        self.ecs.ec_distributions_ACRE = dict()
        v_pcts = [.01,.05,.25,.50,.75,.95,1]
        for acp in list(self.ecs.ecs_ACRE['ACP'].unique()):
            if acp not in self.ecs.ec_distributions_ACRE:
                self.ecs.ec_distributions_ACRE[acp] = dict()
            for fip in [int(fip) for fip in list(self.spatialdata.domain_counties['GEOID'].unique())]:
                if fip not in self.ecs.ec_distributions_ACRE[acp]:
                    self.ecs.ec_distributions_ACRE[acp][fip] = dict()
                for parm in list(self.ecs.ecs_ACRE['Parm'].unique()):
                    vs = [self.ecs.ecs_ACRE[(self.ecs.ecs_ACRE['ACP']==acp)&
                                            (self.ecs.ecs_ACRE['FIPS']==fip)&
                                            (self.ecs.ecs_ACRE['Parm']==parm)][pct].item() for pct in ['pct_1','pct_5','pct_25','pct_50','pct_75','pct_95','pct_100']]
                    self.ecs.ec_distributions_ACRE[acp][fip][parm] = ACPDistributionType(percentiles=v_pcts,values=vs)
        v_pcts = [.01,.50,1]
        for acp in list(self.ecs.ecs_ACRE_simple['ACP'].unique()):
            if acp not in self.ecs.ec_distributions_ACRE:
                self.ecs.ec_distributions_ACRE[acp] = dict()
            dt = dict()
            for parm in list(self.ecs.ecs_ACRE_simple['Parm'].unique()):
                vs = [self.ecs.ecs_ACRE_simple[(self.ecs.ecs_ACRE_simple['ACP']==acp)&
                                               (self.ecs.ecs_ACRE_simple['Parm']==parm)][pct].item() for pct in ['pct_1','pct_50','pct_100']]
                dt[parm] = ACPDistributionType(percentiles=v_pcts,values=vs)
            for fip in [int(fip) for fip in list(self.spatialdata.domain_counties['GEOID'].unique())]:
                if fip not in self.ecs.ec_distributions_ACRE[acp]:
                    self.ecs.ec_distributions_ACRE[acp][fip] = dict()
                self.ecs.ec_distributions_ACRE[acp][fip] = dt

    def _get_baseline_yields_swatpdf(self):
        """Get baseline nutrient yields for field using SWAT-based PDFs for agricultural HRUs"""
        dt_baseline_kgha_key_parm = {parm:self.nutrientdata.baseline_distributions_SWAT[parm].rvs() for parm in self.swatoutvars.n_vars + self.swatoutvars.p_vars}
        return dt_baseline_kgha_key_parm

    def _get_baseline_yields_swatraw(self,field_id):
        """Get baseline nutrient yields for field using raw SWAT outputs (i.e. from HRU intersecting field)"""
        if field_id in self.miscstats.dtHRUGISKeyFieldID: HRUGIS = self.miscstats.dtHRUGISKeyFieldID[field_id]
        else: HRUGIS = self.miscstats.dtHRUGISKeyFieldID[random.choice(list(self.miscstats.dtHRUGISKeyFieldID.keys()))] 
        dt_baseline_kgha_key_parm = {parm:self.nutrientdata.swatsims.iloc[self.nutrientdata.swatsims.index==HRUGIS][parm].item() for parm in self.swatoutvars.n_vars + self.swatoutvars.p_vars}
        return dt_baseline_kgha_key_parm
    
    def _get_baseline_yields_acre(self,field_id):
        """Get baseline nutrient yields for field from ACRE simulations"""
        # need to do - _get_baseline_yields_acre takes too long; find a better way
        df = self.nutrientdata.baseline_yields_kgha_ACRE
        FIP = int(self.miscstats.dt_fips_key_fieldid[field_id])
        Soil_Series_counts = df.loc[(df['FIPS']==FIP),'Soil_Series'].value_counts()
        if 'ALL' in list(Soil_Series_counts.index):  Soil_Series_counts.drop(['ALL'],axis='index')
        Soil_Series = random.choices(population=list(Soil_Series_counts.index),weights=list(Soil_Series_counts),k=1)[0]
        Crop_system_counts = df.loc[(df['FIPS']==FIP) & (df['Soil_Series']==Soil_Series),'Crop_system'].value_counts()
        Crop_system = random.choices(population=list(Crop_system_counts.index),weights=list(Crop_system_counts),k=1)[0]
        Irrigated_counts = df.loc[(df['FIPS']==FIP) & (df['Soil_Series']==Soil_Series) & (df['Crop_system']==Crop_system),'Irrigated'].value_counts()
        Irrigated = random.choices(population=list(Irrigated_counts.index),weights=list(Irrigated_counts),k=1)[0]
        Tiled_counts = df.loc[(df['FIPS']==FIP) & (df['Soil_Series']==Soil_Series) & (df['Crop_system']==Crop_system) & (df['Irrigated']==Irrigated),'Tiled'].value_counts()
        Tiled = random.choices(population=list(Tiled_counts.index),weights=list(Tiled_counts),k=1)[0]
        df_subset = df[(df['FIPS']==FIP) & (df['Soil_Series']==Soil_Series) & (df['Crop_system']==Crop_system) & (df['Irrigated']==Irrigated) & (df['Tiled']==Tiled)]
        dtbasekgha = {Parm:df_subset.loc[df['Parm']==Parm,'pct_50'].iloc[[0]].item() for Parm in self.swatoutvars.n_vars + self.swatoutvars.p_vars}
        dtbasekgha = {Parm:float(dtbasekgha[Parm]) for Parm in dtbasekgha}
        return dtbasekgha
    
    def get_baseline_yields(self,field_id,method='swat pdf'):
        """Get baseline nutrient yields for field"""
        if method.upper().replace(' ','').find('SWATRAW') != -1:
            dtbasekgha = self._get_baseline_yields_swatraw(field_id)
        elif method.upper().replace(' ','').find('SWATPDF') != -1:
            dtbasekgha = self._get_baseline_yields_swatpdf(field_id)
        elif method.upper().replace(' ','').find('ACRE') != -1:
            dtbasekgha = self._get_baseline_yields_acre(field_id)
        else:
            sys.exit('ERROR unknown method')
        return dtbasekgha

    def _get_acp_ecs_acremedian(self,field_id,acp_type):
        """Get ECs for each nutrient form / loss pathway - use median value from ACRE distribution"""
        dt_ec_parm = dict()
        FIP = int(self.miscstats.dt_fips_key_fieldid[field_id])
        if acp_type in set(self.ecs.ecs_ACRE['ACP'].unique()):
            dt_ec_parm = {parm:self.ecs.ecs_ACRE[(self.ecs.ecs_ACRE['FIPS']==FIP)&(self.ecs.ecs_ACRE['ACP']==acp_type)&
                                                 (self.ecs.ecs_ACRE['Parm']==parm)]['pct_50'].item() for parm in self.swatoutvars.n_vars + self.swatoutvars.p_vars}
        elif acp_type in set(self.ecs.ecs_ACRE_simple['ACP'].unique()):
            dt_ec_parm = {parm:self.ecs.ecs_ACRE_simple[(self.ecs.ecs_ACRE_simple['ACP']==acp_type)&
                                                        (self.ecs.ecs_ACRE_simple['Parm']==parm)]['pct_50'].item() for parm in self.swatoutvars.n_vars + self.swatoutvars.p_vars}
        else:
            sys.exit('ERROR unrecognized ACP type ',acp_type)
        return dt_ec_parm
        
    def _get_acp_ecs_acre(self,field_id,acp_type):
        """Get ECs for each nutrient form / loss pathway - use random value from ACRE-specified distribution"""
        FIP = int(self.miscstats.dt_fips_key_fieldid[field_id])
        dt_ec_parm = {parm:self.ecs.ec_distributions_ACRE[acp_type][FIP][parm].rvs().item() for parm in self.swatoutvars.n_vars + self.swatoutvars.p_vars}
        return dt_ec_parm
    
    def get_acp_ecs(self,field_id,acp_type,method='acre'):
        """Get ECs for each nutrient form / loss pathway"""
        if method.upper().replace(' ','').find('ACREMEDIAN') != -1:
            dt_ecs_parm = self._get_acp_ecs_acremedian(field_id=field_id,acp_type=acp_type)
        elif method.upper().replace(' ','').find('ACRE') != -1:
            dt_ecs_parm = self._get_acp_ecs_acre(field_id=field_id,acp_type=acp_type)
        else:
            sys.exit('ERROR unrecognized method type')
        return dt_ecs_parm
    
    def get_folium_map(self):
        domain_centroid = self.spatialdata.domain.to_crs('+proj=cea').centroid.to_crs(self.spatialdata.domain.crs)
        domain_map = folium.Map([domain_centroid.y.iloc[0],domain_centroid.x.iloc[0]])
        domainfg = folium.FeatureGroup(name='Domain')
        for _, r in self.spatialdata.domain.iterrows(): 
            folium.GeoJson(data=geopandas.GeoSeries(r['geometry']).to_json(),style_function=lambda x: {"color":"black","fillColor": "none"}).add_to(domainfg)
        countiesfg = folium.FeatureGroup(name='County Boundaries')
        for _, r in self.spatialdata.domain_counties.iterrows(): 
            folium.GeoJson(data=geopandas.GeoSeries(r['geometry']).to_json(),style_function=lambda x: {"color":"grey","fill": False}).add_to(countiesfg)
        fieldsfg = folium.FeatureGroup(name='Field Boundaries')
        for _, r in self.spatialdata.domain_fields.iterrows(): 
            folium.GeoJson(data=geopandas.GeoSeries(r['geometry']).to_json(),style_function=lambda x: {"color":"black","fill": False}).add_to(fieldsfg)
        #hrusfg = folium.FeatureGroup(name='SWAT HRU Boundaries')
        #for _, r in domain_hrus.iterrows(): 
        #   folium.GeoJson(data=geopandas.GeoSeries(r['geometry']).to_json(),style_function=lambda x: {"color":"grey","fill": False}).add_to(hrusfg)
        domainfg.add_to(domain_map)
        countiesfg.add_to(domain_map)
        fieldsfg.add_to(domain_map)
        #hrusfg.add_to(domain_map)
        folium.LayerControl().add_to(domain_map)
        return domain_map