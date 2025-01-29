import os,sys,scipy,numpy,pandas,geopandas,pygeohydro,ee,geemap,us,random

class ACPNamelist:

    def __init__(self,fname_namelist:str):

        self._init_vars(fname_namelist)
        self._read_namelist()
        self._set_subdirectory_and_file_names()
    
    def _init_vars(self,fname_namelist:str):

        #save namelist file name
        self.fname_namelist = fname_namelist
        
        #
        self.name_hucids = 'huc12_ids'
        self.name_project_dir = 'project_directory'
        self.name_acre_dir = 'acre_directory'
        self.name_field_boundary_dir = 'field_boundary_directory'
        self.name_swatoutputhru_files = 'swat_output_hru_files'
        self.name_swathruboundary_files = 'swat_hru_spatial_boundary_files'
        self.name_overwrite_itermediate_flag = 'intermediate_overwrite'
        self.fieldname_areaha               = 'Area_ha'

        #vars
        self.list_huc12_ids    = numpy.nan
        self.dir_project = numpy.nan
        self.dir_acre_data = numpy.nan
        self.dir_field_boundaries = numpy.nan
        self.list_swat_output_fnames = numpy.nan
        self.list_swat_hrushp_fnames = numpy.nan
        self.bool_overwrite_flag = False

        #subdirectory names
        self.dirCountyBoundaries = numpy.nan
        self.dirFieldBoundaries = numpy.nan
        self.dirHUCBoundaries = numpy.nan
        self.dirSWATSimulations = numpy.nan
        self.dirOutput = numpy.nan

        #intermediate output names
        self.fname_domain = numpy.nan
        self.fname_huc12s = numpy.nan
        self.fname_domain_counties = numpy.nan
        self.fname_domain_field_boundaries = numpy.nan
        self.fname_domain_hru_boundaries = numpy.nan
        self.fname_field_hrugis_key = numpy.nan
        self.fname_swatsims_processed = numpy.nan
    
    def _read_namelist(self):

        #read vars
        if not os.path.isfile(self.fname_namelist):
            sys.exit('ERROR could not find namelist file '+self.fname_namelist)
        dt = dict()
        namelist_lines = list(open(self.fname_namelist,'r'))
        for l in namelist_lines:
            try:
                l0 = l.replace('\n','').replace(' ','')
                if len(l0) > 0:
                    if str(l[0:1]).find('#') == -1:
                        l1 = l0.split('=')
                        var_name = str(l1[0])
                        var_vals = str(l1[1])
                        var_vals = var_vals.split(',')
                        for i in range(len(var_vals)):
                            val = var_vals[i]
                            if val.startswith("'") and val.endswith("'"):
                                val = val[1:len(val)-1]
                            elif val.startswith('"') and val.endswith('"'):
                                val = val[1:len(val)-1]
                            var_vals[i] = val
                        dt[var_name] = var_vals
            except:
                sys.exit('ERROR could not read namelist.txt line: ',l)
        for var_name in dt:
            if isinstance(dt[var_name],list) and len(dt[var_name]) == 1:
                dt[var_name] = dt[var_name][0]

        #check for required vars
        list_required_vars = [self.name_hucids, self.name_project_dir, self.name_acre_dir, self.name_field_boundary_dir, 
                              self.name_swatoutputhru_files, self.name_swathruboundary_files]
        for var_name in list_required_vars:
            if var_name not in dt:
                sys.exit('ERROR required variable '+var_name+' not found in namelist.txt')

        #assign required vars
        self.list_huc12_ids = dt[self.name_hucids]
        self.dir_project = dt[self.name_project_dir]
        self.dir_acre_data = dt[self.name_acre_dir]
        self.dir_field_boundaries = dt[self.name_field_boundary_dir]
        self.list_swat_output_fnames = dt[self.name_swatoutputhru_files]
        self.list_swat_hrushp_fnames = dt[self.name_swathruboundary_files]

        #optional vars
        if self.name_overwrite_itermediate_flag in dt and dt[self.name_overwrite_itermediate_flag].upper().find('TRUE') != -1:
            self.bool_overwrite_flag = True
        else:
            self.bool_overwrite_flag = False

    def _set_subdirectory_and_file_names(self):

        self.dirCountyBoundaries = os.path.join(self.dir_project,'CountyBoundaries')
        self.dirFieldBoundaries = os.path.join(self.dir_project,'FieldBoundaries')
        self.dirHUCBoundaries = os.path.join(self.dir_project,'HUCBoundaries')
        self.dirSWATSimulations = os.path.join(self.dir_project,'SWAT')
        self.dirOutput = os.path.join(self.dir_project,'Output')

        self.fname_domain = os.path.join(self.dirHUCBoundaries,'domain.gpkg')
        self.fname_huc12s = os.path.join(self.dirHUCBoundaries,'domain_huc12s.gpkg')
        self.fname_domain_counties = os.path.join(self.dirCountyBoundaries,'domain_counties.gpkg')
        self.fname_domain_field_boundaries = os.path.join(self.dirFieldBoundaries,'domain_field_boundaries.gpkg')
        self.fname_domain_hru_boundaries = os.path.join(self.dirSWATSimulations,'domain_SWATHRU_boundaries.gpkg')
        self.fname_field_hrugis_key = os.path.join(self.dirFieldBoundaries,'field_hrugis_key.csv')
        self.fname_swatsims_processed = os.path.join(self.dirSWATSimulations,'swat_sims_processed.csv')
        
class ACPDistributionType(scipy.stats.rv_continuous):
    """
    Generic customizable distribution class - to be used to define nutrient yield and/or effectiveness coefficient distributions
    """
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
    
    class SpatialData:
        domain                         = numpy.nan
        domain_huc12s                  = numpy.nan # empty unless overwrite = True
        domain_fields                  = numpy.nan
        domain_counties                = numpy.nan
        
    class MiscStats:
        domain_area_ha                 = numpy.nan
        dt_fips_key_fieldid            = numpy.nan
        dt_areaha_key_fieldid          = numpy.nan
        acp_types                      = numpy.nan
        hrugiss_agonly                 = numpy.nan # used?
        field_ids                      = numpy.nan
        
    class NutrientData:
        baseline_yields_kgha_ACRE      = numpy.nan
        baseline_distributions_SWAT    = numpy.nan
        swatsims                       = numpy.nan

    class EfficiencyCoefficientData:
        ecs_ACRE                       = numpy.nan
        ecs_ACRE_simple                = numpy.nan
        ec_distributions_ACRE          = numpy.nan
        
    def __init__(self,acpnamelist:ACPNamelist):
        """
        init function
            Args:
                INPUT_OVERWRITE_FLAG (boolean): if true, redo all intermediate data processessing tasks
        """

        self._init_vars()
        self._init_dirs(acpnamelist)
        self._read_spatial_data(acpnamelist)
        self._read_nutrient_and_ec_data(acpnamelist)

    def _init_vars(self):
        
        self.miscstats = ACPData.MiscStats()
        self.spatialdata = ACPData.SpatialData()
        self.nutrientdata = ACPData.NutrientData()
        self.ecs = ACPData.EfficiencyCoefficientData()
        
    def _init_dirs(self,acpnamelist:ACPNamelist):

        if not os.path.isdir(acpnamelist.dir_project): 
            os.mkdir(acpnamelist.dir_project)
        if not os.path.isdir(acpnamelist.dirCountyBoundaries): 
            os.mkdir(acpnamelist.dirCountyBoundaries)
        if not os.path.isdir(acpnamelist.dirFieldBoundaries): 
            os.mkdir(acpnamelist.dirFieldBoundaries)
        if not os.path.isdir(acpnamelist.dirHUCBoundaries): 
            os.mkdir(acpnamelist.dirHUCBoundaries)
        if not os.path.isdir(acpnamelist.dirSWATSimulations): 
            os.mkdir(acpnamelist.dirSWATSimulations)
        if not os.path.isdir(acpnamelist.dirOutput): 
            os.mkdir(acpnamelist.dirOutput)
        
    def _read_spatial_data(self,acpnamelist:ACPNamelist):
        
        self._setDomain(acpnamelist)
        self._setCounties(acpnamelist)
        self._setFieldBoundaries(acpnamelist)
        self._setHRUBoundaries(acpnamelist)
        self._evalHRUFieldOverlap(acpnamelist)

    def _read_nutrient_and_ec_data(self,acpnamelist:ACPNamelist):

        self._readBaselineNutrientYields_ACRE(acpnamelist)
        self._readSWATSimulations(acpnamelist)
        self._setBaselineNutrientYieldDistributions_SWAT()
        self._readECs_ACRE(acpnamelist)
        self._setECDistributions_ACRE()
        
    def _setDomain(self,acpnamelist:ACPNamelist):
        """
        Set/read domain boundary info gdf
        """
        if os.path.isfile(acpnamelist.fname_domain) and not acpnamelist.bool_overwrite_flag:
            self.spatialdata.domain = geopandas.read_file(acpnamelist.fname_domain)
            self.spatialdata.domain[acpnamelist.fieldname_areaha] = self.spatialdata.domain.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.miscstats.domain_area_ha = float(self.spatialdata.domain[acpnamelist.fieldname_areaha].item())
        else:
            self.spatialdata.domain_huc12s = pygeohydro.WBD("huc12").byids("huc12", acpnamelist.list_huc12ids)
            self.spatialdata.domain_huc12s.to_file(acpnamelist.fname_huc12s,driver='GPKG')
            self.spatialdata.domain_huc12s = geopandas.read_file(acpnamelist.fname_huc12s)
            self.spatialdata.domain = self.spatialdata.domain_huc12s.dissolve()
            self.spatialdata.domain.to_file(acpnamelist.fname_domain, driver="GPKG")
            self.spatialdata.domain[acpnamelist.fieldname_areaha] = self.spatialdata.domain.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.miscstats.domain_area_ha = float(self.spatialdata.domain[acpnamelist.fieldname_areaha].item())
            
    def _setCounties(self,acpnamelist:ACPNamelist):
        """
        Set/read US county boundaries into gdf
        """
        if os.path.isfile(acpnamelist.fname_domain_counties) and not acpnamelist.bool_overwrite_flag:
            self.spatialdata.domain_counties = geopandas.read_file(acpnamelist.fname_domain_counties)
        else:
            ee.Authenticate()
            ee.Initialize(project='simmgmttrajectory')
            us_counties_fc  = ee.FeatureCollection("TIGER/2016/Counties")
            us_counties = geemap.ee_to_gdf(us_counties_fc)
            us_counties = us_counties[us_counties.STATEFP != '11'] # remove wash dc
            us_counties['STATE_NAME'] = [str(us.states.lookup(r['STATEFP']).name).title() for _,r in us_counties.iterrows()]
            domain_county_intersect = geopandas.overlay(self.spatialdata.domain, us_counties, how='intersection')
            self.spatialdata.domain_counties = us_counties[us_counties["GEOID"].isin(list(set(domain_county_intersect['GEOID'])))]
            self.spatialdata.domain_counties.to_file(acpnamelist.fname_domain_counties)
            
    def _setFieldBoundaries(self,acpnamelist:ACPNamelist):
        """
        Set/read field boundaries into gdf
        """
        if os.path.isfile(acpnamelist.fname_domain_field_boundaries) and not acpnamelist.bool_overwrite_flag:
            self.spatialdata.domain_fields = geopandas.read_file(acpnamelist.fname_domain_field_boundaries)
            self.miscstats.dt_fips_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields['FIPS']))
            self.miscstats.dt_areaha_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields[acpnamelist.fieldname_areaha]))
        else:
            fields_dfs = list()
            for _,r in self.spatialdata.domain_counties.iterrows():
                dirstate = os.path.join(acpnamelist.fname_fieldboundaries,str(us.states.lookup(str(r['STATE_NAME']).title()).abbr).lower())
                fname = os.path.join(dirstate,str(r['NAME']).title()+'.shp')
                fields_dfs.append(geopandas.read_file(fname))
            fields = pandas.concat(fields_dfs)
            self.spatialdata.domain_fields = geopandas.overlay(self.spatialdata.domain, fields, how='intersection')
            self.spatialdata.domain_fields[self.names.fieldname_areaha] = self.spatialdata.domain_fields.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.spatialdata.domain_fields.rename(columns={"name": "county_name"},inplace=True)
            self.spatialdata.domain_fields = geopandas.sjoin(left_df=self.spatialdata.domain_fields,right_df=self.spatialdata.domain_counties,how='inner')
            self.spatialdata.domain_fields['field_id'] = self.spatialdata.domain_fields.index
            self.miscstats.dt_fips_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields['FIPS']))
            self.miscstats.dt_areaha_key_fieldid = dict(zip(self.spatialdata.domain_fields['field_id'],self.spatialdata.domain_fields[acpnamelist.fieldname_areaha]))
            self.spatialdata.domain_fields.to_file(acpnamelist.fname_domain_field_boundaries,driver='GPKG')
        self.miscstats.field_ids = tuple(self.miscstats.dt_areaha_key_fieldid.keys())
    
    def _setHRUBoundaries(self,acpnamelist:ACPNamelist):
        """
        Set/read HRU boundaries into gdf
        """
        if os.path.isfile(acpnamelist.fname_domain_hru_boundaries) and not acpnamelist.bool_overwrite_flag:
            self.spatialdata.domain_hrus = geopandas.read_file(acpnamelist.fname_domain_hru_boundaries)
        else:
            hru_dfs = list()
            for fname_swathru in acpnamelist.list_swat_hrushp_fnames:
                df = geopandas.read_file(fname_swathru)
                modelid = "ERROR"
                if fname_swathru.find('LEFW') != -1: modelid = 'LEFW'
                if fname_swathru.find('UEFW') != -1: modelid = 'UEFW'
                df['modelid'] = modelid
                hru_dfs.append(df)
            hrus = pandas.concat(hru_dfs)
            self.spatialdata.domain_hrus = geopandas.overlay(self.spatialdata.domain,hrus.to_crs(self.spatialdata.domain.crs.to_wkt()),how='intersection')
            self.spatialdata.domain_hrus[acpnamelistfieldname_areaha] = self.spatialdata.domain_hrus.to_crs({'proj':'cea'})['geometry'].area * 0.0001 #confirm units are ha
            self.spatialdata.domain_hrus.drop(['objectid','shape_Length','Shape_Area'],axis=1,inplace=True)
            self.spatialdata.domain_hrus['HRUGIS'] = self.spatialdata.domain_hrus['modelid'] + self.spatialdata.domain_hrus['HRUGIS']
            self.spatialdata.domain_hrus.to_file(acpnamelist.fname_domain_hru_boundaries,driver='GPKG')
        
    def _evalHRUFieldOverlap(self,acpnamelist:ACPNamelist):
        """
        Determine the areal-dominant HRU for each field
        """
        if os.path.isfile(acpnamelist.fname_field_hrugis_key) and not acpnamelist.bool_overwrite_flag:
            ls = list(open(acpnamelist.fname_field_hrugis_key,'r'))
            ls = [l.replace('\n','').split(',') for l in ls] 
            self.miscstats.dtHRUGISKeyFieldID = {int(ls[i][0]):str(ls[i][1]) for i in range(1,len(ls))}
        else:
            domain_fields_overlay = geopandas.overlay(self.spatialdata.domain_fields, self.spatialdata.domain_hrus.to_crs(self.spatialdata.domain_fields.crs), how='intersection')
            max_area_rows = domain_fields_overlay.groupby("field_id")["Area"].idxmax()
            domain_fields_overlay = domain_fields_overlay.loc[max_area_rows]
            self.miscstats.dtHRUGISKeyFieldID = dict(zip(domain_fields_overlay['field_id'],domain_fields_overlay['HRUGIS']))
            with open(acpnamelist.fname_field_hrugis_key,'w') as openfile:
                openfile.write('field_id,hrugis\n')
                for field_id in self.miscstats.dtHRUGISKeyFieldID:
                    openfile.write(str(field_id)+','+str(self.miscstats.dtHRUGISKeyFieldID[field_id])+'\n')
    
    def _readBaselineNutrientYields_ACRE(self,acpnamelist:ACPNamelist):
        """
        Read mean annual simulations from ACRE inputs
        """
        csv_fips = os.path.join(acpnamelist.dir_acre_data,'FIPS')+'.csv'
        _, file_ext = os.path.splitext(csv_fips)
        if file_ext.find('csv') != -1: df_fips = pandas.read_csv(filepath_or_buffer=csv_fips)
        dfs = list()
        for i,r in self.spatialdata.domain_counties.iterrows():
            dir_state = os.path.join(acpnamelist.dir_acre_data,str(r['STATE_NAME']).title())
            csv_file_county = os.path.join(dir_state,str(r['NAME']).title())+'.csv'
            _, file_ext = os.path.splitext(os.path.join(acpnamelist.dir_acre_data,dir_state,csv_file_county))
            if file_ext.find('csv') != -1:
                df = pandas.read_csv(filepath_or_buffer=csv_file_county)
                df['state'] = str(r['STATE_NAME']).title()
                df['county'] = str(r['NAME']).title()
                dfs.append(df)
        self.nutrientdata.baseline_yields_kgha_ACRE = pandas.concat(dfs)
        self.nutrientdata.baseline_yields_kgha_ACRE = pandas.merge(left=self.nutrientdata.baseline_yields_kgha_ACRE,right=df_fips,left_on=['state','county'],right_on=['STATE_NAME','NAME'])
        self.nutrientdata.baseline_yields_kgha_ACRE = self.nutrientdata.baseline_yields_kgha_ACRE.drop(['STATE_NAME', 'NAME'], axis=1)

    def _readSWATSimulations(self,acpnamelist:ACPNamelist): 
        """
        Read mean annual simulations from input swat output.hru files
        """
        if os.path.isfile(acpnamelist.fname_swatsims_processed):
            self.nutrientdata.swatsims = pandas.read_csv(acpnamelist.fname_swatsims_processed)
            self.nutrientdata.swatsims.set_index('HRUGIS',inplace=True)
        else:
            dfs = list()
            for fname in acpnamelist.list_swat_output_fnames:
                if fname.find('notill') == -1:
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
            self.nutrientdata.swatsims.to_csv(acpnamelist.fname_swatsims_processed,index=True)
    
    def _setBaselineNutrientYieldDistributions_SWAT(self):
        """
        Set/create statistical distribution objects using SWAT simulated nutrient yields form agricultural HRUs
        """
        self.nutrientdata.baseline_distributions_SWAT = dict()
        if 'LU_CODE' not in list(self.nutrientdata.swatsims.columns):
            if self.nutrientdata.swatsims.index.name != 'HRUGIS': self.nutrientdata.swatsims.set_index('HRUGIS',inplace=True)
            if self.spatialdata.domain_hrus.index.name != 'HRUGIS': self.spatialdata.domain_hrus.set_index('HRUGIS',inplace=True)
            self.nutrientdata.swatsims = self.nutrientdata.swatsims.join(other=self.spatialdata.domain_hrus[['LU_CODE']],on='HRUGIS')
        for name in ['ORGN','ORGP','NSURQ','NLAT','GWQN','SOLP']:
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

    def _readECs_ACRE(self,acpnamelist:ACPNamelist):
        """
        Read effectiveness coefficient distributions from ACRE directory
        """
        dfs = list()
        for acp in ['CONSERVATION','Contour Farming','Fert 75%','Fert 90%','Filterstrip','MIN_TILL','Ponds 50%','Ponds 75%','Terraces And Waterway','Terraces Only','Waterway Only']:
            csv_acp = os.path.join(acpnamelist.dir_acre_data,acp)+'.csv'
            _, file_ext = os.path.splitext(csv_acp)
            if file_ext.find('csv') != -1:
                df = pandas.read_csv(filepath_or_buffer=csv_acp)
                df['ACP'] = acp
                dfs.append(df)
        self.ecs.ecs_ACRE = pandas.concat(dfs)
        self.ecs.ecs_ACRE.loc[self.ecs.ecs_ACRE["Parm"]=="DP", "Parm"] = "SOLP" # assume 'DP' = 'SOLP'
        csv_acp = os.path.join(acpnamelist.dir_acre_data,'Simple_con')+'.csv'
        _, file_ext = os.path.splitext(csv_acp)
        if file_ext.find('csv') != -1:
            self.ecs.ecs_ACRE_simple = pandas.read_csv(filepath_or_buffer=csv_acp)
            self.ecs.ecs_ACRE_simple = self.ecs.ecs_ACRE_simple.rename(columns={"Practice": "ACP","Constituent": "Parm","min": "pct_1","median": "pct_50","max": "pct_100"})
            self.ecs.ecs_ACRE_simple['pct_1'] = self.ecs.ecs_ACRE_simple['pct_1']/100
            self.ecs.ecs_ACRE_simple['pct_50'] = self.ecs.ecs_ACRE_simple['pct_50']/100
            self.ecs.ecs_ACRE_simple['pct_100'] = self.ecs.ecs_ACRE_simple['pct_100']/100
        self.miscstats.acp_types = list(set(list(self.ecs.ecs_ACRE['ACP'].unique())+list(self.ecs.ecs_ACRE_simple['ACP'].unique())))

    def _setECDistributions_ACRE(self):
        """
        Create effectiveness coefficient distribution objects per ACP-type, county, and nutrient form and loss pathway
        """
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

    def _get_baseline_yields_swatpdf(self,field_id):
        """
        Get baseline nutrient yields for field using SWAT-based PDFs for agricultural HRUs
            Args:
                field_id (string): field unique identifier
            Returns:
                dictionary: baseline yield per nutrient form and loss pathway 
        """
        Parms = ['NSURQ','ORGN','NLAT','GWQN','ORGP','SOLP']
        dtbasekgha = {Parm:self.nutrientdata.baseline_distributions_SWAT[Parm].rvs() for Parm in Parms}
        return dtbasekgha

    def _get_baseline_yields_swatraw(self,field_id):
        """
        Get baseline nutrient yields for field using raw SWAT outputs (i.e. from HRU intersecting field)
            Args:
                field_id (string): field unique identifier
            Returns:
                dictionary: baseline yield per nutrient form and loss pathway 
        """
        Parms = ['NSURQ','ORGN','NLAT','GWQN','ORGP','SOLP']
        HRUGIS = ''
        if field_id in self.miscstats.dtHRUGISKeyFieldID: 
            HRUGIS = self.miscstats.dtHRUGISKeyFieldID[field_id]
        else: 
            HRUGIS = self.miscstats.dtHRUGISKeyFieldID[random.choice(list(self.miscstats.dtHRUGISKeyFieldID.keys()))] 
        dtbasekgha = {Parm:self.nutrientdata.swatsims.iloc[self.nutrientdata.swatsims.index==HRUGIS][Parm].item() for Parm in Parms}
        return dtbasekgha
    
    def _get_baseline_yields_acre(self,field_id):
        """
        Get baseline nutrient yields for field from ACRE simulations
            Args:
                field_id (string): field unique identifier
            Returns:
                dictionary: baseline yield per nutrient form and loss pathway 
        """
        Parms = ['NSURQ','ORGN','NLAT','GWQN','ORGP','SOLP']
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
        dtbasekgha = {Parm:df_subset.loc[df['Parm']==Parm,'pct_50'].iloc[[0]].item() for Parm in Parms}
        dtbasekgha = {Parm:float(dtbasekgha[Parm]) for Parm in dtbasekgha}
        return dtbasekgha
    
    def get_baseline_yields(self,field_id,method='swat pdf'):
        """
        Get baseline nutrient yields for field. Default method is to use SWAT-based PDFs.
            Args:
                field_id (string): field unique identifier
                method (string): 'swat pdf', 'swat raw', or 'acre'
            Returns:
                dictionary: baseline yield per nutrient form and loss pathway 
        """
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
        """
        Get ECs for each nutrient form / loss pathway - use median value from ACRE distribution
            Args:
                field_id (string): field unique identifier
                acp_type (string): type of ACP
            Returns:
                dictionary: ECs per nutrient form and loss pathway 
        """
        dt_ec_parm = dict()
        Parms = ['NSURQ','ORGN','NLAT','GWQN','ORGP','SOLP']
        FIP = int(self.miscstats.dt_fips_key_fieldid[field_id])
        if acp in set(self.ecs.ecs_ACRE['ACP'].unique()):
            dt_ec_parm = {Parm:self.ecs.ecs_ACRE[(self.ecs.ecs_ACRE['FIPS']==FIP)&(self.ecs.ecs_ACRE['ACP']==acp)&
                                                 (self.ecs.ecs_ACRE['Parm']==Parm)]['pct_50'].item() for Parm in Parms}
        elif acp in set(self.ecs.ecs_ACRE_simple['ACP'].unique()):
            dt_ec_parm = {Parm:self.ecs.ecs_ACRE_simple[(self.ecs.ecs_ACRE_simple['ACP']==acp)&
                                                        (self.ecs.ecs_ACRE_simple['Parm']==Parm)]['pct_50'].item() for Parm in Parms}
        else:
            sys.exit('ERROR unrecognized ACP type: ',acp)
        return dt_ec_parm
        
    def _get_acp_ecs_acre(self,field_id,acp_type):
        """
        Get ECs for each nutrient form / loss pathway - use random value from ACRE-specified distribution
            Args:
                field_id (string): field unique identifier
                acp_type (string): type of ACP
            Returns:
                dictionary: ECs per nutrient form and loss pathway 
        """
        Parms = ['NSURQ','ORGN','NLAT','GWQN','ORGP','SOLP']
        FIP = int(self.miscstats.dt_fips_key_fieldid[field_id])
        dt_ec_parm = {Parm:self.ecs.ec_distributions_ACRE[acp_type][FIP][Parm].rvs().item() for Parm in Parms}
        return dt_ec_parm
    
    def get_acp_ecs(self,field_id,acp_type,method='acre'):
        """
        Get ECs for each nutrient form / loss pathway
            Args:
                field_id (string): field unique identifier
                acp_type (string): type of ACP
            Returns:
                dictionary: ECs per nutrient form and loss pathway 
        """
        if method.upper().replace(' ','').find('ACREMEDIAN') != -1:
            dt_ecs_parm = self._get_acp_ecs_acremedian(field_id=field_id,acp_type=acp_type)
        elif method.upper().replace(' ','').find('ACRE') != -1:
            dt_ecs_parm = self._get_acp_ecs_acre(field_id=field_id,acp_type=acp_type)
        else:
            sys.exit('ERROR unrecognized method type')
        return dt_ecs_parm