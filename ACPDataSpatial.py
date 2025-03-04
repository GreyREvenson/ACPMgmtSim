import os,numpy,pandas,geopandas,pygeohydro,ee,us,ACPNamelist,folium
        
class ACPDataSpatial:
    """Class to hold the data"""
    
    class SpatialData:
        """Spatial data"""
        domain                         = numpy.nan
        domain_huc12s                  = numpy.nan # empty unless overwrite = True
        domain_fields                  = numpy.nan
        domain_counties                = numpy.nan
        
    class MiscStats:
        """Miscellaneous"""
        ## need to do - better places /classifications for these things? Do we even need?
        domain_area_ha                 = numpy.nan
        dt_fips_key_fieldid            = numpy.nan
        dt_areaha_key_fieldid          = numpy.nan
        hrugiss_agonly                 = numpy.nan # used?
        field_ids                      = numpy.nan

    def __init__(self,acpnamelist:ACPNamelist.ACPNamelist):
        """Initialize ACPData"""
        if acpnamelist.vars.verbose: print('Initializing ACP data class')
        self._init_vars(acpnamelist)
        self._init_dirs(acpnamelist)
        self._read_spatial_data(acpnamelist)

    def _init_vars(self,acpnamelist:ACPNamelist.ACPNamelist):
        """Initialize class variables"""
        if acpnamelist.vars.verbose: print('Initializing variables')
        self.miscstats    = ACPDataSpatial.MiscStats()
        self.spatialdata  = ACPDataSpatial.SpatialData()
        
    def _init_dirs(self,acpnamelist:ACPNamelist.ACPNamelist):
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
        if not os.path.isdir(acpnamelist.dirnames.swat): 
            os.mkdir(acpnamelist.dirnames.swat)
        if not os.path.isdir(acpnamelist.dirnames.output): 
            os.mkdir(acpnamelist.dirnames.output)
        
    def _read_spatial_data(self,acpnamelist:ACPNamelist.ACPNamelist):
        """Make or read spatial data"""
        if acpnamelist.vars.verbose: print('Making or reading project spatial data')
        self._setDomain(acpnamelist)
        self._setCounties(acpnamelist)
        self._setFieldBoundaries(acpnamelist)
        self._setHRUBoundaries(acpnamelist)
        self._evalHRUFieldOverlap(acpnamelist)
        
    def _setDomain(self,acpnamelist:ACPNamelist.ACPNamelist):
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
            
    def _setCounties(self,acpnamelist:ACPNamelist.ACPNamelist):
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
            us_counties = geopandas.GeoDataFrame.from_features(ee.FeatureCollection(us_counties_fc).getInfo())
            us_counties = us_counties[us_counties.STATEFP != '11'] # remove wash dc
            us_counties['STATE_NAME'] = [str(us.states.lookup(r['STATEFP']).name).title() for _,r in us_counties.iterrows()]
            domain_county_intersect = geopandas.overlay(self.spatialdata.domain, us_counties, how='intersection')
            self.spatialdata.domain_counties = us_counties[us_counties["GEOID"].isin(list(set(domain_county_intersect['GEOID'])))]
            self.spatialdata.domain_counties.to_file(acpnamelist.fnames.counties)
            
    def _setFieldBoundaries(self,acpnamelist:ACPNamelist.ACPNamelist):
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
    
    def _setHRUBoundaries(self,acpnamelist:ACPNamelist.ACPNamelist):
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
        
    def _evalHRUFieldOverlap(self,acpnamelist:ACPNamelist.ACPNamelist):
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