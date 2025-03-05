import os,sys,scipy,numpy,pandas,random
import ACPNamelist,ACPDataSpatial

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
        
class ACPDataNutrientYields:
    """Class to hold nutrient yield data"""
        
    class ACRE:
        baseline_yield_kgha_ACRE       = numpy.nan

    class SWAT:
        baseline_kgha_distributions    = numpy.nan
        swat_output_hru                = numpy.nan
        swat_output_rch                = numpy.nan
        mean_annual_TN_basin_load_kgyr = numpy.nan
        std_annual_TN_basin_load_kgyr  = numpy.nan
        mean_annual_TP_basin_load_kgyr = numpy.nan
        std_annual_TP_basin_load_kgyr  = numpy.nan
        n_vars                         = ('NSURQ','ORGN','NLAT','GWQN')
        p_vars                         = ('ORGP','SOLP')
        
    def __init__(self,acpnamelist:ACPNamelist.ACPNamelist,acpspatial:ACPDataSpatial.ACPDataSpatial,method:str='ACRE'):
        """Initialize ACPData"""
        if acpnamelist.vars.verbose: print('Initializing ACP data class')
        self._init_vars(acpnamelist,method)
        self._read_data(acpnamelist,acpspatial)

    def _init_vars(self,acpnamelist:ACPNamelist.ACPNamelist,method:str):
        """Initialize class variables"""
        if acpnamelist.vars.verbose: print('Initializing variables')
        if   method.upper().find('ACRE') != -1: self.nutrientdata = ACPDataNutrientYields.ACRE()
        elif method.upper().find('SWAT') != -1: self.nutrientdata = ACPDataNutrientYields.SWAT()

    def _read_data(self,acpnamelist:ACPNamelist.ACPNamelist,acpspatial:ACPDataSpatial.ACPDataSpatial):
        if acpnamelist.vars.verbose: print('Reading data')
        if   isinstance(self.nutrientdata,ACPDataNutrientYields.ACRE): self._read_data_ACRE(acpnamelist,acpspatial)
        elif isinstance(self.nutrientdata,ACPDataNutrientYields.SWAT): self._read_data_SWAT(acpnamelist)

    def _read_data_ACRE(self,acpnamelist:ACPNamelist.ACPNamelist,acpspatial:ACPDataSpatial.ACPDataSpatial):
        """Read mean annual simulations from ACRE inputs"""
        if acpnamelist.vars.verbose: print('Reading baseline nutrient yields from ACRE data')
        csv_fips = os.path.join(acpnamelist.dirnames.acre_data,'FIPS')+'.csv'
        _, file_ext = os.path.splitext(csv_fips)
        if file_ext.find('csv') != -1: df_fips = pandas.read_csv(filepath_or_buffer=csv_fips)
        dfs = list()
        for i,r in acpspatial.spatialdata.domain_counties.iterrows():
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
        self.nutrientdata.baseline_yields_kgha_ACRE = pandas.merge(left=self.nutrientdata.baseline_yields_kgha_ACRE,
                                                                   right=df_fips,left_on=['state','county'],right_on=['STATE_NAME','NAME'])
        self.nutrientdata.baseline_yields_kgha_ACRE = self.nutrientdata.baseline_yields_kgha_ACRE.drop(['STATE_NAME', 'NAME'], axis=1)

    def _read_data_SWAT(self,acpnamelist:ACPNamelist.ACPNamelist):
        """Read swat output data"""
        self._read_data_SWAT_outputhru(acpnamelist=acpnamelist)
        self._read_data_SWAT_outputrch(acpnamelist=acpnamelist)
        self._set_nutrient_yield_distributions_SWAT(acpnamelist=acpnamelist)

    def get_swat_mean_annual_loads(self,subbasin_num:int):
        """Get mean annual TN and TP loads for subbasin"""
        mean_annual_TN_basin_load_kgyr = self.nutrientdata.swat_output_rch.loc[self.nutrientdata.swat_output_rch['rch']==subbasin_num,'TNkgyr'].mean()
        mean_annual_TP_basin_load_kgyr = self.nutrientdata.swat_output_rch.loc[self.nutrientdata.swat_output_rch['rch']==subbasin_num,'TPkgyr'].mean()
        return mean_annual_TN_basin_load_kgyr, mean_annual_TP_basin_load_kgyr
        
    def _read_data_SWAT_outputrch(self,acpnamelist:ACPNamelist.ACPNamelist): 
        if acpnamelist.vars.verbose: print('Reading SWAT baseline reach nutrient data')
        if hasattr(acpnamelist.vars, 'input_swat_output_rch'):
            fname = acpnamelist.vars.input_swat_output_rch[0]
            if acpnamelist.vars.verbose: print('    Reading: '+fname)
            ls = list(open(fname,'r'))
            names = ['rch','MON','TNkgyr','TPkgyr']
            vals = {name:list() for name in names}
            for i in range(9,len(ls)):
                vs = ls[i].replace('\n','').split()
                if float(vs[3]) < 1900:
                    vals['rch'].append(int(vs[1]))
                    vals['MON'].append(int(vs[3]))
                    vals['TNkgyr'].append(float(vs[8]))
                    vals['TPkgyr'].append(float(vs[9]))
            self.nutrientdata.swat_output_rch = pandas.DataFrame(vals)

    def _read_data_SWAT_outputhru(self,acpnamelist:ACPNamelist.ACPNamelist): 
        """Read mean annual simulations from input swat output.hru files"""
        ## TO DO - refactor to exclude study area specific code - i.e., LEFW, UEFW, etc.
        if acpnamelist.vars.verbose: print('Reading SWAT baseline HRU nutrient data')
        dfs = list()
        for fname in acpnamelist.vars.input_swat_output_hru:
            if fname.find('notill') == -1:
                if acpnamelist.vars.verbose: print('    Reading: '+fname)
                ls = list(open(fname,'r'))
                dt = {'LULC':[0,4],
                      'HRUGIS':[10,20],
                      'MON':[30,34],
                      'AREAkm2':[34,45],
                      'ORGNkg/ha':[114,125],
                      'ORGPkg/ha':[124,135],
                      'NSURQkg/ha':[144,155],
                      'NLATQkg/ha':[154,165],
                      'NO3Lkg/ha':[164,175],
                      'NO3GWkg/ha':[174,185],
                      'SOLPkg/ha':[184,195]}
                vs = {name:list() for name in dt}
                for i in range(9,len(ls)):
                    if float(ls[i][dt['MON'][0]:dt['MON'][1]]) < 1900:
                        for name in dt:
                            if name.find('HRUGIS') != -1: v = str(ls[i][dt[name][0]:dt[name][1]]).replace(' ','')
                            elif name.find('LULC') != -1: v = str(ls[i][dt[name][0]:dt[name][1]]).replace(' ','')
                            else: v = float(ls[i][dt[name][0]:dt[name][1]])
                            vs[name].append(v)
                df = pandas.DataFrame(vs)
                modelid = -1
                if fname.find('LEFW') != -1: modelid = 'LEFW'
                elif fname.find('UEFW') != -1: modelid = 'UEFW'
                df['HRUGIS'] = modelid + df['HRUGIS']
                dfs.append(df)
        self.nutrientdata.swat_output_hru = pandas.concat(dfs)
        self.nutrientdata.swat_output_hru.set_index('HRUGIS',inplace=True)
        rn_cols = {'ORGNkg/ha':'ORGN','ORGPkg/ha':'ORGP','NSURQkg/ha':'NSURQ','NLATQkg/ha':'NLAT','NO3GWkg/ha':'GWQN','SOLPkg/ha':'SOLP'}
        self.nutrientdata.swat_output_hru.rename(columns=rn_cols,inplace=True)

    def _set_nutrient_yield_distributions_SWAT(self,acpnamelist:ACPNamelist.ACPNamelist):
        """Set/create statistical distribution objects using SWAT simulated nutrient yields form agricultural HRUs"""
        if acpnamelist.vars.verbose: print('Making statistical distribution objects for baseline nutrient yields from SWAT simulations')
        self.nutrientdata.baseline_kgha_distributions = dict()
        for name in self.nutrientdata.n_vars + self.nutrientdata.p_vars:
            iqr = self.nutrientdata.swat_output_hru[name].quantile(0.75)-self.nutrientdata.swat_output_hru[name].quantile(0.25)
            upb = self.nutrientdata.swat_output_hru[name].quantile(0.75) + 3 * iqr
            lwb = self.nutrientdata.swat_output_hru[name].quantile(0.25) - 3 * iqr
            val = self.nutrientdata.swat_output_hru[(self.nutrientdata.swat_output_hru[name] >= lwb)&(self.nutrientdata.swat_output_hru[name] <= upb)&(~self.nutrientdata.swat_output_hru['LULC'].isin(['SEPT','WETN','WATR','URHD','URLD','FRSD']))][name] #['SEPT','WETN','WATR','URHD','URLD','FRSD'] #exclude these lu codes
            kde = scipy.stats.gaussian_kde(val)
            xk = numpy.arange(min(val),max(val),(max(val)-min(val))*0.01)
            pk = [kde.pdf(v).item() for v in xk]
            pk_norm = tuple(p/sum(pk) for p in pk)
            self.nutrientdata.baseline_kgha_distributions[name] =  scipy.stats.rv_discrete(a=0, b=max(val), name=name, values=(xk, pk_norm))

    def get_baseline_yields(self,field_id):
        """Get baseline nutrient yields for field"""
        if isinstance(self.nutrientdata,ACPDataNutrientYields.ACRE): dtbasekgha = self._get_baseline_yields_acre(field_id)
        if isinstance(self.nutrientdata,ACPDataNutrientYields.SWAT): dtbasekgha = self._get_baseline_yields_swatpdf()
        else: sys.exit('ERROR unknown method')
        return dtbasekgha
    
    def _get_baseline_yields_swatpdf(self):
        """Get baseline nutrient yields for field using SWAT-based PDFs for agricultural HRUs"""
        dt_baseline_kgha_key_parm = {parm:self.nutrientdata.baseline_kgha_distributions[parm].rvs() for parm in self.nutrientdata.n_vars + self.nutrientdata.p_vars}
        return dt_baseline_kgha_key_parm

    def _get_baseline_yields_swatraw(self,field_id):
        """Get baseline nutrient yields for field using raw SWAT outputs (i.e. from HRU intersecting field)"""
        if field_id in self.miscstats.dtHRUGISKeyFieldID: HRUGIS = self.miscstats.dtHRUGISKeyFieldID[field_id]
        else: HRUGIS = self.miscstats.dtHRUGISKeyFieldID[random.choice(list(self.miscstats.dtHRUGISKeyFieldID.keys()))] 
        dt_baseline_kgha_key_parm = {parm:self.nutrientdata.swat_output_hru.iloc[self.nutrientdata.swat_output_hru.index==HRUGIS][parm].item() for parm in self.swatoutvars.n_vars + self.swatoutvars.p_vars}
        return dt_baseline_kgha_key_parm
    
    def _get_baseline_yields_acre(self,field_id):
        """Get baseline nutrient yields for field from ACRE simulations"""
        # need to do - this function takes too long; find a better way
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