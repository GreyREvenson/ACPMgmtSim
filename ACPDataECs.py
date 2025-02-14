import os,sys,scipy,numpy,pandas,ACPTrajectoryNamelist

class ECDistribution(scipy.stats.rv_continuous):
    """Generic customizable distribution class"""
    def __init__(self, percentiles, values):
        super().__init__()
        self.percentiles = percentiles
        self.values = values
        self.interp_func = scipy.interpolate.interp1d(percentiles, values, kind='linear', fill_value="extrapolate")
    def _cdf(self, x):
        return numpy.interp(x, self.values, self.percentiles)
    def rvs(self, size=1, random_state=None):
        return self.interp_func(numpy.random.rand(size))
        
class ACPDataECs:
    """Class to hold EC data"""
    
    class ECDataACRE:
        ecs_ACRE         = numpy.nan
        ecs_ACRE_simple  = numpy.nan
        ecs_ACRE_dists   = numpy.nan

    class ECDataSWAT:
        ecs_SWAT         = numpy.nan
        ecs_SWAT_dists   = numpy.nan

    def __init__(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist,method:str='ACRE'):
        """Initialize ACPData class"""
        if acpnamelist.vars.verbose: print('Initializing ACPDataECs instance')
        self._init_vars(method)

    def _init_vars(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist,method:str):
        """Initialize class variables"""
        if acpnamelist.vars.verbose: print('Initializing variables')
        if   method.upper().find('ACRE') != -1: self.ecs = ACPDataECs.ECDataACRE()
        elif method.upper().find('SWAT') != -1: self.ecs = ACPDataECs.ECDataSWAT()

    def _read_data(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist,method:str):
        """Read EC data"""
        if acpnamelist.vars.verbose: print('Initializing variables')
        if   method.upper().find('ACRE') != -1: self._read_data_ACRE()
        elif method.upper().find('SWAT') != -1: self._read_data_SWAT()
        else: sys.exit('ERROR unrecognized method of reading effectiveness efficiency data')
        self._set_acp_types()

    def _read_data_SWAT(self):
        """placeholder function to read ecs from swat outputs?"""
        pass

    def _read_data_ACRE(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Read effectiveness coefficient distributions from ACRE directory"""
        if acpnamelist.vars.verbose: print('Reading effectiveness coefficients from ACRE data')
        acps_ACRE = ['CONSERVATION','Contour Farming','Fert 75%','Fert 90%','Filterstrip','MIN_TILL',
                     'Ponds 50%','Ponds 75%','Terraces And Waterway','Terraces Only','Waterway Only']
        dfs = list()
        for acp in acps_ACRE:
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
            rn_cols = {"Practice": "ACP","Constituent": "Parm","min": "pct_1","median": "pct_50","max": "pct_100"}
            self.ecs.ecs_ACRE_simple = self.ecs.ecs_ACRE_simple.rename(columns=rn_cols)
            self.ecs.ecs_ACRE_simple['pct_1'] = self.ecs.ecs_ACRE_simple['pct_1']/100
            self.ecs.ecs_ACRE_simple['pct_50'] = self.ecs.ecs_ACRE_simple['pct_50']/100
            self.ecs.ecs_ACRE_simple['pct_100'] = self.ecs.ecs_ACRE_simple['pct_100']/100

    def _set_ec_distributions(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Make effectiveness coefficient distribution objects per ACP-type, county, and nutrient form and loss pathway"""
        if acpnamelist.vars.verbose: print('Making statistical distribution objects for effectiveness coefficients from ACRE data')
        self.ecs.ecs_ACRE_dists = dict()
        v_pcts = [.01,.05,.25,.50,.75,.95,1]
        for acp in list(self.ecs.ecs_ACRE['ACP'].unique()):
            if acp not in self.ecs.ecs_ACRE_dists:
                self.ecs.ecs_ACRE_dists[acp] = dict()
            for fip in [int(fip) for fip in list(self.spatialdata.domain_counties['GEOID'].unique())]:
                if fip not in self.ecs.ecs_ACRE_dists[acp]:
                    self.ecs.ecs_ACRE_dists[acp][fip] = dict()
                for parm in list(self.ecs.ecs_ACRE['Parm'].unique()):
                    vs = [self.ecs.ecs_ACRE[(self.ecs.ecs_ACRE['ACP']==acp)&
                                            (self.ecs.ecs_ACRE['FIPS']==fip)&
                                            (self.ecs.ecs_ACRE['Parm']==parm)][pct].item() for pct in ['pct_1','pct_5','pct_25','pct_50','pct_75','pct_95','pct_100']]
                    self.ecs.ecs_ACRE_dists[acp][fip][parm] = ECDistribution(percentiles=v_pcts,values=vs)
        v_pcts = [.01,.50,1]
        for acp in list(self.ecs.ecs_ACRE_simple['ACP'].unique()):
            if acp not in self.ecs.ecs_ACRE_dists:
                self.ecs.ecs_ACRE_dists[acp] = dict()
            dt = dict()
            for parm in list(self.ecs.ecs_ACRE_simple['Parm'].unique()):
                vs = [self.ecs.ecs_ACRE_simple[(self.ecs.ecs_ACRE_simple['ACP']==acp)&
                                               (self.ecs.ecs_ACRE_simple['Parm']==parm)][pct].item() for pct in ['pct_1','pct_50','pct_100']]
                dt[parm] = ECDistribution(percentiles=v_pcts,values=vs)
            for fip in [int(fip) for fip in list(self.spatialdata.domain_counties['GEOID'].unique())]:
                if fip not in self.ecs.ecs_ACRE_dists[acp]:
                    self.ecs.ecs_ACRE_dists[acp][fip] = dict()
                self.ecs.ecs_ACRE_dists[acp][fip] = dt

    def get_acp_ecs(self,field_id,acp_type:str,method:str='acre'):
        """Get ECs for each nutrient form / loss pathway"""
        if   method.upper().replace(' ','').find('ACREMEDIAN') != -1: dt_ecs_parm = self._get_acp_ecs_acremedian(field_id=field_id,acp_type=acp_type)
        elif method.upper().replace(' ','').find('ACRE')       != -1: dt_ecs_parm = self._get_acp_ecs_acre(field_id=field_id,acp_type=acp_type)
        else: sys.exit('ERROR unrecognized method type')
        return dt_ecs_parm

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
    
    def _set_acp_types(self):
        """Determine types of ACPs available in the data"""
        if   isinstance(self.ecs,type(ACPDataECs.ECDataACRE())):
            self.acp_types = list(set(list(self.ecs.ecs_ACRE['ACP'].unique())+list(self.ecs.ecs_ACRE_simple['ACP'].unique())))
        elif isinstance(self.ecs,type(ACPDataECs.ECDataSWAT())):
            pass
