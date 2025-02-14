import os,sys,numpy,pandas,ACPTrajectoryNamelist    

class ACPDataCost:
    """Class to hold the data"""

    class EQIP:
        raw = numpy.nan

    class OtherPriceDataSource:
        df = numpy.nan

    def __init__(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist,method:str='eqip'):
        """Initialize"""
        if acpnamelist.vars.verbose: print('Initializing instance of ACPDataCost')
        self._init_vars(acpnamelist)
        self._init_dirs(acpnamelist)
        self._read_spatial_data(acpnamelist)
        self._read_nutrient_and_ec_data(acpnamelist)

    def _init_vars(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist,method:str):
        """Initialize variables"""
        if acpnamelist.vars.verbose: print('Initializing variables')
        self.costs = ACPDataCost.EQIP()

    def _read_data(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist,method:str):
        if method.upper().find('EQIP') != -1:
            self._read_data_eqip(acpnamelist=acpnamelist)
             
    def _read_data_eqip(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
         
        df = pandas.read_csv(acpnamelist.fnames.eqip)
#huc_level	huc_name	program	obligation_fy	practice_code	practice_name	practice_count	dollars_obligated	suppressed

    def get_cost(self,acp_type,acp_area_ha):
        """Get cost"""
        dt = {'Constructed_Wetland':656,
              'Cover_Crop':340,
              'MIN_TILL':345,
              'Ponds 50%':378,
              'Ponds 75%':378,
              'Waterway Only':412,
              'Filterstrip':393
              }
        if acp_type not in dt: sys.exit('ERROR get_cost unrecognized ACP type '+acp_type)
        rws = self.costs.costdata.loc[(self.costs.costdata['Practice_code']==dt[acp_type])]
        rws = rws.sample(1)
        return rws

    def _readCostInfo(self,acpnamelist:ACPTrajectoryNamelist.ACPNamelist):
        """Read cost information"""
        if acpnamelist.vars.verbose: print('Reading cost information')
        if os.path.isfile(acpnamelist.fnames.cost):
            if acpnamelist.vars.verbose: print('    Reading: '+acpnamelist.fnames.cost)
            self.costs.costdata = pandas.read_csv(acpnamelist.fnames.cost)
        else:
            print('ERROR unable to read cost file '+acpnamelist.fnames.cost)
        