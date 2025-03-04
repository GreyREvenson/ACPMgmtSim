import os,sys,numpy,pandas,random
from src import ACPNamelist,ACPDataSpatial

class ACPDataCosts:
    """Class to hold the data"""

    class CustomPriceList:
        """User supplied price list"""
        raw = numpy.nan

    class EQIP:
        """placeholder until raw eqip data can be found/processed"""
        raw = numpy.nan

    class NRCSPracticeCodeMap:
        """"""
        dt = numpy.nan

    def __init__(self,acpnamelist:ACPNamelist.ACPNamelist,acpspatial:ACPDataSpatial.ACPDataSpatial,method:str='custom'):
        """Initialize"""
        if acpnamelist.vars.verbose: print('Initializing instance of ACPDataCost')
        self._init_vars(acpnamelist,method)
        self._read_data(acpnamelist,acpspatial,method)

    def _init_vars(self,acpnamelist:ACPNamelist.ACPNamelist,method:str):
        """Initialize variables"""
        if acpnamelist.vars.verbose: print('Initializing variables')
        if method.upper().find('CUSTOM') != -1: self.costs = ACPDataCosts.CustomPriceList()
        elif method.upper().find('EQIP') != -1: self.costs = ACPDataCosts.EQIP() 
        self.nrcs_practice_code_map = ACPDataCosts.NRCSPracticeCodeMap()

    def _read_data(self,acpnamelist:ACPNamelist.ACPNamelist,acpspatial:ACPDataSpatial.ACPDataSpatial,method:str):
        """Read the data"""
        if acpnamelist.vars.verbose: print('Reading practice cost data')
        if method.upper().find('CUSTOM') != -1: self._read_data_custom(acpnamelist)
        else: sys.exit('ERROR ACPDataCost._read_data unknown data type')

    def _read_data_custom(self,acpnamelist:ACPNamelist.ACPNamelist):
        """Read custom cost info"""
        if acpnamelist.vars.verbose: print('    Reading: '+acpnamelist.fnames.cost)
        if not os.path.isfile(acpnamelist.fnames.cost): sys.exit('ERROR could not find cost file '+acpnamelist.fnames.cost)
        self.costs.raw = pandas.read_csv(acpnamelist.fnames.cost,low_memory=False)
        self._prep_custom_price_info()

    def _prep_custom_price_info(self):
        """Function to process custom price info"""

        # estimate missing cost values
        nan_indices = numpy.argwhere(numpy.isnan(self.costs.raw['Actual_usd_ac']))
        price_eqip = numpy.delete(self.costs.raw['EQIP_usd_ac'],nan_indices)
        price_actual = numpy.delete(self.costs.raw['Actual_usd_ac'],nan_indices)
        slope = numpy.polyfit(price_eqip,price_actual,1)[0].item()
        self.costs.raw['Actual_usd_ac'] = numpy.where(numpy.isnan(self.costs.raw['Actual_usd_ac']),slope*self.costs.raw['EQIP_usd_ac'],self.costs.raw['Actual_usd_ac'])
        
        # convert ac to ha
        constant_ac_to_ha = 0.404686
        self.costs.raw['EQIP_usd_ac'] = self.costs.raw['EQIP_usd_ac'] * constant_ac_to_ha
        self.costs.raw['Actual_usd_ac'] = self.costs.raw['Actual_usd_ac'] * constant_ac_to_ha
        self.costs.raw['Foregone_income_usd_ac'] = self.costs.raw['Foregone_income_usd_ac'] * constant_ac_to_ha
        self.costs.raw.rename(columns={'EQIP_usd_ac': 'EQIP_usd_ha','Actual_usd_ac': 'Actual_usd_ha','Foregone_income_usd_ac': 'Foregone_income_usd_ha'},inplace=True)

    def _read_data_eqip(self):
        """Placeholder for when raw eqip data can be ingested"""
        pass

    def set_name_to_practicecode_map(self,dt:dict=None):
        """Set mapping between ACP type and NRCS practice codes"""
        if dt is not None:
            self.nrcs_practice_code_map.dt = dt
        else:
            self.nrcs_practice_code_map.dt = {'Constructed_Wetland':656,
                                              'Cover_Crop':340,
                                              'MIN_TILL':345,
                                              'Ponds 50%':378,
                                              'Ponds 75%':378,
                                              'Waterway Only':412,
                                              'Filterstrip':393}

    def get_cost(self,acp_type:str):
        """Get cost"""
        if isinstance(self.costs,type(ACPDataCosts.CustomPriceList())): return self._get_cost_custom(acp_type)
        elif isinstance(self.costs,type(ACPDataCosts.EQIP())):          return self._get_cost_eqip(acp_type)

    def _get_cost_custom(self,acp_type):
        """get build and opportunity cost"""
        if not isinstance(self.nrcs_practice_code_map.dt,dict): self.set_name_to_practicecode_map()
        practice_code = self.nrcs_practice_code_map.dt[acp_type]
        practice_rows = self.costs.raw.loc[(self.costs.raw['Practice_code']==practice_code)]
        cost_actual_usd_ha = random.choice(list(practice_rows['EQIP_usd_ha']))
        cost_opportunity_usd_ha = random.choice(list(practice_rows['Actual_usd_ha']))
        return cost_actual_usd_ha, cost_opportunity_usd_ha

    def _get_cost_eqip(self,actp_type):
        """placeholder"""
        pass

        