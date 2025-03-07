import os,sys,numpy

class ACPNamelist:
    """Class to hold names and user inputs"""

    class DirectoryNames:
        """Directory names"""
        project         = numpy.nan
        acre_data       = numpy.nan
        spatial         = numpy.nan
        swat            = numpy.nan 
        output          = numpy.nan 

    class FileNames:
        """File names including full path"""
        namelist           = numpy.nan
        hucs               = numpy.nan
        domain             = numpy.nan
        counties           = numpy.nan
        fields             = numpy.nan
        hrus               = numpy.nan
        field_hrugis       = numpy.nan
        swat_output_hru    = numpy.nan 
        swat_output_rch    = numpy.nan
        cost               = numpy.nan
        
    class VariableNames:
        """Variable names"""
        area_ha = 'Area_ha'

    class Variables:
        """Variables"""
        overwrite_flag        = False
        verbose               = False
        input_huc12_ids       = numpy.nan
        input_hru_shps        = numpy.nan
        input_swat_output_hru = numpy.nan
        input_swat_output_rch = numpy.nan
        input_fields          = numpy.nan
        file_inputs           = numpy.nan
        acp_types             = numpy.nan

    def __init__(self,fname_namelist:str):
        """Initialize namelist"""
        self._init_vars()
        self._read_namelist(fname_namelist)
        self._set_user_inputs()
        self._set_names()

    def _init_vars(self):
        """Initialize variables"""
        self.dirnames = ACPNamelist.DirectoryNames()
        self.fnames   = ACPNamelist.FileNames()
        self.varnames = ACPNamelist.VariableNames()
        self.vars     = ACPNamelist.Variables()
        
    def _set_names(self):
        """Set static directory and file names"""
        self._set_subdirectory_and_file_names()
        self._set_static_file_names()

    def _set_subdirectory_and_file_names(self):
        """Set project subdirectory names"""
        self.dirnames.spatial = os.path.join(self.dirnames.project,'spatial')
        self.dirnames.swat    = os.path.join(self.dirnames.project,'swat')
        self.dirnames.output  = os.path.join(self.dirnames.project,'outputs')

    def _set_static_file_names(self):
        """Set static files names for intermediate output files"""
        self.fnames.domain       = os.path.join(self.dirnames.spatial,   'domain.gpkg')
        self.fnames.hucs         = os.path.join(self.dirnames.spatial,   'hucs.gpkg')
        self.fnames.counties     = os.path.join(self.dirnames.spatial,   'counties.gpkg')
        self.fnames.fields       = os.path.join(self.dirnames.spatial,   'fields.gpkg')
        self.fnames.hrus         = os.path.join(self.dirnames.swat,      'hrus.gpkg')

    def _remove_whitespace_outside_quotes(self,line:str):
        result = []
        in_quote = False
        quote_char = None
        for char in line:
            if char in ('"', "'"):
                if in_quote and quote_char == char:
                    in_quote = False
                else:
                    in_quote = True
                    quote_char = char
                result.append(char)
            elif not in_quote and char.isspace():
                continue
            else:
                result.append(char)
        return ''.join(result)

    def _read_namelist(self,fname_namelist:str):
        """Read namelist file into generic dictionary"""
        self.fnames.namelist = fname_namelist
        if not os.path.isfile(self.fnames.namelist):
            sys.exit('ERROR could not find namelist file '+self.fnames.namelist)
        self.vars.file_inputs = dict()
        namelist_lines = list(open(self.fnames.namelist,'r'))
        for l in namelist_lines:
            try:
                l0 = self._remove_whitespace_outside_quotes(line=l)
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
                        self.vars.file_inputs[var_name] = var_vals
            except:
                sys.exit('ERROR could not read namelist.txt line: '+l)
        for var_name in self.vars.file_inputs:
            if isinstance(self.vars.file_inputs[var_name],list) and len(self.vars.file_inputs[var_name]) == 1:
                self.vars.file_inputs[var_name] = self.vars.file_inputs[var_name][0]

    def _set_user_inputs(self):
        """Set variables using read-in values"""
        name_project_dir      = 'project_directory'
        name_huc12s           = 'huc12_ids'
        name_acre_dir         = 'acre_directory'
        name_swat_output_hru  = 'swat_output_hru_files'
        name_field_bounds     = 'field_boundary_directory'
        name_overwrite        = 'intermediate_overwrite'
        name_verbose          = 'verbose'
        name_cost_file        = 'cost_file'
        name_acp_types        = 'acp_types'
        name_swat_output_hru  = 'swat_output_hru_files'
        name_swat_output_rch  = 'swat_output_rch_files'
        req = [name_huc12s,name_project_dir,name_acre_dir,name_swat_output_hru,name_field_bounds,name_cost_file]
        for name in req:
            if name not in self.vars.file_inputs: sys.exit('ERROR required variable '+name+' not found in namelist file')
        self.vars.input_huc12_ids = self.vars.file_inputs[name_huc12s]
        if name_project_dir in self.vars.file_inputs:
            self.dirnames.project = os.path.abspath(self.vars.file_inputs[name_project_dir])
        if name_acre_dir in self.vars.file_inputs:  
            self.dirnames.acre_data = os.path.abspath(self.vars.file_inputs[name_acre_dir])
        if name_field_bounds in self.vars.file_inputs:
            self.vars.input_fields = os.path.abspath(self.vars.file_inputs[name_field_bounds])
        if name_cost_file in self.vars.file_inputs:
            if isinstance(self.vars.file_inputs[name_cost_file],str): self.fnames.cost = os.path.abspath(self.vars.file_inputs[name_cost_file])
            else: self.fnames.cost = [os.path.abspath(fname) for fname in self.vars.file_inputs[name_cost_file]]
        if name_acp_types in self.vars.file_inputs:
            if isinstance(self.vars.file_inputs[name_acp_types],str): self.vars.acp_types = [self.vars.file_inputs[name_acp_types]]
            else: self.vars.acp_types = self.vars.file_inputs[name_acp_types]
        if name_swat_output_hru in self.vars.file_inputs:
            if isinstance(self.vars.file_inputs[name_swat_output_hru],str): self.vars.input_swat_output_rch = [os.path.abspath(self.vars.file_inputs[name_swat_output_hru])]
            else: self.vars.input_swat_output_hru = [os.path.abspath(fname) for fname in self.vars.file_inputs[name_swat_output_hru]]
        if name_swat_output_rch in self.vars.file_inputs:
            if isinstance(self.vars.file_inputs[name_swat_output_rch],str): self.vars.input_swat_output_rch = [os.path.abspath(self.vars.file_inputs[name_swat_output_rch])]
            else: self.vars.input_swat_output_rch = [os.path.abspath(fname) for fname in self.vars.file_inputs[name_swat_output_rch]]
        if name_overwrite in self.vars.file_inputs and self.vars.file_inputs[name_overwrite].upper().find('TRUE') != -1:
            self.vars.overwrite_flag = True
        if name_verbose in self.vars.file_inputs and self.vars.file_inputs[name_verbose].upper().find('TRUE') != -1:
            self.vars.verbose = True