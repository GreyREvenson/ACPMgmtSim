import copy,random,sys,numpy,matplotlib
import ACPNamelist,ACPDataSpatial,ACPDataNutrientYields,ACPDataEffectivenessCoefficients,ACPDataCosts

class ACPScenario:
    """Class representing a single ACP management trajectory"""

    class ScenarioRecords:
        """Individual ACP attributes. List index = ACP number. Note lists will later be converted to numpy array"""

        def __init__(self):
            """init records"""
            self.acp_type                   = list()
            self.county_fip                 = list()
            self.field_id                   = list()
            self.field_area_ha              = list()
            self.treated_area_ha            = list()
            self.treated_area_ha_cumulative = list()
            self.baseline_TN_kg             = list()
            self.delta_TN_kg                = list()
            self.delta_TN_kg_cumulative     = list()
            self.baseline_TP_kg             = list()
            self.delta_TP_kg                = list()
            self.delta_TP_kg_cumulative     = list()
            self.cost_usd                   = list()

    class NutrientFormPathwayRecords:
        """Nutrient form and pathway load reduction calculation records - acp_number is related to the list index in ScenarioRecords"""

        def __init__(self):
            """init records"""
            self.acp_number                 = list()
            self.nutrient_form_pathway      = list()
            self.baseline_yield_kgha        = list()
            self.baseline_load_kg           = list()
            self.ec                         = list()
            self.delta_load_kg              = list()

    class CostRecords:
        """ACP cost records - acp_number is related to the list index in ScenarioRecords"""

        def __init__(self):
            """init records"""
            self.acp_number                 = list()
            self.cost_actual_usd            = list()
            self.cost_opp_usd               = list()
            self.cost_constructed_area_ha   = list()
            self.cost_actual_usd_cumulative = list()
            self.cost_opp_usd_cumulative    = list()
            self.cost_constructed_area_ha_cumulative = list()

    def __init__(self,name:str,
                 acpnamelist:ACPNamelist.ACPNamelist,
                 acpspatial=ACPDataSpatial.ACPDataSpatial,
                 baseline_data=ACPDataNutrientYields.ACPDataNutrientYields,
                 ec_data=ACPDataEffectivenessCoefficients.ACPDataECs,
                 cost_data=ACPDataCosts.ACPDataCosts,
                 n_acps:int=100):
        """Initialize a ACP trajectory"""
        self.name = name
        self._init_vars() 
        self._evaluate(acpnamelist=acpnamelist,
                       acpspatial=acpspatial,
                       baseline_data=baseline_data,
                       ec_data=ec_data,
                       cost_data=cost_data,
                       n_acps=n_acps)
        self._post_processing()

    def _init_vars(self):
        """init vars"""
        self.scenario_records              = ACPScenario.ScenarioRecords()
        self.nutrient_form_pathway_records = ACPScenario.NutrientFormPathwayRecords()
        self.cost_records                  = ACPScenario.CostRecords()
    
    def _evaluate(self,acpnamelist:ACPNamelist.ACPNamelist,
                  acpspatial:ACPDataSpatial.ACPDataSpatial,
                  baseline_data:ACPDataNutrientYields.ACPDataNutrientYields,
                  ec_data:ACPDataEffectivenessCoefficients.ACPDataECs,
                  cost_data:ACPDataCosts.ACPDataCosts,
                  n_acps:int=100):
        """Evaluate the ACP trajectory"""
        field_ids = copy.deepcopy(acpspatial.miscstats.field_ids)
        for i_project in range(n_acps):
            field_id = random.choice(field_ids)
            acp_type = random.choice(acpnamelist.vars.acp_types)
            self._calc_load_reduction(acp_number=i_project,field_id=field_id,acp_type=acp_type,baseline_nutrient_yields=baseline_data,acp_ecs=ec_data,acp_spatial=acpspatial)
            self._calc_cost(acp_number=i_project,field_id=field_id,acp_type=acp_type,acp_data_spatial=acpspatial,acp_data_costs=cost_data)

    def _post_processing(self):
        """Post processing"""
        self._convert_arrays()
        self._calc_cumulative_vars()

    def _convert_arrays(self):
        """Convert trajectory lists to numpy arrays after all append calls - numpy.append is time expensive"""
        self.scenario_records.acp_type                           = numpy.array(self.scenario_records.acp_type)
        self.scenario_records.county_fip                         = numpy.array(self.scenario_records.county_fip)
        self.scenario_records.field_id                           = numpy.array(self.scenario_records.field_id)
        self.scenario_records.field_area_ha                      = numpy.array(self.scenario_records.field_area_ha)
        self.scenario_records.treated_area_ha                    = numpy.array(self.scenario_records.treated_area_ha)
        self.scenario_records.baseline_TN_kg                     = numpy.array(self.scenario_records.baseline_TN_kg)
        self.scenario_records.delta_TN_kg                        = numpy.array(self.scenario_records.delta_TN_kg)
        self.scenario_records.baseline_TP_kg                     = numpy.array(self.scenario_records.baseline_TP_kg)
        self.scenario_records.delta_TP_kg                        = numpy.array(self.scenario_records.delta_TP_kg)
        self.scenario_records.cost_usd                           = numpy.array(self.scenario_records.cost_usd)
        self.nutrient_form_pathway_records.acp_number            = numpy.array(self.nutrient_form_pathway_records.acp_number)
        self.nutrient_form_pathway_records.nutrient_form_pathway = numpy.array(self.nutrient_form_pathway_records.nutrient_form_pathway)
        self.nutrient_form_pathway_records.baseline_yield_kgha   = numpy.array(self.nutrient_form_pathway_records.baseline_yield_kgha)
        self.nutrient_form_pathway_records.baseline_load_kg      = numpy.array(self.nutrient_form_pathway_records.baseline_load_kg)
        self.nutrient_form_pathway_records.ec                    = numpy.array(self.nutrient_form_pathway_records.ec)
        self.nutrient_form_pathway_records.delta_load_kg         = numpy.array(self.nutrient_form_pathway_records.delta_load_kg)
        self.cost_records.acp_number                             = numpy.array(self.cost_records.acp_number)
        self.cost_records.cost_actual_usd                        = numpy.array(self.cost_records.cost_actual_usd)
        self.cost_records.cost_opp_usd                           = numpy.array(self.cost_records.cost_opp_usd)
        self.cost_records.cost_constructed_area_ha               = numpy.array(self.cost_records.cost_constructed_area_ha)
    
    def _calc_cumulative_vars(self):
        """Calculate cumulative vars after trajectory is simulated"""
        self.scenario_records.treated_area_ha_cumulative      = numpy.cumsum(self.scenario_records.treated_area_ha)
        self.scenario_records.delta_TN_kg_cumulative          = numpy.cumsum(self.scenario_records.delta_TN_kg)
        self.scenario_records.delta_TP_kg_cumulative          = numpy.cumsum(self.scenario_records.delta_TP_kg)
        self.cost_records.cost_actual_usd_cumulative          = numpy.cumsum(self.cost_records.cost_actual_usd)
        self.cost_records.cost_opp_usd_cumulative             = numpy.cumsum(self.cost_records.cost_opp_usd)
        self.cost_records.cost_constructed_area_ha_cumulative = numpy.cumsum(self.cost_records.cost_constructed_area_ha)

    def _calc_cost(self,acp_number:int,field_id:str,acp_type:str,acp_data_spatial:ACPDataSpatial.ACPDataSpatial,acp_data_costs:ACPDataCosts.ACPDataCosts):
        """Calculate ACP cost"""
        fieldareaha = acp_data_spatial.miscstats.dt_areaha_key_fieldid[field_id]
        cost_actual_usd_ha, cost_opportunity_usd_ha = acp_data_costs.get_cost(acp_type=acp_type)
        self.cost_records.acp_number.append(acp_number)
        self.cost_records.cost_actual_usd.append(cost_actual_usd_ha*fieldareaha)          
        self.cost_records.cost_opp_usd.append(cost_opportunity_usd_ha*fieldareaha)
        self.cost_records.cost_constructed_area_ha.append(fieldareaha)
            
    def _calc_load_reduction(self,acp_number:int,field_id:str,acp_type:str,baseline_nutrient_yields:ACPDataNutrientYields.ACPDataNutrientYields,acp_ecs:ACPDataEffectivenessCoefficients.ACPDataECs,acp_spatial:ACPDataSpatial.ACPDataSpatial):
        """Calculate baseline and ACP affected change in TN and TP"""
        baseline_tn_kg, baseline_tp_kg = 0., 0.
        delta_tn_kg, delta_tp_kg = 0.,0.
        dtbasekgha = baseline_nutrient_yields.get_baseline_yields(field_id=field_id)
        dfacpec = acp_ecs.get_acp_ecs(field_id=field_id,acp_type=acp_type) 
        fieldareaha = acp_spatial.miscstats.dt_areaha_key_fieldid[field_id]
        for Parm in baseline_nutrient_yields.nutrientdata.n_vars:
            ibaseline_tn_kg = dtbasekgha[Parm]*fieldareaha
            idelta_tn_kg = ibaseline_tn_kg*dfacpec[Parm]
            baseline_tn_kg += ibaseline_tn_kg
            delta_tn_kg += idelta_tn_kg
            self.nutrient_form_pathway_records.acp_number.append(acp_number)
            self.nutrient_form_pathway_records.nutrient_form_pathway.append(Parm)
            self.nutrient_form_pathway_records.baseline_yield_kgha.append(dtbasekgha[Parm])
            self.nutrient_form_pathway_records.baseline_load_kg.append(ibaseline_tn_kg)
            self.nutrient_form_pathway_records.ec.append(dfacpec[Parm])
            self.nutrient_form_pathway_records.delta_load_kg.append(idelta_tn_kg)
        for Parm in baseline_nutrient_yields.nutrientdata.p_vars:
            ibaseline_tp_kg = dtbasekgha[Parm]*fieldareaha
            idelta_tp_kg = ibaseline_tp_kg*dfacpec[Parm]
            baseline_tp_kg += ibaseline_tp_kg
            delta_tp_kg += idelta_tp_kg
            self.nutrient_form_pathway_records.acp_number.append(acp_number)
            self.nutrient_form_pathway_records.nutrient_form_pathway.append(Parm)
            self.nutrient_form_pathway_records.baseline_yield_kgha.append(dtbasekgha[Parm])
            self.nutrient_form_pathway_records.baseline_load_kg.append(ibaseline_tp_kg)
            self.nutrient_form_pathway_records.ec.append(dfacpec[Parm])
            self.nutrient_form_pathway_records.delta_load_kg.append(idelta_tp_kg)
        self.scenario_records.acp_type.append(acp_type)
        self.scenario_records.county_fip.append(acp_spatial.miscstats.dt_fips_key_fieldid[field_id])
        self.scenario_records.field_id.append(field_id)
        self.scenario_records.field_area_ha.append(acp_spatial.miscstats.dt_areaha_key_fieldid[field_id])
        self.scenario_records.treated_area_ha.append(acp_spatial.miscstats.dt_areaha_key_fieldid[field_id])
        self.scenario_records.baseline_TN_kg.append(baseline_tn_kg)
        self.scenario_records.delta_TN_kg.append(delta_tn_kg)
        self.scenario_records.baseline_TP_kg.append(baseline_tp_kg)
        self.scenario_records.delta_TP_kg.append(delta_tp_kg)
        
class ACPScenarioSet:
    """Class representing a set of multiple ACP management scenarios"""
    def __init__(self,name:str,acpnamelist:ACPNamelist.ACPNamelist,
                 acpspatial:ACPDataSpatial.ACPDataSpatial,
                 baseline_data:ACPDataNutrientYields.ACPDataNutrientYields,
                 ec_data:ACPDataEffectivenessCoefficients.ACPDataECs,
                 cost_data:ACPDataCosts.ACPDataCosts,
                 n_scenarios:int=1,n_acps:int=100):
        """Initialize set of ACP scenarios"""
        self.name = name
        self.domain_areaha = acpspatial.miscstats.domain_area_ha
        self.bucket = list()
        for i in range(n_scenarios):
            self.bucket.append(ACPScenario(name=name+str(i+1),acpnamelist=acpnamelist,acpspatial=acpspatial,baseline_data=baseline_data,ec_data=ec_data,cost_data=cost_data,n_acps=n_acps))

    def graph_delta_tn_tp(self,acpspatial:ACPDataSpatial.ACPDataSpatial,baseline_data:ACPDataNutrientYields.ACPDataNutrientYields,x_var:str='count',y_var:str='kg yr',basin_out_rch_num:int=0):
        """Make scenario line graph"""
        xs, xtitle = list(), ''
        if x_var == 'count': 
            for acp_scenario in self.bucket:
                _xs = numpy.arange(start=0,stop=len(acp_scenario.scenario_records.delta_TN_kg_cumulative)+1,dtype=acp_scenario.scenario_records.delta_TN_kg_cumulative.dtype)
                #_xs = numpy.insert(_xs,0,0,0)
                xs.append(_xs)
            xtitle = 'Number of ACPs implemented'
        elif x_var == 'cost':
            for acp_scenario in self.bucket:
                _xs = acp_scenario.cost_records.cost_actual_usd_cumulative + acp_scenario.cost_records.cost_opp_usd_cumulative
                _xs = _xs / 1000000.
                _xs = numpy.insert(_xs,0,0,0)
                xs.append(_xs)
            xtitle = 'Cumulative cost (millions of usd)'
        elif x_var == 'treated area':
            for acp_scenario in self.bucket:
                _xs = (acp_scenario.scenario_records.treated_area_ha_cumulative/acpspatial.miscstats.domain_area_ha)*100.
                _xs = numpy.insert(_xs,0,0,0)
                xs.append(_xs)
            xtitle = 'Cumulative treated area (% of Watershed Area)'
        else:
            sys.exit('ERROR unrecognized x axis variable type')
        ys_TN, ys_TP, ytitle_TN, ytitle_TP = list(), list(), '', ''
        if y_var == '% mean annual load':
            if not hasattr(baseline_data.nutrientdata,'mean_annual_TN_basin_load_kgyr') or not hasattr(baseline_data.nutrientdata,'mean_annual_TP_basin_load_kgyr'):
                sys.exit('ERROR missing mean annual load data')
            if basin_out_rch_num == 0:
                sys.exit('ERROR missing basin outlet reach number')
            mean_annual_TN_basin_load_kgyr, mean_annual_TP_basin_load_kgyr = baseline_data.get_swat_mean_annual_loads(basin_out_rch_num)
            baseline_data.get_swat_mean_annual_loads(basin_out_rch_num)
            for acp_scenario in self.bucket:
                _ys_TN = (acp_scenario.scenario_records.delta_TN_kg_cumulative/mean_annual_TN_basin_load_kgyr)*100.
                _ys_TN = numpy.insert(_ys_TN,0,0,0)
                ys_TN.append(_ys_TN)
                _ys_TP = (acp_scenario.scenario_records.delta_TP_kg_cumulative/mean_annual_TP_basin_load_kgyr)*100.
                _ys_TP = numpy.insert(_ys_TP,0,0,0)
                ys_TP.append(_ys_TP)
            ytitle_TN = 'Cumulative TN load reduction\n(% of SWAT simulated mean annual load)'
            ytitle_TP = 'Cumulative TP load reduction\n(% of SWAT simulated mean annual load)'
        elif y_var == 'kg yr':
            for acp_scenario in self.bucket:
                _ys_TN = acp_scenario.scenario_records.delta_TN_kg_cumulative
                _ys_TN = numpy.insert(_ys_TN,0,0,0)
                ys_TN.append(_ys_TN)
                _ys_TP = acp_scenario.scenario_records.delta_TP_kg_cumulative
                _ys_TP = numpy.insert(_ys_TP,0,0,0)
                ys_TP.append(_ys_TP)
            ytitle_TN = 'Cumulative TN load reduction (kg yr)'
            ytitle_TP = 'Cumulative TP load reduction (kg yr)'
        else:
            sys.exit('ERROR unrecognized y axis variable type')
        fig, (ax1,ax2) = matplotlib.pyplot.subplots(1, 2)
        fig.set_size_inches(8.5, 5)
        for i in range(len(xs)):
            ax1.plot(xs[i],ys_TN[i]*-1.)
            ax2.plot(xs[i],ys_TP[i]*-1.)
        ax1.set_xlabel(xtitle)
        ax1.set_ylabel(ytitle_TN)
        if y_var == 'kg yr':
            ax1.axhline(y=baseline_data.nutrientdata.mean_annual_TN_basin_load_kgyr*-1.) #color='r', linestyle='--', linewidth=2, label='y = 0.5')
        ax2.set_xlabel(xtitle)
        ax2.set_ylabel(ytitle_TP)
        if y_var == 'kg yr':
            ax2.axhline(y=baseline_data.nutrientdata.mean_annual_TP_basin_load_kgyr*-1.)
        fig.suptitle('Simulated ACP impacts on TN and TP loads at East Fork watershed outlet')
        fig.tight_layout()
        return fig

    def graph_baseline_boxplots(self):
        """Make baseline yield boxplots"""
        fig, (ax1,ax2) = matplotlib.pyplot.subplots(1, 2)
        fig.set_size_inches(8.5, 2.5)
        xs_tn_kgha, xs_tp_kgha = list(), list()
        for acp_scenario in self.bucket:
            _xs_tn_kgha = acp_scenario.scenario_records.baseline_TN_kg/acp_scenario.scenario_records.field_area_ha
            _xs_tp_kgha = acp_scenario.scenario_records.baseline_TP_kg/acp_scenario.scenario_records.field_area_ha
            xs_tn_kgha = numpy.append(xs_tn_kgha,_xs_tn_kgha)
            xs_tp_kgha = numpy.append(xs_tp_kgha,_xs_tp_kgha)
        ax1.boxplot(xs_tn_kgha,vert=False)
        ax1.set_title('Baseline TN (kg/ha/yr)')
        ax1.set_yticklabels('')
        ax2.boxplot(xs_tp_kgha,vert=False)
        ax2.set_title('Baseline TP (kg/ha/yr)')
        ax2.set_yticklabels('')
        fig.suptitle('Field-scale baseline nutrient yields across all management trajectories')
        fig.tight_layout()
        return fig
    
    def graph_effective_efficiency_boxplots(self):
        """Make baseline yield boxplots"""
        fig, (ax1,ax2) = matplotlib.pyplot.subplots(1, 2)
        fig.set_size_inches(8.5, 2.5)
        xs_tn_efec, xs_tp_efec = list(), list()
        for acp_scenario in self.bucket:
            base_tn_kg = numpy.where(acp_scenario.scenario_records.baseline_TN_kg==0,1,acp_scenario.scenario_records.baseline_TN_kg)
            base_tp_kg = numpy.where(acp_scenario.scenario_records.baseline_TP_kg==0,1,acp_scenario.scenario_records.baseline_TP_kg)
            _xs_tn_efec = (acp_scenario.scenario_records.delta_TN_kg/base_tn_kg)*-100.
            _xs_tp_efec = (acp_scenario.scenario_records.delta_TP_kg/base_tp_kg)*-100.
            xs_tn_efec = numpy.append(xs_tn_efec,_xs_tn_efec)
            xs_tp_efec = numpy.append(xs_tp_efec,_xs_tp_efec)
        ax1.boxplot(xs_tn_efec,vert=False)
        ax1.set_title('Change in TN (% annual baseline load)')
        ax1.set_yticklabels('')
        ax2.boxplot(xs_tp_efec,vert=False)
        ax2.set_title('Change TP (% annual baseline load)')
        ax2.set_yticklabels('')
        fig.suptitle('Field-scale ACP-affected load reductions across all management trajectories')
        fig.tight_layout()
        return fig
    
    def graph_acp_piechart(self):
        """Create pie chart of ACP selections"""
        dt = dict()
        for acp_scenarioectory in self.bucket:
            acps, counts = numpy.unique(acp_scenarioectory.scenario_records.acp_type, return_counts=True)
            for i in range(len(acps)):
                if acps[i] not in dt:
                    dt[acps[i]] = list()
                dt[acps[i]].append(counts[i])
        _len = max([len(dt[acp]) for acp in dt])
        for acp in dt:
            if len(dt[acp]) != _len:
                for _ in range(_len-len(dt[acp])):
                    dt[acp].append(0)
            dt[acp] = numpy.mean(dt[acp])
        fig, ax = matplotlib.pyplot.subplots(1, 1)
        fig.set_size_inches(8.5, 5)
        name_chg = [['Constructed_Wetland','Constructed wetland'],
                    ['Cover_Crop','Cover crops'],
                    ['Drain_WM','Drainage water management'],
                    ['Fert 75%','Fertilizer reduction (75%)'],
                    ['Fert 90%','Fertilizer reduction (90%)'],
                    ['MIN_TILL','Minimum tillage'],
                    ['Riparian_buffer','Riparian buffer'],
                    ['Ponds 50%','Ponds (50%)'],
                    ['Ponds 75%','Ponds (75%)'],
                    ['Contour Farming','Contour farming'],
                    ['Conservation_Crop_Rotation','Conservation crop rotation'],
                    ['CONSERVATION','General conservation'],
                    ['Blind_Inlet','Blind inlet'],
                    ['two_state_ditch','Two stage ditch'],
                    ['Waterway Only','Grassed waterway'],
                    ['Terraces Only','Terraces'],
                    ['Terraces And Waterway','Terraces and grassed waterways']]
        for name_pair in name_chg:
            if name_pair[1] in dt:
                dt[name_pair[1]] = dt[name_pair[0]]
                del dt[name_pair[0]]
        ax.pie(dt.values(),labels=[key+' ('+"{:.2f}".format(dt[key])+')' for key in dt])
        ax.set_title('Average number of ACPs per trajectory')
        fig.tight_layout()
        return fig
