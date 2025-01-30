import ACPTrajectoryData,copy,random

class ACPTrajectory:
    """Class representing a single ACP management trajectory"""

    class TrajectoryRecords:
        """Individual ACP attributes. List index = ACP number"""
        acp_type                   = list()
        county_fip                 = list()
        field_id                   = list()
        field_area_ha              = list()
        treated_area_ha            = list()
        treated_area_ha_cumulative = list()
        baseline_TN_kg             = list()
        delta_TN_kg                = list()
        delta_TN_kg_cumulative     = list()
        baseline_TP_kg             = list()
        delta_TP_kg                = list()
        delta_TP_kg_cumulative     = list()
        cost_usd                   = list()

    class NutrientFormPathwayRecords:
        """Nutrient form and pathway load reduction calculation records - acp_number is related to the list index in TrajectoryRecords"""
        acp_number                 = list()
        nutrient_form_pathway      = list()
        baseline_yield_kgha        = list()
        baseline_load_kg           = list()
        ec                         = list()
        delta_load_kg              = list()
        
    def __init__(self,name:str,acpdata:ACPTrajectoryData.ACPData,baseline_method:str='acre pdf',ec_method:str='acre',n_acps:int=100):
        """Initialize a ACP trajectory"""
        self.name = name 
        self._evaluate(acpdata,baseline_method,ec_method,n_acps)
        
    def _evaluate(self,acpdata:ACPTrajectoryData.ACPData,baseline_method:str,ec_method:str,n_acps:int):
        """Evaluate the ACP trajectory"""
        self.trajectory_records = ACPTrajectory.TrajectoryRecords()
        self.nutrient_form_pathway_records = ACPTrajectory.NutrientFormPathwayRecords()
        field_ids = copy.deepcopy(acpdata.miscstats.field_ids)
        for i_project in range(n_acps):
            field_id = random.choice(field_ids)
            acp_type = random.choice(acpdata.miscstats.acp_types)
            baseline_tn_kg,delta_tn_kg,baseline_tp_kg,delta_tp_kg = self._calc_load_readuction(acp_number=i_project,acpdata=acpdata,field_id=field_id,
                                                                                               acp_type=acp_type,baseline_method=baseline_method,ec_method=ec_method)
            self.trajectory_records.acp_type.append(acp_type)
            self.trajectory_records.county_fip.append(acpdata.miscstats.dt_fips_key_fieldid[field_id])
            self.trajectory_records.field_id.append(field_id)
            self.trajectory_records.field_area_ha.append(acpdata.miscstats.dt_areaha_key_fieldid[field_id])
            self.trajectory_records.treated_area_ha.append(acpdata.miscstats.dt_areaha_key_fieldid[field_id])
            self.trajectory_records.baseline_TN_kg.append(baseline_tn_kg)
            self.trajectory_records.delta_TN_kg.append(delta_tn_kg)
            self.trajectory_records.baseline_TP_kg.append(baseline_tp_kg)
            self.trajectory_records.delta_TP_kg.append(delta_tp_kg)

    def _calc_cost(self,acpdata:ACPTrajectoryData.ACPData,field_id:str,acp_type:str):
        """Calculate ACP cost"""
        pass
            
    def _calc_load_readuction(self,acp_number:int,acpdata:ACPTrajectoryData.ACPData,field_id:str,acp_type:str,baseline_method:str='acre pdf',ec_method:str='acre'):
        """Calculate baseline and ACP affected change in TN and TP"""
        delta_tn_kg = 0.
        delta_tp_kg = 0.
        baseline_tn_kg = 0.
        baseline_tp_kg = 0.
        dtbasekgha = acpdata.get_baseline_yields(field_id=field_id,method=baseline_method)
        dfacpec = acpdata.get_acp_ecs(field_id=field_id,acp_type=acp_type,method=ec_method) 
        fieldareaha = acpdata.miscstats.dt_areaha_key_fieldid[field_id]
        for Parm in acpdata.swatoutvars.n_vars:
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
        for Parm in acpdata.swatoutvars.p_vars:
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
        return baseline_tn_kg,delta_tn_kg,baseline_tp_kg,delta_tp_kg
        
class ACPTrajectorySet:
    """Class representing a set of multiple ACP management trajectories"""
    def __init__(self,name:str,acpdata:ACPTrajectoryData.ACPData,baseline_method:str='acre pdf',ec_method:str='acre',n_trajectories=1,n_acps:int=100):
        """Initialize set of ACP trajectories"""
        self.name = name
        self.traj_buck = list()
        for i in range(len(n_trajectories)):
            self.traj_buck.append(ACPTrajectory(name=name+str(i+1),acpdata=acpdata,baseline_method=baseline_method,ec_method=ec_method,n_acps=n_acps))