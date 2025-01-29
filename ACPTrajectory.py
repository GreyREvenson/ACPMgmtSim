import os,sys,numpy,copy,random
from ACPTrajectoryData import ACPData

class ACPTrajectory:

    class TrajectoryRecords:
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
        acp_number                 = list()
        nutrient_form_pathway      = list()
        baseline_yield_kgha        = list()
        baseline_load_kg           = list()
        ec                         = list()
        delta_load_kg              = list()
        
    def __init__(self,name:str,acpdata:ACPData,baseline_method:str='acre pdf',ec_method:str='acre',n_acps:int=100):

        self.name = name 
        self._evaluate(acpdata,baseline_method,ec_method,n_acps)
        
    def _evaluate(self,acpdata:ACPData,baseline_method:str,ec_method:str,n_acps:int):
        """
        Calculate baseline and ACP-affected change in TN and TP for given field and type of ACP
            Args:
                acpdata (ACPTrajectoryData:ACPData)
                baseline_method (string): 'swat pdf', 'swat raw', or 'acre'
                ec_method (string): 'acre' or 'acre median'
                n_acps (int): numpber of ACPs in trajectory
            Returns:
                baseline_tn_kg (float): baseline TN load (kg) (+)
                delta_tn_kg (float): ACP-affected change in TN load (kg) (+ or -)
                baseline_tp_kg (float): baseline TP load (kg) (+)
                delta_tp_kg (float): ACP-affected change in TP load (kg) (+ or -)
                elines (list of lists): calculation record keeping lines
        """
        
        #init empty record arrays
        self.trajectory_records = ACPTrajectory.TrajectoryRecords()
        self.nutrient_form_pathway_records = ACPTrajectory.NutrientFormPathwayRecords()

        #create list of field ids
        field_ids = copy.deepcopy(acpdata.miscstats.field_ids)

        #iterate over acp project count
        for i_project in range(n_acps):

            #randomly select field and ACP type
            field_id = random.choice(field_ids)
            acp_type = random.choice(acpdata.miscstats.acp_types)

            #estimate ACP impacts
            baseline_tn_kg,delta_tn_kg,baseline_tp_kg,delta_tp_kg = self._calc_load_readuction(acp_number=i_project,acpdata=acpdata,field_id=field_id,
                                                                                               acp_type=acp_type,baseline_method=baseline_method,ec_method=ec_method)

            #add to records
            self.trajectory_records.acp_type.append(acp_type)
            self.trajectory_records.county_fip.append(acpdata.miscstats.dt_fips_key_fieldid[field_id])
            self.trajectory_records.field_id.append(field_id)
            self.trajectory_records.field_area_ha.append(acpdata.miscstats.dt_areaha_key_fieldid[field_id])
            self.trajectory_records.treated_area_ha.append(acpdata.miscstats.dt_areaha_key_fieldid[field_id])
            self.trajectory_records.baseline_TN_kg.append(baseline_tn_kg)
            self.trajectory_records.delta_TN_kg.append(delta_tn_kg)
            self.trajectory_records.baseline_TP_kg.append(baseline_tp_kg)
            self.trajectory_records.delta_TP_kg.append(delta_tp_kg)

    def _calc_cost(self,acpdata:ACPData,field_id:str,acp_type:str):

        pass
            
    def _calc_load_readuction(self,acp_number:int,acpdata:ACPData,field_id:str,acp_type:str,baseline_method:str='acre pdf',ec_method:str='acre'):
        """
        Calculate baseline and ACP-affected change in TN and TP for given field and type of ACP
            Args:
                acpdata (ACPTrajectoryData:ACPData)
                field_id (string): field unique identifier
                acp_type (string): type of ACP
                baseline_method (string): 'swat pdf', 'swat raw', or 'acre'
                ec_method (string): 'acre' or 'acre median'
            Returns:
                baseline_tn_kg (float): baseline TN load (kg) (+)
                delta_tn_kg (float): ACP-affected change in TN load (kg) (+ or -)
                baseline_tp_kg (float): baseline TP load (kg) (+)
                delta_tp_kg (float): ACP-affected change in TP load (kg) (+ or -)
                elines (list of lists): calculation record keeping lines
                    header = ['type of acp', 'county fip', 'field id', 'field area (ha)', 
                              'nutrient form/pathway', 'baseline yield (kg/ha)', baseline load (kg)', 
                              'effectiveness coefficient', 'affected change in load (kg)'] 
        """

        #init field-level variables to calculate
        delta_tn_kg = 0.
        delta_tp_kg = 0.
        baseline_tn_kg = 0.
        baseline_tp_kg = 0.

        #get baseline nutrient yields and effectiveness coefficients
        dtbasekgha = acpdata.get_baseline_yields(field_id=field_id,method=baseline_method)
        dfacpec = acpdata.get_acp_ecs(field_id=field_id,acp_type=acp_type,method=ec_method) 

        #get field area (ha)
        fieldareaha = acpdata.miscstats.dt_areaha_key_fieldid[field_id]

        #iterate over N form/pathways
        for Parm in ('NSURQ','ORGN','NLAT','GWQN'):

            #calculate baseline load for form/pathways
            ibaseline_tn_kg = dtbasekgha[Parm]*fieldareaha

            #calculate ACP affected change in load for form/pathway
            idelta_tn_kg = ibaseline_tn_kg*dfacpec[Parm]

            #add to field-level totals
            baseline_tn_kg += ibaseline_tn_kg
            delta_tn_kg += idelta_tn_kg

            #add to form/pathway records
            self.nutrient_form_pathway_records.acp_number.append(acp_number)
            self.nutrient_form_pathway_records.nutrient_form_pathway.append(Parm)
            self.nutrient_form_pathway_records.baseline_yield_kgha.append(dtbasekgha[Parm])
            self.nutrient_form_pathway_records.baseline_load_kg.append(ibaseline_tn_kg)
            self.nutrient_form_pathway_records.ec.append(dfacpec[Parm])
            self.nutrient_form_pathway_records.delta_load_kg.append(idelta_tn_kg)

        #iterate over P form/pathways
        for Parm in ('ORGP','SOLP'):

            #calculate baseline load for form/pathway
            ibaseline_tp_kg = dtbasekgha[Parm]*fieldareaha

            #calculate ACP affected change in load for form/pathway
            idelta_tp_kg = ibaseline_tp_kg*dfacpec[Parm]

            #add to field-level totals
            baseline_tp_kg += ibaseline_tp_kg
            delta_tp_kg += idelta_tp_kg

            #add to form/pathway records
            self.nutrient_form_pathway_records.acp_number.append(acp_number)
            self.nutrient_form_pathway_records.nutrient_form_pathway.append(Parm)
            self.nutrient_form_pathway_records.baseline_yield_kgha.append(dtbasekgha[Parm])
            self.nutrient_form_pathway_records.baseline_load_kg.append(ibaseline_tp_kg)
            self.nutrient_form_pathway_records.ec.append(dfacpec[Parm])
            self.nutrient_form_pathway_records.delta_load_kg.append(idelta_tp_kg)

        #return field-level baseline loads and ACP affected change
        return baseline_tn_kg,delta_tn_kg,baseline_tp_kg,delta_tp_kg
        