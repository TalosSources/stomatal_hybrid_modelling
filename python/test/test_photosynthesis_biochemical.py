from python.photosynthesis_biochemical import photosynthesis_biochemical
import torch

import inspect

pb_params = set(inspect.signature(photosynthesis_biochemical).parameters)

def compare_with_matlab():
    
    test_data = [
        {'Cc': torch.tensor([277.0525]), 'IPAR': torch.tensor([282.8877]), 'Csl': torch.tensor([384.4200]), 'ra': torch.tensor([10.4263]), 'rb': torch.tensor([10.1525]), 'Ts': torch.tensor([18.5000]), 'Pre': torch.tensor([843.3400]), 'Ds': torch.tensor([1360.6628]), 'Psi_L': torch.tensor([-0.0255]), 'Psi_sto_50': torch.tensor([-2.5000]), 'Psi_sto_00': torch.tensor([-0.5000]), 'CT': 3, 'Vmax': torch.tensor([42.]), 'Ha': torch.tensor([72.]), 'FI': torch.tensor([0.0810]), 'Oa': torch.tensor([210000.]), 'Do': torch.tensor([800.]), 'a1': torch.tensor([5.]), 'go': torch.tensor([0.0100]), 'gmes': torch.tensor([torch.inf]), 'rjv': torch.tensor([2.1000]), 'DS': torch.tensor([0.5000])},#, 'Rn': torch.tensor([739.9410]), 'QG': torch.tensor([84.1242]),},
        {'Cc': torch.tensor([329.6417]), 'IPAR': torch.tensor([29.9987]), 'Csl': torch.tensor([385.2900]), 'ra': torch.tensor([28.1468]), 'rb': torch.tensor([16.4111]), 'Ts': torch.tensor([4.2335]), 'Pre': torch.tensor([829.6200]), 'Ds': torch.tensor([96.8961]), 'Psi_L': torch.tensor([-0.0168]), 'Psi_sto_50': torch.tensor([-2.5000]), 'Psi_sto_00': torch.tensor([-0.5000]), 'CT': 3, 'Vmax': torch.tensor([42.]), 'Ha': torch.tensor([72.]), 'FI': torch.tensor([0.0810]), 'Oa': torch.tensor([210000.]), 'Do': torch.tensor([800.]), 'a1': torch.tensor([5.]), 'go': torch.tensor([0.0100]), 'gmes': torch.tensor([torch.inf]), 'rjv': torch.tensor([2.1000]), 'DS': torch.tensor([0.5000])},#, 'Rn': torch.tensor([147.1678]), 'QG': torch.tensor([74.7663]),},
        {'Cc': torch.tensor([318.8641]), 'IPAR': torch.tensor([1.0586]), 'Csl': torch.tensor([364.3700]), 'ra': torch.tensor([23.6161]), 'rb': torch.tensor([16.5921]), 'Ts': torch.tensor([8.8375]), 'Pre': torch.tensor([837.3350]), 'Ds': torch.tensor([279.8077]), 'Psi_L': torch.tensor([-0.0235]), 'Psi_sto_50': torch.tensor([-2.5000]), 'Psi_sto_00': torch.tensor([-0.5000]), 'CT': 3, 'Vmax': torch.tensor([42.]), 'Ha': torch.tensor([72.]), 'FI': torch.tensor([0.0810]), 'Oa': torch.tensor([210000.]), 'Do': torch.tensor([800.]), 'a1': torch.tensor([5.]), 'go': torch.tensor([0.0100]), 'gmes': torch.tensor([torch.inf]), 'rjv': torch.tensor([2.1000]), 'DS': torch.tensor([0.5000])},#, 'Rn': torch.tensor([-17.7087]), 'QG': torch.tensor([-18.7501]),},
        {'Cc': torch.tensor([351.1793]), 'IPAR': torch.tensor([42.5301]), 'Csl': torch.tensor([397.3400]), 'ra': torch.tensor([15.4748]), 'rb': torch.tensor([13.6120]), 'Ts': torch.tensor([-3.4835]), 'Pre': torch.tensor([828.4850]), 'Ds': torch.tensor([112.5857]), 'Psi_L': torch.tensor([-0.3036]), 'Psi_sto_50': torch.tensor([-2.5000]), 'Psi_sto_00': torch.tensor([-0.5000]), 'CT': 3, 'Vmax': torch.tensor([42.]), 'Ha': torch.tensor([72.]), 'FI': torch.tensor([0.0810]), 'Oa': torch.tensor([210000.]), 'Do': torch.tensor([800.]), 'a1': torch.tensor([5.]), 'go': torch.tensor([0.0100]), 'gmes': torch.tensor([torch.inf]), 'rjv': torch.tensor([2.1000]), 'DS': torch.tensor([0.5000])},#, 'Rn': torch.tensor([11.1133]), 'QG': torch.tensor([-2.5979]),},
        {'Cc': torch.tensor([285.0780]), 'IPAR': torch.tensor([221.2713]), 'Csl': torch.tensor([395.0300]), 'ra': torch.tensor([13.4421]), 'rb': torch.tensor([12.2295]), 'Ts': torch.tensor([19.1835]), 'Pre': torch.tensor([840.4300]), 'Ds': torch.tensor([1319.9600]), 'Psi_L': torch.tensor([-0.0402]), 'Psi_sto_50': torch.tensor([-2.5000]), 'Psi_sto_00': torch.tensor([-0.5000]), 'CT': 3, 'Vmax': torch.tensor([42.]), 'Ha': torch.tensor([72.]), 'FI': torch.tensor([0.0810]), 'Oa': torch.tensor([210000.]), 'Do': torch.tensor([800.]), 'a1': torch.tensor([5.]), 'go': torch.tensor([0.0100]), 'gmes': torch.tensor([torch.inf]), 'rjv': torch.tensor([2.1000]), 'DS': torch.tensor([0.5000])},#, 'Rn': torch.tensor([582.7670]), 'QG': torch.tensor([-2.5979])}
    ]

    for x in test_data:
        CcF,An,rs,Rdark,F755nm,GAM,gsCO2 = photosynthesis_biochemical(**x)
        print(f"obtained rs = {rs}")

    expected_output = {
        'rs': torch.tensor([466.5876]),
        'rs': torch.tensor([563.7118]),
        'rs': torch.tensor([992.6278]),
        'rs': torch.tensor([951.7378]),
        'rs': torch.tensor([453.3037]),
    }

    expected_LE = {
        torch.tensor([188.5540]),
        torch.tensor([59.3475]),
        torch.tensor([6.8946]),
        torch.tensor([39.4273]),
        torch.tensor([192.7630]),
    }

def compare_rs(predictors):
    for x in predictors:
        pb_predictors = {k: v for k, v in x.items() if k in pb_params}
        CcF,An,rs,Rdark,F755nm,GAM,gsCO2 = photosynthesis_biochemical(**pb_predictors)

        #print(f"using pb_predictors:")
        #print(pb_predictors)

        pred_rs = x['rs']

         # compare that the input rs and the obtained rs are the same
        print(f"initial rs: {pred_rs}, output_rs: {rs} (ratio={pred_rs / rs})")