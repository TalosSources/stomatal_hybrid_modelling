"""
This file contains implementations of the differentiable pipelines
relating rs and Q_LE.
To include: 
* the normal pipeline, with differentiable PM Equation
* variations of the normal pipeline, including the function mapping the output of the FCN to a positive number [DROPPED IN SOME WAY]
* the inv pipeline, a FCN computes rs, trained to match the inverted PM rs_hat reconstruction from Q_LE
Actually, we won't separate pipelines as distinct entity, but rather let the config specify individual components freely.
"""
import inspect
import torch

from photosynthesis_biochemical import photosynthesis_biochemical as pb
import differentiable_relations



def make_pipeline(rs_model, gsco2_model, Vmax_model, output_rs=False):

    # this pipeline calls the pb module with the 2 models inserted, using the constants and predictors as inputs, 
    # and computes Q_LE from the outputs 
    pb_params = set(inspect.signature(pb).parameters)
    qle_params = set(inspect.signature(differentiable_relations.Q_LE).parameters)

    # Actual call to the pb and PM modules, called for each ctx
    def subpipeline(predictors):
        pb_predictors = {k: v for k, v in predictors.items() if k in pb_params}
        rs = pb(**pb_predictors, rs_model=rs_model, gsCO2_model=gsco2_model, Vmax_model=Vmax_model)

        qle_predictors = {k: v for k, v in predictors.items() if k in qle_params}
        qle_predictors.__delitem__("rs") # TODO: temporary debug fix. will need to remove rs from predictors one day

        # COMPUTE Q_LE USING THE ESTIMATED RS
        Q_LE = differentiable_relations.Q_LE(rs=rs, **qle_predictors)

        return (Q_LE, rs) if output_rs else Q_LE
    
    # NOTE: this 2 methods as well as the global calculation could be moved to differentiable_relations.py
    def W_on_m2_to_mm_on_H(Q):
        # For now simply divide by 684 NOTE: Not 100% sure it's the way to go
        return Q / 684
    
    def mm_on_H_to_W_on_m2(Q):
        return Q * 684
    
    """
    This is the differentiable pipeline, taking as input an object containing
    values for the 4 different contexts, computing the Q_LE term for each of 
    them, and aggregating them to obtain a Q_LE comparable with observations
    Assumes it receives predictors p, on this format:
    * single value tensors for general EG, ELitter, etc., and also EIn_H and EIn_L 
    * For key ctx, a map of predictors for pb and Q_LE, with values for this ctx.
    Args: 
        x: A pair (single, global), where single is a map ctx -> predictor_dict, 
            and global is a dict of context-independant values. 
    """
    def pipeline(x):

        # ctx_to_predictors maps context to predictor dicts, and g is the global variables dict 
        ctx_to_predictors, g = x

        Q_LE_wm2 = 0

        # For each context, find the associated Q_LE value using the subpipeline. 
        # Sum all the relevant results
        rs_values = []
        for _, preds in ctx_to_predictors.items():
            # TODO: All ctxs as a batch
            ctx_res = subpipeline(preds)
            if output_rs:
                ctx_res, rs = ctx_res
                rs_values.append(rs)

            # put nan values to 0
            #ctx_res[torch.isnan(ctx_res)] = 0.
            #ctx_res = torch.where(torch.isnan(ctx_res), torch.tensor(0.0), ctx_res).detach()
            #ctx_res.requires_grad_()  # Re-enable gradient tracking for the modified tensor
            Q_LE_wm2 = Q_LE_wm2 + ctx_res 

        # Convert the result to mm/H, to perform the computation below
        Q_LE_mmH = W_on_m2_to_mm_on_H(Q_LE_wm2)

        # Use the global variables to find the global Q_LE, which can be compared with observations
        HLTotal = Q_LE_mmH + g["EIn_H"] + g["EIn_L"]
        ET = HLTotal + g["EG"] + g["ELitter"] + g["ESN"] + g["ESN_In"] + g["EWAT"] + g["EICE"] + g["EIn_urb"] + g["EIn_rock"]

        final_Q_LE = mm_on_H_to_W_on_m2(ET)

        return (final_Q_LE, rs_values) if output_rs else final_Q_LE

    return pipeline