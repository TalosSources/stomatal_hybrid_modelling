"""
This file contains implementations of the differentiable pipelines
relating rs and Q_LE.
To include: 
* the normal pipeline, with differentiable PM Equation
* variations of the normal pipeline, including the function mapping the output of the FCN
* the inv pipeline, a FCN computes rs, trained to match the inverted PM rs_hat reconstruction from Q_LE
"""
import inspect
import torch

from photosynthesis_biochemical import photosynthesis_biochemical as pb
import differentiable_relations



def make_pipeline(gsCO2_model, Vmax_model, output_rs=False):

    print(f"outputrs = {output_rs}")
    # define constants. Doing our function structure that way allows to store the constants cleanly
    # constants = ... # NOTE: Might not be necessary

    # the pipeline calls the pb module with the 2 models inserted, using the constants and predictors as inputs, and converts the outputs to Q_LE 
    pb_params = set(inspect.signature(pb).parameters)
    qle_params = set(inspect.signature(differentiable_relations.Q_LE).parameters)

    def subpipeline(predictors):
        pb_predictors = {k: v for k, v in predictors.items() if k in pb_params}
        _,_,rs,_,_,_,_ = pb(**pb_predictors, gsCO2_model=gsCO2_model, Vmax_model=Vmax_model)
        #print(f"input_rs={input_rs} while output_rs = {rs}")
        qle_predictors = {k: v for k, v in predictors.items() if k in qle_params}
        qle_predictors.__delitem__("rs") # TODO: temporary debug fix

        #pred_rs = predictors["rs"]
        #print(f"given_rs: {pred_rs}, computed_rs: {rs}")
        #if not torch.isnan(pred_rs) and not torch.isinf(pred_rs) and not torch.isnan(rs) and not torch.isinf(rs):
        #    if torch.abs(rs - pred_rs) > 0.5 and rs*pred_rs > 0.001:
        #        print(f"DIFFERENT!")
        #        print(f"and pb_predictors are: {pb_predictors}")

        # COMPUTE Q_LE USING THE PREDICTOR RS
        #Q_LE = differentiable_relations.Q_LE(rs=pred_rs, **qle_predictors)

        # COMPUTE Q_LE USING THE ESTIMATED RS
        Q_LE = differentiable_relations.Q_LE(rs=rs, **qle_predictors)

        return (Q_LE, rs) if output_rs else Q_LE
        #print(f"\n\nsubpipeline got {predictors}\n\n and output {Q_LE} (shape={Q_LE.size()})\n\n")
        #return Q_LE
    
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
        single_pipeline_values, g = x
        # predictors contains a map between context and predictor maps for this context.
        Q_LE_wm2 = torch.zeros(32) # TODO: Need to replace by the actual batch_size
        rs_values = []
        #print(f"the acc is {Q_LE_wm2}, size={Q_LE_wm2.size()}")
        for ctx, preds in single_pipeline_values.items():
            #print(f"[ctx={ctx}] ", end='')
            ctx_res = subpipeline(preds)
            if output_rs:
                ctx_res, rs = ctx_res
                rs_values.append(rs)
            # put nan values to 0
            ctx_res[torch.isnan(ctx_res)] = 0.
            Q_LE_wm2 += ctx_res
        Q_LE_mmH = W_on_m2_to_mm_on_H(Q_LE_wm2)
        #print(f"predicted Q_LE_mmH: {Q_LE_mmH}")

        HLTotal = Q_LE_mmH + g["EIn_H"] + g["EIn_L"]
        ET = HLTotal + g["EG"] + g["ELitter"] + g["ESN"] + g["ESN_In"] + g["EWAT"] + g["EICE"] + g["EIn_urb"] + g["EIn_rock"]
        #print(f"and then q_le_mmh={Q_LE_mmH}, ET=q_le_mmh+{g['EIn_H'] + g['EIn_L']+g['EG'] + g['ELitter'] + g['ESN'] + g['ESN_In'] + g['EWAT'] + g['EICE'] + g['EIn_urb'] + g['EIn_rock']}")

        final_Q_LE = mm_on_H_to_W_on_m2(ET)
        #print(f"yielding final_Q_LE = {final_Q_LE}")

        return (final_Q_LE, rs_values) if output_rs else final_Q_LE

    return pipeline