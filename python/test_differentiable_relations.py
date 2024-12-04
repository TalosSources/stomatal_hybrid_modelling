import inspect

import differentiable_relations

"""
Takes as input an array of predictors, each line containing named parameters:
rs, ra, Rn, QG, Ds, Ts
"""
qle_params = set(inspect.signature(differentiable_relations.Q_LE).parameters)
compute_rs_params = set(inspect.signature(differentiable_relations.compute_rs).parameters)
def test_Q_LE_invertible(predictors):
    print(f"===================0Test Q_LE_invertible===================0")
    for x in predictors:
        #print(f"using predictor = {x}")
        # extract parameters used in Q_LE
        qle_predictors = {k: v for k, v in x.items() if k in qle_params}

        # compute Q_LE with those parameters
        LE = differentiable_relations.Q_LE(**qle_predictors)

        # extract parameters used in compute_rs
        compute_rs_predictors = {k: v for k, v in x.items() if k in compute_rs_params}

        # compute rs with Q_LE and those parameters
        output_rs = differentiable_relations.compute_rs(Q_LE=LE, **compute_rs_predictors)

        # compare that the initial rs and the obtained rs are the same
        print(f"initial rs: {qle_predictors['rs']}, output_rs: {output_rs}")

"""
Takes as input an array of predictors, each line containing a tuple with named parameters:
rs, ra, Rn, QG, Ds, Ts
and:
Q_LE
"""
def test_Q_LE_using_inverted_PM(predictors):
    print(f"===================0Test Q_LE_using_inverted_PM===================0")
    for x, y in predictors:
        #print(f"using predictor = {x}, output = {y}")

        input_LE = y

        # extract parameters used in compute_rs
        compute_rs_predictors = {k: v for k, v in x.items() if k in compute_rs_params}

        # compute rs with Q_LE and those parameters
        output_rs = differentiable_relations.compute_rs(Q_LE=input_LE, **compute_rs_predictors)

        # compare the obtained rs with the initial rs
        print(f"initial rs: {x['rs']}, output_rs: {output_rs}")

    