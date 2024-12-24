# TODO: kind of optional this one
import python.train as train


# TODO: Use something like that in unit testing
test_data = [
    (({"ctxa":{"k1":1.0,"k2":2.0}, "ctxb":{"k1":3.0,"k2":4.0}}, {"g1":11.0, "g2":12.0}), 5.0),
    (({"ctxa":{"k1":21.0,"k2":22.0}, "ctxb":{"k1":23.0,"k2":24.0}}, {"g1":13.0, "g2":14.0}), 6.0),
    (({"ctxa":{"k1":31.0,"k2":32.0}, "ctxb":{"k1":33.0,"k2":34.0}}, {"g1":15.0, "g2":16.0}), 7.0),
]
batch = train.ctx_dict_batch(test_data, 2)
print(f"batch: {batch}")
