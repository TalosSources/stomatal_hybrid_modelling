

def printGradInfo(x, name):
    print(f"{name} value: {x}")
    print(f"{name} type: {type(x)}")
    print(f"{name} requireGrad: {x.requires_grad}")
    print(f"{name} grad: {x.grad}")