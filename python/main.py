from photosynthesis_biochemical import photosynthesis_biochemical

import numpy as np
import torch

CT = 3 # control flag

Cc = torch.tensor([1])
IPAR = torch.tensor([1])
Csl = torch.tensor([1])
ra = torch.tensor([1])
rb = torch.tensor([1])
Ts = torch.tensor([1])
Pre = torch.tensor([1])
Ds = torch.tensor([1])
Psi_L = torch.tensor([1])
Psi_sto_50 = torch.tensor([2])
Psi_sto_00 = torch.tensor([1])
Vmax = torch.tensor([1])
DS = torch.tensor([1])
Ha = torch.tensor([1])
FI = torch.tensor([1])
Oa = torch.tensor([1])
Do = torch.tensor([1])
a1 = torch.tensor([1])
go = torch.tensor([1])
gmes = torch.tensor([1])
rjv = torch.tensor([1])

result = photosynthesis_biochemical(Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds, Psi_L,Psi_sto_50,Psi_sto_00, CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv)

print(f"result = {result}")