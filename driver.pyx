LatticeRegistry = {}
ModelRegistry = {}

# try to make this more general. perhaps move later?
from lattice import linear
from models import FA

LatticeRegistry[linear.lattice_name] = linear
ModelRegistry[FA.model_name] = FA
