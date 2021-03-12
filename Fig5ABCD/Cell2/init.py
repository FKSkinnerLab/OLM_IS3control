from mpi4py import MPI
from neuron import h

h.load_file("SynParamSearch.hoc")
execfile("RunSims.py")
h.quit()
sys.exit(0)