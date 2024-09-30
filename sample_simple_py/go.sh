#mpiexec -np 2 compA : -np 1 compB : -np 2 compC : -np 3 python3 compD.py
source /work/jh210022a/c26011/pyenv/odyssey/bin/activate
module load waitio
module load python
mpiexec -np 1 compA : -np 1 compB : -np 1 compC : -np 1 python3 compD.py 
#mpiexec -np 1 compA : -np 1 compB : -np 1 compC 
#mpiexec -np 1 compA : -np 1 python3 compD.py 
