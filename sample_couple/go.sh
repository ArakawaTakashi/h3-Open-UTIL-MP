#mpiexec -np 2 compA : -np 1 compB : -np 2 compC : -np 3 python3 compD.py
#mpiexec -np 1 compA : -np 1 compB : -np 1 compC : -np 1 python3 compD.py 
mpiexec -np 1 compA : -np 1 python3 compD.py 
