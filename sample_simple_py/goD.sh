#!/bin/bash
#PJM -L rscgrp=debug-a
#PJM -L node=1
#PJM --mpi proc=1
###PJM --omp thread=6
#PJM -L elapse=00:20:00
###PJM -g jh180023o
#PJM -g jh210022a
####PJM -j
####PJM -e err
#PJM -N waitio_test

module purge
module load intel impi
module load waitio

#rm *log*
#rm *out*
#rm fort*
#rm *err

##export LD_LIBRARY_PATH=/work/jh210022a/share/waitio/opt/waitio/lib/:$LD_LIBRARY_PATH
export COUPLER_JOB_ID=1
export WAITIO_HYBRID_THRESHOULD=524288
export WAITIO_SCRATCH_DIR="."

###export WAITIO_COMM_TYPE=hybrid
export WAITIO_COMM_TYPE=file
###export WAITIO_COMM_TYPE=socket
###export WAITIO_COMM_TYPE=verbs
###export WAITIO_MAXFILE_SIZE_PER_PROC=104857
###export WAITIO_MAXFILE_SIZE_PER_PROC=1048576
export WAITIO_MAXFILE_SIZE_PER_PROC=8048576
###export WAITIO_MAXFILE_SIZE_PER_PROC=67108864
###export WAITIO_MAXFILE_SIZE_PER_PROC=67108864000

echo "waitio-serv START"
export WAITIO_MASTER_HOST=`/work/share/waitio/bin/waitio-serv -c`
echo "waitio-serv END"
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=2
export WAITIO_NPB=3


module unload intel impi
module load gcc ompi

###source /work/jh210022a/c26011/h3opp/bin/activate

###mpiexec --oversubscribe -np 128 ./compC.aquarius
mpiexec -np 1 python3 ./compD.py
