#!/bin/bash
#PJM -L rscgrp=debug-o
#PJM -L node=1x1x1:torus
#PJM --mpi proc=1
#PJM -L elapse=00:20:00
####PJM -g jh180023a
#PJM -g jh210022o
##3#PJM -j
####PJM -e err
#PJM -N waitio_test

module load fj
module load fjmpi
module load waitio

##export LD_LIBRARY_PATH=/work/01/jh210022a/share/waitio/opt/waitio/lib/:$LD_LIBRARY_PATH
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
export WAITIO_MASTER_HOST=`/work/share/waitio/bin/waitio-serv-a64fx -c`
#export WAITIO_MASTER_HOST=wo0001
echo "waitio-serv END"
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=1
export WAITIO_NPB=3

mpiexec -n 1 ./compB.odyssey

