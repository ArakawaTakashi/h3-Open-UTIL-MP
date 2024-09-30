#!/bin/bash
#PJM -L rscgrp=debug-o
#PJM -L node=1x1x1:torus
##PJM -L node=1
#PJM --mpi proc=1
###PJM --omp thread=6
#PJM -L elapse=00:20:00
###PJM -g jh180023o
#PJM -g jh210022o
####PJM -j
####PJM -e err
#PJM -N waitio_test

module load fj
module load fjmpi
module load waitio

#rm *log*
#rm *out*
#rm fort*
#rm *err

##export LD_LIBRARY_PATH=/work/jh210022a/share/waitio/opt/waitio/lib/:$LD_LIBRARY_PATH
export COUPLER_JOB_ID=1
export WAITIO_HYBRID_THRESHOULD=524288
export WAITIO_SCRATCH_DIR="."

rm -rf ${WAITIO_SCRATCH_DIR}/.waitio-file-${COUPLER_JOB_ID}/
###export WAITIO_COMM_TYPE=hybrid
export WAITIO_COMM_TYPE=file
###export WAITIO_COMM_TYPE=socket
###export WAITIO_COMM_TYPE=verbs
###export WAITIO_MAXFILE_SIZE_PER_PROC=104857
###export WAITIO_MAXFILE_SIZE_PER_PROC=1048576
export WAITIO_MAXFILE_SIZE_PER_PROC=8048576
###export WAITIO_MAXFILE_SIZE_PER_PROC=67108864
###export WAITIO_MAXFILE_SIZE_PER_PROC=67108864000

export WAITIO_MASTER_HOST=`hostname`
export WAITIO_MASTER_PORT=7100
export WAITIO_PBID=0
export WAITIO_NPB=3

echo "waitio-serv START"
hostname
/work/share/waitio/bin/waitio-serv-a64fx -d -m $WAITIO_MASTER_HOST
echo "waitio-serv END"

mpiexec -np 1 ./compA.odyssey
