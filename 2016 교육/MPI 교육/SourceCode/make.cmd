#!/bin/bash
# @ job_type = MPICH
# @ step_name = reuslt
# @ output = $(step_name).out
# @ error = $(step_name).err
# @ resources = ConsumableCpus(1) ConsumableMemory(1000mb)
# @ wall_clock_limit = 00:05:00
# @ class = normal
# @ node =1
# @ tasks_per_node = 4
# @ queue

export RSHCOMMAND=/opt/ibmll/LoadL/full/bin/llspawn
time mpirun -np $LOADL_TOTAL_TASKS -machinefile $LOADL_HOSTFILE test.x
