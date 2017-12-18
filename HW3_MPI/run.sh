#!/bin/bash


make

jobid=$(qsub jobfile_part1_node1| cut -d '.' -f 1)
#jobid=$(qsub -W depend=afterok:${jobid} jobfile_part1_node1| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part1_node2| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part1_node4| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part1_node8| cut -d '.' -f 1)

jobid=$(qsub -W depend=afterok:${jobid} jobfile_part2_node1| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part2_node2| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part2_node4| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part2_node8| cut -d '.' -f 1)

jobid=$(qsub -W depend=afterok:${jobid} jobfile_part3_node1| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part3_node2| cut -d '.' -f 1)
jobid=$(qsub -W depend=afterok:${jobid} jobfile_part3_node4| cut -d '.' -f 1)
#jobid=$(qsub -W depend=afterok:${jobid} jobfile_part3_node8| cut -d '.' -f 1)

qsub -W depend=afterok:${jobid} jobfile_part3_node8


