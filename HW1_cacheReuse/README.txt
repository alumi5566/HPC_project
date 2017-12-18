


- The makefile cover all the need object file

- dgemm_2_1_0.c: Problem1, Part#2: comparison of dgemm2, dgemm1 and dgemm0
  obj file: dgemm_2_1_0
  job file: jobfile_dgemm_2_1_0

- dgemm_3_2_1_0.c: Problem1, Part#3: comparison of dgemm3, dgemm2, dgemm1 and dgemm0
  obj file: dgemm_3_2_1_0
  job file: jobfile_dgemm_3_2_1_0

- block_1_2_ijk_jik.c/ block_1_2_jki_kji.c/ block_1_2_kij_ikj.c:
  Problem2, Part#3: comparison of non-blocking and blocking algorithm
  obj file: block_1_2_ijk_jik/ block_1_2_jki_kji/ block_1_2_kij_ikj
  job file: jobfile_block_1_2_A/ jobfile_block_1_2_B/ jobfile_block_1_2_C

- block_both.c : Problem2, Part#4: use both register and cache blocking
  obj file: block_both/ block_both_O0/ block_both_O1/ block_both_O2/ block_both_O3
	depends on the argument when you compile
  job file:jobfile_block_both/ jobfile_block_both_O0/ jobfile_block_both_O1/ jobfile_block_both_O2/ jobfile_block_both_O3

- make clean to delete all the object file

- ./rmJobFile.sh will delete all the result file

- The performance result in my report is in /Result_1014