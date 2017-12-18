
- The makefile cover all the need object file

- callMypack.c: Problem 1, comparison of LAPCKE and mydgetrf
  obj file: ./callMypack [matrix_size] (default is 1600)
  job file: jobfile_callMypack

- callpack.c: Problem 2, comparison of mydgetrf and blocked version
  obj file: ./callpack [matrix_size] (default is 1600)
  job file: jobfile_calllapack

- make clean to delete all the object file

- ./rmJobFile.sh will delete all the result file