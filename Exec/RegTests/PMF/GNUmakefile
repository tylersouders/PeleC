# AMReX
DIM = 3
COMP = gnu
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
<<<<<<< HEAD
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = TRUE
=======
USE_MPI = FALSE
USE_OMP = FALSE
USE_CUDA = FALSE
>>>>>>> development
USE_HIP = FALSE
USE_DPCPP = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
USE_REACT = TRUE
USE_SUNDIALS_PP = FALSE
USE_EB = FALSE
USE_MASA = FALSE
Eos_dir := Fuego
<<<<<<< HEAD
Chemistry_Model := drm19
Transport_dir := Simple
USE_SINGE = TRUE
=======
Chemistry_Model := LiDryer
Transport_dir := Simple
>>>>>>> development

# GNU Make
Bpack := ./Make.package
Blocs := .
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
DEFINES += -DAMREX_GPUS_PER_NODE=2
