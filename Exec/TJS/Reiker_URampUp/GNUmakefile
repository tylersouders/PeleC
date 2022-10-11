# AMReX
DIM = 2
COMP = intel
PRECISION = DOUBLE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_DPCPP = FALSE
# USE_SUNDIALS_PP = TRUE
# USE_CUDA_CVODE = TRUE

# NO_MPI_CHECKING = TRUE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleC
USE_REACT = TRUE
USE_MASA = FALSE
USE_EB = TRUE
Eos_Model := Fuego
Chemistry_Model := JL4
#Chemistry_Model := Ghani_C3H8
#Chemistry_Model := propane_cerfacs_2step
Transport_Model := Simple

# GNU Make
Bpack := ./Make.package
Blocs := .
PELEC_HOME := ../../..
include $(PELEC_HOME)/Exec/Make.PeleC
