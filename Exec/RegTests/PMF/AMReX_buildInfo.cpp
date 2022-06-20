
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2022-06-13 15:33:55.559986";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/work2/08181/tsouders/frontera/PeleBuild/2022_0216_Turb/PeleC/Exec/RegTests/PMF";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux login3.frontera.tacc.utexas.edu 3.10.0-1160.45.1.el7.x86_64 #1 SMP Wed Oct 13 17:20:51 UTC 2021 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "/work2/08181/tsouders/frontera/PeleBuild/2022_0216_Turb/PeleC/Submodules/AMReX";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "intel";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "8.3.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = "";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = "";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = "";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  int num_modules = 1;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";
  static const char AUX1[] = "CHEMISTRY";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return AUX1;

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";
  static const char AUX1[] = "LiDryer";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return AUX1;

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "v0.2.1-144-g248d366a-dirty";
  static const char HASH2[] = "22.02-7-g9558c8a16";
  static const char HASH3[] = "v0.1-961-gae9ca56d-dirty";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;
    case 3: return HASH3;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

#ifdef AMREX_USE_CUDA
const char* buildInfoGetCUDAVersion() {

  static const char CUDA_VERSION[] = "";
  return CUDA_VERSION;
}
#endif

}
