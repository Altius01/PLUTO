#define PHYSICS MHD
#define DIMENSIONS 2
#define GEOMETRY CARTESIAN
#define BODY_FORCE NO
#define COOLING NO
#define RECONSTRUCTION LINEAR
#define TIME_STEPPING RK3
#define NTRACER 0
#define PARTICLES NO
#define USER_DEF_PARAMETERS 0

/* -- physics dependent declarations -- */

#define EOS IDEAL
#define ENTROPY_SWITCH NO
#define DIVB_CONTROL DIV_CLEANING
#define BACKGROUND_FIELD NO
#define AMBIPOLAR_DIFFUSION NO
#define RESISTIVITY NO
#define HALL_MHD NO
#define THERMAL_CONDUCTION NO
#define VISCOSITY EXPLICIT
#define LES EXPLICIT
#define RADIATION NO
#define ROTATING_FRAME NO

/* -- user-defined parameters (labels) -- */

/* [Beg] user-defined constants (do not change this line) */

#define CHECK_DIVB_CONDITION YES

/* [End] user-defined constants (do not change this line) */
