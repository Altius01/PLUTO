#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            4

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  INCLUDE_LES                    NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  CS_WIND                        0
#define  RHO_AMB                        1
#define  CS_AMB                         2
#define  V_CSM                          3

/* [Beg] user-defined constants (do not change this line) */

#define  INCLUDE_JDIR                   NO
#define  INTERNAL_BOUNDARY              YES
#define  INTERNAL_BOUNDARY_CFL          YES
#define  SHOCK_FLATTENING               MULTID
#define  CHAR_LIMITING                  YES
#define  LIMITER                        VANLEER_LIM

/* [End] user-defined constants (do not change this line) */
