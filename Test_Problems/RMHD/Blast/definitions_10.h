#define  PHYSICS                        RMHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            8

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  RADIATION                      NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT

/* -- user-defined parameters (labels) -- */

#define  RHO_IN                         0
#define  PRS_IN                         1
#define  RHO_OUT                        2
#define  PRS_OUT                        3
#define  BMAG                           4
#define  THETA                          5
#define  PHI                            6
#define  RADIUS                         7

/* [Beg] user-defined constants (do not change this line) */

#define  CT_EMF_AVERAGE                 UCT_GFORCE
#define  CT_EN_CORRECTION               YES
#define  LIMITER                        MINMOD_LIM
#define  ASSIGN_VECTOR_POTENTIAL        YES

/* [End] user-defined constants (do not change this line) */
