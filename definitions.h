#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 MP5
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   DIV_CLEANING
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    EXPLICIT
#define  HALL_MHD                       EXPLICIT
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  RADIATION                      NO
#define  ROTATING_FRAME                 NO
#define  LES                            EXPLICIT
#define  MHD_LES                        EXPLICIT


#define INIT_CS                         0.1
#define INIT_YS                         0.1
#define INIT_DS                         0.1
#define DYNAMIC_PROCEDURE               NO



/* -- user-defined parameters (labels) -- */


/* [Beg] user-defined constants (do not change this line) */

#define  ROTATE                         -1
#define  SHOCK_FLATTENING               NO
#define  WARNING_MESSAGES               YES
#define  CHECK_DIVB_CONDITION           YES

/* [End] user-defined constants (do not change this line) */
