#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            0

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   EIGHT_WAVES
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */


/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   1.0e-6
#define  UNIT_LENGTH                    1.0e5
#define  UNIT_VELOCITY                  1.0e5
#define  INTERNAL_BOUNDARY              YES
#define  MULTIPLE_LOG_FILES             NO
#define  UNIT_PRESSURE                  sqrt(4*CONST_PI*UNIT_DENSITY*pow(UNIT_VELOCITY,2))

/* [End] user-defined constants (do not change this line) */
