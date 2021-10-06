#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
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
#define  RING_AVERAGE                   8
#define  MULTIPLE_LOG_FILES             YES
#define  LIMITER                        MC_LIM

/* [End] user-defined constants (do not change this line) */
