#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK2
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            8

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  RSTAR                          0
#define  RHOSTAR                        1
#define  GAMMA                          2
#define  CS_VESC                        3
#define  VA_VESC                        4
#define  VROT_VESC                      5
#define  RSS                            6
#define  POLYIDX                        7

/* [Beg] user-defined constants (do not change this line) */

#define  LIMITER                        MINMOD_LIM

/* [End] user-defined constants (do not change this line) */
