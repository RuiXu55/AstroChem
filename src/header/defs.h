#ifndef DEFINITIONS_H
#define DEFINITIONS_H 

/*----------------------------------------------------------------------------*/
/*  macros which define physics and algorithm
 *  (user modified via configure) */
#define CHEMISTRY

#define MAXLEN 256

typedef double Real;

/*----------------------------------------------------------------------------*/
/* general purpose macros (never modified) */
#ifndef MIN
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif
#define SIGN(a) ( ((a) < 0.) ? -1. : 1. )
#define SQR(x) ( (x)*(x) )
#define STR(x) #x
#define SQRT2 1.4142135623730951
#define ONE_OVER_SQRT2 0.7071067811865475
#define PI       3.14159265358979323846
#define ONE_3RD  0.3333333333333333
#define TWO_3RDS 0.6666666666666667
#define FOUR_3RDS 1.333333333333333
#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+20
#define AU 1.49598e13

#endif /* DEFINITIONS_H */
