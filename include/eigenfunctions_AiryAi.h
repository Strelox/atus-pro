#ifndef __eigenfunctions_AiryAi_h__
#define __eigenfunctions_AiryAi_h__

#include "gsl/gsl_sf.h"
#include "math.h"

typedef double (*TEF)(const double, const double);

double AiryEF_1( const double a, const double x )
{
return 1.426104629e+00*sqrt(a)*gsl_sf_airy_Ai(a*x-2.338107e+00, 0);
}

double AiryEF_2( const double a, const double x )
{
return 1.245157319e+00*sqrt(a)*gsl_sf_airy_Ai(a*x-4.087949e+00, 0);
}

double AiryEF_3( const double a, const double x )
{
return 1.155796749e+00*sqrt(a)*gsl_sf_airy_Ai(a*x-5.520560e+00, 0);
}

double AiryEF_4( const double a, const double x )
{
return 1.097874722e+00*sqrt(a)*gsl_sf_airy_Ai(a*x-6.786708e+00, 0);
}

double AiryEF_5( const double a, const double x )
{
return 1.055592004e+00*sqrt(a)*gsl_sf_airy_Ai(a*x-7.944134e+00, 0);
}

double AiryEF_6( const double a, const double x )
{
return 1.022575597e+00*sqrt(a)*gsl_sf_airy_Ai(a*x-9.022651e+00, 0);
}

double AiryEF_7( const double a, const double x )
{
return 9.956488919e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.004017e+01, 0);
}

double AiryEF_8( const double a, const double x )
{
return 9.730099788e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.100852e+01, 0);
}

double AiryEF_9( const double a, const double x )
{
return 9.535427771e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.193602e+01, 0);
}

double AiryEF_10( const double a, const double x )
{
return 9.365103494e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.282878e+01, 0);
}

double AiryEF_11( const double a, const double x )
{
return 9.214018166e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.369149e+01, 0);
}

double AiryEF_12( const double a, const double x )
{
return 9.078491613e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.452783e+01, 0);
}

double AiryEF_13( const double a, const double x )
{
return 8.955789207e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.534076e+01, 0);
}

double AiryEF_14( const double a, const double x )
{
return 8.843826175e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.613269e+01, 0);
}

double AiryEF_15( const double a, const double x )
{
return 8.740978533e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.690563e+01, 0);
}

double AiryEF_16( const double a, const double x )
{
return 8.645957847e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.766130e+01, 0);
}

double AiryEF_17( const double a, const double x )
{
return 8.557725562e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.840113e+01, 0);
}

double AiryEF_18( const double a, const double x )
{
return 8.475432936e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.912638e+01, 0);
}

double AiryEF_19( const double a, const double x )
{
return 8.398377853e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-1.983813e+01, 0);
}

double AiryEF_20( const double a, const double x )
{
return 8.325973228e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.053733e+01, 0);
}

double AiryEF_21( const double a, const double x )
{
return 8.257723450e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.122483e+01, 0);
}

double AiryEF_22( const double a, const double x )
{
return 8.193206466e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.190137e+01, 0);
}

double AiryEF_23( const double a, const double x )
{
return 8.132060095e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.256761e+01, 0);
}

double AiryEF_24( const double a, const double x )
{
return 8.073971303e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.322416e+01, 0);
}

double AiryEF_25( const double a, const double x )
{
return 8.018667811e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.387156e+01, 0);
}

double AiryEF_26( const double a, const double x )
{
return 7.965911380e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.451030e+01, 0);
}

double AiryEF_27( const double a, const double x )
{
return 7.915492411e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.514082e+01, 0);
}

double AiryEF_28( const double a, const double x )
{
return 7.867225614e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.576353e+01, 0);
}

double AiryEF_29( const double a, const double x )
{
return 7.820946384e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.637881e+01, 0);
}

double AiryEF_30( const double a, const double x )
{
return 7.776507909e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.698699e+01, 0);
}

double AiryEF_31( const double a, const double x )
{
return 7.733778694e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.758839e+01, 0);
}

double AiryEF_32( const double a, const double x )
{
return 7.692640529e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.818331e+01, 0);
}

double AiryEF_33( const double a, const double x )
{
return 7.652986805e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.877201e+01, 0);
}

double AiryEF_34( const double a, const double x )
{
return 7.614721019e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.935475e+01, 0);
}

double AiryEF_35( const double a, const double x )
{
return 7.577755572e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-2.993176e+01, 0);
}

double AiryEF_36( const double a, const double x )
{
return 7.542010734e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.050327e+01, 0);
}

double AiryEF_37( const double a, const double x )
{
return 7.507413704e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.106947e+01, 0);
}

double AiryEF_38( const double a, const double x )
{
return 7.473897871e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.163056e+01, 0);
}

double AiryEF_39( const double a, const double x )
{
return 7.441402114e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.218671e+01, 0);
}

double AiryEF_40( const double a, const double x )
{
return 7.409870245e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.273810e+01, 0);
}

double AiryEF_41( const double a, const double x )
{
return 7.379250478e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.328488e+01, 0);
}

double AiryEF_42( const double a, const double x )
{
return 7.349495004e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.382721e+01, 0);
}

double AiryEF_43( const double a, const double x )
{
return 7.320559585e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.436523e+01, 0);
}

double AiryEF_44( const double a, const double x )
{
return 7.292403216e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.489907e+01, 0);
}

double AiryEF_45( const double a, const double x )
{
return 7.264987820e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.542886e+01, 0);
}

double AiryEF_46( const double a, const double x )
{
return 7.238277986e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.595471e+01, 0);
}

double AiryEF_47( const double a, const double x )
{
return 7.212240712e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.647675e+01, 0);
}

double AiryEF_48( const double a, const double x )
{
return 7.186845220e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.699507e+01, 0);
}

double AiryEF_49( const double a, const double x )
{
return 7.162062731e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.750980e+01, 0);
}

double AiryEF_50( const double a, const double x )
{
return 7.137866326e-01*sqrt(a)*gsl_sf_airy_Ai(a*x-3.802101e+01, 0);
}

TEF AIRYEF[] = {&AiryEF_1, &AiryEF_2, &AiryEF_3, &AiryEF_4, &AiryEF_5, &AiryEF_6, &AiryEF_7, &AiryEF_8, &AiryEF_9, &AiryEF_10, &AiryEF_11, &AiryEF_12, &AiryEF_13, &AiryEF_14, &AiryEF_15, &AiryEF_16, &AiryEF_17, &AiryEF_18, &AiryEF_19, &AiryEF_20, &AiryEF_21, &AiryEF_22, &AiryEF_23, &AiryEF_24, &AiryEF_25, &AiryEF_26, &AiryEF_27, &AiryEF_28, &AiryEF_29, &AiryEF_30, &AiryEF_31, &AiryEF_32, &AiryEF_33, &AiryEF_34, &AiryEF_35, &AiryEF_36, &AiryEF_37, &AiryEF_38, &AiryEF_39, &AiryEF_40, &AiryEF_41, &AiryEF_42, &AiryEF_43, &AiryEF_44, &AiryEF_45, &AiryEF_46, &AiryEF_47, &AiryEF_48, &AiryEF_49, &AiryEF_50};
#endif
