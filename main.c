#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>
#include <assert.h>
#include <stdlib.h>

//control if floating-point arithmetic is run
#define floating_arithmetic 0

//general debugging
#define debug_matrices 0
#define debug_iterations 0

//individual debugging
#define Mkk 0
#define num_denom 0
#define quot 0
#define theta_sumdiff 0
#define theta_LR 0
#define wrap2pi_debug 0
#define lookupindex_debug 0
#define cos_sin_LR 0
#define lookup_cos_sin 0

#define FIXED_FRACTIONAL_PART 16

#define invfact2 0.5f
#define invfact3 0.16666666666666666f
#define invfact4 0.04166666666666666f
#define invfact5 0.00833333333333333f
#define invfact6 0.00138888888888888f
#define invfact7 0.00019841269841269f
#define invfact8 2.48015873015873e-05f

#define invfact2fix 32768
#define invfact3fix 10922
#define invfact4fix 2730
#define invfact5fix 546
//#define inv3 0.3333333333333333f


float lookup_cos_float[128] = {

    1.000000f,     0.998795f,     0.995185f,     0.989177f,     0.980785f,
    0.970031f,     0.956940f,     0.941544f,     0.923880f,     0.903989f,
    0.881921f,     0.857729f,     0.831470f,     0.803208f,     0.773010f,
    0.740951f,     0.707107f,     0.671559f,     0.634393f,     0.595699f,
    0.555570f,     0.514103f,     0.471397f,     0.427555f,     0.382683f,
    0.336890f,     0.290285f,     0.242980f,     0.195090f,     0.146730f,
    0.098017f,     0.049068f,    -0.000000f,    -0.049068f,    -0.098017f,
   -0.146730f,    -0.195090f,    -0.242980f,    -0.290285f,    -0.336890f,
   -0.382683f,    -0.427555f,    -0.471397f,    -0.514103f,    -0.555570f,
   -0.595699f,    -0.634393f,    -0.671559f,    -0.707107f,    -0.740951f,
   -0.773010f,    -0.803208f,    -0.831470f,    -0.857729f,    -0.881921f,
   -0.903989f,    -0.923880f,    -0.941544f,    -0.956940f,    -0.970031f,
   -0.980785f,    -0.989177f,    -0.995185f,    -0.998795f,    -1.000000f,
   -0.998795f,    -0.995185f,    -0.989177f,    -0.980785f,    -0.970031f,
   -0.956940f,    -0.941544f,    -0.923880f,    -0.903989f,    -0.881921f,
   -0.857729f,    -0.831470f,    -0.803208f,    -0.773010f,    -0.740951f,
   -0.707107f,    -0.671559f,    -0.634393f,    -0.595699f,    -0.555570f,
   -0.514103f,    -0.471397f,    -0.427555f,    -0.382683f,    -0.336890f,
   -0.290285f,    -0.242980f,    -0.195090f,    -0.146730f,    -0.098017f,
   -0.049068f,    -0.000000f,     0.049068f,     0.098017f,     0.146730f,
    0.195090f,     0.242980f,     0.290285f,     0.336890f,     0.382683f,
    0.427555f,     0.471397f,     0.514103f,     0.555570f,     0.595699f,
    0.634393f,     0.671559f,     0.707107f,     0.740951f,     0.773010f,
    0.803208f,     0.831470f,     0.857729f,     0.881921f,     0.903989f,
    0.923880f,     0.941544f,     0.956940f,     0.970031f,     0.980785f,
    0.989177f,     0.995185f,     0.998795f
};

float lookup_sin_float[128] = {

    0.000000f,     0.049068f,     0.098017f,     0.146730f,     0.195090f,
    0.242980f,     0.290285f,     0.336890f,     0.382683f,     0.427555f,
    0.471397f,     0.514103f,     0.555570f,     0.595699f,     0.634393f,
    0.671559f,     0.707107f,     0.740951f,     0.773010f,     0.803208f,
    0.831470f,     0.857729f,     0.881921f,     0.903989f,     0.923880f,
    0.941544f,     0.956940f,     0.970031f,     0.980785f,     0.989177f,
    0.995185f,     0.998795f,     1.000000f,     0.998795f,     0.995185f,
    0.989177f,     0.980785f,     0.970031f,     0.956940f,     0.941544f,
    0.923880f,     0.903989f,     0.881921f,     0.857729f,     0.831470f,
    0.803208f,     0.773010f,     0.740951f,     0.707107f,     0.671559f,
    0.634393f,     0.595699f,     0.555570f,     0.514103f,     0.471397f,
    0.427555f,     0.382683f,     0.336890f,     0.290285f,     0.242980f,
    0.195090f,     0.146730f,     0.098017f,     0.049068f,     0.000000f,
   -0.049068f,    -0.098017f,    -0.146730f,    -0.195090f,    -0.242980f,
   -0.290285f,    -0.336890f,    -0.382683f,    -0.427555f,    -0.471397f,
   -0.514103f,    -0.555570f,    -0.595699f,    -0.634393f,    -0.671559f,
   -0.707107f,    -0.740951f,    -0.773010f,    -0.803208f,    -0.831470f,
   -0.857729f,    -0.881921f,    -0.903989f,    -0.923880f,    -0.941544f,
   -0.956940f,    -0.970031f,    -0.980785f,    -0.989177f,    -0.995185f,
   -0.998795f,    -1.000000f,    -0.998795f,    -0.995185f,    -0.989177f,
   -0.980785f,    -0.970031f,    -0.956940f,    -0.941544f,    -0.923880f,
   -0.903989f,    -0.881921f,    -0.857729f,    -0.831470f,    -0.803208f,
   -0.773010f,    -0.740951f,    -0.707107f,    -0.671559f,    -0.634393f,
   -0.595699f,    -0.555570f,    -0.514103f,    -0.471397f,    -0.427555f,
   -0.382683f,    -0.336890f,    -0.290285f,    -0.242980f,    -0.195090f,
   -0.146730f,    -0.098017f,    -0.049068f
};

//scale factor: 2^20
int lookup_cos_32_inter[128] = {

     1048576,      1047313,      1043527,      1037227,      1028428,
     1017151,      1003425,       987281,       968758,       947901,
      924761,       899394,       871859,       842224,       810560,
      776944,       741455,       704181,       665210,       624636,
      582558,       539076,       494295,       448324,       401273,
      353255,       304386,       254783,       204567,       153858,
      102778,        51451,           -0,       -51451,      -102778,
     -153858,      -204567,      -254783,      -304386,      -353255,
     -401273,      -448324,      -494295,      -539076,      -582558,
     -624636,      -665210,      -704181,      -741455,      -776944,
     -810560,      -842224,      -871859,      -899394,      -924761,
     -947901,      -968758,      -987281,     -1003425,     -1017151,
    -1028428,     -1037227,     -1043527,     -1047313,     -1048576,
    -1047313,     -1043527,     -1037227,     -1028428,     -1017151,
    -1003425,      -987281,      -968758,      -947901,      -924761,
     -899394,      -871859,      -842224,      -810560,      -776944,
     -741455,      -704181,      -665210,      -624636,      -582558,
     -539076,      -494295,      -448324,      -401273,      -353255,
     -304386,      -254783,      -204567,      -153858,      -102778,
      -51451,           -0,        51451,       102778,       153858,
      204567,       254783,       304386,       353255,       401273,
      448324,       494295,       539076,       582558,       624636,
      665210,       704181,       741455,       776944,       810560,
      842224,       871859,       899394,       924761,       947901,
      968758,       987281,      1003425,      1017151,      1028428,
     1037227,      1043527,      1047313
};

//scale factor: 2^20
int lookup_sin_32_inter[128] = {

           0,        51451,       102778,       153858,       204567,
      254783,       304386,       353255,       401273,       448324,
      494295,       539076,       582558,       624636,       665210,
      704181,       741455,       776944,       810560,       842224,
      871859,       899394,       924761,       947901,       968758,
      987281,      1003425,      1017151,      1028428,      1037227,
     1043527,      1047313,      1048576,      1047313,      1043527,
     1037227,      1028428,      1017151,      1003425,       987281,
      968758,       947901,       924761,       899394,       871859,
      842224,       810560,       776944,       741455,       704181,
      665210,       624636,       582558,       539076,       494295,
      448324,       401273,       353255,       304386,       254783,
      204567,       153858,       102778,        51451,            0,
      -51451,      -102778,      -153858,      -204567,      -254783,
     -304386,      -353255,      -401273,      -448324,      -494295,
     -539076,      -582558,      -624636,      -665210,      -704181,
     -741455,      -776944,      -810560,      -842224,      -871859,
     -899394,      -924761,      -947901,      -968758,      -987281,
    -1003425,     -1017151,     -1028428,     -1037227,     -1043527,
    -1047313,     -1048576,     -1047313,     -1043527,     -1037227,
    -1028428,     -1017151,     -1003425,      -987281,      -968758,
     -947901,      -924761,      -899394,      -871859,      -842224,
     -810560,      -776944,      -741455,      -704181,      -665210,
     -624636,      -582558,      -539076,      -494295,      -448324,
     -401273,      -353255,      -304386,      -254783,      -204567,
     -153858,      -102778,       -51451
};


typedef int32_t fixed_t;
typedef int64_t dw_fixed_t;

float wrap2pi_float(float angle) {
    float twopi = 2.0f * 3.141592f;
    return angle - twopi * floor ( angle / twopi);
}

fixed_t wrap2pi_32(fixed_t angle)
{
    fixed_t twopi_32 = 6588397;

    if (angle < 0)
    {
        return angle - twopi_32 * ((angle / twopi_32) - 1);
    }
    else
    {
        return angle - twopi_32 * (angle / twopi_32);
    }
}

int lookupindex(float angle) {
    float a = wrap2pi_float(angle);
    float b = 3.141592f / 64.0f;
    int index = (int)(a / b);
    return index;
}

fixed_t lookupindex_32(fixed_t angle)
{
    fixed_t a = wrap2pi_32(angle);

    fixed_t b = 51472;
    int index = (a / b);
    return index;
}

float lookup_cos(float angle) {
    return lookup_cos_float[lookupindex(angle)];
}

float lookup_sin(float angle) {
    return lookup_sin_float[lookupindex(angle)];
}


fixed_t lookup_cos_32(fixed_t angle) {
    return lookup_cos_32_inter[lookupindex_32(angle)];
}

fixed_t lookup_sin_32(fixed_t angle) {
    return lookup_sin_32_inter[lookupindex_32(angle)];
}

float fix2float(fixed_t in) {
    return ((float) in / (float) (1 << FIXED_FRACTIONAL_PART));
}

fixed_t float2fix(float in) {
    return (fixed_t) (in * (1 << FIXED_FRACTIONAL_PART));
}

fixed_t fixed_multiply(fixed_t a, fixed_t b) {
    return ((dw_fixed_t) a * (dw_fixed_t) b) / (1 << 16);
}

float aprx_atan_float(float a)
{
    //printf("1/x: %f\n", a);

    if(a > 0.5f && a <= 1.0f)
    {
        return 0.644f * a + 0.142f;
    }
    else if (a >= -0.5f && a <= 0.5f)
    {
        //printf("0.928f * a: %f\n", 0.928f * a);
        //printf("(0.928f * a) * 2^20: %f\n", (0.928f * a) * 1048576);
        return 0.928f * a;
    }
    else
    {
        return 0.644f * a - 0.142f;
    }
}

float aprx_atan2_float(float x)
{
    //printf("\nx: %f\n", x);

    if(abs(x) <= 1.0f)
    {
        return aprx_atan_float(x);
    }
    else if (x > 1.0f)
    {
        return 3.141592f / 2.0f  -  aprx_atan_float(1 / x);
    }
    else
    {
        return -3.141592f / 2.0f  -  aprx_atan_float(1 / x);
    }
}

fixed_t aprx_atan(fixed_t a)
{
    if(a > 524288 && a <= 1048576)
    {
        return (fixed_t)(((675283 * (dw_fixed_t)(a)) >> 20) + 148898);
    }
    else if (a >= -524288 && a <= 524288)
    {
        //printf("((973079 * a) >> 20): %d\n", ((973079 * a) >> 20));
        return (fixed_t)((973079 * (dw_fixed_t)(a)) >> 20);
    }
    else
    {
        return (fixed_t)(((675283 * (dw_fixed_t)(a)) >> 20) - 148898);
    }
}

fixed_t aprx_atan2(fixed_t x)
{
    //printf("\nx: %d\n", x);

    //printf("((dw_fixed_t)(1) << 63) / x: %lld\n", ((dw_fixed_t)(1) << 62) / x);

    //int x_reciprocal = ((1 << 31) / x) << 9;
    //fixed_t x_reciprocal = (fixed_t)((((dw_fixed_t)(1) << 63) / x) >> 23);
        //couldn't use cuz sign was being lost, not sure how

    fixed_t x_reciprocal = (fixed_t)((((dw_fixed_t)(1) << 62) / x) >> 22);

    //printf("1/x: %d\n", x_reciprocal);

    if(abs(x) <= 1048576)
    {
        return aprx_atan(x);
    }
    else if (x > 1048576)
    {

        return 1647099 - aprx_atan(x_reciprocal);
    }
    else
    {
        return -1647099 - aprx_atan(x_reciprocal);
    }
}

/*float aprx_cot_float(float x) {
    if(x > 1 && x <= 1.77) {
        return -0.347 * x + 1.12;
    }
    else if (x > 1.77 && x <= 3.307) {
        return -0.138 * x + 0.75;
    }
    else if (x >


    }
}*/

float taylor_sin(float x) {
    return x - (x * x * x) * invfact3 + (x * x * x * x * x) * invfact5;
}

float taylor_cos(float x) {
    return 1.0f - (x * x) * invfact2 + (x * x * x * x) * invfact4;
}

fixed_t taylor_sin_fix(fixed_t x) {
    fixed_t x2 = fixed_multiply(x, x);
    fixed_t x3 = fixed_multiply(x2, x);
    fixed_t x5 = fixed_multiply(x3, x2);
    return x - x3 * invfact3fix + x5 * invfact5fix;
}

fixed_t taylor_cos_fix(fixed_t x) {
    fixed_t x2 = fixed_multiply(x, x);
    fixed_t x4 = fixed_multiply(x2, x2);
    return float2fix(1.0f) - x2 * invfact2fix + x4 * invfact4fix;
}

/*float taylor_atan(float x) {
    return x - (x * x * x)
}*/

//prints a 2x2 matrix with floating-point numbers
void print2f(float temp[2][2])
{
    for (unsigned int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < 2; j++)
        {
            printf("%12f ", temp[i][j]);
        }

        printf("\n");
    }

    printf("\n");
}

//prints a 4x4 matrix with floating-point numbers
void print4f(float temp[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            printf("%12f ", temp[i][j]);
        }

        printf("\n");
    }

    printf("\n");
}

//prints a 4x4 matrix with integers
void print4i(int temp[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            printf("%12i ", temp[i][j]);
        }

        printf("\n");
    }

    printf("\n");
}

//transposes a matrix
void transpose(float* restrict p, float mat[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        *(p + i * 4 + 0) = mat[0][i];
        *(p + i * 4 + 1) = mat[1][i];
        *(p + i * 4 + 2) = mat[2][i];
        *(p + i * 4 + 3) = mat[3][i];
    }
}

void transpose_32(fixed_t* restrict p, fixed_t mat[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        *(p + i * 4 + 0) = mat[0][i];
        *(p + i * 4 + 1) = mat[1][i];
        *(p + i * 4 + 2) = mat[2][i];
        *(p + i * 4 + 3) = mat[3][i];
    }
}

//performs matrix multiplication
void multiply(float* restrict p, float temp1[4][4], float temp2[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            for (unsigned int k = 0; k < 4; k++)
            {
                *(p + i * 4 + j) += temp1[i][k] * temp2[k][j];
            }
        }
    }
}

void multiply_32(fixed_t* restrict p, fixed_t temp1[4][4], fixed_t temp2[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            for (unsigned int k = 0; k < 4; k++)
            {
                *(p + i * 4 + j) += (fixed_t)((dw_fixed_t)(temp1[i][k] * (dw_fixed_t)(temp2[k][j])) >> 20);
            }
        }
    }
}

//sets all matrix indices to 0
void zero(float* restrict p)
{
    for (unsigned int i = 0; i < 4; i++)
    {
        *(p + i * 4 + 0) = 0;
        *(p + i * 4 + 1) = 0;
        *(p + i * 4 + 2) = 0;
        *(p + i * 4 + 3) = 0;
    }
}

void zero_32(fixed_t* restrict p)
{
    for (unsigned int i = 0; i < 4; i++)
    {
        *(p + i * 4 + 0) = 0;
        *(p + i * 4 + 1) = 0;
        *(p + i * 4 + 2) = 0;
        *(p + i * 4 + 3) = 0;
    }
}

//copies a matrix into another matrix
void copy(float* restrict p, float mat[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            *(p + i * 4 + j) = mat[i][j];
        }
    }
}

void copy_32(fixed_t* restrict p, fixed_t mat[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            *(p + i * 4 + j) = mat[i][j];
        }
    }
}

int testfunc() {
    fixed_t inv2 = float2fix(invfact2);
    fixed_t inv3 = float2fix(invfact3);
    fixed_t inv4 = float2fix(invfact4);
    fixed_t inv5 = float2fix(invfact5);
    fixed_t inv6 = float2fix(invfact6);
    fixed_t inv7 = float2fix(invfact7);
    fixed_t inv8 = float2fix(invfact8);
    fixed_t thrity = float2fix(30.0f);
    printf( " Fix: %d\n" , inv2 );
    printf( " Fix: %d\n" , inv3);
    printf( " Fix: %d\n" , inv4 );
    printf( " Fix: %d\n" , inv5 );
    printf( " Fix: %d\n" , inv6 );
    printf( " Fix: %d\n" , inv7 );
    printf( " Fix: %d\n" , inv8 );
    printf( " Fix: %d\n" , thrity );
    //float flt = fix2float(fix);
    //printf(" Float: %f\n", flt);
}

int main()
{
    struct timeval start, end;
    gettimeofday(&start, NULL);

    //initialize matrices
    float M[4][4] = {{31, 77, -11, 26}, {-42, 14, 79, -53}, {-68, -10, 45, 90}, {34, 16, 38, -19}};
    float U[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    float V[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

    short int M_12[4][4] = {{31, 77, -11, 26}, {-42, 14, 79, -53}, {-68, -10, 45, 90}, {34, 16, 38, -19}};
    fixed_t M_32[4][4];                                                                                                     //scale factor: 2^31 / 2^11 = 2^20
    //fixed_t U_32[4][4] = {{0x7FFFFFFF, 0, 0, 0}, {0, 0x7FFFFFFF, 0, 0}, {0, 0, 0x7FFFFFFF, 0}, {0, 0, 0, 0x7FFFFFFF}};    //scale factor: 2^31 / 2^0 = 2^31
    //fixed_t V_32[4][4] = {{0x7FFFFFFF, 0, 0, 0}, {0, 0x7FFFFFFF, 0, 0}, {0, 0, 0x7FFFFFFF, 0}, {0, 0, 0, 0x7FFFFFFF}};	//scale factor: 2^31
    fixed_t U_32[4][4] = {{1048576, 0, 0, 0}, {0, 1048576, 0, 0}, {0, 0, 1048576, 0}, {0, 0, 0, 1048576}};                  //scale factor: 2^20
    fixed_t V_32[4][4] = {{1048576, 0, 0, 0}, {0, 1048576, 0, 0}, {0, 0, 1048576, 0}, {0, 0, 0, 1048576}};	                //scale factor: 2^20


    //modified U and V matrices
    float U_mod[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    float V_mod[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

    //fixed_t U_mod_32[4][4] = {{0x7FFFFFFF, 0, 0, 0}, {0, 0x7FFFFFFF, 0, 0}, {0, 0, 0x7FFFFFFF, 0}, {0, 0, 0, 0x7FFFFFFF}};	//scale factor: 2^31
    //fixed_t V_mod_32[4][4] = {{0x7FFFFFFF, 0, 0, 0}, {0, 0x7FFFFFFF, 0, 0}, {0, 0, 0x7FFFFFFF, 0}, {0, 0, 0, 0x7FFFFFFF}};	//scale factor: 2^31
    fixed_t U_mod_32[4][4] = {{1048576, 0, 0, 0}, {0, 1048576, 0, 0}, {0, 0, 1048576, 0}, {0, 0, 0, 1048576}};	                //scale factor: 2^20
    fixed_t V_mod_32[4][4] = {{1048576, 0, 0, 0}, {0, 1048576, 0, 0}, {0, 0, 1048576, 0}, {0, 0, 0, 1048576}};	                //scale factor: 2^20

    //transposes of modified U and V matrices
    float U_mod_T[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    float V_mod_T[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    fixed_t U_mod_T_32[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};	//scale factor: 2^20
    fixed_t V_mod_T_32[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};	//scale factor: 2^20

    //temporary matrix
    float temp[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    fixed_t temp_32[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    //angles
    float theta_sum;
    float theta_diff;
    float theta_l;
    float theta_r;

    fixed_t theta_sum_32;
    fixed_t theta_diff_32;
    fixed_t theta_l_32;
    fixed_t theta_r_32;

    //matrix indices (use register int later)
    float Mii;
    float Mij;
    float Mji;
    float Mjj;

    fixed_t Mii_32;
    fixed_t Mij_32;
    fixed_t Mji_32;
    fixed_t Mjj_32;

    //matrix index operations
    float sum_num;
    float sum_denom;
    float diff_num;
    float diff_denom;
    float sum_quot;
    float diff_quot;

    fixed_t sum_num_32;
    fixed_t sum_denom_32;
    fixed_t diff_num_32;
    fixed_t diff_denom_32;
    fixed_t sum_quot_32;
    fixed_t diff_quot_32;

    fixed_t cos_l;
    fixed_t sin_l;
    fixed_t cos_r;
    fixed_t sin_r;

    printf("Start\n\nM:\n");
    print4f(M);
    printf("\n\n");

    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            M_32[i][j] = (fixed_t)round(M_12[i][j] * 1048576);
        }
    }

    printf("Modified, M_32 = round(M_12 * 2^20)\n\nM_32:\n");
    print4i(M_32);
    printf("\n\n");

    //main program
    //for (unsigned int sweep = 1; sweep < 2; sweep++)
    for (unsigned int sweep = 1; sweep < 30; sweep++)
    {
        //double for loop for matrix indexing
        for (unsigned int i = 0; i < 3; i++)
        {
            for (unsigned int j = i+1; j < 4; j++)
            {
                if (debug_iterations)
                {
                    printf("Sweep %d, Pair (%d - %d)\n\n", sweep, i+1, j+1);
                }

                //floating-point operations
                if (floating_arithmetic)
                {
                    Mii = M[i][i];
                    Mij = M[i][j];
                    Mji = M[j][i];
                    Mjj = M[j][j];

                    sum_num = Mji + Mij;
                    sum_denom = Mjj - Mii;
                    diff_num = Mji - Mij;
                    diff_denom = Mjj + Mii;

                    sum_quot = sum_num / sum_denom;
                    diff_quot = diff_num / diff_denom;

                    theta_sum = aprx_atan2_float(sum_quot);
                    theta_diff = aprx_atan2_float(diff_quot);

                    theta_l = (theta_sum - theta_diff) / 2;
                    theta_r = (theta_sum + theta_diff) / 2;

                    U_mod[i][i] = cos(theta_l);
                    U_mod[i][j] = -sin(theta_l);
                    U_mod[j][i] = sin(theta_l);
                    U_mod[j][j] = cos(theta_l);

                    V_mod[i][i] = cos(theta_r);
                    V_mod[i][j] = -sin(theta_r);
                    V_mod[j][i] = sin(theta_r);
                    V_mod[j][j] = cos(theta_r);

                    /*
                    U_mod[i][i] = lookup_cos(theta_l);
                    U_mod[i][j] = -lookup_sin(theta_l);
                    U_mod[j][i] = lookup_sin(theta_l);
                    U_mod[j][j] = lookup_cos(theta_l);

                    V_mod[i][i] = lookup_cos(theta_r);
                    V_mod[i][j] = -lookup_sin(theta_r);
                    V_mod[j][i] = lookup_sin(theta_r);
                    V_mod[j][j] = lookup_cos(theta_r);
                    */

                    transpose(*U_mod_T, U_mod);
                    transpose(*V_mod_T, V_mod);

                    copy(*temp, U);
                    zero(*U);
                    multiply(*U, temp, U_mod_T);

                    copy(*temp, V);
                    zero(*V);
                    multiply(*V, temp, V_mod_T);

                    zero(*temp);
                    multiply(*temp, U_mod, M);

                    zero(*M);
                    multiply(*M, temp, V_mod_T);

                    zero(*U_mod);
                    zero(*V_mod);

                    U_mod[0][0] = 1;
                    U_mod[1][1] = 1;
                    U_mod[2][2] = 1;
                    U_mod[3][3] = 1;

                    V_mod[0][0] = 1;
                    V_mod[1][1] = 1;
                    V_mod[2][2] = 1;
                    V_mod[3][3] = 1;
                }


                //fixed-point operations
                Mii_32 = M_32[i][i];
                Mij_32 = M_32[i][j];
                Mji_32 = M_32[j][i];
                Mjj_32 = M_32[j][j];

                sum_num_32 = Mji_32 + Mij_32;
                sum_denom_32 = Mjj_32 - Mii_32;
                diff_num_32 = Mji_32 - Mij_32;
                diff_denom_32 = Mjj_32 + Mii_32;

                sum_quot_32 = (fixed_t)((((dw_fixed_t)(sum_num_32) << 31) / sum_denom_32) >> 11);
                diff_quot_32 = (fixed_t)((((dw_fixed_t)(diff_num_32) << 31) / diff_denom_32) >> 11);

                theta_sum_32 = aprx_atan2(sum_quot_32);
                theta_diff_32 = aprx_atan2(diff_quot_32);

                theta_l_32 = (theta_sum_32 - theta_diff_32) >> 1;
                theta_r_32 = (theta_sum_32 + theta_diff_32) >> 1;

                cos_l = lookup_cos_32(theta_l_32);
                sin_l = lookup_sin_32(theta_l_32);
                cos_r = lookup_cos_32(theta_r_32);
                sin_r = lookup_sin_32(theta_r_32);

                U_mod_32[i][i] = cos_l;
                U_mod_32[i][j] = -sin_l;
                U_mod_32[j][i] = sin_l;
                U_mod_32[j][j] = cos_l;

                V_mod_32[i][i] = cos_r;
                V_mod_32[i][j] = -sin_r;
                V_mod_32[j][i] = sin_r;
                V_mod_32[j][j] = cos_r;

                transpose_32(*U_mod_T_32, U_mod_32);
                transpose_32(*V_mod_T_32, V_mod_32);

                copy_32(*temp_32, U_32);
                zero_32(*U_32);
                multiply_32(*U_32, temp_32, U_mod_T_32);

                copy_32(*temp_32, V_32);
                zero_32(*V_32);
                multiply_32(*V_32, temp_32, V_mod_T_32);

                zero_32(*temp_32);
                multiply_32(*temp_32, U_mod_32, M_32);

                zero_32(*M_32);
                multiply_32(*M_32, temp_32, V_mod_T_32);

                zero_32(*U_mod_32);
                zero_32(*V_mod_32);

                U_mod_32[0][0] = 1048576;
                U_mod_32[1][1] = 1048576;
                U_mod_32[2][2] = 1048576;
                U_mod_32[3][3] = 1048576;

                V_mod_32[0][0] = 1048576;
                V_mod_32[1][1] = 1048576;
                V_mod_32[2][2] = 1048576;
                V_mod_32[3][3] = 1048576;


                //debugging printf statements
                if (Mkk)
                {
                    printf("M[i][i] = %f\n", Mii);
                    printf("M[i][j] = %f\n", Mij);
                    printf("M[j][i] = %f\n", Mji);
                    printf("M[j][j] = %f\n", Mjj);
                    printf("\n");

                    printf("M_32[i][i] = %d\n", Mii_32);
                    printf("M_32[i][j] = %d\n", Mij_32);
                    printf("M_32[j][i] = %d\n", Mji_32);
                    printf("M_32[j][j] = %d\n", Mjj_32);
                    printf("\n");
                }

                if (num_denom)
                {
                    printf("sum_num =       M[j][i] + M[i][j] =     %f\n", sum_num);
                    printf("sum_denom =     M[j][j] - M[i][i] =     %f\n", sum_denom);
                    printf("diff_num =      M[j][i] - M[i][j] =     %f\n", diff_num);
                    printf("diff_denom =    M[j][j] + M[i][i] =     %f\n", diff_denom);
                    printf("\n");

                    printf("sum_num_32 =       M[j][i]_32 + M[i][j]_32 =     %d\n", sum_num_32);
                    printf("sum_denom_32 =     M[j][j]_32 - M[i][i]_32 =     %d\n", sum_denom_32);
                    printf("diff_num_32 =      M[j][i]_32 - M[i][j]_32 =     %d\n", diff_num_32);
                    printf("diff_denom_32 =    M[j][j]_32 + M[i][i]_32 =     %d\n", diff_denom_32);
                    printf("\n");
                }

                if (quot)
                {
                    printf("sum_quot =      sum_num / sum_denom =       %f\n", sum_quot);
                    printf("diff_quot =     diff_num / diff_denom =     %f\n", diff_quot);
                    printf("\n");

                    printf("sum_quot_32 =      sum_num_32 / (sum_denom_32 >> 20) =       %d\n", sum_quot_32);
                    printf("diff_quot_32 =     diff_num_32 / (diff_denom_32 >> 20) =     %d\n", diff_quot_32);
                    printf("\n");
                }

                if (theta_sumdiff)
                {
                    printf("theta_sum =     atan(sum_quot) =        %f\n", atan(sum_quot));
                    printf("theta_diff =    atan(diff_quot) =       %f\n", atan(diff_quot));
                    printf("\n");

                    printf("theta_sum =     aprx_atan2_float(sum_quot) =        %f\n", theta_sum);
                    printf("theta_diff =    aprx_atan2_float(diff_quot) =       %f\n", theta_diff);
                    printf("\n");

                    printf("theta_sum * 2^20 =     aprx_atan2_float(sum_quot) =        %f\n", round(theta_sum * 1048576));
                    printf("theta_diff * 2^20 =    aprx_atan2_float(diff_quot) =       %f\n", round(theta_diff * 1048576));
                    printf("\n");

                    printf("theta_sum_32 =     aprx_atan(sum_quot_32) =        %d\n", theta_sum_32);
                    printf("theta_diff_32 =    aprx_atan(diff_quot_32) =       %d\n", theta_diff_32);
                    printf("\n");

                }

                if (theta_LR)
                {
                    printf("theta_l =       (theta_sum - theta_diff) / 2 =    %f\n", theta_l);
                    printf("theta_r =       (theta_sum + theta_diff) / 2 =    %f\n", theta_r);
                    printf("\n");

                    printf("theta_l_32 =       (theta_sum_32 - theta_diff_32) / 2 =    %d\n", theta_l_32);
                    printf("theta_r_32 =       (theta_sum_32 + theta_diff_32) / 2 =    %d\n", theta_r_32);
                    printf("\n");

                }

                if (wrap2pi_debug)
                {
                    printf("wrap2pi_float(theta_l): %f\n", wrap2pi_float(theta_l));
                    printf("wrap2pi_float(theta_r): %f\n", wrap2pi_float(theta_r));
                    printf("\n");

                    printf("wrap2pi_32(theta_l_32): %d\n", wrap2pi_32(theta_l_32));
                    printf("wrap2pi_32(theta_r_32): %d\n", wrap2pi_32(theta_r_32));
                    printf("\n");
                }

                if (lookupindex_debug)
                {
                    printf("lookupindex(theta_l): %d\n", lookupindex(theta_l));
                    printf("lookupindex(theta_r): %d\n", lookupindex(theta_r));
                    printf("\n");

                    printf("lookupindex_32(theta_l_32): %d\n", lookupindex_32(theta_l_32));
                    printf("lookupindex_32(theta_r_32): %d\n", lookupindex_32(theta_r_32));
                    printf("\n");
                }

                if (cos_sin_LR)
                {
                    printf("cos(theta_l): %f\n", cos(theta_l));
                    printf("sin(theta_l): %f\n", sin(theta_l));
                    printf("cos(theta_r): %f\n", cos(theta_r));
                    printf("sin(theta_r): %f\n", sin(theta_r));
                    printf("\n");

                    printf("cos(theta_l) * 2^20: %f\n", cos(theta_l) * 1048576);
                    printf("sin(theta_l) * 2^20: %f\n", sin(theta_l) * 1048576);
                    printf("cos(theta_r) * 2^20: %f\n", cos(theta_r) * 1048576);
                    printf("sin(theta_r) * 2^20: %f\n", sin(theta_r) * 1048576);
                    printf("\n");
                }

                if (lookup_cos_sin)
                {
                    printf("lookup_cos(theta_l): %f\n", lookup_cos(theta_l));
                    printf("lookup_sin(theta_l): %f\n", lookup_sin(theta_l));
                    printf("lookup_cos(theta_r): %f\n", lookup_cos(theta_r));
                    printf("lookup_sin(theta_r): %f\n", lookup_sin(theta_r));
                    printf("\n");

                    printf("lookup_cos(theta_l) * 2^20: %f\n", lookup_cos(theta_l) * 1048576);
                    printf("lookup_sin(theta_l) * 2^20: %f\n", lookup_sin(theta_l) * 1048576);
                    printf("lookup_cos(theta_r) * 2^20: %f\n", lookup_cos(theta_r) * 1048576);
                    printf("lookup_sin(theta_r) * 2^20: %f\n", lookup_sin(theta_r) * 1048576);
                    printf("\n");

                    printf("cos_l: %d\n", cos_l);
                    printf("sin_l: %d\n", sin_l);
                    printf("cos_r: %d\n", cos_r);
                    printf("sin_r: %d\n", sin_r);
                    printf("\n");
                }

                if (debug_matrices)
                {
                    printf("U:\n");
                    print4f(U);

                    printf("V:\n");
                    print4f(V);

                    printf("U_32:\n");
                    print4i(U_32);

                    printf("V_32:\n");
                    print4i(V_32);
                }

                if (debug_iterations)
                {
                    printf("M:\n");
                    print4f(M);
                    printf("\n\n");

                    printf("M_32:\n");
                    print4i(M_32);
                    printf("\n\n");
                }
            }
        }
    }

    printf("End\n\n");

    if (floating_arithmetic)
    {
        printf("M:\n");
        print4f(M);
        printf("\n\n");
    }

    printf("M_32:\n");
    print4i(M_32);
    printf("\n\n");

    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            M[i][j] = (float)(M_32[i][j]) / 1048576 * 1.0;
        }
    }

    printf("M (M_32 scaled down):\n");
    print4f(M);
    printf("\n\n");

    gettimeofday(&end, NULL);
    double time_taken = (end.tv_sec * 1 + (end.tv_usec) / 1 ) - (start.tv_sec * 1 + (start.tv_usec) / 1 ); // in seconds
    printf("time program took %d microseconds to execute\n", (int)time_taken);


    return 0;
}