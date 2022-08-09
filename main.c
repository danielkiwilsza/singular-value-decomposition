#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/time.h>
#include <assert.h>
#include <stdlib.h>

#define debug 0
#define debug_matrices 0
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

float lookup_sin_float[128] = {
        0.00000f,
        0.04907f,
        0.09802f,
        0.14673f,
        0.19509f,
        0.24298f,
        0.29028f,
        0.33689f,
        0.38268f,
        0.42756f,
        0.47140f,
        0.51410f,
        0.55557f,
        0.59570f,
        0.63439f,
        0.67156f,
        0.70711f,
        0.74095f,
        0.77301f,
        0.80321f,
        0.83147f,
        0.85773f,
        0.88192f,
        0.90399f,
        0.92388f,
        0.94154f,
        0.95694f,
        0.97003f,
        0.98079f,
        0.98918f,
        0.99518f,
        0.99880f,
        1.00000f,
        0.99880f,
        0.99518f,
        0.98918f,
        0.98079f,
        0.97003f,
        0.95694f,
        0.94154f,
        0.92388f,
        0.90399f,
        0.88192f,
        0.85773f,
        0.83147f,
        0.80321f,
        0.77301f,
        0.74095f,
        0.70711f,
        0.67156f,
        0.63439f,
        0.59570f,
        0.55557f,
        0.51410f,
        0.47140f,
        0.42756f,
        0.38268f,
        0.33689f,
        0.29028f,
        0.24298f,
        0.19509f,
        0.14673f,
        0.09802f,
        0.04907f,
        0.00000f,
        -0.04907f,
        -0.09802f,
        -0.14673f,
        -0.19509f,
        -0.24298f,
        -0.29028f,
        -0.33689f,
        -0.38268f,
        -0.42756f,
        -0.47140f,
        -0.51410f,
        -0.55557f,
        -0.59570f,
        -0.63439f,
        -0.67156f,
        -0.70711f,
        -0.74095f,
        -0.77301f,
        -0.80321f,
        -0.83147f,
        -0.85773f,
        -0.88192f,
        -0.90399f,
        -0.92388f,
        -0.94154f,
        -0.95694f,
        -0.97003f,
        -0.98079f,
        -0.98918f,
        -0.99518f,
        -0.99880f,
        -1.00000f,
        -0.99880f,
        -0.99518f,
        -0.98918f,
        -0.98079f,
        -0.97003f,
        -0.95694f,
        -0.94154f,
        -0.92388f,
        -0.90399f,
        -0.88192f,
        -0.85773f,
        -0.83147f,
        -0.80321f,
        -0.77301f,
        -0.74095f,
        -0.70711f,
        -0.67156f,
        -0.63439f,
        -0.59570f,
        -0.55557f,
        -0.51410f,
        -0.47140f,
        -0.42756f,
        -0.38268f,
        -0.33689f,
        -0.29028f,
        -0.24298f,
        -0.19509f,
        -0.14673f,
        -0.09802f,
        -0.04907f
};

float lookup_cos_float[128] = {
        1.00000f,
        0.99880f,
        0.99518f,
        0.98918f,
        0.98079f,
        0.97003f,
        0.95694f,
        0.94154f,
        0.92388f,
        0.90399f,
        0.88192f,
        0.85773f,
        0.83147f,
        0.80321f,
        0.77301f,
        0.74095f,
        0.70711f,
        0.67156f,
        0.63439f,
        0.59570f,
        0.55557f,
        0.51410f,
        0.47140f,
        0.42756f,
        0.38268f,
        0.33689f,
        0.29028f,
        0.24298f,
        0.19509f,
        0.14673f,
        0.09802f,
        0.04907f,
        -0.00000f,
        -0.04907f,
        -0.09802f,
        -0.14673f,
        -0.19509f,
        -0.24298f,
        -0.29028f,
        -0.33689f,
        -0.38268f,
        -0.42756f,
        -0.47140f,
        -0.51410f,
        -0.55557f,
        -0.59570f,
        -0.63439f,
        -0.67156f,
        -0.70711f,
        -0.74095f,
        -0.77301f,
        -0.80321f,
        -0.83147f,
        -0.85773f,
        -0.88192f,
        -0.90399f,
        -0.92388f,
        -0.94154f,
        -0.95694f,
        -0.97003f,
        -0.98079f,
        -0.98918f,
        -0.99518f,
        -0.99880f,
        -1.00000f,
        -0.99880f,
        -0.99518f,
        -0.98918f,
        -0.98079f,
        -0.97003f,
        -0.95694f,
        -0.94154f,
        -0.92388f,
        -0.90399f,
        -0.88192f,
        -0.85773f,
        -0.83147f,
        -0.80321f,
        -0.77301f,
        -0.74095f,
        -0.70711f,
        -0.67156f,
        -0.63439f,
        -0.59570f,
        -0.55557f,
        -0.51410f,
        -0.47140f,
        -0.42756f,
        -0.38268f,
        -0.33689f,
        -0.29028f,
        -0.24298f,
        -0.19509f,
        -0.14673f,
        -0.09802f,
        -0.04907f,
        -0.00000f,
        0.04907f,
        0.09802f,
        0.14673f,
        0.19509f,
        0.24298f,
        0.29028f,
        0.33689f,
        0.38268f,
        0.42756f,
        0.47140f,
        0.51410f,
        0.55557f,
        0.59570f,
        0.63439f,
        0.67156f,
        0.70711f,
        0.74095f,
        0.77301f,
        0.80321f,
        0.83147f,
        0.85773f,
        0.88192f,
        0.90399f,
        0.92388f,
        0.94154f,
        0.95694f,
        0.97003f,
        0.98079f,
        0.98918f,
        0.99518f,
        0.99880f
};

typedef int32_t fixed_t;
typedef int64_t dw_fixed_t;

float wrap2pi_float(float angle) {
    float twopi = 2.0f * 3.141592f;
    return angle - twopi * floor ( angle / twopi);
}

int lookupindex(float angle) {
    float a = wrap2pi_float(angle);
    float b = 3.141592f / 64.0f;
    int index = (int)(a / b);
    return index;
}

float lookup_sin(float angle) {
    printf("Actual: %.4f\n", sin(angle));
    float ret = lookup_sin_float[lookupindex(angle)];
    printf("Approx: %.4f\n", ret);
    return ret;
}

float lookup_cos(float angle) {
    return lookup_cos_float[lookupindex(angle)];
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

float aprx_atan_float(float a) {

    if(a > 0.5f && a <= 1.0f)
        return 0.644f * a + 0.142f;
    else if(a >= -0.5f && a <= 0.5f)
        return 0.928f * a;
    else
        return 0.644f * a - 0.142f;
}

float aprx_atan2_float(float x) {
    if(abs(x) <= 1.0f) 
        return aprx_atan_float(x);
    else if(x > 1.0f)
        return 3.141592f / 2.0f  -  aprx_atan_float(1 / x);
    else 
        return -3.141592f / 2.0f  -  aprx_atan_float(1 / x);
}

float aprx_sin(float x) {
    assert(x >= 0);
    if(x >= 0.0f && x <= 0.5f) {
        return 0.958f * x;
    }
    else if(x > 0.5f && x <= 1.005f) {
        return 0.724f * x + 0.117f;
    }
    else if(x > 1.005f && x <= 1.275f) {
        return 0.418f * x + 0.424f;
    }
    else {
        return 0.147f * x + 0.769f;
    }
}

float aprx_sin2(float x) {
    printf("Gave sin wrap to %.4f\n", x);
    float a = wrap2pi_float(x);
    printf("Sin wrap to %.4f\n", a);
    printf("Actual: %.4f\n", sin(a));
    if(a < 3.141592f / 2.0f) {
        printf("Approx: %.4f\n", aprx_sin(a));
        return aprx_sin(a);
    }
    else if(a < 3.141592f) {
        printf("Approx: %.4f\n", aprx_sin(3.141592f - a));
        return aprx_sin(3.141592f - a);
    }
    else if(a < 1.5f * 3.141592f) {
        printf("Approx: %.4f\n", -aprx_sin(a - 3.141592f));
        return -aprx_sin(a - 3.141592f);
    }
    else {
        printf("Approx: %.4f\n", -aprx_sin(2.0f * 3.141592f - a));
        return -aprx_sin(2.0f * 3.141592f - a);
    }

}

float aprx_cos2(float x) {
    return aprx_sin2(3.141592f / 2.0f - x);
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

    //modified U and V matrices
    float U_mod[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    float V_mod[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

    //transposes of modified U and V matrices
    float U_mod_T[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    float V_mod_T[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    //temporary matrix
    float temp[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};

    //angles
    float theta_sum;
    float theta_diff;
    float theta_l;
    float theta_r;

    //matrix indices (use register int later)
    float Mii;
    float Mij;
    float Mji;
    float Mjj;

    //matrix index operations
    float sum_num;
    float sum_denom;
    float diff_num;
    float diff_denom;
    float sum_quot;
    float diff_quot;


    printf("Start\n\nM:\n");
    print4f(M);
    printf("\n\n");


    //main program
    for (unsigned int sweep = 1; sweep < 10; sweep++)
    {
        //double for loop for matrix indexing
        for (unsigned int i = 0; i < 3; i++)
        {
            for (unsigned int j = i+1; j < 4; j++)
            {
                printf("Sweep %d, Pair (%d - %d)\n\n", sweep, i+1, j+1);

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

                /*
                theta_sum = atan(sum_quot);
                theta_diff = atan(diff_quot);*/
                
                
                theta_sum = aprx_atan2_float(sum_quot);
                theta_diff = aprx_atan2_float(diff_quot);

                theta_l = (theta_sum - theta_diff) / 2;
                theta_r = (theta_sum + theta_diff) / 2;
                
               /* U_mod[i][i] = cos(theta_l);
                U_mod[i][j] = -sin(theta_l);
                U_mod[j][i] = sin(theta_l);
                U_mod[j][j] = cos(theta_l);

                V_mod[i][i] = cos(theta_r);
                V_mod[i][j] = -sin(theta_r);
                V_mod[j][i] = sin(theta_r);
                V_mod[j][j] = cos(theta_r);*/
                
                
                U_mod[i][i] = aprx_cos2(theta_l);
                U_mod[i][j] = -aprx_sin2(theta_l);
                U_mod[j][i] = aprx_sin2(theta_l);
                U_mod[j][j] = aprx_cos2(theta_l);

                V_mod[i][i] = aprx_cos2(theta_r);
                V_mod[i][j] = -aprx_sin2(theta_r);
                V_mod[j][i] = aprx_sin2(theta_r);
                V_mod[j][j] = aprx_cos2(theta_r);

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

                //reset U_mod and V_mod to identity matrices
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

                if (debug)
                {
                    printf("M[i][i] = %f\n", Mii);
                    printf("M[i][j] = %f\n", Mij);
                    printf("M[j][i] = %f\n", Mji);
                    printf("M[j][j] = %f\n", Mjj);
                    printf("\n");

                    printf("sum_num =       M[j][i] + M[i][j] =     %f\n", sum_num);
                    printf("sum_denom =     M[j][j] - M[i][i] =     %f\n", sum_denom);
                    printf("diff_num =      M[j][i] - M[i][j] =     %f\n", diff_num);
                    printf("diff_denom =    M[j][j] + M[i][i] =     %f\n", diff_denom);
                    printf("\n");

                    printf("sum_quot =      sum_num / sum_denom =       %f\n", sum_quot);
                    printf("diff_quot =     diff_num / diff_denom =     %f\n", diff_quot);
                    printf("\n");

                    printf("theta_sum =     atan(sum_quot) =        %f\n", theta_sum);
                    printf("theta_diff =    atan(diff_quot) =       %f\n", theta_diff);
                    printf("\n");

                    printf("theta_l =       (theta_sum - theta_diff) / 2 =    %f\n", theta_l);
                    printf("theta_r =       (theta_sum + theta_diff) / 2 =    %f\n", theta_r);
                    printf("\n");

                    printf("wrap2pi_float(theta_l): %f\n", wrap2pi_float(theta_l));
                    printf("wrap2pi_float(theta_r): %f\n", wrap2pi_float(theta_r));
                    printf("lookupindex(theta_l): %d\n", lookupindex(theta_l));
                    printf("lookupindex(theta_r): %d\n", lookupindex(theta_r));
                    printf("\n");

                    printf("cos(theta_l): %f\n", cos(theta_l));
                    printf("sin(theta_l): %f\n", sin(theta_l));
                    printf("cos(theta_r): %f\n", cos(theta_r));
                    printf("sin(theta_r): %f\n", sin(theta_r));
                    printf("\n");

                    printf("lookup_cos(theta_l): %f\n", lookup_cos(theta_l));
                    printf("lookup_sin(theta_l): %f\n", lookup_sin(theta_l));
                    printf("lookup_cos(theta_r): %f\n", lookup_cos(theta_r));
                    printf("lookup_sin(theta_r): %f\n", lookup_sin(theta_r));
                    printf("\n");
                }

                if (debug_matrices)
                {
                    printf("U:\n");
                    print4f(U);

                    printf("V:\n");
                    print4f(V);
                }

                printf("M:\n");
                print4f(M);
                printf("\n\n");

            }
        }
    }

    gettimeofday(&end, NULL);
    double time_taken = (end.tv_sec * 1000 + (end.tv_usec) / 1000 ) - (start.tv_sec * 1000 + (start.tv_usec) / 1000 ); // in seconds
    printf("time program took %f milliseconds to execute\n", time_taken);

    return 0;
}