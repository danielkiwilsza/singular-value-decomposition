#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>

#define debug 1
#define debug_matrices 1

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
    for (unsigned int sweep = 1; sweep < 5; sweep++)
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

                theta_sum = atan(sum_quot);
                theta_diff = atan(diff_quot);

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

                    printf("cos(theta_l): %f\n", cos(theta_l));
                    printf("sin(theta_l): %f\n", sin(theta_l));
                    printf("cos(theta_r): %f\n", cos(theta_r));
                    printf("sin(theta_r): %f\n", sin(theta_r));
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