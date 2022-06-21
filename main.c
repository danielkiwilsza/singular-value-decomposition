#include <stdio.h>
#include <math.h>

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
void transpose(float* p, float mat[4][4])
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            *(p + i * 4 + j) = mat[j][i];
        }
    }
}

//performs matrix multiplication
void multiply(float* p, float temp1[4][4], float temp2[4][4])
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
void zero(float* p)
{
    for (unsigned int i = 0; i < 4; i++)
    {
        for (unsigned int j = 0; j < 4; j++)
        {
            *(p + i * 4 + j) = 0;
        }
    }
}

//copies a matrix into another matrix
void copy(float* p, float mat[4][4])
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

                theta_sum = atan( (M[j][i] + M[i][j]) / (M[j][j] - M[i][i]) );
                theta_diff = atan( (M[j][i] - M[i][j]) / (M[j][j] + M[i][i]) );

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

                /*
                printf("\ncos(theta_l): %f\n", cos(theta_l));
                printf("sin(theta_l): %f\n", sin(theta_l));
                printf("cos(theta_r): %f\n", cos(theta_r));
                printf("sin(theta_r): %f\n\n", sin(theta_r));
                */

                printf("M:\n");
                print4f(M);
                printf("\n\n");

                /*
                printf("U:\n");
                print4f(U);

                printf("V:\n");
                print4f(V);
                */

            }
        }
    }

    return 0;
}