// Copyright 2022 Rogelio_Kiera

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>

#define random

bool input_matrix(double** matrix, int n, int m);
void output_matrix(double** a, int n, int m);
bool input_number(int* n);

void sle(double** matrix, int n, int m, double* roots);
void direct_Gauss(double** mat, int row, int col);
void inverse_Gauss(double** mat, int row, int col, double* roots);
double det(double** m, int n);
double** invert_matrix(double** mat, int n);
void check(double **mat, double **invmat, int n);

int main(void) {
    int n;
    double t1, t2;
    printf("input size of matrix n = ");
    if (input_number(&n)) {
        double** arr = new double* [n];
        for (int i = 0; i < n; i++)
            arr[i] = new double[n];
#ifdef random
        srand(time(NULL));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                arr[i][j] = (double)rand() / RAND_MAX;
        output_matrix(arr, n, n);
        printf("\n\n");
        if (!(det(arr, n) <= 1e-12 && det(arr, n) >= -1e-12)) {
            double** temp = new double* [n];
            for (int i = 0; i < n; i++)
                temp[i] = new double[n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    temp[i][j] = arr[i][j];
            t1 = omp_get_wtime();
            arr = invert_matrix(arr, n);
            t2 = omp_get_wtime();
            output_matrix(arr, n, n);
            check(temp, arr, n);
            printf("Time in sec: %0.10lf, in min: %0.10lf\n", t2 - t1, (t2 - t1)/60.0);
            for (int i = 0; i < n; i++)
                delete[]temp[i];
            delete[]temp;
        } else {
            printf("n/a\n");
        }
#else
        if (input_matrix(arr, n, n)) {
            if (!(det(arr, n) <= 1e-12 && det(arr, n) >= -1e-12)) {
                double** mat = new double* [n];
                for (int i = 0; i < n; i++)
                    mat[i] = new double[n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        mat[i][j] = arr[i][j];
                t1 = omp_get_wtime();
                arr = invert_matrix(arr, n);
                t2 = omp_get_wtime();
                output_matrix(arr, n, n);
                check(mat, arr, n);
                printf("Time in sec: %0.10lf, in min: %0.10lf\n", t2 - t1, (t2 - t1) / 60.0);
                for (int i = 0; i < n; i++)
                    delete[]mat[i];
                delete[]mat;
            }
            else {
                printf("n/a\n");
            }
        }
#endif
        for (int i = 0; i < n; i++)
            delete []arr[i];
        delete []arr;
    }
    printf("\n");
    system("pause");
    return 0;
}

void output_matrix(double** a, int n, int m) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (j == m - 1 && i != n - 1) {
                printf("%.6lf\n", a[i][j]);
            }
            else if (i == n - 1 && j == m - 1) {
                printf("%.6lf", a[i][j]);
            }
            else {
                printf("%.6lf ", a[i][j]);
            }
        }
    }
    printf("\n\n");
}

bool input_matrix(double** matrix, int n, int m) {
    int flag = true;
    char c;
    for (int i = 0; i < n && flag; i++) {
        for (int j = 0; j < m && flag; j++) {
            if (scanf_s("%lf%c", &matrix[i][j], &c) == 2 && (c == ' ' || c == '\n')) {}
            else {
                printf("n/a\n");
                flag = false;
            }
        }
    }
    return flag;
}

bool input_number(int* n) {
    bool flag = true;
    char c;
    if (scanf_s("%d%c", n, &c) == 2 && (c == ' ' || c == '\n')) {}
    else {
        printf("n/a\n");
        flag = false;
    }
    return flag;
}

void sle(double** matrix, int n, int m, double* roots) {
    direct_Gauss(matrix, n, m);
    inverse_Gauss(matrix, n, m, roots);
}

void direct_Gauss(double** mat, int row, int col) {
    double d = 0.0;
    double* c = (double*)malloc(col * sizeof(double));
    if (c != NULL) {
        for (int f = 0; f < row; f++) {
            for (int k = f; k < row - 1; k++) {
                for (int j = f; j < row - 1; j++) {
                    if (mat[j][f] == 0.0) {
                        for (int i = 0; i < col; i++)
                            c[i] = mat[j][i];
                        for (int i = 0; i < col; i++)
                            mat[j][i] = mat[j + 1][i];
                        for (int i = 0; i < col; i++)
                            mat[j + 1][i] = c[i];
                    }
                }
            }
            for (int i = f + 1; i < row; i++) {
                d = mat[i][f];
                for (int j = f; j < col; j++)
                    mat[i][j] = mat[i][j] - mat[f][j] * d / mat[f][f];
            }
        }
    } else {
        printf("error memory allocation\n");
    }
    free(c);
}

void inverse_Gauss(double** mat, int row, int col, double* roots) {
    double s;
    roots[row - 1] = mat[row - 1][col - 1] / mat[row - 1][col - 2];
    for (int i = row - 2; i >= 0; i--) {
        s = 0.0;
        roots[i] = 0.0;
        for (int j = col - 2; j != i; j--)
            s += mat[i][j] * roots[j];
        roots[i] = (1.0 / mat[i][i]) * (mat[i][col - 1] - s);
    }
}

double** invert_matrix(double** mat, int n) {
    double* b = (double*)malloc(n * sizeof(double));
    double** temp = new double* [n];
    for (int i = 0; i < n; i++)
        temp[i] = new double[n];
    double** g = (double**)malloc(n * (n + 1) * sizeof(double) + n * sizeof(double*));
    if (g != NULL && b != NULL) {
        double* ptr = (double*)(g + n);
        for (int i = 0; i < n; i++)
            g[i] = ptr + (n + 1) * i;
        for (int i = 0; i < n; i++) {
            for (int i1 = 0; i1 < n; ++i1) {
                for (int j = 0; j < n; ++j)
                    g[i1][j] = mat[i1][j];
                i1 == i ? g[i1][n] = 1.0 : g[i1][n] = 0.0;
            }
            sle(g, n, n + 1, b);
            for (int i1 = 0; i1 < n; ++i1)
                temp[i1][i] = b[i1];
        }
    } else {
        printf("error memory allocation");
    }
    free(g);
    free(b);
    for (int i = 0; i < n; i++)
        delete []mat[i];
    delete []mat;
    return temp;
}

double det(double** m, int n) {
    double det = 1.0;
    double** temp = (double**)malloc(n * n * sizeof(double) + n * sizeof(double*));
    if (temp != NULL) {
        double* ptr = (double*)(temp + n);
        for (int i = 0; i < n; i++)
            temp[i] = ptr + n * i;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                temp[i][j] = m[i][j];
        direct_Gauss(temp, n, n);
        for (int i = 0; i < n; i++)
            det *= temp[i][i];
    } else {
        printf("error memory allocation\n");
    }
    free(temp);
    return det;
}

void check(double** mat, double** invmat, int n) {
    double norma = 0.0;
    double** temp = (double**)malloc(n * n * sizeof(double) + n * sizeof(double*));
    if (temp != NULL) {
        double* ptr = (double*)(temp + n);
        for (int i = 0; i < n; i++)
            temp[i] = ptr + n * i;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp[i][j] = 0.0;
                for (int k = 0; k < n; k++)
                    temp[i][j] += mat[i][k] * invmat[k][j];
            }
        }
        output_matrix(temp, n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                i != j ? temp[i][j] = temp[i][j] * temp[i][j] : temp[i][j] = (temp[i][j] - 1.0) * (temp[i][j] - 1.0);
                norma += temp[i][j];
            }
        }
        norma = sqrt(norma);
        printf("norma: %10.10E\n\n", norma);
    } else {
        printf("error memory allocation\n");
    }
    free(temp);
}
