// Copyright 2022 Rogeliok

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#define random

int input_matrix(double** matrix, int* n, int* m);
void output_matrix(double** a, int n, int m);
bool input_number(int* n);

void sle(double** matrix, int n, int m, double* roots);
int direct_Gauss(double** mat, int row, int col);
int iverse_Gauss(double** mat, int row, int col, double* roots);
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
            printf("\n\n");
            check(temp, arr, n);
            printf("\n\n");
            printf("Time in sec: %0.10lf, in min: %0.10lf", t2 - t1, (t2 - t1)/60.0);
            for (int i = 0; i < n; i++)
                delete[]temp[i];
            delete[]temp;
        } else {
            printf("n/a");
        }
#else
        if (input_matrix(arr, &n, &n)) {
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
                printf("\n\n");
                check(mat, arr, n);
                printf("\n\n");
                printf("Time in sec: %0.10lf, in min: %0.10lf", t2 - t1, (t2 - t1) / 60.0);
                for (int i = 0; i < n; i++)
                    delete[]mat[i];
                delete[]mat;
            }
            else {
                printf("n/a");
            }
        }
#endif
        for (int i = 0; i < n; i++)
            delete []arr[i];
        delete []arr;
    }
    return 0;
}
#pragma region functions
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
}

int input_matrix(double** matrix, int* n, int* m) {
    int flag = 1;
    char c;
    for (int i = 0; i < *n && flag; i++) {
        for (int j = 0; j < *m && flag; j++) {
            if (scanf_s("%lf%c", &matrix[i][j], &c) == 2 && (c == ' ' || c == '\n')) {}
            else {
                flag = 0;
            }
        }
    }
    if (!flag)
        printf("n/a");
    return flag;
}

bool input_number(int* n) {
    bool flag = true;
    char c;
    if (scanf_s("%d%c", n, &c) == 2 && (c == ' ' || c == '\n')) {}
    else {
        flag = false;
    }
    if (!flag)
        printf("n/a");
    return flag;
}

void sle(double** matrix, int n, int m, double* roots) {
    direct_Gauss(matrix, n, m);
    iverse_Gauss(matrix, n, m, roots);
}

int direct_Gauss(double** mat, int row, int col) {
    double d = 0.0;
    double* c = new double[col];
    for (int f = 0; f < row; f++) {  // свидение к ступенчатому виду
        for (int k = f; k < row - 1; k++) {  // сортировка строк
            for (int j = f; j < row - 1; j++) {
                if (mat[j][f] == 0.0) {  // swap строк
                    for (int i = 0; i < col; i++)
                        c[i] = mat[j][i];
                    for (int i = 0; i < col; i++)
                        mat[j][i] = mat[j + 1][i];
                    for (int i = 0; i < col; i++)
                        mat[j + 1][i] = c[i];
                }
            }
        }
        for (int i = f + 1; i < row; i++) {  // сердце кода
            d = mat[i][f];
            for (int j = f; j < col; j++)
                mat[i][j] = mat[i][j] - mat[f][j] * d / mat[f][f];
        }
    }
    delete []c;
    return 0;
}

int iverse_Gauss(double** mat, int row, int col, double* roots) {
    double s = 0.0;
    roots[row - 1] = mat[row - 1][col - 1] / mat[row - 1][col - 2];
    for (int i = row - 2; i >= 0; i -= 1) {  // нахождение корней
        s = 0.0;
        for (int j = col - 2; j > 0; j -= 1)
            s += mat[i][j] * roots[j];
        roots[i] = (1.0 / mat[i][i]) * (mat[i][col - 1] - s);
    }
    return 0;
}

double** invert_matrix(double** mat, int n) {
    //double* b = malloc(n * sizeof(double));
    double* b = new double[n];
    /*double** g = malloc(n * (n + 1) * sizeof(double) + n * sizeof(double*));
    double* ptr1 = (double*)(g + n);
    for (int i = 0; i < n; i++)
        g[i] = ptr1 + (n + 1) * i;*/
    double **g = new double* [n];
    for (int i = 0; i < n; i++)
        g[i] = new double[n + 1];
    double** m = new double* [n];
    for (int i = 0; i < n; i++)
        m[i] = new double[n];
    for (int i = 0; i < n; i++) {  // основной цикл
        for (int i1 = 0; i1 < n; ++i1) {  // заполнене матрицы g для гаусса
            for (int j = 0; j < n; ++j)
                g[i1][j] = mat[i1][j];
        }
        for (int j = 0; j < n; j++) {
            if (j == i) g[j][n] = 1.0;
            else g[j][n] = 0.0;
        }
        for (int j = 0; j < n; j++)
            b[j] = 0.0;
        sle(g, n, n + 1, b);
        for (int i1 = 0; i1 < n; ++i1)  // зануление g от греха подальше
            for (int j = 0; j < n; ++j)
                g[i1][j] = 0.0;
        for (int j = 0; j < n; j++)  // передача значений элементов обратной матрицы из цикла
            m[j][i] = b[j];
    }
    //free(g);
    //free(b);
    for (int i = 0; i < n; i++) {
        delete []g[i];
        delete []mat[i];
    }
    delete []g;
    delete []mat;
    delete []b;
    return m;
}

double det(double** m, int n) {
    double** temp = new double* [n];
    for (int i = 0; i < n; i++)
        temp[i] = new double[n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            temp[i][j] = m[i][j];
    double det = 1.0;
    direct_Gauss(temp, n, n);
    for (int i = 0; i < n; i++)
        det *= temp[i][i];
    for (int i = 0; i < n; i++)
        delete []temp[i];
    delete []temp;
    return det;
}

void check(double** mat, double** invmat, int n) {
    double** temp = new double* [n];
    double norma = 0.0;
    for (int i = 0; i < n; i++)
        temp[i] = new double[n];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            temp[i][j] = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                temp[i][j] += mat[i][k] * invmat[k][j];
    output_matrix(temp, n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            temp[i][j] = temp[i][j] * temp[i][j];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            norma += temp[i][j];
    norma = norma - n;
    printf("\n\n");
    printf("norma: %0.50lf", norma);
    for (int i = 0; i < n; i++)
        delete[]temp[i];
    delete[]temp;
}
#pragma endregion
