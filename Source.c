#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void leastSquaresMethod(double**, const double *, const double *, const int, const int);
void gaussMethod(double**, double *, const int);
void displayPolinom(double*, const int);
void displayArrays(const double *, const double*, const int);
void calculatePolynoms(double *X, double *Y, const int N, const int maxPowerPolynom);
double **allocationMemoryMatrix(const int countRows, const int countColumns);
void inputElements(double *X, double *Y, const int countElements);
void deleteMatrix(double **Matrix, const int countColumns);

int main()
{
	int maxPowerPolynom;
	int countElements;
	printf("Enter countElements : ");
	scanf_s("%d", &countElements);
	printf("Enter maxPowerPolynom : ");
	scanf_s("%d", &maxPowerPolynom);
	double *X = (double*)calloc(countElements, sizeof(double));
	double *Y = (double*)calloc(countElements, sizeof(double));
	inputElements(X, Y, countElements);
	displayArrays(X, Y, countElements);
	calculatePolynoms(X, Y, countElements, maxPowerPolynom);
	free(Y);
	free(X);
	putchar('\n');
	system("pause");
	return 0;
}

void inputElements(double *X, double *Y, const int countElements)
{
	for (int i = 0; i < countElements; ++i) {
		printf("Enter X[%d] : ", i + 1);
		scanf_s("%lf", &X[i]);
		printf("Enter Y[%d] : ", i + 1);
		scanf_s("%lf", &Y[i]);
	}
	return;
}

void calculatePolynoms(double *X,double *Y,const int N,const int maxPowerPolynom)
{
	double **Matrix = allocationMemoryMatrix(maxPowerPolynom + 1, maxPowerPolynom + 2);
	double *multiplierArray = (double*)calloc(maxPowerPolynom + 1, sizeof(double));
	for (int i = 1; i <= maxPowerPolynom; ++i)
	{
		leastSquaresMethod(Matrix, X, Y, N, i);
		gaussMethod(Matrix, multiplierArray, i + 1);
		displayPolinom(multiplierArray, i);
		double f0 = -multiplierArray[0] / multiplierArray[1];
		printf("\nM = %lf\n", f0);
		
	}
	deleteMatrix(Matrix, maxPowerPolynom + 2);
	free(multiplierArray); multiplierArray = NULL;
	return;
}

double **allocationMemoryMatrix(const int countRows,const int countColumns)
{
	double** Matrix = (double**)calloc(countRows, sizeof(double));
	if (Matrix == NULL)
		return NULL;
	for (int i = 0; i < countRows; ++i){
		Matrix[i] = (double*)calloc(countColumns, sizeof(double));
		if (Matrix == NULL)
			return NULL;
	}
	return Matrix;
}

void deleteMatrix(double **Matrix,const int countColumns)
{
	for (int i = 0; i < countColumns; ++i)
		free(Matrix[i]);
	free(Matrix); Matrix = NULL;
	return;
}

void leastSquaresMethod(double** matrix, const double *X, const double *Y, const int N, const int M)
{
	double dSum;
	for (int k = 0; k <= M; ++k)
	{
		for (int j = 0; j <= M; ++j)
		{
			dSum = 0;
			for (int i = 0; i < N; ++i)
			{
				dSum += pow(X[i], j + k);
			}
			matrix[k][j] = dSum;
		}
		dSum = 0;
		for (int i = 0; i < N; ++i)
		{
			dSum += Y[i] * pow(X[i], k);
		}
		matrix[k][M + 1] = dSum;
	}
	return;
}

void swapEduations(double** matrix, int N, int k, int i)
{
	if (k == i)
		return;
	double *Temp = matrix[k];
	matrix[k] = matrix[i];
	matrix[i] = Temp;
	return;
}
void gaussMethod(double** Matrix, double *vector_X, const int N)
{
	double aMax;
	double mainElement;
	for (int k = 0; k < N; ++k)
	{
		aMax = fabs(Matrix[k][k]);
		int pos = k;
		for (int i = k + 1; i < N; ++i)
		{
			if (fabs(Matrix[i][k]) > aMax) {
				aMax = fabs(Matrix[i][k]);
				pos = i;
			}
		}
		swapEduations(Matrix, N, k, pos);
		mainElement = Matrix[k][k];
		for (int j = 0; j < N + 1; ++j)
		{
			Matrix[k][j] = Matrix[k][j] / mainElement;
		}
		for (int i = k + 1; i < N; ++i)
		{
			double M = Matrix[i][k];
			for (int j = k; j < N + 1; ++j)
			{
				Matrix[i][j] -= M * Matrix[k][j];
			}
		}
	}
	int k = N - 1;
	vector_X[k] = Matrix[k][N];
	for (int i = N - 2; i >= 0; --i)
	{
		double dSum = 0;
		for (int j = k; j > i; --j)
		{
			dSum += Matrix[i][j] * vector_X[j];
		}
		vector_X[i] = (Matrix[i][N] - dSum) / Matrix[i][i];
	}
	return;
}

void delimitation(const int w)
{
	for (int i = 0; i < w; ++i)
	{
		putchar('-');
	}
	return;
}
void displayPolinom(double* A, const int M)
{
	const char ch = 'x';
	printf("\nPolinom %d degree :\n\n", M);
	printf("P%d(%c) = %+lf", M, ch, A[0]);
	for (int i = 1; i <= M; ++i)
	{
		printf("%+lf*%c^%d", A[i], ch, i);
	}
	putchar('\n');
	delimitation(80);
	return;
}

void displayArrays(const double* X, const double* Y, const int N)
{
	const int delim = 16 * N;
	printf("\t\tInitial data\n");
	delimitation(delim + 2);
	putchar('\n');
	printf("#: ");
	for (int i = 0; i < N; ++i)
	{
		(i + 1) != N ? printf("%2d\t | ", i + 1) : printf("%2d\t |\n", i + 1);
	}
	delimitation(delim + 2);
	printf("\nX: ");
	for (int i = 0; i < N; ++i)
	{
		(i + 1) != N ? printf("%lf | ", X[i]) : printf("%lf |\n", X[i]);
	}
	delimitation(delim + 2);
	printf("\nY: ");
	for (int i = 0; i < N; ++i)
	{
		(i + 1) != N ? printf("%lf | ", Y[i]) : printf("%lf |\n", Y[i]);
	}
	delimitation(delim + 2);
	printf("\n");
	return;
}