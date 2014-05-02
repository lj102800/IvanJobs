//**
// * demo for merge sort
// */
//
//#include <cstdio>
//#include <cstdlib>
//
//#define MAX 100
//int c[MAX];
//
//#define M 1000000 // sentinel val used to stand for oo
//
//class MergeSort {
//public:
//	void doIt(char *fileIn, char *fileOut);
//private:
//	void merge(int container[], int p, int q, int r);
//	void mergeSort(int container[] , int p, int r);
//};
//
//void MergeSort::mergeSort(int c[], int p, int r) {
//	if (p < r) {
//		int q = (p + r) / 2;
//		mergeSort(c, p, q);
//		mergeSort(c, q + 1, r);
//		merge(c, p, q, r);
//	}
//}
//
//void MergeSort::merge(int c[], int p, int q, int r) {
//	int n1 = q - p + 1;
//	int n2 = r - q;
//
//	int *L = (int *)malloc((n1 + 1) * sizeof(int));
//	int *R = (int *)malloc((n2 + 1) * sizeof(int));
//
//	int i;
//	for (i = p; i <= q; ++i)
//		L[i - p] = c[i];
//	for (i = q + 1; i <= r; ++i)
//		R[i - q - 1] = c[i];
//
//	L[n1] = M;
//	R[n2] = M;
//
//	int k, j;
//	i = 0;
//	j = 0;
//	for (k = p; k <= r; ++k) {
//		if (L[i] > R[j]) {
//			c[k] = R[j];
//			j++;
//		} else {
//			c[k] = L[i];
//			i++;
//		}
//	}
//
//	free(L);
//	free(R);
//}
//void MergeSort::doIt(char *fileIn, char *fileOut) {
//	FILE *pFin = fopen(fileIn, "r");
//	FILE *pFout = fopen(fileOut, "w");
//
//	int n; 
//	fscanf(pFin, "%d", &n);
//	
//	int cIdx;
//	for (cIdx = 1; cIdx <= n; ++cIdx) {
//		int s;
//		fscanf(pFin, "%d", &s);
//		int i;
//		for (i = 0; i < s; ++i) {
//			fscanf(pFin, "%d", &c[i]);
//		}
//
//		// do merge sort here
//		// here we got s, c 
//		mergeSort(c, 0, s - 1);
//
//		fprintf(pFout, "Case# %d:\n", cIdx);
//
//		for( i = 0; i < s; ++i) {
//			fprintf(pFout, "%d", c[i]);
//			if (i != s - 1) {
//				fprintf(pFout, " ");
//			}
//		}
//		fprintf(pFout, "\n");
//
//	}
//	fclose(pFin);
//	fclose(pFout);
//}
//
