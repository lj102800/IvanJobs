/**
 * demo for QuickSort
 */

#include <cstdio>

#define MAX 100
int c[MAX];

class QuickSort {
public:
	void quickSort(int p, int r);

private:	
	int partition(int p, int r);
};

int QuickSort::partition(int p, int r) {
	int x = c[r];
	int i = p - 1;

	int j;
	for(j = p; j < r; ++j) {
		if (c[j] <= x) {
			i++;
			int tmp = c[i];
			c[i] = c[j];
			c[j] = tmp;
		}
	}

	int tmp = c[i + 1];
	c[i + 1] = c[r];
	c[r] = tmp;

	return i + 1;
}

void QuickSort::quickSort(int p, int r) {
	if (p < r) {
		int q = this->partition(p, r);
		this->quickSort(p, q - 1);
		this->quickSort(q + 1, r);
	}		
}