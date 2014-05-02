#include <cstdio>

class InsertionSort {
public:
		void sort(int a[], int size);
		
};

void InsertionSort::sort(int a[], int size)
{
	int i , j;
	for (i = 1; i < size; ++i) {
			int val = a[i];
			j = i  - 1;
			while (j >= 0 && a[j] > val) {
					a[j + 1] = a[j];	
					j--;
			}
			a[j + 1] = val;
	}
}

int main(int argc, char **argv)
{
	InsertionSort e;
	int container[] = {9, 3, 2, 4, 5, 1, 7, 6};
	e.sort(container, 8);
	
	int i;
	for (i = 0; i < 8; i++)
	{
		printf("%d\n", container[i]);
	}
	
	char c;
	scanf("%c", &c);
	return 0;
}
