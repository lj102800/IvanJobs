/**
 * use array to store a heap.
 * A[1...n], A[1] is root, A[i]'s parent is A[i/2], left sibling is A[2*i], right sibling is A[ 2*i + 1].
 * A[heapsize/2 + 1 ...] all leavies.
 */

#include <cstdio>

#define MAX 100
int c[MAX];

#define L(i) 2 * i
#define R(i) 2 * i + 1
#define P(i) i / 2


class HeapSort {
public:
	void doIt();

private:
	int size; // size of heap
	void buildHeap(int length);
	void maxHeapify(int i);// surpose L(i) and R(i) are all heap already.
	void heapSort();
};

void HeapSort::doIt() {
	this->heapSort();
}
// here we got a heap stored in an array, how can to sort these elements?
// 
void HeapSort::heapSort() {
	while (size > 1) {
		int tmp = c[1];
		c[1] = c[size];
		c[size] = tmp;

		size--;
		maxHeapify(1);
	}
}

/**
 * here we got an array with elements from index 1 to length, 
 * we will not use index 0.
 */
void HeapSort::buildHeap(int length) {
	this->size = length;
	int i;
	for (i = length / 2 + 1; i >=1 ; i--) this->maxHeapify(i);
}

void HeapSort::maxHeapify(int i) {
	int largest = i;
	int l = L(i);
	int r = R(i);

	if (l <= size && c[i] < c[l]) {
		largest = l;
	}
	if (r <= size && c[largest] < c[r]) {
		largest = r;
	}

	if (largest != i) {
		int tmp = c[largest];
		c[largest] = c[i];
		c[i] = tmp;
		maxHeapify(largest);
	}
}
