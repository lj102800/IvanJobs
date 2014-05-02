/**
 * this is about Union Find algorithm.
 */
#include <cstdio>

#define MAX 100
int mat[MAX];

class QuickFindUF {
public:
	QuickFindUF(int N); // initialize union-find data structure with N objects (0 - N-1)
	void uni(int p, int q); // add connection between p and q
	bool connected(int p, int q); // are p and q in the same component?

	int find(int p); // component identifier for p (0 to N - 1)

	int count(); // number of components

private:
	int n;

};

QuickFindUF::QuickFindUF(int N) {
	this->n = N;
	int i;
	for (i = 0; i < N; ++i)
		mat[i] = i;
}

bool QuickFindUF::connected(int p, int q) {
	if (mat[p] == mat[q]) return true;
	else return false;
}

void QuickFindUF::uni(int p, int q) {
	int p_cid = mat[p];
	int q_cid = mat[q];

	if (!this->connected(p, q)) {
		int i;
		for (i = 0; i < this->n; ++i) {
			if (mat[i] == p_cid) mat[i] = q_cid;
		}
	}
}

int QuickFindUF::find(int p) { 
	return mat[p];
}

int QuickFindUF::count() {
	return this->n;
}


int main() {

	int N;
	scanf("%d", &N);
	QuickFindUF uf(N);

	int p, q;
	while (scanf("%d", &p) != EOF) {
		q = scanf("%d", &q);

		if (!uf.connected(p, q)) {
			uf.uni(p, q);
			printf("%d %d\n", p, q);
		}
	}

	
	return 0;
}