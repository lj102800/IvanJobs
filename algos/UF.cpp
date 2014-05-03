/**
 * this is about Union Find algorithm.
 */
#include <cstdio>

#define MAX 100
int P[MAX];

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
		P[i] = i;
}

bool QuickFindUF::connected(int p, int q) {
	if (P[p] == P[q]) return true;
	else return false;
}

void QuickFindUF::uni(int p, int q) {
	int p_cid = P[p];
	int q_cid = P[q];

	if (!this->connected(p, q)) {
		int i;
		for (i = 0; i < this->n; ++i) {
			if (P[i] == p_cid) P[i] = q_cid;
		}
	}
}

int QuickFindUF::find(int p) { 
	return P[p];
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
