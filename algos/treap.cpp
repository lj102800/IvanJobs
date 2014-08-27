#include <cstdio>
#include <cstdlib>
#define REP(i, n) for(int _n = n, i = 0; i < _n; i++)
#define FOR(i, a, b) for(int i = (a), _b = (b); i <= _b; i++)
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))
#define Abs(x) ((x) > 0 ? (x) :(-(x)))
#define L(fmt, ...) do {if(true) printf(fmt"\n", ##__VA_ARGS__);} while(false)

#define MAXN 10000 // max number of nodes
struct Node {
	Node* ch[2]; // ch is short for children
	int r;	// priority
	int v;	// value
	Node() {
		ch[0] = ch[1] = NULL;
		r = v = 0;
	}
	int cmp(int x) const {
		if (x == v) return -1;
		return x < v ? 0 : 1;
	}
};

struct Treap {
	Node DB[MAXN];
	int pos; // current allocate position
	Node* root; // root of this Treap
	Treap() {
		pos = 0; root = NULL;
	}
	// d: 0 left rotate, 1 right rotate
	void rotate(Node*& o, int d) {
		Node* k = o->ch[d^1]; o->ch[d^1] = k->ch[d]; k->ch[d] = o; o = k;
	}

	void insert( Node*& o , int x) {
		if (o == NULL) {
			o = &DB[pos]; pos++;
			o->v = x; o->r = rand();
		} else {
			int d = o->cmp(x);
			insert(o->ch[d], x);
			if (o->ch[d]->r > o->r) rotate(o, d^1);
		}
	}

	void remove( Node*& o, int x ) {
		int d = o->cmp(x);
		if (d == -1) {
			if (o->ch[0] == NULL) o = o->ch[1];
			else if (o->ch[1] == NULL) o = o->ch[0];
			else {
				int d2 = (o->ch[0]->r > o->ch[1]->r ? 1 : 0);
				rotate(o, d2);
				remove(o->ch[d2], x);
			}
		} else remove(o->ch[d], x);
	}

	bool has( Node* o,int x) {
		while(o != NULL) {
			int d = o->cmp(x);
			if (d == -1) return true; // existed
			else o = o->ch[d];
		}
		return false; // not existed;
	}

	int size() {
		return pos;
	}
};


Treap T;
int main() {
	int v;
	while(scanf("%d", &v) != EOF) {
		if (T.has(T.root, v) == false) T.insert(T.root, v);
	}
	printf("Treap size:%d\n", T.size());
	return 0;
}

