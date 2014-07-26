#include <cstring>

#define MAXNODE 1000
#define SIGMA 30

struct Trie {
	int ch[MAXNODE][SIGMA];
	int val[MAXNODE];
	int sz; // node size
	Trie() { sz = 1; memset(ch[0], 0, sizeof(ch[0]));} // initially, only one root node
	int idx(char c) { return c - 'a';} // index of c
	
	void insert(char* s, int v) {
		int u = 0, n = strlen(s);
		for (int i = 0; i < n; i++) {
			int c = idx(s[i]);
			if (!ch[u][c]) {
				memset(ch[sz], 0, sizeof(ch[sz]));
				val[sz] = 0;
				ch[u][c] = sz++;
			}
			u = ch[u][c];
		}
		val[u] = v;
	}
};
