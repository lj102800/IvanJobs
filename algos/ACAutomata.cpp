#include <cstring>
#include <string>
#include <queue>
#include <cstdio>
#include <map>

using namespace std;
const int MAXNODE = 10000;
const int SIGMA_SIZE = 26;
const int MAXS = 200;
struct ACAutomata {
	int ch[MAXNODE][SIGMA_SIZE];
	int f[MAXNODE]; // fail function
	int val[MAXNODE]; // every tail of string has a non-zero val
	int last[MAXNODE]; // output next node of linked list
	int cnt[MAXS];
	int sz;
	
	void init() {
		sz = 1;
		memset(ch[0], 0, sizeof(ch[0]));
		memset(cnt, 0, sizeof(cnt));
		ms.clear();
	}
	
	int idx(char c) { return c - 'a'; }
	
	void insert(char* s, int v){
		int u = 0, n = strlen(s);
		for(int i = 0; i < n; i++) {
			int c = idx(s[i]);
			if (!ch[u][c]) {
				memset(ch[sz], 0, sizeof(ch[sz]));
				val[sz] = 0;
				ch[u][c] = sz++;
			}
			u = ch[u][c];
		}
		val[u] = v;
		ms[string(s)] = v;
	}

	// recursively pring all strings ended with node j
	void print(int j) {
		if (j) {
			cnt[val[j]]++;
			print(last[j]);
		}
	}

	int find(char* T) {
		int n = strlen(T);
		int j = 0;  // curr node idx, initially root node
		for (int i = 0; i < n; i++) { // idx of text string
			int c = idx(T[i]);
			while(j && !ch[j][c]) j = f[j];
			j = ch[j][c];
			if (val[j]) print(j);
			else if (last[j]) print(last[j]); // find
		}
	}
	
	// calculate fail function
	void getFail() {
		queue<int> q;
		f[0] = 0;
		// init queue
		for (int c = 0; c < SIGMA_SIZE; c++) {
			int u = ch[0][c];
			if (u) { f[u] = 0; q.push(u); last[u] = 0;}
		}

		// BFS seq to calculate fail
		while(!q.empty()) {
			int r = q.front(); q.pop();
			for (int c = 0; c < SIGMA_SIZE; c++) {
				int u = ch[r][c];
				if (!u) continue;
				q.push(u);
				int v = f[r];
				while(v && !ch[v][c]) v = f[v];
				f[u] = ch[v][c];
				last[u] = val[f[u]] ? f[u] : last[f[u]];
			}
		}
	}
};

ACAutomata ac;
char text[1000001], P[151][80];
