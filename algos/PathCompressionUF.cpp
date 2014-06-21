// recursive version
int find(int x) {
	if (x != parent[x]) {
		parent[x] = find(parent[x]);
	}
	return parent[x];
}

// iteration version
int find(int x) {
	int k, j, r;
	r = x;
	while(r != parent[r]) 
		r = parent[r];
	k = x;
	while(k != r) {
		j = parent[k];
		parent[k] = r;
		k = j;
	}
	return r;
}
