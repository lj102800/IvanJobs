// Graph Theory

/* ===================================================== *\
 | DAG DFS 
\* ===================================================== */
int edge[V][V], pre[V], post[V], tag;
void dfstag(int cur, int n) {
	// vertex: 0 ~ n-1
	pre[cur] = ++tag;
	for (int i = 0; i < n; ++i) if (edge[cur][i]) {
		if (0 == pre[i]) {
			printf("Tree Edge!\n");
			dfstag(i, n);
		} else {
			if (0 == post[i]) printf("Back Edge!\n");
			else if (pre[i] > pre[cur])
				printf("Down Edge!\n");
			else printf("Cross Edge!\n");
		}
	}
	post[cur] = ++tag;
}
