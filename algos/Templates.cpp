// Graph Theory

/* ===================================================== *\
 | DAG DFS TAG
 | edge[i][j] if i->j has edge? 1->true, 0->false
 | V -> maxV, a little bit bigger than number of vertexs
 | pre and post ? it seems pre is used for tag dfs tree'v level of curr node.
 | tag ?
 | we can use this operation to find Tree Edge, Down Edge, Back Edge and Cross Edge! 
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

/*===========================================================*\
 | find one bridge in a undirected graph
 | INIT: edge, vis, pre, anc, bridge all init to 0
 | CALL: dfs(0, -1, 1, n);
\*===========================================================*/
int bridge, edge[V][V], anc[V], pre[V], vis[V];
void dfs(int cur, int father, int dep, int n) {
	// vertex: 0 ~ n-1
	if (bridge) return ;
	vis[cur] = 1; pre[cur] = anc[cur] = dep;
	for (int i = 0; i < n; ++i) if (edge[cur][i]) {
		if (i != father && 1 == vis[i]) {
			if (pre[i] < anc[cur])
				anc[cur] = pre[i]; // back edge
		}
		if (0 == vis[i]) { // tree edge
			dfs(i, cur, dep + 1, n);
			if (bridge) return ;
			if (anc[i] < anc[cur]) anc[cur] = anc[i];
			if (anc[i] > pre[cur]) {
				printf("find a bridge %d->%d\n", cur, i); 
				bridge = 1; return ;
			}
		}
	}
	vis[cur] = 2;
}

