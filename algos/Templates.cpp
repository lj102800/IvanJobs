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


/*==================================================*\
 | 无向图连通度（割）
 | INIT： edge[][]邻接矩阵；vis[], pre[], anc[], deg[] 置为0
 | CALL: dfs(0, -1, 1, n);
 | k = deg[0], deg[i] + 1 (i = 1...n-1) 为删除该节点后得到的连通图个数。
 | 注意： 0 作为根比较特殊！
 | 没搞明白这个函数的意思：（
\*===================================================*/
int edge[V][V], anc[V], pre[V], vis[V], deg[V];
void dfs(int cur, int father, int dep, int n) {
	// vertex: 0 ~ n-1
	int cnt = 0;
	vis[cur] = 1; pre[cur] = anc[cur] = dep;
	for (int i = 0; i < n; ++i) if (dege[cur][i]) {
		if (i != father && 1 == vis[i]) {
			if (pre[i] < anc[cur])
				anc[cur] = pre[i]; //back edge
		}
		if (0 == vis[i]) { // tree edge
			dfs(i, cur, dep + 1, n);
			++cnt; // branches number
			if (anc[i] < anc[cur]) anc[cur] = anc[i];
			if ((curr == 0 && cnt > 1) ||
				(cnt != 0 && anc[i] >= pre[cur]))
					++deg[cur];
		}
	}
	vis[cur] = 2;
}

/*=======================================================*\
 | 最大团问题 DP + DFS
 | INIT： g[][] 邻接矩阵
 | CALL： res = clique(n);
\*=======================================================*/

int g[V][V], dp[V], stk[V][V], mx;
int dfs(int n, int ns, int dep) {
	if (0 == ns) {
		if (dep > mx) mx = dep;
		return 1;
	}
	int i, j, k, p, cnt;
	for (i = 0; i < ns; ++i) {
		k = stk[dep][i]; cnt = 0;
		if (dep + n - k <= mx) return 0;
		if (dep + dp[k] <= mx) return 0;
		for (j = i + 1; j < ns; ++j) {
			p = stk[dep][j];
			if (g[k][p]) stk[dep + 1][cnt++] = p;
		}
		dfs(n, cnt, dep + 1);
	}
	return 1;
}
int clique(int n) {
	int i, j, ns;
	for (mx = 0, i = n - 1; i >= 0; i--) {
		// vertex: 0 ~ n-1
		for (ns = 0, j = i + 1; j < n; ++j)
			if (g[i][j]) stk[1][ns++] = j;
		dfs(n, ns, 1); dp[i] = mx;
	}
	return mx;
}

