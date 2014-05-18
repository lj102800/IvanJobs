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


/*=============================================*\
 | 欧拉路径O(E): a path that use all edges each just once! 
 | INIT: adj[][]置为图的邻接表； cnt[a] 为a点的临界点个数；
 | CALL： elpath(0); 注意：不要有自向边。
 | 简单测试了一下，感觉不work？
\*===============================================*/
int adj[V][V], idx[V][V], cnt[V], stk[V], top;
int path(int v) {
	for (int w; cnt[v] > 0 ; v = w) {
		stk[top++] = v;
		w = adj[v][ --cnt[v] ];
		adj[w][idx[w][v]] = adj[w][--cnt[w]];
	}
	return v;
}
void elpath(int b, int n) { // euler path begin from vertex b, vertex 0 ~ n-1
	int i, j;
	for (i = 0; i < n; ++i) 
		for (j = 0; j < cnt[i]; ++j)
			idx[i][adj[i][j]] = j;
	printf("%d", b);
	for (top = 0; path(b) == b && top != 0;) {
		b = stk[--top];
		printf("-%d", b);
	}
	printf("\n");
}

/*===================================================*\
 | Dijkstra Array Impl O(N^2)
 | Dijkstra -- 数组实现（在此基础上可直接改成STL的Queue实现）
 | lowcost[] -- beg到其他点的最近距离
 | path[] -- beg 为根展开的树， 记录父亲节点
\*===================================================*/
#define INF 0x03F3F3F3F
const int N;
int path[N], vis[N];
void Dijkstra(int cost[][N], int lowcost[N], int n, int beg) {
	int i, j, min;
	memset(vis, 0, sizeof(vis));
	vis[beg] = 1;
	for (i = 0; i < n; ++i) {
		lowcost[i] = cost[beg][i]; path[i] = beg;
	}
	lowcost[beg] = 0;
	path[beg] = -1;
	int pre = beg;
	for (i = 1; i < n; ++i) {
		min = INF;
		for (j = 0; j < n; ++j) 
			if (vis[j] == 0 &&
			lowcost[pre] + cost[pre][j] < lowcost[j]) {
				lowcost[j] = lowcost[pre] + cost[pre][j];
				path[j] = pre;
			}
		for (j = 0; j < n; ++j) 
			if (vis[j] == 0 && lowcost[j] < min) {
				min = lowcost[j]; pre = j;
			}
		vis[pre] = 1;
	}
}
