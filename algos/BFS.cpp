//#include <iostream>
//#include <vector>
//#include <queue>
//
//using namespace std;
//
//
//#define MAX 100
//
//bool matrix[MAX][MAX];
//char color[MAX]; // 'W' white, 'G' gray, 'B' black
//int d[MAX];
//int p[MAX]; //root's parent is 0
//// index starting from 1
//
//class BFS {
//public :
//	void bfs();
//
//};
//
//// n stands number of elements
//// s stands for start element
//void BFS::bfs(int n, int s) {
//	int i;
//	for (i = 1; i <= n; ++i) {
//		if (i != s) {
//			color[i] = 'W';
//			p[i] = 0;
//			d[i] = MAX;
//		}
//	}
//
//	color[s] = 'G';
//	d[s] = 0;
//	p[s] = 0;
//
//	queue<int> q;
//	q.push(s);
//
//	while (!q.empty()) {
//		int shit = q.front();
//		int i;
//		for ( i = 1; i < n; ++i) {
//			if ()
//		}
//		
//	}
//}