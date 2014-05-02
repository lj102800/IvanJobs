//#include <cstdio>
//#include <queue>
//
//using namespace std;
//
//int mat[] = {153, 370, 371, 407};
//
//int main() {
//	int m, n;
//	while (scanf("%d %d", &m, &n) == 2) {
//		queue<int> q;
//		int i;
//		for (i = 0; i < 4; ++i) {
//			if (mat[i] >= m && mat[i] <= n) {
//				q.push(mat[i]);
//			}
//		}
//
//		if (q.empty())
//			printf("no\n");
//		else {
//			while (!q.empty()) {
//				int tmp = q.front();
//				printf("%d", tmp);
//				q.pop();
//				if (q.empty())
//					printf("\n");
//				else
//					printf(" ");
//			}
//
//		}
//	}
//	return 0;
//}