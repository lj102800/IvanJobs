#include <cstdio>
#include <cstring>

struct node {
	node(){
		cnt = 0;
		int i;
		for (i = 0; i < 26; ++i)
			next[i] = NULL;	
	}
	int cnt;
	node* next[26];
};


#define MAX 1000
node db[MAX];

int curr = 0;
node* root = &db[curr++];


int main() {
	int n;
	char str[100];
	while (scanf("%d", &n) && n != 0){
		int i, j;
		curr = 1;
		for (i = 0; i < n; ++i){
			scanf("%s", str);
			// here we got str, and we need insert it into the trie tree.
			node* p = root;
			for (j = 0; j < strlen(str); ++j){
				int ith = str[j] - 'a';
				if (p->next[ith] != NULL){ //ok, we find it already exists, so add p->cnt++
					p = p->next[ith];
					p->cnt++;
				} else {
					p->next[ith] = &db[curr++];
					p = p->next[ith];
					p->cnt++;
				}
			}
		}
		scanf("%s", str); // this is prefix
		//we will do searching in trie for prefix str
		node* needle = root;
		int res;
		for (j = 0; j < strlen(str); ++j){
			int slot = str[j] - 'a';
			if (needle->next[slot] != NULL){
				needle = needle->next[slot];
				if ( j == strlen(str) - 1){
					res = needle->cnt;
					break;
				}
			} else {
				res = -1;
				break;
			}
		}

		printf("res:%d\n", res);
	}
	return 0;
}