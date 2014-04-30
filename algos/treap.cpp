#include <iostream>

typedef struct _Node {
	_Node* childs[2];
	int p;	// priority for Heap
	int k;  // key for BST
	int cmp(int x) const {
		if (x == k) return -1;
		return x < k ? 0 : 1;
	}
}Node;

void rotate(Node*& node, int d) {
	Node* k = node->childs[d^1];
	node->childs[d^1] = k->childs[d];
	k->childs[d] = node;
	node = k;
}

void insert(Node*& o, int x) {
	if (o == NULL) {
		o = new Node();
		o->childs[0] = o->childs[1] = NULL;
		o->k = x;
		o->p = rand();
	} else {
		int d = o->cmp(x);
		insert(o->childs[d], x);
		if (o->childs[d]->p > o->p) rotate(o, d^1);
	}
}

void remove(Node*& o, int x) {
	int d = o->cmp(x);
	if (d == -1) {
		if (o->childs[0] == NULL) o = o->childs[1];
		else if (o->childs[1] == NULL) o = o->childs[0];
		else {
			int d2 = (o->childs[0]->p > o->childs[1]->p ? 1 : 0);
			rotate(o, d2);
			remove(o->childs[d2], x);
		}
	} else remove(o->childs[d], x);
}

int find(Node* o, int x) {
	while(o != NULL) {
		int d = o->cmp(x);
		if (d == -1) return 1; // true
		else o = o->childs[d];
	}
	return 0; // false
}


