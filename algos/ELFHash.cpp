#include <iostream>

using namespace std;
const unsigned long M = 100;
int ELFHash(char* key) {
	unsigned long h = 0;
	while(*key) {
		h = (h << 4) + *key++;
		unsigned long g = h & 0Xf0000000L;
		if (g) h^= g >> 24;
		h &= ~g;
	}
	return h % M;
}
int main() {
	char line[1024];
	while(cin>>line) {
		cout<<ELFHash(line)<<endl;
	}
	return 0;
}
