#include <iostream>
#include <string>

using namespace std;

int main() {
	string str;
	while (cin>>str) {
		int i;
		char max = str[0];
		for (i = 1; i < str.size(); ++i) {
			if (str[i] > max)
				max = str[i];
		}
		for (i = 0; i < str.size(); ++i) {
			if (str[i] == max) {
				cout<<str[i]<<"(max)";
			} else {
				cout<<str[i];
			}
		}
		cout<<endl;
	}	
	return 0;
}
