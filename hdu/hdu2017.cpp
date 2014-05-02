#include <iostream>
#include <string>
#include <cctype>

using namespace std;

int main() {
	int n;
	string str;
	cin>>n;
	while (n--) {
		cin>>str;
		int i;
		int count = 0;
		for (i = 0; i < str.size(); ++i) {
			if (isdigit(str[i]))
				count++;
		}		
		cout<<count<<endl;	
	}
	return 0;
}
