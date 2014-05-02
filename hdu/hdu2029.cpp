#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

int main() {
	int n;
	string str;
	cin>>n;
	while (n--) {
		cin>>str;
		string rstr = str;
		reverse(rstr.begin(), rstr.end());
		if (str == rstr)
			cout<<"yes"<<endl;
		else 
			cout<<"no"<<endl;		
	}
	
	return 0;
}
