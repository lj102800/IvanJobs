#include <cstdio>
#include <map>
#include <string>


using namespace std;

#define MAX 100
char line[MAX];
char en[50];
char fat[50];

map<string, string> m; // fat -> en

int main(){
	while (gets(line) && line[0] != '\0'){
		sscanf(line, "%s %s", en, fat);
		string fatstr(fat);
		string enstr(en);
		m.insert(pair<string, string>(fatstr, enstr));
	}

	while (gets(line) && line[0] != '\0'){
		string needle(line);
		map<string, string>::iterator it;
		it = m.find(needle);
		if (it != m.end()){
			string tmp = it->second;
			printf("%s\n", tmp.c_str());
		} else {
			printf("eh\n");
		}
	}

	return 0;
}
