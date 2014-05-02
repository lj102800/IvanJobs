#include <cstdio>
#include <cstring>
#include <cctype>

char id[100];
char keyIds[][20] = {   "auto",
		        "break",
		    	"case",
			"char",
			"const",
			"continue",
			"default",
			"do",
			"double",
			"else",
			"enum",
			"extern",
			"float",
			"for",
			"goto",
			"if",
			"int","long","register",
			"return", "short", "signed",
			"sizeof", "static", "struct",
			"switch", "typedef", "union",
			"unsigned", "void", "volatile",
			"while", "asm", "cdecl_cs_ds_es",
			"far", "huge", "interrupt", "near",
			"pascal_ss"
		    };


bool isId() {
	int len = strlen(id);
	int i;
	for (i = 0; i < len; ++i) {
		if (!(isdigit(id[i]) || isalpha(id[i]) || id[i] == '_'))
			return false;
	}

	if (isdigit(id[0]))
		return false;
	int keylen = sizeof(keyIds) / 20;
	for (i = 0; i < keylen; i++) {
		if (strcmp(keyIds[i], id) == 0)
			return false;
	}
	return true;
}	



int main() {
	int n;
	scanf("%d", &n);
	getchar();
	while(n--) {
		gets(id);
	        if (isId())
			printf("yes\n");
		else
			printf("no\n");			
	}	
	return 0;
}
