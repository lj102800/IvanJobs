#include <cstdio>

char word[100];

void diplay(int n){
}    

int main() {
    int C, c;
    scanf("%d", &C);
    while (C--){
        scanf("%d", &c);
        getchar();
        while (c--){
            char ch = getchar();
            int i = 0;
            while(ch != '\n'){
                if (ch != ' '){
                    word[i++] = ch;
                } else {
                    int j;
                    for(j = i - 1; j >= 0; j--){
                        printf("%c", word[j]);
                    }
                    printf(" ");
                    i = 0;        
                }
                ch = getchar();        
            }
            int j;
            for(j = i - 1; j >= 0; j--){
                printf("%c", word[j]);
            }
            printf("\n");        
        }
        if(C != 0){
            printf("\n");
            
        }    
               
    }    
     
    return 0;   
}    
