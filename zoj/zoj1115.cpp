#include <iostream>
#include <string>
#include <cstdio>

using namespace std;
int make(string s){
    int len = s.length();
    int sum = 0;
    int i;
    for(i = 0; i < len; ++i){
        sum += (s[i] - '0');            
    }
    int n;
    while(!(sum > 0 && sum < 10)){
        n = sum;
        sum = 0;
        while(n > 0){
         sum += n % 10;
         n /= 10;   
        }    
    }
    return sum;            
    
}    
bool zero(string s){
    int len = s.length();
    int sum = 0;
    int i;
    for (i=0; i < len; ++i)
        sum += s[i] - '0';
    if (sum == 0)
        return true;
    else
        return false;   
}    

int main() {
    string s;
    while(cin>>s && !zero(s)){
         printf("%d\n", make(s));       
    }        
    return 0;   
}    