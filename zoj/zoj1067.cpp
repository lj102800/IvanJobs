#include <cstdio>
#include <cmath>
#include <vector>
using namespace std;

struct RGB {
    int r;
    int g;
    int b;
};

vector<RGB> target;
vector<RGB> mapped;

float dis(RGB a, RGB b) {
    return sqrt((float)((a.r - b.r)*(a.r - b.r) + (a.g - b.g)*(a.g - b.g) + (a.b - b.b)*(a.b - b.b)));
}    

int main()
{
  int r, g, b;
  int count = 0;
  while (scanf("%d%d%d", &r, &g, &b) && !(r==-1 && g==-1 && b==-1)){
    count++;
    RGB a;
    a.r = r;
    a.g = g;
    a.b = b;
    if (count <= 16){    
        target.push_back(a);    
    } else  {
        mapped.push_back(a);
    }                  
  }
  
  int size = mapped.size();
  int i, j;
  
  for(i = 0; i < size; ++i){
    int min = 1 << 30;
    int idx = 0;
    for (j = 0; j < 16; ++j){
        int tmp = dis(mapped[i], target[j]);
        if (min > tmp){
            min = tmp;
            idx = j;
        }        
    }
    printf("(%d,%d,%d) maps to (%d,%d,%d)\n", mapped[i].r, mapped[i].g, mapped[i].b, target[idx].r, target[idx].g, target[idx].b);            
  }    
            
  
  return 0;
}
