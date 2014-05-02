/**
 *
 */

#include <cstdio>

#define MAX 100
int con[MAX];

class InsertionSort {
public:
    void doIt(char *fileIn, char *fileOut);
};

void InsertionSort::doIt(char *fileIn, char *fileOut) {
    FILE *pFin = fopen(fileIn, "r");
    FILE *pFout = fopen(fileOut, "w");

    int n; // number of cases
    fscanf(pFin, "%d", &n);
    int c_idx; //case index, starting from 1
    for (c_idx = 1; c_idx <= n; ++c_idx) {
        int s; // number of elements
        fscanf(pFin, "%d", &s);
        int i; // looping idx
        for (i = 0; i < s; ++i) {
            fscanf(pFin, "%d", &con[i]);
        }

        // do insertion sort here
        int j;
        for (i = 1; i < s; ++i) {
            // insert to con[0..i-1]
            j = i - 1;
            int val = con[i];
            while (con[j] > val && j >= 0) {
                con[j + 1] = con[j];
                j--;
            }
            con[j] = val;
        }

        fprintf(pFout, "Case# %d:\n", c_idx);
        for (i = 0; i < s; ++i) {
            fprintf(pFout, "%d", con[i]);
            if (i != s - 1) {
                fprintf(pFout, " ");
            }
        }
        fprintf(pFout, "\n");
    }

    fclose(pFin);
    fclose(pFout);
}

int main()
{
    InsertionSort i;
    i.doIt("in.txt", "out.txt");
    return 0;
}
