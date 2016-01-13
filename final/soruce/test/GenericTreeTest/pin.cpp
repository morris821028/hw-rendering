#include <bits/stdc++.h>
using namespace std;

int size = 0;
void dfs(int l, int r) {
	size++;
	if (l > r) {
		printf("1 0 ");
		return ;
	}
	printf("1 ");
	int base = l;
	int gap = (r - l + 1) / (rand()%16+1);
	if (gap < 1)	gap = 1;
	for (int i = 0; l <= r; i++) {
		base += rand()%gap;
		base = min(base, r);
		dfs(l, base-1);
		base++;
		l = base;
	}
	printf("0 ");
}
int main() {
	freopen("in.txt", "w", stdout);
    srand(time(NULL));
    int testcase = 1;
    while (testcase--) {
    	int n = 100000;
    	size = 0;
		dfs(0, n);
		puts("0");
//		printf("size %d\n", size);
	}
	return 0;
}

