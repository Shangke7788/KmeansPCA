#include <stdio.h>
#include <algorithm>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
using namespace std;
typedef long long lint;

struct Node {
	int id;
	double x, y;
	bool operator < (const Node & o) const {
		if (id != o.id) {
			return id < o.id;
		} else {
			return fabs(x * y) < fabs(o.x * o.y);
		}
	}
	Node() {
	}
	Node(double x, double y) {
		id = rand(), this->x = x, this->y = y;
	}
};
const double PI = cos(-1.0);
const int K = 4;
const int N = 110;
Node p[N];

double a[4][3] = {
	{10, 5, 3},
	{12, 20, 7}, 
	{0, -4, 6},
	{-20, 10, 8}
};

int main() {
	srand(7);
	int pn = 0;
	for (int i = 0; i < K; i++) {
		int n = 25 + i + 1;
		for (int j = 0; j < n; j++) {
			double t = a[i][2] * rand() / RAND_MAX;
			double d = 2.0 * PI * rand() / RAND_MAX;
			p[pn++] = Node(t * cos(d) + a[i][0], t * sin(d) + a[i][1]);
		}
	}
	printf("%d %d %d\n", pn, 2, K);
	sort(p, p + pn);
	for (int i = 0; i < pn; i++) {
		printf("%.3lf %.3lf\n", p[i].x, p[i].y);
	}
	printf("// answer: 4 clusters with number of 26, 27, 28, 29\n");
	return 0;
}
