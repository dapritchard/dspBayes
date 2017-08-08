#include <iostream>

int main() {

    int x[4];

    int n = 1 << 4;

    for (int i = 0; i < n; ++i) {

	for (int k = 0; k < 4; ++k) {

	    x[k] = (i & (1<<k)) ? 1 : 0;

	    std::cout << x[k] <<  "  ";
	}
	std::cout << "\n";

    }

    return 0;
}
