#include "rungeKutt.h"

int main()
{
	int size = 13;
	double h = 0.1, e = 0.00001, x = 0, u = 1, u1 = 1;
	RungeKutt a(size, h, e, x, u,2);
	
	a.rk4();
	a.printRK4();
	a.printRKwhs4();
	a.printRKwcs4();
	return 0;
}

