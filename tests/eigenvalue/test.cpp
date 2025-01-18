#include <stdio.h>
int main()
{
	long double x = 1.0;
#if 0
	/* mha: the 'll' modifier is for integral types */
	printf("%llf %llf %llf %llf\n", x, x, x, x);
#endif
	printf("%Lf %Lf %Lf %Lf\n", x, x, x, x);
	/* mha: removed silly %LF specifier */
	printf("Le %Le\n", x);
	printf("LE %LE\n", x);
	printf("Lf %Lf\n", x);
	printf("Lg %Lg\n", x);
	printf("LG %LG\n", x);
	return 0;
}
