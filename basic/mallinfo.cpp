
#if defined(__unix__) ||  defined(__linux__)

#include <malloc.h>
#include <stdio.h>

void display_mallinfo()
{
    struct mallinfo2 mi;

    mi = mallinfo2();

    printf("Usage statistics heap:\n");

    printf("Total non-mmapped bytes (arena):       %ld\n", mi.arena);
    printf("# of free chunks (ordblks):            %ld\n", mi.ordblks);
    printf("# of free fastbin blocks (smblks):     %ld\n", mi.smblks);
    printf("# of mapped regions (hblks):           %ld\n", mi.hblks);
    printf("Bytes in mapped regions (hblkhd):      %ld\n", mi.hblkhd);
    printf("Max. total allocated space (usmblks):  %ld\n", mi.usmblks);
    printf("Free bytes held in fastbins (fsmblks): %ld\n", mi.fsmblks);
    printf("Total allocated space (uordblks):      %ld\n", mi.uordblks);
    printf("Total free space (fordblks):           %ld\n", mi.fordblks);
    printf("Topmost releasable block (keepcost):   %ld\n", mi.keepcost);
}

#else 

void display_mallinfo(){}

#endif 