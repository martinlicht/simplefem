
#if false and ( defined(__unix__) ||  defined(__linux__) )

#include <malloc.h>
#include <stdio.h>

void display_mallinfo()
{
    mallinfo2 mi;

    mi = mallinfo2();

    printf("Heap usage statistics:\n");

    printf("Total non-mmapped bytes (arena):       %zu\t\t%.3f GiB\n", (size_t)mi.arena,    mi.arena    / 1073741824. );
    printf("# of free chunks (ordblks):            %zu\n", (size_t)mi.ordblks );
    printf("# of free fastbin blocks (smblks):     %zu\n", (size_t)mi.smblks );
    printf("# of mapped regions (hblks):           %zu\n", (size_t)mi.hblks );
    printf("Bytes in mapped regions (hblkhd):      %zu\t\t%.3f GiB\n", (size_t)mi.hblkhd,   mi.hblkhd   / 1073741824. );
    printf("Max. total allocated space (usmblks):  %zu\t\t%.3f GiB\n", (size_t)mi.usmblks,  mi.usmblks  / 1073741824. );
    printf("Free bytes held in fastbins (fsmblks): %zu\t\t%.3f GiB\n", (size_t)mi.fsmblks,  mi.fsmblks  / 1073741824. );
    printf("Total allocated space (uordblks):      %zu\t\t%.3f GiB\n", (size_t)mi.uordblks, mi.uordblks / 1073741824. );
    printf("Total free space (fordblks):           %zu\t\t%.3f GiB\n", (size_t)mi.fordblks, mi.fordblks / 1073741824. );
    printf("Topmost releasable block (keepcost):   %zu\t\t%.3f GiB\n", (size_t)mi.keepcost, mi.keepcost / 1073741824. );

    malloc_stats();
}

#else 

#include <stdio.h>
void display_mallinfo(){
    printf("Heap usage statistics: requested but not available in this implementation.\n");
}

#endif 
