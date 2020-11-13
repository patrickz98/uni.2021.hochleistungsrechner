/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

int *
mistake1 (void)
{
    // int buf[] = {1, 1, 2, 3, 4, 5};
    // User static to prevent memory allocation (may lead to memory leaks)
    static int buf[] = {1, 1, 2, 3, 4, 5};
    return buf;
}

int *
mistake2 (void)
{
    // Char is the wrong type here.
    // Use int instead of char to prevent buffer override.
    // int *buf = malloc (sizeof (char) * 4);
    int *buf = malloc(sizeof(int) * 4);

    // Write at wrong position
    // buf[2] = 2;
    buf[1] = 2;
    return buf;
}

int *mistake3(void)
{
    /* In dieser Funktion darf kein Speicher direkt allokiert werden. */

    // Remove unused variable.
    // int mistake2_ = 0;

    // This lead to an override of the function mistake2
    // Call the function instead!
    // int *buf = (int *) &mistake2;
    int *buf = mistake2();

    buf[0] = 3;
    return buf;
}

int *
mistake4(void)
{
    // Four chars are not necessarily an integer.
    // Use sizeof(int) to get system specific size.
    // int *buf = malloc (sizeof (char) * 4);
    int *buf = malloc(sizeof(int));

    // Write integer to the first position.
    // This leads to an exact override of the buffer.
    // buf[4] = 4;
    buf[0] = 4;

    // Don't free the buffer here!
    // This leads to unpredictable results
    // free(buf);

    return buf;
}

int
main(void)
{
    /* Modifizieren Sie diese Zeile nicht! */
    int *p[4] = {
        &mistake1()[1],
        &mistake2()[1],
        mistake3(),
        mistake4()
    };

    printf("1: %d\n", *p[0]);
    printf("2: %d\n", *p[1]);
    printf("3: %d\n", *p[2]);
    printf("4: %d\n", *p[3]);

    /* mhh muss hier noch etwas gefreed werden? */
    /* Fügen sie hier die korrekten aufrufe von free() ein */
    free(p[3]);            /* welcher Pointer war das doch gleich?, TODO: Fixme... :-) */
    free(p[2]);

    // Get first position (pointer) of the array for free operation
    free(p[1] - 1);

    return 0;
}
