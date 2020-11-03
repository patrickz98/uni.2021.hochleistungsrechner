#include <stdio.h>

// Definieren Sie ein enum cardd
typedef enum : unsigned int
{
    N = 0x0001u,
    W = 0x0010u,
    E = 0x0100u,
    S = 0x1000u,
} cardd;

// Definieren Sie ein 3x3-Array Namens map, das Werte vom Typ cardd enthält
char *map[3][3] = {
        {"0", " ", "0"},
        {" ", "0", " "},
        {"0", " ", "0"},
};

// Die Funktion set_dir soll an Position [x, y] den Wert dir in das Array map eintragen
// Überprüfen Sie x und y, um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir(int x, int y, cardd dir)
{
    printf("x=%d, y=%d, dir=%4x\n", x, y, dir);

    // Check for invalid coordinates
    // 0: false, >=1: true
    int invalid = 0;
    invalid += ((x > 2) || (x < 0));
    invalid += ((y > 2) || (y < 0));

    if (invalid)
    {
        return;
    }

    switch (dir)
    {
        case W:
            map[ y ][ x ] = "W";
            break;

        case N:
            map[ y ][ x ] = "N";
            break;

        case E:
            map[ y ][ x ] = "E";
            break;

        case S:
            map[ y ][ x ] = "S";
            break;

        default:
            // Suppress warning
            break;
    }
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map(void)
{
    for (int inx = 0; inx < 3; inx++)
    {
        for (int iny = 0; iny < 3; iny++)
        {
            printf("%c", *map[inx][iny]);
        }

        printf("\n");
    }
}

// In dieser Funktion darf nichts verändert werden!
int main(void)
{
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	show_map();

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);
	set_dir(1, 3, N|S|E);
	set_dir(1, 1, N|S|E|W);

	show_map();

	return 0;
}
