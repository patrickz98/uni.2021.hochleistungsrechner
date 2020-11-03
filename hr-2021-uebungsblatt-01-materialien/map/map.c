#include <stdio.h>

// Definieren Sie ein enum cardd
enum cardd
{
    N,
    W,
    E,
    S,
};

// Definieren Sie ein 3x3-Array Namens map, das Werte vom Typ cardd enthält


// Die Funktion set_dir soll an Position [x, y] den Wert dir in das Array map eintragen
// Überprüfen Sie x und y, um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir(int x, int y, cardd dir)
{

}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map(void)
{

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
