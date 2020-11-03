#include <stdio.h>

// Nach der Korrektur der Funktion call_by_reference löschen
int TODO;

void basic_pointer(int x)
{
	int* adresse_von_x;

	adresse_von_x = &x;

	printf("Der Wert von x ist: %d\n", 0 /* TODO */);
	printf("Die Adresse von x ist %p\n", NULL /* TODO */);
	printf("Adresse von x mittels adresse_von_x %p\n", NULL /* TODO */);
	printf("Wert von x mittels adresse_von_x: %d\n", 0 /* TODO */);
}

void basic_pointer2(int x)
{
	int* adresse_von_x = NULL /* TODO */;
	// Eine andere Variable y erhaelt den Wert von x
	int y = 0 /* TODO */;

	printf("Der Wert von y ist %d\n", 0 /* TODO */);

	// Zuweisung über die Adresse
	x = 10;
	y = *adresse_von_x;

	printf("Der Wert von y ist %d\n", 0 /* TODO */);
}

void basic_pointer_changeValue(int x)
{
	int* adresse_von_x = NULL /* TODO */;

	// Ändern Sie den Wert von x zu 10
	TODO = 10;
	printf("x = %d\n", 0 /* TODO */);

	// Ändern Sie den Wert von x über seine ADRESSE
	TODO = 20;
	printf("x = %d\n", x);
}

void call_by_reference(int* x)
{
	// Ändern Sie den Wert, der an der Adresse steht, die im Wert x gespeichert ist
	TODO = 200;
}

int main(void)
{
	int x = 5;

	basic_pointer(x);
	basic_pointer2(x);
	basic_pointer_changeValue(x);

	printf("Der Wert von x vor der Funktion call_by_reference: %d\n",x);
	call_by_reference(&x);
	printf("Der Wert von x nach der Funktion call_by_reference: %d\n",x);

	return 0;
}
