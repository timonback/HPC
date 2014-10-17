#include <stdio.h>

// Definieren Sie ein enum cardd

typedef enum
{
    N = 0x0001,
    S = 0x0010,
    E = 0x0100,
    W = 0x1000,
} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd map[3][3];

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit

void set_dir(int x, int y, cardd dir)
{
    if (0 <= x && x < 3 && 0 <= y && y < 3)
    {
        //Beinhaltet dir eine Himmelsrichtung?
        if (dir & N || dir & E || dir & S || dir & W)
        {
            //Nur Kombinationen, welche nicht gleichzeitig Nord und Sued, als auch Ost und West sind, sind zugelassen.
            if (!(dir & N && dir & S) && !(dir & E && dir & W))
            {
                map[x][y] = dir;
            }
        }
    }
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben

//Aufgabe 2.0
void show_map(void)
{
    int i,j;
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++)
        {
            if(j!= 0)
            {
                printf("   ");
            }
            
            cardd d = map[i][j];

            switch (d)
            {
                case N:
                    printf("%c", 'N');
                    break;

                case E:
                    printf("%c", 'E');
                    break;

                case S:
                    printf("%c", 'S');
                    break;

                case W:
                    printf("%c", 'W');
                    break;

                default:
                    printf("%c", '0');
            }
        }
        
        printf("\n");
    }
}

//Aufgabe 2.1
void show_map2(void)
{
    int i, j;
    //Durch die globale map iterieren
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            cardd d = map[i][j];

            //Abstandshalter (Space)
            if (j != 0)
            {
                //Immer 2 Leerzeichen
                printf("  ");
            }

            //Ausgeben der Himmelsrichtungen
            if (d & N)
            {
                //nord
                printf("N");

                //Zweite Richtung?
                if (d & E)
                {
                    printf("E");
                }
                else if (d & W)
                {
                    printf("W");
                }
            }
            else if (d & S)
            {
                //Sued
                printf("S");

                //Zweite Richtung?
                if (d & E)
                {
                    printf("E");
                }
                else if (d & W)
                {
                    printf("W");
                }
            }
            else if (d & E)
            {
                //ost
                printf("E");
            }
            else if (d & W)
            {   
                //West
                printf("W");
            }
            else
            {
                //Keine Angegeben
                printf(" 0 ");
            }
        }

        printf("\n");
    }
}

int main(void)
{
    // In dieser Funktion darf nichts verändert werden!
    set_dir(0, 1, N);
    set_dir(1, 0, W);
    set_dir(1, 4, W);
    set_dir(1, 2, E);
    set_dir(2, 1, S);

    set_dir(0, 0, N|W);
    set_dir(0, 2, N|E);
    set_dir(0, 2, N|S);
    set_dir(2, 0, S|W);
    set_dir(2, 2, S|E);
    set_dir(2, 2, E|W);

    show_map2();

    return 0;
}

