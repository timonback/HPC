#include <stdio.h>

// Definieren Sie ein enum cardd

typedef enum
{
    N = 1<<0,
    S = 1<<1,
    E = 1<<2,
    W = 1<<3,
} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd map[3][3];

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit

void set_dir(int x, int y, cardd dir)
{
    if (0 <= x && x < sizeof(map)/sizeof(map[0]) && 0 <= y && y < sizeof(map[0])/sizeof(map[0][0]))
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
    for(i = 0; i < sizeof(map)/sizeof(map[0]); i++)
    {
        for(j = 0; j < sizeof(map[0])/sizeof(map[0][0]); j++)
        {            
            cardd d = map[i][j];

            switch (d)
            {
                case N:
                    printf(" N ");
                    break;
                case N|E:
                    printf(" NE");
                    break;
                case N|W:
                    printf("NW ");
                    break;
                case S:
                    printf(" S ");
                    break;
                case S|E:
                    printf(" SE");
                    break;
                case S|W:
                    printf("SW ");
                    break;
                case E:
                    printf("  E");
                    break;
                case W:
                    printf("W  ");
                    break;
                default:
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

    show_map();

    return 0;
}

