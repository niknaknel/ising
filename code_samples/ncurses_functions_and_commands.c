
//  ncurses: for on-screen display of lattices
//  Commands you will need to call ncurses

#include <ncurses.h>
//#include "curses_koen.c"    

// In the calling program:
//initCurses();
//printCurses(latticename,N,L,10000);
//endwin();

/*********************************************************************/
/* 
   Funksies om ncurses vir die Isingmodel te gebruik. 
   (gebaseer op kode van Etienne Koen)

1. Roep die funksie initCurses heel in die begin van jou Ising main
   program nadat die veranderlikes verklaar is

2. Roep die funksie printCurses binne die timestep loop wat die sweep
   omsluit, na' sweep. Die argumente is:
   lattice  - die skikking met die rooster self 
   nsize    - die aantal elemente in die skikking
   lsize    - die sylengte L van die rooster
   nsleep   - aantal mikrosekondes wat die beeld vertraag moet word

3. Om die ncurses libraries te link, moet jy (benewens die include
   <ncurses.h> wat hieronder gedoen word) in jou kompilering ook die
   opsie -lncurses invoeg, maw byvoorbeeld 

     gcc -(allerhande opsies) dieprogram.c -lm -lncurses;

4. Die grootte van die roosters wat op die skerm vertoon kan word is
   afhanklik van die font wat die shell gebruik. Om groter roosters te
   hanteer, kan jy die shell font kleiner stel deur hulle onder die
   Settings menu bo-aan die shell te verander. Roostergroottes sal
   wissel tussen 10x10 tot 130x130 (vir baie klein fonts).

5. Om die vertoning van die rooster te verlangsaam, kan jy of elke
   10de (of 1000de) sweep vertoon; andersins kan jy ook die funksie
   usleep(getal); na' die roep van printCurses invoeg (waar getal 'n
   integer is wat mikrosekondes tel).

6. Die program kan op enige stadium met control-C gestop word.

*/

/*********************************************************************/
#include<stdio.h>
#include<ncurses.h>
#include<unistd.h>   // needed for usleep
#define TIME_OUT 5

void initCurses(){
	initscr();
	cbreak();
	start_color();
	timeout(TIME_OUT);
	/////////noecho();
	clear();
	refresh();
	init_pair(1, COLOR_GREEN, COLOR_GREEN);
	init_pair(2,COLOR_BLUE,COLOR_BLUE);
	return;
}

void printCurses(int lattice[], int nsize, int lsize, unsigned int nsleep){
	int i;
	move(2,10);

	for (i = 0; i < nsize; i++){
		if (((i+1)%lsize == 0)&&(i != 0)){
			if (lattice[i] > 0){
				attron(COLOR_PAIR(1));
				printw("  \n");
				attroff(COLOR_PAIR(1));
			}
			else{
				attron(COLOR_PAIR(2));
				printw("  \n");
				attroff(COLOR_PAIR(2));
			}
			move(2+(i/lsize),10);
		}
		else{
			if (lattice[i] > 0){
				attron(COLOR_PAIR(1));
				printw("  ");
				attroff(COLOR_PAIR(1));
			}
			else{
				attron(COLOR_PAIR(2));
				printw("  ");
				attroff(COLOR_PAIR(2));
			}
		}
	}
	
	move(2,10);
	refresh();
	usleep(nsleep);
	return;
}

void printCurseshalf(int lattice[], int nsize, int lsize, unsigned int nsleep){
	int i;
	move(2,10);

	for (i = 0; i < nsize/2; i++){
		if (((i+1)%lsize == 0)&&(i != 0)){
			if (lattice[i] > 0){
				attron(COLOR_PAIR(1));
				printw("  \n");
				attroff(COLOR_PAIR(1));
			}
			else{
				attron(COLOR_PAIR(2));
				printw("  \n");
				attroff(COLOR_PAIR(2));
			}
			move(2+(i/lsize),10);
		}
		else{
			if (lattice[i] > 0){
				attron(COLOR_PAIR(1));
				printw("  ");
				attroff(COLOR_PAIR(1));
			}
			else{
				attron(COLOR_PAIR(2));
				printw("  ");
				attroff(COLOR_PAIR(2));
			}
		}
	}
	
	move(2,10);
	refresh();
	usleep(nsleep);
	return;
}
