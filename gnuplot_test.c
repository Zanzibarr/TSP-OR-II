/*
**  grafico.c
**
**  Calcola per punti la funzione y=f(x), salva sul file
**  "dati.txt" le coordinate dei punti e visualizza il grafico
**  utilizzando GNUplot.
**
**  Marco Liverani (liverani@mat.uniroma3.it) - Marzo 2006
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*
 *  Si definisce una funzione f(x) qualsiasi
 */

float f(float x) {
  return(sin(x)*15 + x*x);
}

/*
 *  Funzione principale.
 */

int main(void) {
  FILE *file;
  float x, x0, x1, deltax, y;
  int n;
  
  printf("Inserisci gli estremi dell'intervallo: ");
  scanf("%f %f", &x0, &x1);
  printf("Inserisci il numero di passi: ");
  scanf("%d", &n);
  deltax = (x1-x0)/n;
  x = x0;
  /*
   *   apro il file "dati.txt" in scrittura per registrare le coordinate
   *   dei punti della funzione
   */
  file = fopen("dati.txt", "wt");

  while (x <= x1) {
    /*
     *  calcolo il valore della funzione nel punto x
     */
    y = f(x);

    /*
     *  scrivo sul file le coordinate del punto (x,y)
     */
    fprintf(file, "%f %f\n", x, y);

    /*
     *  incremento la variabile x
     */
    x = x + deltax;

  }
  /*
   *  chiudo il file con le coordinate dei punti da visualizzare
   */
  fclose(file); 

  /*
   *  apro in scrittura il file "comando.txt" per registrarci il comando
   *  che dovra` essere eseguito da GNUplot
   */
  file = fopen("comando.txt", "wt"); 

  /*
   *  scrivo sul file il comando da eseguire 
   */
  fprintf(file, "plot \"dati.txt\" with lines\n");

  /*
   *  chiudo il file su cui ho scritto il comando da eseguire 
   */
  fclose(file); 

  /*
   *  eseguo il programma GNUplot passandogli il nome del file che 
   *  contiene il comando da eseguire
   */
  system("gnuplot -persistent comando.txt");

  remove("dati.txt");
  remove("comando.txt");

  return(0);
}