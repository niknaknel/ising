
/*   
     Example of using a pipe to gnuplot directly from the c program
*/
#include <stdio.h>

int main() 
{
  FILE *plotcommandf;              // gnuplot command file
  FILE *pipef;                     // pipe to gnuplot
  FILE *dataf;                     // your numbers in a data file
  char *cmdname = "command_file.output";  // gnuplot command file
  char *datname = "data_file.output";     // your numbers in a data file
  char *some_name_example = "NAME EXAMPLE";
  char *some_x_label = "MY x LABEL";
  char *some_y_label = "MY Favourite Y LABEL";
  int some_y_value[5] = {9, 4, 5, 2, 4, -2};
  int some_x_value[5] = {1, 2, 5, 9, 14, 16};
  int i;

  // Create the data file
  plotcommandf = fopen(cmdname,"w+");
  dataf = fopen(datname,"w+");
  for (i=0; i < 5; i++) {
    fprintf(dataf,"%d\t%d\n",some_x_value[i],some_y_value[i]);
  }
  fclose(dataf);


  /* Create the gnuplot command file
     
     Here you can specify as many gnuplot commands and options as you
     like. Note that you can substitute variable names from C into gnuplot
     with the %s format as shown below. The escaped quotes \" mean that
     C inserts a quote into the output as needed in the gnuplot command file.

     For example, you can set other options such as
     "set mouse; "
     "set size ratio -1; "
     "set key off"

     You can also automatically create png printable output with

     "set term png"
     "set output "out/ME_vs_t_L050T1.900.png"
     "set term X11"

  */ 
  fprintf(plotcommandf,
	  "set title \"Put your title here with variable name %s\" offset 0,-1\n"
	  "set yrange [-2:*]\n"
	  "set xlabel \"x-label, directly or via a variable %s\"   offset 0,0.5\n"
	  "set ylabel \"y-label, directly or via a variable %s\"   offset 0.5,0\n"
	  "plot '%s' using 1:2 with linespoints\n"
	  ,some_name_example, some_x_label, some_y_label,datname);
  fclose(plotcommandf);


  //initialize pipe, give gnuplot initial commands through cdf.output
  pipef = popen("/usr/bin/gnuplot -persist command_file.output  -", "w");
  fflush(pipef);          // force it through the pipe
  fclose(pipef);

  // Close the gnuplot window simply by typing "q" on the keyboard 

  return 0;
}
