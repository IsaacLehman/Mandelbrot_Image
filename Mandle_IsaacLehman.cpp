// Mandelbrot code from Mr. Ounanes
// https://github.com/ChinksofLight/mandelbrot_cpp?source=post_page-----c7ad6a1bf2d9----------------------
// credit above

/* Isaac Lehman ~ Comp 233 ~ Mandlebrot */

/*
	Dr. Valentines description,
	Ounanes’ code takes each pixel in a 2D matrix and counts how many iterations it survives in a calculating loop.
	Pixels (row, col) values that “converge” in the loop are considered members of the Mandelbrot Set. 
	We define “converge” as “survives a set number of iterations without escaping”.  
	Points that “escape” the loop before the maximum iteration count are not members of Mandelbrot.  
	They are used to generate a spectrum of colors to visualize the Mandelbrot Set.

*/

#include <iostream>
#include <fstream>
#include <complex>
#include <time.h>
#include <omp.h>
#include <string>

using namespace std;


// the heigth and width of the ppm image created
#define WIDTH 3000			// width of picture
#define HEIGHT 3000			// height of picture
#define NUM_COLORS 3		// 3 for r,g,b
#define MAX_ITERS 300		// how many iterations needed to pass to be in mandlebrot set
#define MAX_RGB 255			// Max RGB value
#define NUM_COLORS_DIST 9	// how many colors there are in the color distribution

/* COLORS */
// format {red, green, blue}
int RED[3] = { MAX_RGB, 0, 0 };
int CYAN[3] = { 0, MAX_RGB, MAX_RGB };
int GREEN[3] = { 0, MAX_RGB, 0 };
int BLUE[3] = { 0, MAX_RGB, 0 };
int PURPLE[3] = { MAX_RGB, 0, MAX_RGB };
int YELLOW[3] = { MAX_RGB, MAX_RGB, 0 };
int WHITE[3] = { MAX_RGB, MAX_RGB, MAX_RGB };
int ORANGE[3] = { MAX_RGB, 110, 0 };
int BLACK[3] = { 0, 0, 0 };


// the order of the colors for mandlebrot, when you change, remember to change NUM_COLORS_DIST
int* colorDist[NUM_COLORS_DIST] = { PURPLE, BLACK, GREEN, BLACK, YELLOW, BLACK, CYAN, BLACK, ORANGE};



/*
	does the math for pointer access
*/
int offset2D(int row, int col) {
	return (row * HEIGHT) + col;
}


/*
	Creates the rgb colors for a pixel at (row, col) that compose the image.

	if the point survives a set number of iterations, it is condidered a 
	mandlebrot point, return 0.  Otherwise return the number of 
	iterations the point survived.
*/
int value(int row, int col) {
	complex<float> point((float)row / WIDTH - 1.5, (float)col / HEIGHT - 0.5);
	complex<float> z(0, 0);
	int nb_iter = 0;
	while (abs(z) < 2 && nb_iter <= MAX_ITERS) {
		z = z * z + point;
		nb_iter++;
	}
	if (nb_iter < MAX_ITERS){// if didn't survive
		return nb_iter;
	}else // if in mandlebrot set
		return 0;
}


/*
	Sets the rgb values for index of the color map
	- just for cleaner code
*/
void setRGB_colorMap(int colorMap[MAX_ITERS][NUM_COLORS], int iter, int color[3]) {
	colorMap[iter][0] = color[0]; // red
	colorMap[iter][1] = color[1]; // green
	colorMap[iter][2] = color[2]; // blue
}




/*
	How to calculate a color for a given u value
	--------------------------------------------
	UC = current u
	CC = current color
	UL = low u
	CL = low color
	UH = high u
	CH = high color

	do for each: r, g, b
	CC = (CH-CL) * ((UC-UL) / (UH-UL)) + CL;


	- credits to Condron’s -
*/
int calc_CC(double UC, double UL, double CL, double UH, double CH) {
	double CC = (CH - CL) * ((UC - UL) / (UH - UL)) + CL;
	return (int)CC;
}


/*
	calculate the color for a given pixel
	- maps an iteration to a specific RGB color
	  these colors are then to set a specific 
	  pixel in the mandle picture based on the 
	  number of iterations it survived.

	example for 9 colors...
	distribution		(R,G,B)
	u = 0.00  = black	(0,0,0)
	u > 0.000 = white	(255,255,255)
	u = 0.125 = black	(0,0,0)
	u = 0.25  = purple	(0,0,255)
	u = 0.375 = black	(0,0,0)
	u = 0.50  = red		(255,0,0)
	u = 0.625 = black	(0,0,0)
	u = 0.75  = green	(0,255,0)
	u = 0.875 = black	(0,0,0)
	u = 1.00  = Yellow	(0,255,255)

	Set the colors above in colorDist
	and remember to set NUM_COLORS_DIST to the
	number of colors you used
*/
void fill_color_map(int colorMap[MAX_ITERS][NUM_COLORS], double* u) {
	// variables
	double UC, UL, UH; // as defined above
	int spotH, spotL; // the locations in the colorDist array
	const double interval = 1.0 / (NUM_COLORS_DIST - 1); // the interval between colors
	int color[MAX_RGB]; // the curent color - i.e. CC


	/* Mandle points = black, u = 0*/
	setRGB_colorMap(colorMap, 0, BLACK);

	// not worth parallelising... loop of only 300ish
	for (int iter = 1; iter < MAX_ITERS; iter++)
	{
		UC = u[iter]; 

		double i = interval;
		while (i < UC) { // find the next interval above UC
			i += interval;
		}
		UH = i;
		UL = UH - interval;
		spotH = (int)(i * (NUM_COLORS_DIST - 1)); // zero indexed spot for colorDist
		spotL = spotH - 1;

		// calculate the current color
		for (int i = 0; i < MAX_RGB; i++)
		{
			color[i] = calc_CC(UC, UL, colorDist[spotL][i], UH, colorDist[spotH][i]);
		}

		// update color map with current color
		setRGB_colorMap(colorMap, iter, color);
	}
}



int main() {
	// header
	cout << "Isaac Lehman ~ Comp 233 ~ Mandlebrot" << endl;
	cout << "Width: " << WIDTH << "\tHeight: " << HEIGHT << "\n------------------------------------------" << endl;


	ofstream my_Image("mandelbrot.ppm"); // create ppm file that will store the picture


	if (my_Image.is_open()) { // if file was made succesfully
		// print some documentation, # represents a comment in ppm
		// p3 declares it's a ppm file, the dimensions of the file, and the rgb range.  
		my_Image << "P3\n# Comp 233 ~ Isaac Lehman ~ Mandlebrot image in ASCII PPM\n# Special thanks to Ounane for his code\n";
		my_Image << WIDTH << " " << HEIGHT << "\n255\n";

		/* VARIABLES */
		int counts[MAX_ITERS] = { 0 };		// number of pixels that survived the given number of iterations
		int histogram[MAX_ITERS] = { 0 };	// reports a running total of counts
		double u[MAX_ITERS] = { 0.0 };		// uniform distrobution of points outside mandlebrot set
		int colorMap[MAX_ITERS][NUM_COLORS];// the rgb colors for a given number of iterations survived
		clock_t startT, stopT;				// wallclock timer
		double time;						// time to runProgram
		int* pixelIterCount;				// the number of iterations a pixel at (row, col) has survived


		pixelIterCount = (int*)malloc(WIDTH * HEIGHT * sizeof(int));

		omp_set_num_threads(4); // set our requested number of threads

		/*
			Calculate the pixels for the image
		*/
		startT = clock();	//start stopwatch
		{
#pragma omp parallel for schedule(dynamic)//parallise the outer for loop
			for (int r = 0; r < WIDTH; r++) {
				for (int c = 0; c < HEIGHT; c++) {
					int val = value(r, c);						// get the # of iters (row, col) survived
#pragma omp atomic
					counts[val]++;							// incriment the count, critical to avoid race condition
					*(pixelIterCount + offset2D(r, c)) = val;	// keep track of the number of iters for a pixel
				}
			}

			/*
				Make the histogram
			*/
			int sum = 0;
			for (int i = 1; i < MAX_ITERS; i++) // start at 1, only do points not in mandle set
			{
				sum += counts[i];
				histogram[i] = sum;
			}


			/*
				make the u distrobution
			*/
			int total_nonMandle = histogram[MAX_ITERS - 1]; // # of non-Mandelbrot Set points in our image
			for (int i = 1; i < MAX_ITERS; i++) // start at 1, only do points not in mandle set
			{
				u[i] = (double)histogram[i] / total_nonMandle; // percent of total
			}


			/*
				make the color map
			*/
			fill_color_map(colorMap, u);
		}
		stopT = clock();	//stop stopwatch
		time = (double)(stopT - startT) / CLOCKS_PER_SEC;
		cout << "image created,\tTIME:\t\t" << time << " sec" << endl;



		/*
			Print the image to the ppm file
		*/
		startT = clock();	//start stopwatch
		{
			int size = WIDTH * HEIGHT;
			for (int r = 0; r < WIDTH; r++) {
				int c = 0;
				while (c < HEIGHT) {
					string line = "";
					int z = 0;
					// write 5 rgb value sets to each line in the ppm file
					while (z < 5 && c < HEIGHT) {

						// get the number of iterations a pixel survived
						int iters = *(pixelIterCount + offset2D(r, c));
						// find the coresponding rgb values
						int r = colorMap[iters][0];
						int g = colorMap[iters][1];
						int b = colorMap[iters][2];

						// append to string
						line += to_string(r);
						line += " ";
						line += to_string(g);
						line += " ";
						line += to_string(b);
						line += " ";
						z++;
						c++;
					}
					// write line to file
					my_Image << line << "\n";
				}
			}
		}
		stopT = clock();	//stop stopwatch
		time = (double)(stopT - startT) / CLOCKS_PER_SEC;
		cout << "image printed,\tTIME:\t\t" << time << " sec" << endl;


		// collect garbage
		free(pixelIterCount);
		pixelIterCount = NULL;

		my_Image.close(); // close the file
	}
	else // if file couldn't be opened
		cout << "Could not open the file";


	cout << "\n\n\t\t<Normal Termination>\n\n";
	return 0;
}