#include "cs1037lib-window.h" // for basic drawing and keyboard/mouse input operations 
#include "cs1037lib-button.h" // for basic buttons and other GUI controls 
#include "Cstr.h"       // convenient macros for generating C-style text strings
#include <iostream>     // for cout
#include <vector>
#include <stack>
#include "Basics2D.h"
#include "Table2D.h"
#include "Math2D.h"
#include <cmath>
#include "Image2D.h"
#include "wire.h" // declarations of global variables and functions shared with wire.cpp
#include <stdio.h>
#include "graph.h"//header file for graphcut

using namespace std;



// declarations of global variables used for GUI controls/buttons/dropLists
const char* image_names[] = { "Tsukuba_L","Tsukuba_R" }; // an array of image file names
int im_index = 0;    // index of currently opened image (inital value)
int im_index1 = 1;
const char* mode_names[]  = { "Windows", "Scan-Line","SLGraph","Global","Dijkstra"}; // an array of mode names 
enum Mode  {Windows=0, ScanLine=1,SLGraph=2,Global=3,Dijkstra=4}; // see main.cpp
Mode mode = Windows; // index of the current mode (in the array 'mode_names')
bool view = false; // flag set by the check box
const int cp_height = 34; // height of "control panel" (area for buttons)
const int pad = 10; // width of extra "padding" around image (inside window)
int w = 3; // region growing threshold
int T_box; // handle for T-box (for setting threshold parameter T)

// declarations of global functions used in GUI (see code below "main")

void image_load(int index); // call-back function for dropList selecting image file
void image_load1(int index);
void mode_set(int index);   // call-back function for dropList selecting mode 
void clear();  // call-back function for button "Clear"
void T_set(const char* T_string);
void Caculate();
//void Integral();
//void ComputeDijstra(Point seed);

int main()
{
	  // initializing buttons/dropLists and the window (using cs1037utils methods)
	int blank = CreateTextLabel(""); // adds grey "control panel" for buttons/dropLists, see "cs1037utils.h"
    SetControlPosition(blank,0,0); SetControlSize(blank,1280,cp_height); // see "cs1037utils.h"
	int dropList_files = CreateDropList(2, image_names, im_index, image_load); // the last argument specifies the call-back function, see "cs1037utils.h"
	int dropList_files1 = CreateDropList(2, image_names, im_index1, image_load1);
	int label = CreateTextLabel("Mode:"); // see "cs1037utils.h"
	int dropList_modes = CreateDropList(5, mode_names, mode, mode_set); // the last argument specifies the call-back function, see "cs1037utils.h"
	T_box = CreateTextBox(to_Cstr("w=" << w), T_set);
	int button_clear = CreateButton("Clear",clear); // the last argument specifies the call-back function, see "cs1037utils.h"
	int button_Caculate = CreateButton("Calculate", Caculate);
	
	SetWindowTitle("Stereo Matching");      // see "cs1037utils.h"
    SetDrawAxis(pad,cp_height+pad,false); // sets window's "coordinate center" for GetMouseInput(x,y) and for all DrawXXX(x,y) functions in "cs1037utils" 
	                                      // we set it in the left corner below the "control panel" with buttons

	  // initializing the application
	image_load(im_index);
	image_load1(im_index1);
	SetWindowVisible(true); // see "cs1037utils.h"

	  // while-loop processing keys/mouse interactions 
	while (!WasWindowClosed()) // WasWindowClosed() returns true when 'X'-box is clicked
	{
		char c;
		if (GetKeyboardInput(&c)) // check keyboard
		{ 
			if (c == 'd') clear(); // 
		}
	}

	  // deleting the controls
	DeleteControl(button_clear);    // see "cs1037utils.h"
	DeleteControl(dropList_files);
	DeleteControl(dropList_modes);     
	DeleteControl(label);
	return 0;
}

// call-back function for the 'mode' selection dropList 
// 'int' argument specifies the index of the 'mode' selected by the user among 'mode_names'
void mode_set(int index)
{
	mode = (Mode) index;
	cout << "drawing mode is set to " << mode_names[index] << endl;
	reset_segm();
	draw();
}


void T_set(const char* T_string) {
	sscanf_s(T_string, "w=%d", &w);
	cout << "parameter w is set to " << w << endl;
}


// call-back function for the 'image file' selection dropList
// 'int' argument specifies the index of the 'file' selected by the user among 'image_names'
void image_load(int index) 
{
	im_index = index;
	cout << "loading image file " << image_names[index] << ".bmp" << endl;
	image = loadImage<RGB>(to_Cstr(image_names[index] << ".bmp")); // global function defined in Image2D.h
	int width  = max(400,(int)image.getWidth())*2 + 2*pad + 80;
	int height = max(100,(int)image.getHeight())*2+ 2*pad + cp_height;
	SetWindowSize(width,height); // window height includes control panel ("cp")
   // SetControlPosition(    T_box,     image.getWidth()+pad+5, cp_height+pad);
    //SetControlPosition(check_box_view,image.getWidth()+pad+5, cp_height+pad+25);
	reset_segm();  // clears current "contour" and "region" objects - function in a4.cpp
	draw();
}

void image_load1(int index)
{
	im_index = index;
	cout << "loading image file " << image_names[index] << ".bmp" << endl;
	image1 = loadImage<RGB>(to_Cstr(image_names[index] << ".bmp")); // global function defined in Image2D.h
	int width1 = max(400, (int)image1.getWidth())*2 + 2 * pad + 80;
	int height1 = max(100, (int)image1.getHeight())*2 + 2 * pad + cp_height;
	SetWindowSize(width1, height1); // window height includes control panel ("cp")
	//SetControlPosition(T_box, image.getWidth() + pad + 5, cp_height + pad);
	//SetControlPosition(check_box_view, image.getWidth() + pad + 5, cp_height + pad + 25);
	reset_segm();  // clears current "contour" and "region" objects - function in a4.cpp
	draw();
}

// call-back function for button "Clear"
void clear() { 
	reset_segm(); // clears current "contour" and "region" objects - function in wire.cpp
	image_load1(im_index1);
	draw();
}

void draw(Point mouse)
{ 
	unsigned i;
	// Clear the window to white
	SetDrawColour(255,255,255); DrawRectangleFilled(-pad,-pad,1280,1024);
	
	if (!image.isEmpty()) drawImage(image); // draws image (object defined in wire.cpp) using global function in Image2D.h (1st draw method there)
	else {SetDrawColour(255, 0, 0); DrawText(2,2,"image was not found"); return;}

	if (!image1.isEmpty()) drawImage(image1,Point(image.getWidth()+5*pad,0)); // draws image (object defined in wire.cpp) using global function in Image2D.h (1st draw method there)
	else { SetDrawColour(255, 0, 0); DrawText(2, 2, "image was not found"); return; }
}

void Integral(int d, int n)
{
	unsigned int width = image.getWidth();
	unsigned int height = image.getHeight();
	for (int i = d; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			double temp= double(image[Point(i, j)]) - double(image1[Point(i - d, j)]);
			Initialimage[Point(i, j)] = temp*temp;
		}
	}
	double *columnSum = new double[height]; 
	for (int i = 0; i<width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			if(j==0)
			columnSum[0] = Initialimage[Point(i,0)];
			if(j>0)
			columnSum[j] = columnSum[j-1] + Initialimage[Point(i,j)];
			if(i==0)
			SSDimage[Point(0,j)] = columnSum[j];
			if (i > 0)
			SSDimage[Point(i,j)] = SSDimage[Point(i-1,j)] + columnSum[j];

		}
	}
	/*
	image = Initialimage;
	image1 = SSDimage;
	draw();
	*/
	//geting the windows.
	double twssd, wssd;

	if (n == 0)
	{
		for (int i = w; i < width - w; i++)
		{
			for (int j = w; j < height - w; j++)
			{
				wssd = SSDimage[Point(i + w, j + w)] - SSDimage[Point(i - w, j + w)] - SSDimage[Point(i + w, j - w)] + SSDimage[Point(i - w, j - w)];
				tempimage[Point(i, j)] = wssd;
			}
		}
	}

	if (n > 0)
	{
		for (int i = w; i < width - w; i++)
		{
			for (int j = w; j < height - w; j++)
			{
				wssd = SSDimage[Point(i + w, j + w)] - SSDimage[Point(i - w, j + w)] - SSDimage[Point(i + w, j - w)] + SSDimage[Point(i - w, j - w)];
			    twssd = tempimage[Point(i, j)];
				if (wssd < twssd)
				{
					Depthimage[Point(i, j)] = d * 20;
					tempimage[Point(i, j)] = wssd;
					//cout << Depthimage[Point(i, j)] << endl;
				}
				//cout << tempimage[Point(i, j)] << endl;
			}
		}
	}
}


void DP_move(int j,int dmax)
{
	int n = image.getWidth();
	int cur = 0;
	vector<vector<double>> energy;
	vector<Point> path;
	energy.clear();
	energy.resize(n);
	for (int ln = 0; ln < n; ln++)
		for (int i = 0; i <= dmax; i++)
			energy[ln].push_back(0.0);

	Point current, next;
	double unary = 0;
	double pair = 0;
	double min, cost, tempe;
	cost = 5;
	double cons = 20, sigma = 100;
	while (1)
	{
		current = Point(cur, j);
		next = Point(cur + 1, j);
		int rd = -1;
		int rs = -1;
		for (int d = 1; d <= dmax; d++)
		{
			min = INFINITY;
			rd = -1;
			rs = -1;
			for (int s = 1; s <= dmax; s++)
			{
				Point match = current - Point(s, 0);
				//cout << match.x << "  " << match.y << endl;
				//unary = abs(double(image[current]) - double(image1[match]));
				unary = (double(image[current]) - double(image1[match]))*(double(image[current]) - double(image1[match]));
				unary = abs(double(image[current]) - double(image1[match]));
				if (unary >100) 	unary = 100;
				pair = cost*abs(s - d);
				//double dif = double(image[current]) - double(image[next]);
				//pair = cons * exp(- dif*dif / (2 * sigma * sigma))*abs(s-d);
				tempe = unary + pair;
				if (tempe + energy[cur][s] < min)
				{
					min = tempe + energy[cur][s];
					rs = s;
					rd = d;
					//cout << "i=" << i << ",j=" << j << ",min=" << min << endl;
					//¼ÇÂ¼i,j,min
				}
			}

			if (min != INFINITY)
			{
				Point temp(rs, rd);
				path.push_back(temp);
				energy[cur + 1][rd] = min;
			}
		}

		cur++;
		if (cur == n - 1)
		{
			break;
		}
	}
	double finalmin = INFINITY;
	int last = -1;
	for (int k = 1; k <=dmax; k++)
	{
		if (energy[n - 1][k] < finalmin)
		{
			finalmin = energy[n - 1][k];
			last = k;
		}
	}
	int num = n - 1;
	while (1)
	{
		bool s = 0;
		Depthimage[Point(num,j)] = last*20;
		if (num == 0)
			break;
		int tempindex = (num - 1) * dmax + last;
		last = path[tempindex].x;
		num--;
	}

}

void ComputeDijstra(int y)
{

	//if (!image.pointIn(seed)) return;
	//region.reset(0); // resets 2D table "region" for visualization
					 // Reset 2D arrays "dist" and "toParent"  (erazing all current paths)
	vector<Point> active;
	Active.reset(-1);
	Point current(0, 0);
	dist.reset(INFINITY);
	toParent.reset(NONE);
	dist[current] = 0;
	Active[current] = 1;
	double weight = INFINITY;


	//queue<Point> active;
	//active.push(seed);
	// STUDENTS: YOU NEED TO REPLACE THE CODE BELOW (which is a copy of BFS algorithm in region growing)
	// Create a queue (priority_queue) for "active" points/pixels and
	// traverse pixels to improve paths stored in "toParent" (and "dist")

	//Step 1: Calculating the weights.
	//In function set penalty with Gaussian Function.

	//Step2:Dijstra Algorith. For 4&8 connected graph.

	//cout << "(" << seed.x << "," << seed.y << ")" << endl;
	while (1)
	{
		//if (!image.pointIn(seed)) return;
		//Point k(-1, -1);

		int k = -1;
		for (int i = 0; i < 6; i++)
		{
			Point q = current + DJshift[i];

			if ((dist.pointIn(q)))
			{
				if (i == 0)  weight = Dijhimage[current];
				if (i == 1)  weight = Dijvimage[current];
				if (i == 2)  weight = Dijdimage[current];
				if (i == 3)  weight = Dijhimage[q];
				if (i == 4)  weight = Dijvimage[q];
				if (i == 5)  weight = Dijdimage[q];

				if (Active[q] != 1 && weight!=INFINITY)
				{

			    	if (dist[current] + weight < dist[q])
					{
						dist[q] = dist[current] + weight;
						if (i == 0) { toParent[q] = LEFT; }
						if (i == 1) { toParent[q] = TOP; }
						if (i == 2) { toParent[q] = LEFTDIAGNOL; }
						if (i == 3) { toParent[q] = RIGHT; }
						if (i == 4) { toParent[q] = BOTTOM; }
						if (i == 5) { toParent[q] = RIGHTDIAGNO; }
						//cout << "weight= " << weight << endl;
						//cout << "Direction= " << toParent[q];

					}

					if (Active[q] == -1)
					{
						active.push_back(q);
						Active[q] = 0;
					}
				}

			}

		}

		double min = INFINITY;

		for (int x = 0; x < active.size(); x++)
		{
			Point temp = active[x];
			if (dist[temp] < min)
			{
				min = dist[temp];
				k = x;
			}
		}
		if (k != -1)
			Active[active[k]] = 1;
		current = active[k];
		//cout << "next Point.x=" << active[k].x << "next Point.y=" << active[k].y << endl;
		//cout << "before active_size=" << active.size() << endl;
		active.erase(active.begin() + k);
		//cout << "active_size=" << active.size() << endl;

		if (k == -1)
			break;
	}

	Point temp(image.getWidth() - 1, image.getWidth() - 1);
	int num = image.getWidth() - 1;
	Direction n = toParent[temp];
	/*
	for (int i = 0; i < 1; i++)
	{
		for (int j = 0; j < image.getWidth(); j++)
		{
			cout << toParent[Point(i,j)]<<endl;
		}
	}
	*/

//	while (n != NONE)
	int depth = 0;
	while (num >= 0)
	{

		//cout << n << endl;
		//cout << toParent[temp] << endl;
	
	    if (toParent[temp] == 4)
				depth++;
		if (toParent[temp] == 3)
				depth--;
		if (toParent[temp] == 5)
		{
			Depthimage[Point(num, y)] = depth * 20;
			//cout << "(" << num << "    ,    " << y << ") = " << abs(temp.y - temp.x) << endl;
		}
		num--;
		temp = temp + DJshift[n];
		n = toParent[temp];
		//cout << toParent[temp];
		//cout << "(" << temp.x << "    ,    " << temp.y <<")"<< endl;
		//cout << temp.x - temp.y << endl;
		//if(toParent[temp]<3)
			//cout << toParent[temp]<<endl;
	}

	if (y == image.getHeight()-1)
	//if (y == 190)
	{
			image1 = Depthimage;
			draw();
	}
}


void Caculate()
{
	int dmax = 15;
	if (mode == Windows)
	{

		int count = 0;
		for (int d = 1; d <= dmax; d++)
		{
			//Caculate the Integral Map.
			Integral(d, count);
			count++;
		}

		image1 = Depthimage;
		draw();
	}

	if (mode == ScanLine)
	{

		for (int j = 0; j < image.getHeight(); j++)
		{
			DP_move(j, dmax);
		}
		image1 = Depthimage;
		draw();
	}

	if (mode == SLGraph)
	{
		unsigned int width = image.getWidth();
		unsigned int height = image.getHeight();
		double Dp = 0;
		double wpq = 0;
		double sigma = 100;
		double cons = 20;
		//cout << height << "   " << width << endl; 288,384
		for (int j = 0; j < height; j++)
		{
			//cout << "----------------------" <<j<< endl;
			typedef Graph<double, double, double> GraphType;
			unsigned int NumNode = width*dmax;
			unsigned int NumEdge = 2 * width*dmax - width - dmax;
			GraphType *g = new GraphType(NumNode, NumEdge);
			g->add_node(NumNode);
			for (int i = 0; i < width; i++)
			{
				//cout << "-----------------------" << endl;
				for (int d = 0; d < dmax; d++)
				{
					int n = d*width + i;
					if (d == 0)
					{

						g->add_tweights(n, 10, 10000);
					}
					if (d > 0)
					{
						Dp = abs(double(image[Point(i, j)]) - double(image1[Point(i - d, j)]));
						//cout << "Dp=" << Dp << endl;
						//Dp = Dp*Dp;
						g->add_edge(n, n - width, Dp, Dp);
					}
					if (d == dmax - 1)
					{
						g->add_tweights(n, 10000, 10);
					}
					if (i >= 1)
					{
						//cout << n << endl;
						wpq = abs(double(image[Point(i, j)]) - double(image[Point(i - 1, j)]));
						wpq = cons * exp(-wpq*wpq / (2 * sigma * sigma));
						//wpq = 20;

						g->add_edge(n, n - 1, wpq, wpq);
					}
				}
			}

			double flow = g->maxflow();
			//cout << "Flow = " << flow << endl;

			for (int k = 0; k < width; k++)
			{
				//if(k>=80 && k<85)
				//cout << "************************" << endl;
				for (int d = 0; d < dmax - 1; d++)
				{
					int s = width*d + k;
					//if (k>=80 && k<85)
					//cout << g->what_segment(s) << endl;

					if (g->what_segment(s) != g->what_segment(s + width))
					{
						//cout << "d= " << d << endl;
						Depthimage[Point(k, j)] = d * 20;
					}

				}
			}
			//cout << " &&&&&&&" << endl;
			delete g;
		}

		image1 = Depthimage;
		draw();
	}
	if (mode == Global)
	{
		typedef Graph<double, double, double> GraphType;
		unsigned int width = image.getWidth();
		unsigned int height = image.getHeight();
		unsigned int NumNode = width*height;
		unsigned int NumEdge = (2 + dmax) * width*height - width - height;
		double Dp = 0;
		double wpq1 = 50;
		double wpq2 = 50;
		double sigma = 100;
		double cons = 10;
		GraphType *g = new GraphType(NumNode*dmax, NumEdge);
		g->add_node(NumNode*dmax);

			for (int i = 0; i < height; i++)
			{
				for (int j = 0; j < width; j++)
				{
		        	for (int d = 0; d < dmax; d++)
					{
					int pos = i*width + j;
					int n = d*height*width + pos;
					if (j >= 1)
					{
						//cout << "i=" << i << "   ,j=" << j << endl;
						wpq1 = abs(double(image[Point(j, i)]) - double(image[Point(j - 1, i)]));
						//wpq1 = wpq1*wpq1;
						wpq1 = cons * exp(-wpq1*wpq1 / (2 * sigma * sigma));
						
						g->add_edge(n, n - 1, wpq1, wpq1);
					}
					if (i >= 1)
					{
						wpq2 = abs(double(image[Point(j, i)]) - double(image[Point(j, i - 1)]));
						//wpq2 = wpq2*wpq2;
						wpq2 = cons * exp(-wpq2*wpq2 / (2 * sigma * sigma));
						g->add_edge(n, n - width, wpq2, wpq2);
					}
					if (d == 0)
					{
						g->add_tweights(n, 10, 10000);
					}
					if (d > 0)
					{
						Dp = abs(double(image[Point(j, i)]) - double(image1[Point(j - d, i)]));
						g->add_edge(n, n - height*width, Dp, INFINITY);
					}
					if (d == dmax-1)
					{
						g->add_tweights(n, 10000, 10);
					}
				}
			}

		}
		double flow = g->maxflow();
		cout << "Flow = " << flow << endl;

		for (int a = 0; a < height; a++)
		{
			for (int b = 0; b < width; b++)
			{
				int k = a*width + b;
				for (int d =0; d < dmax-1; d++)
				{
					// cout << "d= " << d << "     k=" << k << endl;
					long int s = height*width*d + k;
				//	if (k == 4820)
					//	cout << g->what_segment(s) << endl;
					if (g->what_segment(s) != g->what_segment(s + width*height))
					{
						// cout << "d= " << d << endl;
						Depthimage[Point(b, a)] = d * 20;

					}

				}
			}
		}
		//cout << " &&&&&&&" << endl;
		delete g;
		image1 = Depthimage;
		draw();
	}
	
	if (mode == Dijkstra)
	{
		unsigned int width = image.getWidth();
		unsigned int height = image.getHeight();
		double wh, wv, wd;
		wh = 15000;
		wv = 15000;
		for (int j = 0; j < height; j++)
		//for (int j = 0; j <=150; j++)
		{
			Dijhimage.reset(INFINITY);
			Dijvimage.reset(INFINITY);
			Dijdimage.reset(INFINITY);

			for (int i = 0; i < width-1; i++)
			{
				wh = abs(double(image1[Point(i + 1, j)]) - double(image1[Point(i, j)]));
				//if (i == 0)cout << wh << endl;
				wv = abs(double(image[Point(i + 1, j)]) - double(image[Point(i, j)]));
	

				for (int t = 0; t < width - 1; t++)
				{
					//if(i<=t)
					Dijhimage[Point(i, t)] = wh;
				}
				for (int k = 0; k < width-1; k++)
				{
				   // if(i>=k)
					Dijvimage[Point(k, i)] = wv;

					//if (k <= i )//&& i-k <= dmax)
					{

						wd = abs(double(image[Point(i, j)]) - double(image1[Point(k, j)]));
						wd = wd*wd;
						Dijdimage[Point(i, k)] = wd;
					}
				}
			}
			/*
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					cout << "----------------" << endl;
					cout << j << "   " << i << endl;
					cout << double(Dijhimage[Point(j,i)]) << endl;
					cout << double(Dijvimage[Point(j,i)]) << endl;
					cout << double(Dijdimage[Point(j,i)]) << endl;
				}
			}
			*/
			ComputeDijstra(j);

		}
		cout << "&&&***" << endl;

	}
	
}


