#pragma once

// THIS FILE CONTAINS DECLARATIONS OF SOME GLOBAL VARIABLES AND FUNCTIONS DEFINED EITHER IN "wire.cpp" OR IN "main.cpp"
// THIS FILE MUST BE INCLUDED IN "main.cpp" and "wire.cpp" AS BOTH USE THESE GLOBAL VARIABLES AND FUNCTIONS 
using namespace std;
enum Direction { RIGHT=0, BOTTOM=1, RIGHTDIAGNO=2, LEFT = 3, TOP = 4, LEFTDIAGNOL = 5, NONE = 10 };
static const Point DJshift[6] = { Point(1,0),Point(0,1),Point(1,1),Point(-1,0),Point(0,-1),Point(-1,-1) };


extern Table2D<RGB> image;  // loaded by GUI in main.cpp
extern Table2D<RGB> image1;  // loaded by GUI in main.cpp
extern Table2D<double> SSDimage;
extern Table2D<double> tempimage;
extern Table2D<double> Initialimage;
extern Table2D<double> penalty;


extern Table2D<double> Dijhimage;
extern Table2D<double> Dijvimage;
extern Table2D<double> Dijdimage;
extern Table2D<double>    dist;     // 2D table of "distances" for current paths to the last seed
extern Table2D<Direction> toParent; // 2D table of "directions" along current paths to the last seed.
extern Table2D<double>    Active;



extern Table2D<int> Depthimage;
extern bool closedContour;    // a flag indicating if contour was closed - set in wire.cpp
extern vector<Point> contour;   // a list of 'contour' points - set in wire.cpp
extern Table2D<int> region;   // a binary 2D mask of 'region' points - set in wire.cpp

void addToContour(Point click);
void addToContourLast(Point click);
void regionGrow(Point seed, double T);
void reset_segm();
stack<Point>* liveWire(Point p);
void computePaths(Point seed);
int contourInterior();

extern bool view; // defined in main.cpp (boolean flag set by a check box)
void draw(Point mouse = Point(-1,-1)); // defined in main.cpp, but it is also called in wire.cpp for visualisation 


////////////////////////////////////////////////////////////
// MyPoint objects extend Points with extra variable "m_priority". 
// MyPoint have an overloaded comparison operators based on their priority.
class MyPoint : public Point {
	double m_priority;
public:
	MyPoint(Point p, double priority) : Point(p), m_priority(priority) {}
	double getPriority() const {return m_priority;}	
    bool operator<(const MyPoint& b) const   {return getPriority() > b.getPriority();}
};
