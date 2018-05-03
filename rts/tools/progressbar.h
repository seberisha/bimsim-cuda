#ifndef RTS_PROGRESSBAR_H
#define RTS_PROGRESSBAR_H

#include <iostream>
#include <sstream>
using namespace std;

static void rtsProgressBar(int percent)
{
	stringstream bar;
	static int x = 0;
	string slash[4];
	slash[0] = "\\";
	slash[1] = "-";
	slash[2] = "/";
	slash[3] = "|";
	bar<<"[";
	for(int i=0; i<40; i++)
	{
		if(percent > (float)i/(float)40 *100)
			bar << "*";
		else
			bar << " ";
	}
	bar<<"]";
	cout << "\r"; // carriage return back to beginning of line
	cout << bar.str() << " " << slash[x] << " " << percent << " %"; // print the bars and percentage
	x++; // increment to make the slash appear to rotate
	if(x == 4)
	x = 0; // reset slash animation
}

#endif
