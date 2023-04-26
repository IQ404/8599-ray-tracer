# 8599 Ray Tracer

This repository aims to record a path tracer that I am currently writing for the course CSC8599.

## Objectives:

- Prototype an offine path tracer (reference: [the book by Peter Shirley](https://raytracing.github.io/books/RayTracingInOneWeekend.html)).
- Improve the ray-tracer by adding more advanced features to it (reference: [the next book by Peter Shirley](https://raytracing.github.io/books/RayTracingTheNextWeek.html)).
- Optimize the path tracer (from the CPU side, e.g. multithreading) and integrate its output into a GUI. Try to achieve >= 10 fps for a scene that contains a relatively small number of items without sacrificing the advanced features as much as possible (reference: [this series by The Cherno](https://www.youtube.com/watch?v=gfW1Fhd9u9Q&list=PLlrATfBNZ98edc5GshdBtREv5asFW3yXl&index=1)).
- Try to make the path tracer real-time. E.g. using DirectX RayTracing (DXR) or Vulkan.


---

## Progression

### April 25th 2023

- Hello Graphic World

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  A trial program used to output PPM format image example
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

/*

Note: PPM image can be viewed by **Portable Anymap Viewer** on Windows

*/

#include <iostream>

using namespace std;

int main()
{
	// image size:

	const int image_width = 256;
	const int image_height = 256;
	const int max_color = 255;

	// output data: (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)

	cout << "P3" << endl						// colors are in ASCII
		<< image_width << " " << image_height << endl		// column  row
		<< max_color << endl;					// value for max color

	// RGB triplets: (rendered from left to right, top to bottom)

	for (int row = image_height - 1; row >= 0; row--)
	{
		cerr << '\r' << "Scanlines Remaining: " << row << ' ' << flush;		// ??? Why do we want std::flush here?
											// \r means writing from the head of the current line

		for (int column = 0; column < image_width; column++)
		{
			// rgb value ranging 0.0 - 1.0
			double r = double(column) / (image_width - 1);	// interpolate the width (left == 0; right == 1)
			double g = double(row) / (image_height - 1);	// interpolate the height (top == 0; bottom == 1)
			double b = 0.25;				// for each pixel the portion of blue is constant

			int red = int(max_color * r);
			int green = int(max_color * g);
			int blue = int(max_color * b);			// ??? what's the difference if we use max_color == 255.999?

			cout << red << " " << green << " " << blue << endl;
		}
	}
	cerr << '\n' << "Done" << endl;
}
```

The output of the above program is as follows: (TODO: producing image of other format using [stb_image.h](https://github.com/nothings/stb))

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/HelloGraphicWorld.jpg" width="600" height="600"></a>

- The `vec3` class


