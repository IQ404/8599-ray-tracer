# 8599 Ray Tracer

This repository aims to record a path tracer that I am currently writing for my MSc module CSC8599.

## Objectives:

- Prototype an offine path tracer.
  
  <ins>Possible references</ins>: [the book by Peter Shirley](https://raytracing.github.io/books/RayTracingInOneWeekend.html)
  
  **Estimated date of completion:** no later than <ins>15th May 2023</ins>

- Improve the ray-tracer by adding more advanced features to it.
  
  <ins>Possible references</ins>: [the next book by Peter Shirley](https://raytracing.github.io/books/RayTracingTheNextWeek.html)
  
  **Estimated date of completion:** <ins>early July 2023</ins>

- Optimize the path tracer (from the CPU side, e.g. multithreading) and integrate its output into an interactive GUI. Try to achieve >= 10 fps for a scene that contains a relatively small number of items, a relatively low sampling rate and a relatively low resolution, so that we don't need to sacrifice the advanced features as much as possible.

  <ins>Possible references</ins>: [the series by The Cherno](https://www.youtube.com/watch?v=gfW1Fhd9u9Q&list=PLlrATfBNZ98edc5GshdBtREv5asFW3yXl&index=1), [Walnut GUI by The Cherno](https://github.com/TheCherno/Walnut), [Futhark](https://github.com/athas/raytracinginoneweekendinfuthark), [boksajak](https://github.com/boksajak/raytracingthenextweek)
  
  **Estimated date of completion:** optimization should be ongoing indefinitely. But the first complete version of the ray tracer running in real-time is aiming to be finished by the <ins>end of May 2023</ins>

- Try to accelerate the path tracer by running code on GPU.
  
  <ins>Possible references</ins>: [CUDA](https://developer.nvidia.com/blog/accelerated-ray-tracing-cuda/), [DXR](https://github.com/theroyn/RealTimeRayTracing), [OptiX](https://developer.nvidia.com/rtx/ray-tracing/optix), VK_NV_ray_tracing
  
  **Estimated date of completion:** The deadline of this module (CSC8599) is <ins>22nd August 2023</ins>, efforts on implementing GPU accelerations should cease no earlier than this date.

## Housekeeping:

- When developing this project, the code is running on <ins>11980HK</ins> and <ins>RTX3080 (laptop version)</ins>.
- The code is written in Visual Studio 2022 with C++20.
- This implementation uses RIGHT-handed coordinate systems.

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

- The `Vector3D` class

This is the class that will be used to represent anything that can be encoded as a 3-tuple in $\mathbb R^3$.

```cpp
/*****************************************************************//**
 * \file   Vector3D.h
 * \brief  The class of 3D vector for representing geometry in R^3 and RGB color
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <cmath>		// to use std::sqrt
#include <cassert>
 //#define NDEBUG		// uncomment this if we don't want assertion (e.g. when we want things like inf)

class Vector3D
{
	double v[3];

public:

	// Constructors:

	Vector3D()
		: v{ 0,0,0 }
	{

	}

	Vector3D(double x, double y, double z)
		: v{ x,y,z }
	{

	}

	// Operators:

	Vector3D operator-() const					// additive inverse
	{
		return Vector3D{ -v[0], -v[1], -v[2] };
	}

	double& operator[](int i)
	{
		assert(i == 0 || i == 1 || i == 2);
		return v[i];
	}

	double operator[](int i) const
	{
		assert(i == 0 || i == 1 || i == 2);
		return v[i];
	}

	Vector3D& operator+=(const Vector3D& u)
	{
		v[0] += u.v[0];
		v[1] += u.v[1];
		v[2] += u.v[2];

		return *this;
	}

	Vector3D& operator*=(const double d)
	{
		v[0] *= d;
		v[1] *= d;
		v[2] *= d;

		return *this;
	}

	Vector3D& operator/=(const double d)
	{
		assert(d != 0.0);

		v[0] /= d;
		v[1] /= d;
		v[2] /= d;

		return *this;
	}

	// Methods:

	double x() const
	{
		return v[0];
	}

	double y() const
	{
		return v[1];
	}

	double z() const
	{
		return v[2];
	}

	double squared_length() const
	{
		return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	}

	double length() const
	{
		return std::sqrt(squared_length());
	}

};

// Unitility Functions:

inline std::ostream& operator<<(std::ostream& os, const Vector3D& v)
{
	return os << v.x() << ' ' << v.y() << ' ' << v.z();
}

inline Vector3D operator+(const Vector3D& a, const Vector3D& b)
{
	return Vector3D{ a.x() + b.x(), a.y() + b.y(), a.z() + b.z() };
}

inline Vector3D operator-(const Vector3D& a, const Vector3D& b)
{
	return Vector3D{ a.x() - b.x(), a.y() - b.y(), a.z() - b.z() };
}

inline Vector3D operator*(const Vector3D& a, const Vector3D& b)			// Note that this is NOT dot or cross product!
{
	return Vector3D{ a.x() * b.x(), a.y() * b.y(), a.z() * b.z() };
}

inline Vector3D operator*(double d, const Vector3D& v)
{
	return Vector3D{ d * v.x(), d * v.y(), d * v.z() };
}

inline Vector3D operator*(const Vector3D& v, double d)
{
	return d * v;
}

inline Vector3D operator/(const Vector3D& v, double d)
{
	return (1 / d) * v;
}

inline double dot(const Vector3D& a, const Vector3D& b)
{
	return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

inline Vector3D cross(const Vector3D& a, const Vector3D& b)
{
	return Vector3D
	{
		a.y() * b.z() - a.z() * b.y(),
		a.z() * b.x() - a.x() * b.z(),
		a.x() * b.y() - a.y() * b.x()
	};
}

inline Vector3D unit_vector(const Vector3D& v)
{
	return v / v.length();
}

// For better code readability (as Vector3D will represent things with different physical meanings):

using point3D = Vector3D;
using colorRGB = Vector3D;

#endif // !VECTOR3D_H
```

- The `color.h` header

A function is created for outputting a `colorRGB` represented by a `vector3D` to the standard output as follows:

```cpp
/*****************************************************************//**
 * \file   color.h
 * \brief  Utility Functions for outputting ColorRGB object as pixel
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/
 
#ifndef COLOR_H
#define COLOR_H

#include "Vector3D.h"

#include <iostream>

void write_color(std::ostream& os, colorRGB pixel_color)
{
	// Assume the rgb component values of colorRGB is in range [0.0, 1.0], and the output integer value is in range [0, 255].

	os << int(255.999 * pixel_color.x()) << ' '
	   << int(255.999 * pixel_color.y()) << ' '
	   << int(255.999 * pixel_color.z()) << '\n';

	// ??? Explain why .999 can be necessary.
	// ??? What is the difference if static_cast<int>() is used instead of int()?

}

#endif // !COLOR_H
```

### April 26th 2023

- The `Ray` class

```cpp
/*****************************************************************//**
 * \file   Ray.h
 * \brief  The class representing the ray
 * 
 * \author Xiaoyang
 * \date   April 2023
 *********************************************************************/

#ifndef RAY_H
#define RAY_H

#include "Vector3D.h"

class Ray
{
	Point3D orig;
	Vector3D dir;

public:

	// Constructors:

	Ray()	// both origin and direction are initialized to (0,0,0), and thus the line degenerates to the point (0,0,0).
	{

	}

	Ray(const Point3D& origin, const Vector3D& direction)
		: orig{ origin }, dir{ direction }
	{

	}

	// Getters:

	Point3D origin() const
	{
		return orig;
	}

	Vector3D direction() const
	{
		return dir;
	}

	Point3D at(double t) const		// The point P(t) = A + tB, where A is the origin of the ray and B is the direction of the line.
	{
		return orig + t * dir;
	}
};

#endif // !RAY_H
```

- Create the background scene

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  main file for 8599 ray tracer
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

// Note: PPM image can be viewed by **Portable Anymap Viewer** on Windows

#include <iostream>
#include "Vector3D.h"
#include "color.h"
#include "Ray.h"

// Currently, we want a blue-to-white gradient background.

ColorRGB ray_color(const Ray& ray)		// this currently returns the color of the background in the direction of the ray
{
	Vector3D unit_direction = unit_vector(ray.direction());
	// Now, unit_direction.y() is between [-1,1], we normalize this range to [0,1] as follows:
	double height_weighting = 0.5 * (unit_direction.y() + 1.0);
	// Apply linear interpolation to blend white (at the bottom) and blue (at the top):
	// Note: the interpolation is NOT linear on the viewport.
	return (1 - height_weighting) * ColorRGB { 1.0, 1.0, 1.0 } + height_weighting * ColorRGB{ 0.5,0.7,1.0 };
}

int main()
{
	// Parameters of output image:

	const double aspect_ratio = 16.0 / 9.0;		// x/y
	const int image_width = 400;
	const int image_height = int(image_width / aspect_ratio);	// ??? use static_cast<int>()?

	// Color Settings:

	const int max_color = 255;

	// Camera & Viewport Settings:
	// Note: the point on the viewport plane is assumed to be at the centre of the corresponding pixel on the final image.

	double viewport_height = 2.0;
	double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
	double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane).
	Point3D origin{ 0.0,0.0,0.0 };	// where camera locates.
	Vector3D horizontal{ viewport_width, 0.0,0.0 };		// for calculating the left-to-right offset of the endpoint on the viewport
	Vector3D vertical{ 0.0,viewport_height,0.0 };		// for calculating the bottom-to-top offset of the endpoint on the viewport
	Point3D bottom_left = origin - Vector3D{ 0.0,0.0,focal_length } - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint


	// Output Data:
	// (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)

	std::cout << "P3" << '\n'								// colors are in ASCII		(??? Explain the meaning)
		<< image_width << ' ' << image_height << '\n'		// column  row
		<< max_color << '\n';								// value for max color

	// RGB triplets: (each rendered as a pixel, from left to right, top to bottom)

	for (int row = image_height - 1; row >= 0; row--)
	{
		std::cerr << '\r' << "Scanlines Remaining: " << row << ' ' << std::flush;		// ??? Why do we want std::flush here?
																						// Note: \r means writing from the head of the current line

		for (int column = 0; column < image_width; column++)
		{
			Vector3D horizontal_offset = (double(column) / (image_width - 1)) * horizontal;
			Vector3D vertical_offset = (double(row) / (image_height - 1)) * vertical;
			Ray ray{ origin, bottom_left + horizontal_offset + vertical_offset - origin };
			ColorRGB pixel_color = ray_color(ray);
			write_color(std::cout, pixel_color);
		}
	}
	std::cerr << '\n'
			  << "Done."
			  << '\n';
}
```

The output is as follows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/GradientBackground.jpg" width="600" height="600"></a>

### April 27th 2023

- Render a ray-traced sphere

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  main file for 8599 ray tracer
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

// Note: PPM image can be viewed by **Portable Anymap Viewer** on Windows

#include <iostream>
#include "Vector3D.h"
#include "color.h"
#include "Ray.h"

bool is_hitting_sphere(const Point3D& center, double radius, const Ray& ray)	// center/radius of the sphere
{
	// hard-coded to determine the number of root(s) of the qudratic equation which represents how the ray intersects with the sphere:
	Vector3D sphere_factor = ray.origin() - center;
	double A = dot(ray.direction(), ray.direction());
	double B = 2.0 * dot(ray.direction(), sphere_factor);
	double C = dot(sphere_factor, sphere_factor) - radius * radius;
	double discriminant = B * B - 4 * A * C;

	return (discriminant > 0);		// ??? should we include == 0 case?
}

// 

ColorRGB ray_color(const Ray& ray)		// this currently returns the color of what the ray directly hits (the sphere or the background)
{
	// case: hitting a sphere
	
	if (is_hitting_sphere({ 0.0,0.0,-1.0 }, 0.5, ray))		// we place a sphere with radius == 0.5 on z-axis where its center is on the viewport plane.
	{
		return ColorRGB{ 1.0,1.0,0.0 };		// a sphere in yellow
	}

	// case: hitting the background (currently, we want a blue-to-white gradient background)
	
	Vector3D unit_direction = unit_vector(ray.direction());
	// Now, unit_direction.y() is between [-1,1], we normalize this range to [0,1] as follows:
	double height_weighting = 0.5 * (unit_direction.y() + 1.0);
	// Apply linear interpolation to blend white (at the bottom) and blue (at the top):
	// Note: the interpolation is NOT linear on the viewport.
	return (1 - height_weighting) * ColorRGB { 1.0, 1.0, 1.0 } + height_weighting * ColorRGB{ 0.5,0.7,1.0 };
}

int main()
{
	// Parameters of output image:

	const double aspect_ratio = 16.0 / 9.0;		// x/y
	const int image_width = 400;
	const int image_height = int(image_width / aspect_ratio);	// ??? use static_cast<int>()?

	// Color Settings:

	const int max_color = 255;

	// Camera & Viewport Settings:
	// Note: the point on the viewport plane is assumed to be at the centre of the corresponding pixel on the final image.

	double viewport_height = 2.0;
	double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
	double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane).
	Point3D origin{ 0.0,0.0,0.0 };	// where camera locates.
	Vector3D horizontal{ viewport_width, 0.0,0.0 };		// for calculating the left-to-right offset of the endpoint on the viewport
	Vector3D vertical{ 0.0,viewport_height,0.0 };		// for calculating the bottom-to-top offset of the endpoint on the viewport
	Point3D bottom_left = origin - Vector3D{ 0.0,0.0,focal_length } - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint


	// Output Data:
	// (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)

	std::cout << "P3" << '\n'								// colors are in ASCII		(??? Explain the meaning)
		<< image_width << ' ' << image_height << '\n'		// column  row
		<< max_color << '\n';								// value for max color

	// RGB triplets: (each rendered as a pixel, from left to right, top to bottom)

	for (int row = image_height - 1; row >= 0; row--)
	{
		std::cerr << '\r' << "Scanlines Remaining: " << row << ' ' << std::flush;		// ??? Why do we want std::flush here?
																						// Note: \r means writing from the head of the current line

		for (int column = 0; column < image_width; column++)
		{
			Vector3D horizontal_offset = (double(column) / (image_width - 1)) * horizontal;
			Vector3D vertical_offset = (double(row) / (image_height - 1)) * vertical;
			Ray ray{ origin, bottom_left + horizontal_offset + vertical_offset - origin };
			ColorRGB pixel_color = ray_color(ray);
			write_color(std::cout, pixel_color);
		}
	}
	std::cerr << '\n'
			  << "Done."
			  << '\n';
}
```

The output is as follows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/FirstSphere.jpg" width="700" height="600"></a>

- Shade the sphere with its surface normals

We will map $rgb$ linearly from the unit normal vector $(n_x, n_y, n_z)$.

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  main file for 8599 ray tracer
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

// Note: PPM image can be viewed by **Portable Anymap Viewer** on Windows

#include <iostream>
#include "Vector3D.h"
#include "color.h"
#include "Ray.h"

#include <chrono>

double sphere_hitting_index(const Point3D& center, double radius, const Ray& ray)	// center/radius of the sphere
{
	// hard-coded to determine the number of root(s) of the qudratic equation which represents how the ray intersects with the sphere:
	Vector3D sphere_factor = ray.origin() - center;
	double A = dot(ray.direction(), ray.direction());
	double B = 2.0 * dot(ray.direction(), sphere_factor);
	double C = dot(sphere_factor, sphere_factor) - radius * radius;
	double discriminant = B * B - 4 * A * C;

	if (discriminant < 0.0)
	{
		return -1.0;
	}

	return ((-B) - sqrt(discriminant)) / (2.0 * A);		// Note: -sqrt() will always result the smaller root, which reveals the "front" (here means larger z-value) point of intersection.
}

// 

ColorRGB ray_color(const Ray& ray)		// this currently returns the color of what the ray directly hits (the sphere or the background)
{
	const Point3D sphere_center{ 0.0,0.0,-1.0 };

	double hitting_index = sphere_hitting_index(sphere_center, 0.5, ray);	// we place a sphere with radius == 0.5 on z-axis where its center is on the viewport plane.

	// case: hitting a sphere
	
	if (hitting_index > 0.0)	// Note that this means we only consider cases where the "front" (here means larger z-value) intersection exists AND is strictly in front of the camera
	{
		Vector3D surface_normal = unit_vector(ray.at(hitting_index) - sphere_center);
		
		return 0.5 * ColorRGB{ surface_normal.x() + 1, surface_normal.y() + 1, surface_normal.z() + 1 };
	}

	// case: hitting the background (currently, we want a blue-to-white gradient background)
	// Note: currently, if the sphere is too big so that the camera is on its surface or in it, then this program will not display it.
	
	Vector3D unit_direction = unit_vector(ray.direction());
	// Now, unit_direction.y() is between [-1,1], we normalize this range to [0,1] as follows:
	double height_weighting = 0.5 * (unit_direction.y() + 1.0);
	// Apply linear interpolation to blend white (at the bottom) and blue (at the top):
	return (1 - height_weighting) * ColorRGB { 1.0, 1.0, 1.0 } + height_weighting * ColorRGB{ 0.5,0.7,1.0 };
}

int main()
{
	// Parameters of output image:

	const double aspect_ratio = 16.0 / 9.0;		// x/y
	const int image_width = 400;
	const int image_height = int(image_width / aspect_ratio);	// ??? use static_cast<int>()?

	// Color Settings:

	const int max_color = 255;

	// Camera & Viewport Settings:
	// Note: the point on the viewport plane is assumed to be at the centre of the corresponding pixel on the final image.

	double viewport_height = 2.0;
	double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
	double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane).
	Point3D origin{ 0.0,0.0,0.0 };	// where camera locates.
	Vector3D horizontal{ viewport_width, 0.0,0.0 };		// for calculating the left-to-right offset of the endpoint on the viewport
	Vector3D vertical{ 0.0,viewport_height,0.0 };		// for calculating the bottom-to-top offset of the endpoint on the viewport
	Point3D bottom_left = origin - Vector3D{ 0.0,0.0,focal_length } - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint


	// Output Data:
	// (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)

	std::cout << "P3" << '\n'								// colors are in ASCII		(??? Explain the meaning)
		<< image_width << ' ' << image_height << '\n'		// column  row
		<< max_color << '\n';								// value for max color

	// benchmark
	auto start = std::chrono::high_resolution_clock::now();

	// RGB triplets: (each rendered as a pixel, from left to right, top to bottom)

	for (int row = image_height - 1; row >= 0; row--)
	{
		std::cerr << '\r' << "Scanlines Remaining: " << row << ' ' << std::flush;		// ??? Why do we want std::flush here?
																						// Note: \r means writing from the head of the current line

		for (int column = 0; column < image_width; column++)
		{
			Vector3D horizontal_offset = (double(column) / (image_width - 1)) * horizontal;
			Vector3D vertical_offset = (double(row) / (image_height - 1)) * vertical;
			Ray ray{ origin, bottom_left + horizontal_offset + vertical_offset - origin };
			ColorRGB pixel_color = ray_color(ray);
			write_color(std::cout, pixel_color);
		}
	}

	// benchmark
	auto end = std::chrono::high_resolution_clock::now();

	std::cerr << '\n'
			  << "Done."
			  << '\n';

	// benchmark
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cerr << "It took " << elapsed.count() << " milliseconds." << std::endl;
}
```

The output is as follows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/sphere_with_surface_normal.jpg" width="800" height="600"></a>

### April 28th 2023

- An abstract class (`Hittable`) representing any objects in the world that is hittable:

We add an extra control of "visibility": only hits that are in $[t_{min}, t_{max}]$ are considered. (Recall the ray is represented as $P(t)=A+tB$)

```cpp
/*****************************************************************//**
 * \file   Hittable.h
 * \brief  The abstract class for anything that can be hitted by the ray
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef HITTABLE_H
#define HITTABLE_H

#include "Ray.h"

struct HitRecord		// ??? How is HitRecord useful as it does not record the information about the ray?
{
	Point3D point;
	Vector3D normal;
	double t;
};

class Hittable
{
public:
	virtual bool is_hit_by(const Ray& ray, double t_min, double t_max, HitRecord& record) const = 0;
};

#endif // !HITTABLE_H
```

- Now we can drive `Sphere` class from `Hittable`:

```cpp
/*****************************************************************//**
 * \file   Sphere.h
 * \brief  Class declaration for ray-hittable sphere
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef SPHERE_H
#define SPHERE_H

#include "Hittable.h"

/*
C++ side note:
If we don't provide class access specifier in inheritance,
it would be public if the sub-class is declared as struct and be private if the sub-class is declared as class
*/

class Sphere : public Hittable
{
	Point3D center;
	double radius;

public:

	Sphere()	// currently, radius is uninitialized! what default value should we set?
	{

	}

	Sphere(Point3D centre, double r)
		: center{ centre }, radius{r}
	{

	}

	virtual bool is_hit_by(const Ray& ray, double t_min, double t_max, HitRecord& record) const override;
};

#endif // !SPHERE_H
```

```cpp
/*****************************************************************//**
 * \file   Sphere.cpp
 * \brief  Class definition for ray-hittable sphere
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#include "Sphere.h"

// Side-note: virtual keyword is not needed AND not allowed to be outside a class declaration, and thus for a virtual function defintion outside the class we also can declare it with override

bool Sphere::is_hit_by(const Ray& ray, double t_min, double t_max, HitRecord& record) const
{
	Vector3D sphere_factor = ray.origin() - center;
	double A = ray.direction().squared_length();
	double half_B = dot(ray.direction(), sphere_factor);
	double C = sphere_factor.squared_length() - radius * radius;
	double discriminant = half_B * half_B - A * C;

	if (discriminant < 0.0)		// at this stage, discriminant == 0 case isn't rejected.
	{
		return false;
	}

	double dist = std::sqrt(discriminant);
	double root = ((-half_B) - dist) / A;	// the "near" root
	if (root < t_min || t_max < root)
	{
		root = ((-half_B) + dist) / A;		// the "far" root
		if (root < t_min || t_max < root)	// in-range far root is ACCEPTED when near root is rejected
		{
			return false;
		}
	}
	// record the intersection:
	record.t = root;
	record.point = ray.at(root);
	record.normal = (record.point - center) / radius;

	return true;
}
```

### April 29th 2023

- Add the information of which side of the surface that the ray is coming from to `HitRecord`

```cpp
/*****************************************************************//**
 * \file   Hittable.h
 * \brief  The abstract class for anything that can be hitted by the ray
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef HITTABLE_H
#define HITTABLE_H

#include "Ray.h"

struct HitRecord		// ??? How is HitRecord useful as it does not record the information about the ray?
{
	Point3D point;
	Vector3D normal;
	double t;
	bool is_hitting_front_face;

	// Assume for any geometric entity rendered by this ray tracer, there is a defined front face (and thus a back face) and thus a front (i.e. outward) normal. (??? What about Klein bottle/ring?)
	// We choose to set the normal to always towards the ray:
	inline void set_normal(const Ray& ray, const Vector3D& outward_normal)		// Setting the normal and the front-back bool
	// ??? Isn't in-class method implicitly inlined? Can we delete the explicit inline keyword?
	{
		is_hitting_front_face = (dot(ray.direction(), outward_normal) < 0.0);	// NOTE: Hitting at right angle is counted as hitting from outside
		normal = is_hitting_front_face ? (outward_normal) : (-outward_normal);
	}
};

class Hittable
{
public:
	virtual bool is_hit_by(const Ray& ray, double t_min, double t_max, HitRecord& record) const = 0;
};

#endif // !HITTABLE_H
```

Hence, the `is_hit_by` method of the `Sphere` class becomes:

```cpp
/*****************************************************************//**
 * \file   Sphere.cpp
 * \brief  Class definition for ray-hittable sphere
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#include "Sphere.h"

// Side-note: virtual keyword is not needed AND not allowed to be outside a class declaration, and thus for a virtual function defintion outside the class we also can declare it with override

bool Sphere::is_hit_by(const Ray& ray, double t_min, double t_max, HitRecord& record) const
{
	Vector3D sphere_factor = ray.origin() - center;
	double A = ray.direction().squared_length();
	double half_B = dot(ray.direction(), sphere_factor);
	double C = sphere_factor.squared_length() - radius * radius;
	double discriminant = half_B * half_B - A * C;

	if (discriminant < 0.0)		// at this stage, discriminant == 0 case isn't rejected.
	{
		return false;
	}

	double dist = std::sqrt(discriminant);
	double root = ((-half_B) - dist) / A;	// the "near" root
	if (root < t_min || t_max < root)
	{
		root = ((-half_B) + dist) / A;		// the "far" root
		if (root < t_min || t_max < root)	// in-range far root is ACCEPTED when near root is rejected
		{
			return false;
		}
	}
	// record the intersection:
	record.t = root;
	record.point = ray.at(root);
	record.set_normal(ray, (record.point - center) / radius);

	return true;
}
```

- Create a class `CompositeHittable` to represent an entity that is hittable by the ray and is consist of multiple `Hittable` objects

Since `CompositeHittable` is hittable, it's natural to make it a `Hittable`:

```cpp
/*****************************************************************//**
 * \file   CompositeHittable.h
 * \brief  The declaration of the class representing hittable that is made by hittables
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef COMPOSITEHITTABLE_H
#define COMPOSITEHITTABLE_H

#include <memory>	// to use smart pointers
#include <vector>
#include "Hittable.h"

/*
C++ side notes on shared pointers:
	
	std::shared_ptr<T> is an pointer encapsulated with reference-counting semantics such that it is automatically deleted when the counting is 0.
	
	shared_ptr<T> sptr = make_shared<T>(T_constructor_params ...);	// make_shared allocates the object on the heap

	Moreover, there is no way to manage a stack allocated object with a shared pointer!
*/

class CompositeHittable : public Hittable
{
	std::vector<std::shared_ptr<Hittable>> components;
	
public:

	// Constructors:

	CompositeHittable()
	{

	}

	CompositeHittable(std::shared_ptr<Hittable> hittable_object)
	{
		components.push_back(hittable_object);
	}

	// Methods:

	void clear()
	{
		components.clear();
	}

	void add(std::shared_ptr<Hittable> hittable_object)
	{
		components.push_back(hittable_object);
	}

	virtual bool is_hit_by(const Ray& ray, double t_min, double t_max, HitRecord& record) const override;
};

#endif // !COMPOSITEHITTABLE_H
```

```cpp
/*****************************************************************//**
 * \file   CompositeHittable.cpp
 * \brief  The definitions of the class representing hittable that is made by hittables
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#include "CompositeHittable.h"

bool CompositeHittable::is_hit_by(const Ray& ray, double t_min, double t_max, HitRecord& record) const
{
	HitRecord component_record;
	bool is_hit = false;
	double visibility = t_max;

	for (const auto& component : components)
	{
		if (component->is_hit_by(ray, t_min, visibility, component_record))
		{
			is_hit = true;
			visibility = component_record.t;
			record = component_record;
		}
	}

	return is_hit;
}
```

- Create the main header for common utility

```cpp
/*****************************************************************//**
 * \file   RayTracingToolbox.h
 * \brief  Numeric constants, utility functions and common headers for the ray tracer
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef RAYTRACINGTOOLBOX_H
#define RAYTRACINGTOOLBOX_H

#include <cmath>
#include <limits>
#include <memory>

// Numeric Constants:

const double positive_infinity = std::numeric_limits<double>::infinity();	// (std::numeric_limits<double>::max() < std::numeric_limits<double>::infinity()) == true
const double pi = 3.141592653589793'2385;	// ??? Is the precision of a double able to store the digits after 3.141592653589793?

// Utility Functions:

inline double degrees_to_radians(double degrees)
{
	return (degrees / 180.0) * pi;
}

// Common Headers:

#include "Vector3D.h"
#include "Ray.h"


#endif // !RAYTRACINGTOOLBOX_H
```

- Currently, our `main.cpp` becomes:

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  The renderer for 8599 ray tracer
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

// Note: PPM image can be viewed by **Portable Anymap Viewer** on Windows

#include <iostream>
#include <chrono>

#include "RayTracingToolbox.h"

#include "color.h"
#include "CompositeHittable.h"
#include "Sphere.h"


ColorRGB ray_color(const Ray& ray, const Hittable& world)		// this currently returns the color of what the ray directly hits (the sphere or the background)
{
	HitRecord record;

	// Hitting something: render the normal towards the ray (outwards if the front/back normal is at right angle to the ray)
	if (world.is_hit_by(ray, 0, positive_infinity, record))
	{
		return 0.5 * (record.normal + ColorRGB{ 1.0,1.0,1.0 });		// Reminder: normal is unit vector
	}

	// Not hitting anything: render the sky
	double interpolation_factor = 0.5 * (unit_vector(ray.direction()).y() + 1.0);	// Normalized to [0,1]
	return (1.0 - interpolation_factor) * ColorRGB { 1.0, 1.0, 1.0 } + interpolation_factor * ColorRGB{ 0.5,0.7,1.0 };
}

int main()
{
	// Parameters of output image:
	const double aspect_ratio = 16.0 / 9.0;		// x/y
	const int image_width = 400;
	const int image_height = int(image_width / aspect_ratio);	// ??? use static_cast<int>()?

	// Color Settings:
	const int max_color = 255;

	// Creating the (objects in the) world:
	CompositeHittable world;	// empty world
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, -100.5, -1.0 }, 100.0));	// add the ground
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, 0.0, -1.0 }, 0.5));	// add a ball on the ground

	// Camera & Viewport Settings:
	// Note: the point on the viewport plane is assumed to be at the centre of the corresponding pixel on the final image.
	double viewport_height = 2.0;
	double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
	double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane).
	Point3D origin{ 0.0,0.0,0.0 };	// where camera locates.
	Vector3D horizontal{ viewport_width, 0.0,0.0 };		// for calculating the left-to-right offset of the endpoint on the viewport
	Vector3D vertical{ 0.0,viewport_height,0.0 };		// for calculating the bottom-to-top offset of the endpoint on the viewport
	Point3D bottom_left = origin - Vector3D{ 0.0,0.0,focal_length } - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint


	// Rendering (i.e. output data):
	// (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)
	std::cout << "P3" << '\n'								// colors are in ASCII		(??? Explain the meaning)
		<< image_width << ' ' << image_height << '\n'		// column  row
		<< max_color << '\n';								// value for max color
	// benchmark
	auto start = std::chrono::high_resolution_clock::now();
	// RGB triplets: (each rendered as a pixel, from left to right, top to bottom)

	for (int row = image_height - 1; row >= 0; row--)
	{
		std::cerr << '\r' << "Scanlines Remaining: " << row << ' ' << std::flush;		// ??? Why do we want std::flush here?
																						// Note: \r means writing from the head of the current line

		for (int column = 0; column < image_width; column++)
		{
			Vector3D horizontal_offset = (double(column) / (image_width - 1)) * horizontal;
			Vector3D vertical_offset = (double(row) / (image_height - 1)) * vertical;
			Ray ray{ origin, bottom_left + horizontal_offset + vertical_offset - origin };
			ColorRGB pixel_color = ray_color(ray, world);
			write_color(std::cout, pixel_color);
		}
	}
	// benchmark
	auto end = std::chrono::high_resolution_clock::now();
	std::cerr << '\n'
			  << "Done."
			  << '\n';
	// benchmark
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cerr << "It took " << elapsed.count() << " milliseconds." << std::endl;
}
```

which renders the following image:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/surface_normal_with_ground.jpg" width="640" height="600"></a>

### April 30th 2023

- The first attempt to multi-threading:

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  The renderer for 8599 ray tracer
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

 // Note: PPM image can be viewed by **Portable Anymap Viewer** on Windows

#include <iostream>
#include <vector>
#include <chrono>		// for benchmark
#include <execution>	// for multi-threading

#include "RayTracingToolbox.h"

#include "color.h"
#include "CompositeHittable.h"
#include "Sphere.h"


ColorRGB ray_color(const Ray& ray, const Hittable& world)		// this currently returns the color of what the ray directly hits (the sphere or the background)
{
	HitRecord record;

	// Hitting something: render the normal towards the ray (outwards if the front/back normal is at right angle to the ray)
	if (world.is_hit_by(ray, 0, positive_infinity, record))
	{
		return 0.5 * (record.normal + ColorRGB{ 1.0,1.0,1.0 });		// Reminder: normal is unit vector
	}

	// Not hitting anything: render the sky
	double interpolation_factor = 0.5 * (unit_vector(ray.direction()).y() + 1.0);	// Normalized to [0,1]
	return (1.0 - interpolation_factor) * ColorRGB { 1.0, 1.0, 1.0 } + interpolation_factor * ColorRGB{ 0.5,0.7,1.0 };
}

int main()
{
	// Parameters of output image:
	const double aspect_ratio = 16.0 / 9.0;		// x/y
	const int image_width = 400;
	const int image_height = int(image_width / aspect_ratio);	// ??? use static_cast<int>()?

	// Color Settings:
	const int max_color = 255;

	// Creating the (objects in the) world:
	CompositeHittable world;	// empty world
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, -100.5, -1.0 }, 100.0));	// add the ground
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, 0.0, -1.0 }, 0.5));	// add a ball on the ground

	// Camera & Viewport Settings:
	// Note: the point on the viewport plane is assumed to be at the centre of the corresponding pixel on the final image.
	double viewport_height = 2.0;
	double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
	double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane).
	Point3D origin{ 0.0,0.0,0.0 };	// where camera locates.
	Vector3D horizontal{ viewport_width, 0.0,0.0 };		// for calculating the left-to-right offset of the endpoint on the viewport
	Vector3D vertical{ 0.0,viewport_height,0.0 };		// for calculating the bottom-to-top offset of the endpoint on the viewport
	Point3D bottom_left = origin - Vector3D{ 0.0,0.0,focal_length } - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint


	// Rendering (i.e. output data):
	// (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)
	std::cout << "P3" << '\n'								// colors are in ASCII		(??? Explain the meaning)
		<< image_width << ' ' << image_height << '\n'		// column  row
		<< max_color << '\n';								// value for max color

	// Preparations for multi-threading:
	std::vector<std::vector<ColorRGB>> image;
	image.resize(image_height);
	for (auto& row : image)
	{
		row.resize(image_width);
	}
	std::vector<int> rows(image_height);
	std::vector<int> columns(image_width);
	for (int i = 0; i < image_height; i++)
	{
		rows[i] = image_height - 1 - i;
	}
	for (int j = 0; j < image_width; j++)
	{
		columns[j] = j;
	}

	// benchmark
	auto start = std::chrono::high_resolution_clock::now();
	// RGB triplets: (For PPM format: each rgb triplet is rendered as a pixel, from left to right, top to bottom)
	// Multi-threading:
	std::for_each(std::execution::par, rows.begin(), rows.end(),
		[&](int row)
		{
			std::for_each(std::execution::par, columns.begin(), columns.end(),
			[&](int column)
				{
					Vector3D horizontal_offset = (double(column) / (image_width - 1)) * horizontal;
					Vector3D vertical_offset = (double(row) / (image_height - 1)) * vertical;
					Ray ray{ origin, bottom_left + horizontal_offset + vertical_offset - origin };
					ColorRGB pixel_color = ray_color(ray, world);
					image[image_height - 1 - row][column] = pixel_color;
				}
			);
		}
	);
	//// Single threading:
	//for (int row = image_height - 1; row >= 0; row--)
	//{
	//	std::cerr << '\r' << "Scanlines Remaining: " << row << ' ' << std::flush;		// ??? Why do we want std::flush here?
	//	// Note: \r means writing from the head of the current line
	//
	//	for (int column = 0; column < image_width; column++)
	//	{
	//		Vector3D horizontal_offset = (double(column) / (image_width - 1)) * horizontal;
	//		Vector3D vertical_offset = (double(row) / (image_height - 1)) * vertical;
	//		Ray ray{ origin, bottom_left + horizontal_offset + vertical_offset - origin };
	//		ColorRGB pixel_color = ray_color(ray, world);
	//		image[image_height - 1 - row][column] = pixel_color;
	//	}
	//}
	for (const auto& row : image)
	{
		for (const auto& pixel_color : row)
		{
			write_color(std::cout, pixel_color);
		}
	}
	// benchmark
	auto end = std::chrono::high_resolution_clock::now();
	std::cerr << '\n'
		<< "Done."
		<< '\n';
	// benchmark
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cerr << "It took " << elapsed.count() << " milliseconds." << std::endl;
}
```

The results are as follows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/FirstMultithreadingBenchmark.jpg" width="640" height="240"></a>

- Create utility functions for
  - generating random reals
  - clamping reals

```cpp
/*****************************************************************//**
 * \file   RayTracingToolbox.h
 * \brief  Numeric constants, utility functions and common headers for the ray tracer
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef RAYTRACINGTOOLBOX_H
#define RAYTRACINGTOOLBOX_H

#include <cmath>
#include <limits>
#include <memory>
#include <cstdlib>
//#include <ctime>
#include <cassert>

// Numeric Constants:

const double positive_infinity = std::numeric_limits<double>::infinity();	// (std::numeric_limits<double>::max() < std::numeric_limits<double>::infinity()) == true
const double pi = 3.141592653589793'2385;	// ??? Is the precision of a double able to store the digits after 3.141592653589793?

// Utility Functions:

inline double degrees_to_radians(double degrees)
{
	return (degrees / 180.0) * pi;
}

inline double random_real_number()	// returns value in [0,1)
// Note that we want to exclude 1
{
	// TODO: ??? How to randomize the seed? (We probably don't want to use std::time since we are aiming for a real-time application)

	/*
	Note that the value of RAND_MAX is implementation defined.
	With this in mind, the following expression is designed to prevent incorrect return value due to integer overflow and the limited precision of floating point calculation.
	2147483648 is chosen because, while it is large enough, we have (2147483647 / (2147483647 + 1.0)) < 1 evaluated to true.
	*/
	return (std::rand() % 2147483648) / ((RAND_MAX < 2147483648) ? (RAND_MAX + 1.0) : (2147483647 + 1.0));
}

inline double random_real_number(double min, double max)	// returns value in [min,max)
{
	assert(min <= max);
	return min + (max - min) * random_real_number();
}

inline double clamp(double r, double min, double max)
{
	assert(min <= max);
	if (r < min)
	{
		return min;
	}
	if (r > max)
	{
		return max;
	}
	return r;
}

// Common Headers:

#include "Vector3D.h"
#include "Ray.h"


#endif // !RAYTRACINGTOOLBOX_H
```

### May 1st 2023

- Create the `Camera` class:

```cpp
/*****************************************************************//**
 * \file   Camera.h
 * \brief  The class for the camera
 * 
 * \author Xiaoyang Liu
 * \date   May 2023
 *********************************************************************/

#ifndef CAMERA_H
#define CAMERA_H

#include "RayTracingToolbox.h"

class Camera
{
	Point3D origin;			// where camera locates.
	Vector3D horizontal;	// for calculating the left-to-right offset of the endpoint on the viewport
	Vector3D vertical;		// for calculating the bottom-to-top offset of the endpoint on the viewport
	Point3D bottom_left;	// the bottom-left point on the viewpoint

public:

	// Constructors:

	Camera()
	{
		const double aspect_ratio = 16.0 / 9.0;		// x/y
		double viewport_height = 2.0;
		double viewport_width = viewport_height * aspect_ratio;		// viewport has the same aspect ratio as the image if the pixels on the display is square shaped.
		double focal_length = 1.0;		// this is the distance from the camera to the viewport (projection plane).
		origin = Point3D{ 0.0,0.0,0.0 };	// where camera locates.
		horizontal = Vector3D{ viewport_width, 0.0,0.0 };		// for calculating the left-to-right offset of the endpoint on the viewport
		vertical = Vector3D{ 0.0,viewport_height,0.0 };		// for calculating the bottom-to-top offset of the endpoint on the viewport
		bottom_left = origin - Vector3D{ 0.0,0.0,focal_length } - (horizontal / 2.0) - (vertical / 2.0);		// the bottom-left point on the viewpoint
	}

	// Methods:

	Ray extract_ray(double horizontal_offset_factor, double vertical_offset_factor)
	{
		return Ray{ origin, bottom_left + horizontal_offset_factor * horizontal + vertical_offset_factor * vertical - origin };
	}
};

#endif // !CAMERA_H
```

- Implementing basic MSAA:

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  The renderer for 8599 ray tracer
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

 /*
 Note:
	PPM image can be viewed by** Portable Anymap Viewer** on Windows.
 
 */  

// ---------------------------Control Panel---------------------------
#define Multithread 1
#define Antialiasing 1
// -------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <chrono>		// for benchmark
#include <execution>	// for multi-threading

#include "RayTracingToolbox.h"

#include "color.h"
#include "CompositeHittable.h"
#include "Sphere.h"
#include "Camera.h"


ColorRGB ray_color(const Ray& ray, const Hittable& world)		// this currently returns the color of what the ray directly hits (the sphere or the background)
{
	HitRecord record;

	// Hitting something: render the normal towards the ray (outwards if the front/back normal is at right angle to the ray)
	if (world.is_hit_by(ray, 0, positive_infinity, record))
	{
		return 0.5 * (record.normal + ColorRGB{ 1.0,1.0,1.0 });		// Reminder: normal is unit vector
	}

	// Not hitting anything: render the sky
	double interpolation_factor = 0.5 * (unit_vector(ray.direction()).y() + 1.0);	// Normalized to [0,1]
	return (1.0 - interpolation_factor) * ColorRGB { 1.0, 1.0, 1.0 } + interpolation_factor * ColorRGB{ 0.5,0.7,1.0 };
}

int main()
{
	// Parameters of output image:
	const double aspect_ratio = 16.0 / 9.0;		// x/y
	const int image_width = 400;
	const int image_height = int(image_width / aspect_ratio);	// ??? use static_cast<int>()?
#if Antialiasing
	const int samples_per_pixel = 100;		// for MSAA
#else
	const int samples_per_pixel = 1;		// for MSAA
#endif // Antialiasing


	// Color Settings:
	const int max_color = 255;

	// Creating the (objects in the) world:
	CompositeHittable world;	// empty world
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, -100.5, -1.0 }, 100.0));	// add the ground
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, 0.0, -1.0 }, 0.5));	// add a ball on the ground

	// Camera:
	Camera camera;

	// Rendering (i.e. output data):
	// (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)
	std::cout << "P3" << '\n'								// colors are in ASCII		(??? Explain the meaning)
		<< image_width << ' ' << image_height << '\n'		// column  row
		<< max_color << '\n';								// value for max color

	// Preparations for multi-threading:
	std::vector<std::vector<ColorRGB>> image;
	image.resize(image_height);
	for (auto& row : image)
	{
		row.resize(image_width);
	}
	std::vector<int> rows(image_height);
	std::vector<int> columns(image_width);
	for (int i = 0; i < image_height; i++)
	{
		rows[i] = image_height - 1 - i;
	}
	for (int j = 0; j < image_width; j++)
	{
		columns[j] = j;
	}

	// benchmark
	auto start = std::chrono::high_resolution_clock::now();
	// RGB triplets: (For PPM format: each rgb triplet is rendered as a pixel, from left to right, top to bottom)
#if Multithread
	// Multi-threading:
	std::for_each(std::execution::par, rows.begin(), rows.end(),
		[&](int row)
		{
			std::for_each(std::execution::par, columns.begin(), columns.end(),
			[&](int column)
				{
#if Antialiasing
					ColorRGB pixel_color;	// (0,0,0)
					for (int s = 0; s < samples_per_pixel; s++)
					{
						double horizontal_offset_factor = (column + random_real_number()) / image_width;
						double vertical_offset_factor = (row + random_real_number()) / image_height;
						Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
						pixel_color += ray_color(ray, world);
					}
					image[image_height - 1 - row][column] = pixel_color;
#else
					double horizontal_offset_factor = (column + 0.5) / image_width;
					double vertical_offset_factor = (row + 0.5) / image_height;
					Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
					ColorRGB pixel_color = ray_color(ray, world);
					image[image_height - 1 - row][column] = pixel_color;
#endif // Antialiasing
				}
			);
		}
	);
#else
	// Single threading:
	for (int row = image_height - 1; row >= 0; row--)
	{
		std::cerr << '\r' << "Scanlines Remaining: " << row << ' ' << std::flush;		// ??? Why do we want std::flush here?
		// Note: \r means writing from the head of the current line
	
		for (int column = 0; column < image_width; column++)
		{
#if Antialiasing
			ColorRGB pixel_color;	// (0,0,0)
			for (int s = 0; s < samples_per_pixel; s++)
			{
				double horizontal_offset_factor = (column + random_real_number()) / image_width;
				double vertical_offset_factor = (row + random_real_number()) / image_height;
				Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
				pixel_color += ray_color(ray, world);
			}
			image[image_height - 1 - row][column] = pixel_color;
#else
			double horizontal_offset_factor = (column + 0.5) / image_width;
			double vertical_offset_factor = (row + 0.5) / image_height;
			Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
			ColorRGB pixel_color = ray_color(ray, world);
			image[image_height - 1 - row][column] = pixel_color;
#endif // Antialiasing
		}
	}
#endif // Multithread
	// Output the pixel data:
	for (const auto& row : image)
	{
		for (const auto& pixel_color : row)
		{
			write_color(std::cout, pixel_color, samples_per_pixel);
		}
	}
	// benchmark
	auto end = std::chrono::high_resolution_clock::now();
	std::cerr << '\n'
		<< "Done."
		<< '\n';
	// benchmark
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cerr << "\nIt took " << elapsed.count() << " milliseconds.\n";
}
```
where the `write_color` function needs to account for the sampling rate:

```cpp
/*****************************************************************//**
 * \file   color.h
 * \brief  Utility Functions for outputting ColorRGB object as pixel
 * 
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef COLOR_H
#define COLOR_H

#include "Vector3D.h"

#include <iostream>

void write_color(std::ostream& os, ColorRGB pixel_color, int samples_per_pixel)
{
	// Assume the rgb component values of colorRGB for each sample is in range [0.0, 1.0], and the output integer value is in range [0, 255].
	// Assume the rgb component values of colorRGB is the sum of all the samples

	double r = pixel_color.x() / samples_per_pixel;
	double g = pixel_color.y() / samples_per_pixel;
	double b = pixel_color.z() / samples_per_pixel;

	os << int(256 * clamp(r, 0.0, 0.999)) << ' '
	   << int(256 * clamp(g, 0.0, 0.999)) << ' '
	   << int(256 * clamp(b, 0.0, 0.999)) << '\n';

	// ??? Explain why the 0.999 can be necessary.
	// ??? What is the difference if static_cast<int>() is used instead of int()?

}

#endif // !COLOR_H
```

For  `sampling rate == 1` (i.e. no antialiasing) where the ray shoots at the center of the pixel, the output shows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/MSAAx1.jpg" width="800" height="460"></a>

For  `sampling rate == 4`, the output shows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/MSAAx4.jpg" width="800" height="460"></a>

For  `sampling rate == 100`, the output shows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/MSAAx100.jpg" width="800" height="460"></a>

### May 2nd 2023

- A basic diffuse/matte material model:
  
  - takes incident light's color, reduce its brightness;
  - randomizes reflection direction.

// TODO: write notes here

```cpp
/*****************************************************************//**
 * \file   Vector3D.h
 * \brief  The class of 3D vector for representing geometry in R^3 and RGB color
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <cmath>		// to use std::sqrt
#include <cassert>
//#define NDEBUG		// uncomment this if we don't want assertion (e.g. when we want things like inf)
#include "RayTracingToolbox.h"

class Vector3D
{
	double v[3];

public:

	// Constructors:

	Vector3D()
		: v{ 0,0,0 }
	{

	}

	Vector3D(double x, double y, double z)
		: v{ x,y,z }
	{

	}

	// Operators:

	Vector3D operator-() const					// additive inverse
	{
		return Vector3D{ -v[0], -v[1], -v[2] };
	}

	double& operator[](int i)
	{
		assert(i == 0 || i == 1 || i == 2);
		return v[i];
	}

	double operator[](int i) const
	{
		assert(i == 0 || i == 1 || i == 2);
		return v[i];
	}

	Vector3D& operator+=(const Vector3D& u)
	{
		v[0] += u.v[0];
		v[1] += u.v[1];
		v[2] += u.v[2];

		return *this;
	}

	Vector3D& operator*=(const double d)
	{
		v[0] *= d;
		v[1] *= d;
		v[2] *= d;

		return *this;
	}

	Vector3D& operator/=(const double d)
	{
		assert(d != 0.0);

		v[0] /= d;
		v[1] /= d;
		v[2] /= d;

		return *this;
	}

	// Methods:

	double x() const
	{
		return v[0];
	}

	double y() const
	{
		return v[1];
	}

	double z() const
	{
		return v[2];
	}

	double squared_length() const
	{
		return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	}

	double length() const
	{
		return std::sqrt(squared_length());
	}

	/*
	C++ side notes: static method can be called without instantiation, (and thus) it cannot depend on any non-static data
	*/

	inline static Vector3D random()
	{
		return Vector3D{ random_real_number(), random_real_number(), random_real_number() };
	}

	inline static Vector3D random(double min, double max)
	{
		return Vector3D{ random_real_number(min,max), random_real_number(min,max), random_real_number(min,max) };
	}
};

// Unitility Functions:

inline std::ostream& operator<<(std::ostream& os, const Vector3D& v)
{
	return os << v.x() << ' ' << v.y() << ' ' << v.z();
}

inline Vector3D operator+(const Vector3D& a, const Vector3D& b)
{
	return Vector3D{ a.x() + b.x(), a.y() + b.y(), a.z() + b.z() };
}

inline Vector3D operator-(const Vector3D& a, const Vector3D& b)
{
	return Vector3D{ a.x() - b.x(), a.y() - b.y(), a.z() - b.z() };
}

inline Vector3D operator*(const Vector3D& a, const Vector3D& b)			// Note that this is NOT dot or cross product!
{
	return Vector3D{ a.x() * b.x(), a.y() * b.y(), a.z() * b.z() };
}

inline Vector3D operator*(double d, const Vector3D& v)
{
	return Vector3D{ d * v.x(), d * v.y(), d * v.z() };
}

inline Vector3D operator*(const Vector3D& v, double d)
{
	return d * v;
}

inline Vector3D operator/(const Vector3D& v, double d)
{
	return (1 / d) * v;
}

inline double dot(const Vector3D& a, const Vector3D& b)
{
	return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

inline Vector3D cross(const Vector3D& a, const Vector3D& b)
{
	return Vector3D
	{
		a.y() * b.z() - a.z() * b.y(),
		a.z() * b.x() - a.x() * b.z(),
		a.x() * b.y() - a.y() * b.x()
	};
}

inline Vector3D unit_vector(const Vector3D& v)
{
	return v / v.length();
}

inline Vector3D random_in_unit_sphere()
{
	while (true)
	{
		Vector3D p = Vector3D::random(-1, 1);
		if (p.squared_length() >= 1.0)
		{
			continue;
		}
		return p;
	}
}

// For better code readability (as Vector3D will represent things with different physical meanings):

using Point3D = Vector3D;
using ColorRGB = Vector3D;

#endif // !VECTOR3D_H
```

// TODO: write notes here

```cpp
/*****************************************************************//**
 * \file   main.cpp
 * \brief  The renderer for 8599 ray tracer
 *
 * \author Xiaoyang Liu
 * \date   April 2023
 *********************************************************************/

 /*
 Note:
	PPM image can be viewed by** Portable Anymap Viewer** on Windows.
 
 */  

// ---------------------------Control Panel---------------------------
#define Multithread 1
#define Antialiasing 1
// -------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <chrono>		// for benchmark
#include <execution>	// for multi-threading

#include "RayTracingToolbox.h"

#include "color.h"
#include "CompositeHittable.h"
#include "Sphere.h"
#include "Camera.h"


ColorRGB ray_color(const Ray& ray, const Hittable& world, const int bounce_depth)
{
	if (bounce_depth <= 0)
	{
		return ColorRGB{ 0.0,0.0,0.0 };		// representing no light
	}

	HitRecord record;

	if (world.is_hit_by(ray, 0, positive_infinity, record))
	{
		// Energy reduced by half for each bounce
		return 0.5 * ray_color(Ray{ record.point, (record.point + record.normal + random_in_unit_sphere()) - record.point }, world, bounce_depth - 1);
		// recall that record.normal is computed to be always on the reflection side
	}

	// Not hitting anything: render the sky
	double interpolation_factor = 0.5 * (unit_vector(ray.direction()).y() + 1.0);	// Normalized to [0,1]
	return (1.0 - interpolation_factor) * ColorRGB { 1.0, 1.0, 1.0 } + interpolation_factor * ColorRGB{ 0.5,0.7,1.0 };
}

int main()
{
	// Parameters of output image:
	const double aspect_ratio = 16.0 / 9.0;		// x/y
	const int image_width = 400;
	const int image_height = int(image_width / aspect_ratio);	// ??? use static_cast<int>()?
#if Antialiasing
	const int samples_per_pixel = 100;		// for MSAA
#else
	const int samples_per_pixel = 1;		// for MSAA
#endif // Antialiasing
	const int max_bounce_depth = 50;

	// Color Settings:
	const int max_color = 255;

	// Creating the (objects in the) world:
	CompositeHittable world;	// empty world
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, -100.5, -1.0 }, 100.0));	// add the ground
	world.add(std::make_shared<Sphere>(Point3D{ 0.0, 0.0, -1.0 }, 0.5));	// add a ball on the ground

	// Camera:
	Camera camera;

	// Rendering (i.e. output data):
	// (Note that by using > operator in Windows Command Prompt the contents of std::cout can be redirected to a file while the contents of std::cerr remains in the terminal)
	std::cout << "P3" << '\n'								// colors are in ASCII		(??? Explain the meaning)
		<< image_width << ' ' << image_height << '\n'		// column  row
		<< max_color << '\n';								// value for max color

	// Preparations for multi-threading:
	std::vector<std::vector<ColorRGB>> image;
	image.resize(image_height);
	for (auto& row : image)
	{
		row.resize(image_width);
	}
	std::vector<int> rows(image_height);
	std::vector<int> columns(image_width);
	for (int i = 0; i < image_height; i++)
	{
		rows[i] = image_height - 1 - i;
	}
	for (int j = 0; j < image_width; j++)
	{
		columns[j] = j;
	}

	// benchmark
	auto start = std::chrono::high_resolution_clock::now();
	// RGB triplets: (For PPM format: each rgb triplet is rendered as a pixel, from left to right, top to bottom)
#if Multithread
	// Multi-threading:
	std::for_each(std::execution::par, rows.begin(), rows.end(),
		[&](int row)
		{
			std::for_each(std::execution::par, columns.begin(), columns.end(),
			[&](int column)
				{
#if Antialiasing
					ColorRGB pixel_color;	// (0,0,0)
					for (int s = 0; s < samples_per_pixel; s++)
					{
						double horizontal_offset_factor = (column + random_real_number()) / image_width;
						double vertical_offset_factor = (row + random_real_number()) / image_height;
						Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
						pixel_color += ray_color(ray, world, max_bounce_depth);
					}
					image[image_height - 1 - row][column] = pixel_color;
#else
					double horizontal_offset_factor = (column + 0.5) / image_width;
					double vertical_offset_factor = (row + 0.5) / image_height;
					Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
					ColorRGB pixel_color = ray_color(ray, world, max_bounce_depth);
					image[image_height - 1 - row][column] = pixel_color;
#endif // Antialiasing
				}
			);
		}
	);
#else
	// Single threading:
	for (int row = image_height - 1; row >= 0; row--)
	{
		std::cerr << '\r' << "Scanlines Remaining: " << row << ' ' << std::flush;		// ??? Why do we want std::flush here?
		// Note: \r means writing from the head of the current line
	
		for (int column = 0; column < image_width; column++)
		{
#if Antialiasing
			ColorRGB pixel_color;	// (0,0,0)
			for (int s = 0; s < samples_per_pixel; s++)
			{
				double horizontal_offset_factor = (column + random_real_number()) / image_width;
				double vertical_offset_factor = (row + random_real_number()) / image_height;
				Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
				pixel_color += ray_color(ray, world, max_bounce_depth);
			}
			image[image_height - 1 - row][column] = pixel_color;
#else
			double horizontal_offset_factor = (column + 0.5) / image_width;
			double vertical_offset_factor = (row + 0.5) / image_height;
			Ray ray = camera.extract_ray(horizontal_offset_factor, vertical_offset_factor);
			ColorRGB pixel_color = ray_color(ray, world, max_bounce_depth);
			image[image_height - 1 - row][column] = pixel_color;
#endif // Antialiasing
		}
	}
#endif // Multithread
	// Output the pixel data:
	for (const auto& row : image)
	{
		for (const auto& pixel_color : row)
		{
			write_color(std::cout, pixel_color, samples_per_pixel);
		}
	}
	// benchmark
	auto end = std::chrono::high_resolution_clock::now();
	std::cerr << '\n'
		<< "Done."
		<< '\n';
	// benchmark
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cerr << "\nIt took " << elapsed.count() << " milliseconds.\n";
}
```

The output is as follows:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/G1.jpg" width="600" height="500"></a>

- Apply basic Gamma Correction

For `gamma == 2`:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/G2.jpg" width="540" height="300"></a>

For `gamma == 3`:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/G3.jpg" width="540" height="300"></a>

For `gamma == 4`:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/G4.jpg" width="540" height="300"></a>

### May 3rd 2023

- Fix shadow acne

// TODO: write notes here

If we starts intersection detection at `0.0` (i.e. with shadow ance): (`gamma == 2`)

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/WithShadowAcne.jpg" width="540" height="400"></a>

If we starts intersection detection at `0.001` (i.e. without shadow ance): (`gamma == 2`)

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/WithoutShadowAcne.jpg" width="540" height="400"></a>

Note that, since the ray "eaten" by the surface will be likely to bounce many times inside the sphere, heal shadow acne is also benefical to performance:

<img src="https://github.com/IQ404/8599-ray-tracer/blob/main/Sample%20Images/ShadowAcneImpactOnPerformance.jpg" width="540" height="300"></a>
