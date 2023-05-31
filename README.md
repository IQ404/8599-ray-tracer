# 8599 Ray Tracer

This repository aims to record a path tracer that I am currently writing for my MSc module CSC8599.

## Related Repositories

- [The offline prototype of the 8599 ray tracer](https://github.com/IQ404/8599-ray-tracer-prototype)
- [8599 ray tracer with the Walnut GUI](https://github.com/IQ404/8599-ray-tracer-gui)

## Objectives:

- Prototype an offline ray tracer.
  
  <ins>References</ins>: [the book by Peter Shirley](https://raytracing.github.io/books/RayTracingInOneWeekend.html)
  
  **Estimated date of completion:** Completed. :heavy_check_mark:

- Create an interactive GUI renderer framework, running multithreaded on CPU.

  <ins>References</ins>: [the series by The Cherno](https://www.youtube.com/watch?v=gfW1Fhd9u9Q&list=PLlrATfBNZ98edc5GshdBtREv5asFW3yXl&index=1), [Walnut GUI by The Cherno](https://github.com/TheCherno/Walnut), [Futhark](https://github.com/athas/raytracinginoneweekendinfuthark), [boksajak](https://github.com/boksajak/raytracingthenextweek)
  
  **Estimated date of completion:** Completed. :heavy_check_mark:

- Integrate the offline ray tracer into the interactive GUI renderer framework. Implementing:

  - non-physical path tracing
  - Whitted-Style ray tracing
  - physically-based path tracing (i.e. take rendering equation & BRDF into account)

  <ins>References</ins>: [GAMES101](https://sites.cs.ucsb.edu/~lingqi/teaching/games101.html) (this is the chinese version of [cs180](https://sites.cs.ucsb.edu/~lingqi/teaching/cs180.html))

  **Estimated date of completion:** <ins>15th June 2023</ins>
  
- Improve the ray-tracer by adding more advanced features to it (currently planning to focus on real time denoising).
  
  <ins>Possible references</ins>: [the next book by Peter Shirley](https://raytracing.github.io/books/RayTracingTheNextWeek.html), [GAMES202](https://sites.cs.ucsb.edu/~lingqi/teaching/games202.html) (this is the chinese version of [CS292F](https://sites.cs.ucsb.edu/~lingqi/teaching/cs292f.html))
  
  **Estimated date of completion:** by the <ins>end of July 2023</ins>

- If I still have time, I will try to accelerate the path tracer by running code on GPU.
  
  <ins>Possible references</ins>: [CUDA](https://developer.nvidia.com/blog/accelerated-ray-tracing-cuda/), [DXR](https://github.com/theroyn/RealTimeRayTracing), [OptiX](https://developer.nvidia.com/rtx/ray-tracing/optix), [Vulkan](https://github.com/GPSnoopy/RayTracingInVulkan)
  
  **Estimated date of completion:** <ins>22nd August 2023</ins> which is the deadline of this module (CSC8599).

## Housekeeping:

- When developing this project, the code is running on <ins>11980HK</ins> and <ins>RTX3080 (laptop version)</ins>.
- The code is written in Visual Studio 2022 with C++20.
- This implementation uses RIGHT-handed coordinate systems.

## Progression

- [An offine ray tracer prototype](https://github.com/IQ404/8599-ray-tracer-prototype/blob/main/README.md) (Finished)
- [Integrating the ray tracer with the Walnut GUI](https://github.com/IQ404/8599-ray-tracer-with-walnut-gui/blob/master/README.md) (Ongoing...)
