# 8599 Ray Tracer

This repository aims to record a path tracer that I am currently writing for my MSc module CSC8599.

## Related Repositories

- [8599-ray-tracer-prototype](https://github.com/IQ404/8599-ray-tracer-prototype)
- [8599-ray-tracer-with-walnut-gui](https://github.com/IQ404/8599-ray-tracer-with-walnut-gui)

## Objectives:

- Prototype an offine ray tracer.
  
  <ins>References</ins>: [the book by Peter Shirley](https://raytracing.github.io/books/RayTracingInOneWeekend.html)
  
  **Estimated date of completion:** Completed. :heavy_check_mark:

- Optimize the path tracer (from the CPU side, e.g. multithreading) and integrate its output into an interactive GUI. Try to achieve >= 10 fps for a scene that contains a relatively small number of items, a relatively low sampling rate and a relatively low resolution, so that we don't need to sacrifice the advanced features as much as possible.

  <ins>Possible references</ins>: [the series by The Cherno](https://www.youtube.com/watch?v=gfW1Fhd9u9Q&list=PLlrATfBNZ98edc5GshdBtREv5asFW3yXl&index=1), [Walnut GUI by The Cherno](https://github.com/TheCherno/Walnut), [Futhark](https://github.com/athas/raytracinginoneweekendinfuthark), [boksajak](https://github.com/boksajak/raytracingthenextweek)
  
  **Estimated date of completion:** optimization should be ongoing indefinitely. But the first complete version of the ray tracer running in real-time is aiming to be finished by the <ins>end of May 2023</ins>

- Improve the ray-tracer by adding more advanced features to it.
  
  <ins>Possible references</ins>: [the next book by Peter Shirley](https://raytracing.github.io/books/RayTracingTheNextWeek.html), [GAMES101](https://sites.cs.ucsb.edu/~lingqi/teaching/games101.html) (this is the chinese version of [cs180](https://sites.cs.ucsb.edu/~lingqi/teaching/cs180.html)), [GAMES202](https://sites.cs.ucsb.edu/~lingqi/teaching/games202.html) (this is the chinese version of [CS292F](https://sites.cs.ucsb.edu/~lingqi/teaching/cs292f.html))
  
  **Estimated date of completion:** I'm still thinking about what "advanced" features am I aiming for. Say, if it's features like motion blur, it should not take too long, but if it's modern techniques on denoising, it may then take me relatively a long time to read and implement recent papers. <ins>Currently I prefer the latter because that is more crucial to my PhD proposal. If so, it is possible that I may not have time to do GPU acceleration for this ray tracer and I will postpone learning vulkan until I start my PhD.</ins>

- If time allows, try to accelerate the path tracer by running code on GPU.
  
  <ins>Possible references</ins>: [CUDA](https://developer.nvidia.com/blog/accelerated-ray-tracing-cuda/), [DXR](https://github.com/theroyn/RealTimeRayTracing), [OptiX](https://developer.nvidia.com/rtx/ray-tracing/optix), [Vulkan](https://github.com/GPSnoopy/RayTracingInVulkan)
  
  **Estimated date of completion:** <ins>22nd August 2023</ins> which is the deadline of this module (CSC8599).

## Housekeeping:

- When developing this project, the code is running on <ins>11980HK</ins> and <ins>RTX3080 (laptop version)</ins>.
- The code is written in Visual Studio 2022 with C++20.
- This implementation uses RIGHT-handed coordinate systems.

## Progression

- [An offine ray tracer prototype](https://github.com/IQ404/8599-ray-tracer-prototype/blob/main/README.md) (Finished)
- [Integrating the ray tracer with the Walnut GUI](https://github.com/IQ404/8599-ray-tracer-with-walnut-gui/blob/master/README.md) (Ongoing...)
