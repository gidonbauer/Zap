# Zap - Differentiating Shocks in the two-dimensional Burgers equation

This repository contains a Shock-Tracking Finite Volume solver for the two-dimensional Burgers equation

$$\partial_t u + \partial_x \frac{1}{2} u^2 + \partial_y \frac{1}{2} u^2 = 0$$

$$\Leftrightarrow{} \partial_t u + u \partial_x u + u \partial_y u = 0$$

with given initial- and boundary conditions.

## Quickstart

Building and running the examples:
```console
$ cmake -Bbuild
$ cd bin
$ make -j
$ ./example_quarter_circle
$ ./example_x_ramp
```

The examples can be run with different parameters, to see all options run the example with the `--help` flag.

The solution is written to the directory `output/<executable name>/` in the files `u*.grid` and `t*.mat`.
It can be rendered to a video using the provided renderer, for example the default quarter circle example can be rendered like with the command
```console
$ ./renderer --scale 50 ../output/example_quarter_circle/u_cell_based_25x25.grid ../output/example_quarter_circle/t_cell_based_25x25.mat solution.mp4
```

The renderer depends on the library _FreeType_ and _FFmpeg_.
Those dependencies must be installed to compile the renderer.
If you do not want to install the dependencies for the renderer, you can disable building it by setting the cmake-option `ZAP_BUILD_RENDERER` to `OFF`.
Not installing _FFmpeg_ will only create problems at runtime.
The renderer depends on Unix functions (`pipe`, `write`, ...) and will only work on Unix compliant operating systems, e.g. MacOS or Linux (only tested on MacOS).

## Dependencies

- [AD.hpp](https://www.stce.rwth-aachen.de/): For algorithmic differentiation.
- [Eigen](https://eigen.tuxfamily.org/): For all matrix operations.
- [FFmpeg](https://ffmpeg.org/): For rendering the solution.
- [FreeType](https://freetype.org/): For rendering text.
- [stb_image_write](https://github.com/nothings/stb/): For saving the canvas in jpeg-format.

_FFmpeg_ and _FreeType_ must be installed by the user, _AD.hpp_, _Eigen_, and _stb\_image\_write_ are provided in the `ThirdParty` folder.
