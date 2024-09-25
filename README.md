# Zap - Differentiating Shocks in the two-dimensional Burgers equation

This repository contains a Finite Volume solver for the two-dimensional Burgers equation
$$\partial_t u + \partial_x \frac{1}{2} u^2 + \partial_y \frac{1}{2} u^2 = 0$$
$$\Leftrightarrow{} \partial_t u + u \partial_x u + u \partial_y u = 0$$
with given initial- and boundary conditions.

## Quickstart
TODO: This is very outdated!

Building and running the solver (an appropriate value for nx and ny is 400)
```console
$ cmake -Bbuild
$ cd bin
$ make -j
$ ./run_solver <nx> <ny>
```

The solution is written to the output directory in the files `u.dat` and `t.dat`.
It can be rendered to a video using the provided (recommended) renderer
```console
$ ./renderer ../output/u.dat ../output/t.dat solution.mp4
```
or the (not recommended) legacy renderer
```console
$ ./legacy_renderer ../output/u.dat solution.mp4
```

The legacy renderer (`legacy_renderer.cpp`) does not depend on `FreeType`, but does not have all the functions of the recommended renderer.
If you do not want to install the dependencies for the renderer, you can disable building it by setting the cmake-option `ZAP_BUILD_RENDERER` to `OFF`.
Not installing FFmpeg will only create problems at runtime.
Both renderers depend on Unix functions and will only work on Unix compliant operating systems (`pipe`, `write`, ...), e.g. MacOS or Linux (only tested on MacOS).

## Dependencies

- [Eigen](https://eigen.tuxfamily.org/): For all matrix operations.
- [FFmpeg](https://ffmpeg.org/): For rendering the solution.
- [FreeType](https://freetype.org/): For rendering text.
- [stb_image_write](https://github.com/nothings/stb/): For saving the canvas in .jpeg format.

FFmpeg and FreeType must be installed by the user, Eigen and stb\_image\_write are provided in the `ThirdParty` folder.
