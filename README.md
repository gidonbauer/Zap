# Zap - Differentiating Shocks in the two-dimensional Burgers equation

This repository contains a Finite Volume solver for the two-dimensional Burgers equation
$$\partial_t u + \partial_x \frac{1}{2} u^2 + \partial_y \frac{1}{2} u^2 = 0$$
$$\Leftrightarrow{} \partial_t u + u \partial_x u + u \partial_y u = 0$$
with given initial- and boundary conditions.

## Quickstart

Building and running the solver (an appropriate value for nx and ny is 400)
```console
$ cmake -Bbuild
$ cd bin
$ make -j
$ ./run_solver <nx> <ny>
```

The solution is written to the output directory in the files `u.dat` and `t.dat`.
It can be rendered to a video using the provided renderer

```console
$ ./renderer ../output/u.dat ../output/t.dat solution.mp4
```

## Dependencies

- [Eigen](https://eigen.tuxfamily.org/): For all matrix operations.
- [ffmpeg](https://ffmpeg.org/): For rendering the solution.
- [FreeType](https://freetype.org/): For rendering text.

FFmpeg and FreeType must be installed by the user, Eigen is provided in the `ThirdPary` folder.