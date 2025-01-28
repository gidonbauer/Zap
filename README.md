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

## Custom data formats

The solver save intermediate steps in a custom binary data format.
All time steps and solution at that time are stored in two files.
All iterations are appended to the corresponding file.
The files contain a simple header followed by the data.

### Structured Matrix Data

Matrix-like data is stored in the "incremental matrix format".
This format is used for saving the time steps, which are essentially a #box[$1 times 1$-matrix] and for saving the solution of the Godunov-type solver without shock tracking that is used as a comparison.
The header has the following layout:
1. Six bytes magic string "`$SHOCK`" to identify the format.
2. Three bytes for the data type stored in the file, the representation for each data type can be seen in the table below.
3. Number of rows encoded as a 64-bit unsigned integer.
4. Number of columns encoded as a 64-bit unsigned integer.
5. Single byte indicating if the data is saved in a row-major order; 1 if data is row-major, 0 if data is column-major.
The matrix data for each time step is then simply stored as a binary block of data.
This is done repeatedly and the information from the header is used to reconstruct the data when reading it.
The shape and size of the matrix may not change for any time step, this is obviously the case for both use cases of this format.
This format was inspired by the `NumPy` `npy` data format.

|  | `unsigned integer` | `signed integer` | `floating-point` |
|:-|:------------------:|:----------------:|:----------------:|
|  8-bit | `" u8"` | `" i8"` | - |
| 16-bit | `"u16"` | `"i16"` | - |
| 32-bit | `"u32"` | `"i32"` | `"f32"` |
| 64-bit | `"u64"` | `"i64"` | `"f64"` |

### Unstructured Grid Data

The solution of the shock tracking solver cannot be saved in the "incremental matrix format" because it cannot encode the cut cells.
Therefore, a new "incremental grid format" is introduced that stores the cells directly.
The header has the following layout:
1. Six bytes magic string "`%SHOCK`" to identify the format.
2. Three bytes for the data type stored in the file, the representation for each data type can be seen in the table above.
3. Starting point of the grid in $x$-direction $x_\mathrm{min}$, binary value of data type defined in 2.
4. End point of the grid in $x$-direction $x_\mathrm{max}$, binary value of data type defined in 2.
5. Number of cells in $x$-direction $n_x$, encoded as a 64-bit unsigned integer.
6. Starting point of the grid in $y$-direction $y_\mathrm{min}$, binary value of data type defined in 2.
7. End point of the grid in $y$-direction $y_\mathrm{max}$, binary value of data type defined in 2.
8. Number of cells in $y$-direction $n_y$, encoded as a 64-bit unsigned integer.
9. Number of dimensions of the conserved quantity $u$ (here always 1 because $u \in \mathbb{R}^1$), encoded as a 64-bit unsigned integer.

Then the cells are saved one by one. Cartesian cells are saved in the following way:
1. Single byte with value 0 indicating a cartesian cell
1. Starting point of the cell in $x$-direction $x_{i-\tfrac{1}{2}}$
1. Cell size in $x$-direction $\Delta x$
1. Starting point of the cell in $y$-direction $y_{j-\tfrac{1}{2}}$
1. Cell size in $x$-direction $\Delta y$
1. Cell value $U_{i,j}$

The cut cells are saved as:
1. Single byte with value 1 indicating a cut cell
1. Starting point of the cell in $x$-direction $x_{i-\tfrac{1}{2}}$
1. Cell size in $x$-direction $\Delta x$
1. Starting point of the cell in $y$-direction $y_{j-\tfrac{1}{2}}$
1. Cell size in $x$-direction $\Delta y$
1. Single byte indicating the cut type as in the table below
1. Two numbers for the first cut point
1. Two numbers for the second cut point
1. Left subcell value $U^L_{i,j}$ where the left subcell is the cell containing the bottom left corner.
1. Right subcell value $U^R_{i,j}$

This is repeated for every time step.
The number of whole cell may not change but the number of cut cells and therefore the number of subcell, can change arbitrarily.

| | | | | |
|:-:|:-:|:-:|:-:|:-:|
|          | **`BOTTOM`** | **`RIGHT`** | **`TOP`** | **`LEFT`** |
| **`BOTTOM`** | --- | `BOTTOM_RIGHT` | `MIDDLE_VERT` | `BOTTOM_LEFT` |
| **`RIGHT`**  | `BOTTOM_RIGHT` | --- | `TOP_RIGHT` | `MIDDLE_HORI` |
| **`TOP`**    | `MIDDLE_VERT` | `TOP_RIGHT` | --- | `TOP_LEFT` |
| **`LEFT`**   | `BOTTOM_LEFT` | `MIDDLE_HORI` | `TOP_LEFT` | --- |

## Dependencies

- [AD.hpp](https://www.stce.rwth-aachen.de/): For algorithmic differentiation.
- [Eigen](https://eigen.tuxfamily.org/): For all matrix operations.
- [FFmpeg](https://ffmpeg.org/): For rendering the solution.
- [FreeType](https://freetype.org/): For rendering text.
- [stb_image_write](https://github.com/nothings/stb/): For saving the canvas in jpeg-format.
- [{FMT}](https://fmt.dev/latest/index.html): Text formatting
- [GoogleTest](https://github.com/google/googletest): Unit testing
- [benchmark](https://github.com/google/benchmark): Benchmarking

_FFmpeg_ and _FreeType_ must be installed by the user, _AD.hpp_, _Eigen_, and _stb\_image\_write_ are provided in the `ThirdParty` folder.