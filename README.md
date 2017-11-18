AOLineSegmentSampling
=====================

This source code is built on top of [PBRT-V3](http://pbrt.org). If you happen to use this source code in your research, please cite the code using the following entry:

> @article{singh17variance,<br>
>    author = "Singh, Gurprit and Miller, Bailey and Jarosz, Wojciech",<br>
>    title = "Variance and Convergence Analysis of Monte Carlo Line and Segment Sampling",<br>
>    journal = "Computer Graphics Forum (Proceedings of EGSR)",<br>
>    year = "2017",<br>
>    volume = "36",<br>
>    number = "4",<br>
>    month = "June",<br>
>    doi = "10.1111/cgf.13226",<br>
>    publisher = "The Eurographics Association",<br>
>    keywords = "stochastic sampling, signal processing, Fourier transform, Power spectrum"<br>
>}

Specifics
---------

Ambient Occlusion using points, segments and line samples is implemented as separate
integrators in pbrtv3 (description below).
```
av.cpp (Ambient Occlusion with points)
avls.cpp (Ambient Occlusion using segments)
avl.cpp (Ambient Occlusion using lines)
```

Each pixel is sampled at the center and a ray is shot from the center of the
pixel on the image plane (have a look at the `sampler.cpp` file).
```
Sampler::Sampler(int64_t samplesPerPixel) : samplesPerPixel(samplesPerPixel) {}
CameraSample Sampler::GetCameraSample(const Point2i &pRaster) {
    CameraSample cs;
    ///Ray is shot always from the center of the pixel on the Image plane.
    cs.pFilm = (Point2f)pRaster + Point2f(0.5f,0.5f);// Get2D();
    cs.time = Get1D();
    cs.pLens = Get2D();
    return cs;
}
```
We only vary the number of secondary rays directly by passing it as an
integrator parameter(s).

### New Samplers ###

We added regular (grid) and multi-jittered samplers.
You can also set `jitter [false]` in the stratified sampler to get
regular grid samples.

### Visualize Integrand ###

To visualize an integrand within a pixel, first define a cropwindow that is
one pixel wide, set `visualizeIntegrand[true]`,
choose a `regular` sampler and set the sample count to `262144` (in `points.pbrt
    set "integer numSecondarySamples" [262144]` whereas is `line_segments.pbrt
    set "integer numPhiSamples" [512] "integer numThetaSamples" [512]"`).

### Fourier Transform Integrand ###

To visualize/compute the Fourier power spectrum of an integrand within a pixel,
make sure to choose the `cropwindow` that is one pixel wide (e.g. `[410 411 305 306]`).
All the integrand power spectra added in the paper are directly computed
from PBRT using `regular` sampler.

Building source code
--------------------

To check out code together with all dependencies, be sure to use the
`--recursive` flag when cloning the repository, i.e.
```bash
$ git clone --recursive https://github.com/sinbag/ao-line-segment-sampling.git
```
If you accidentally already cloned pbrt without this flag (or to update an
pbrt source tree after a new submodule has been added, run the following
command to also fetch the dependencies:
```bash
$ git submodule update --init --recursive
```
As mentioned before, we simply add our ambient occlusion using points, segments and line samples as different integrators in PBRT-V3. Please read below for more details about PBRT-V3.

PBRT-V3
-------
pbrt uses [cmake](http://www.cmake.org/) for its build system.  On Linux
and OS X, cmake is available via most package management systems.  For
Windows, or to build it from source, see the [cmake downloads
page](http://www.cmake.org/download/).

* For command-line builds on Linux and OS X, once you have cmake installed,
  create a new directory for the build, change to that directory, and run
  `cmake [path to pbrt-v3]`. A Makefile will be created in that
  current directory.  Run `make -j4`, and pbrt, the obj2pbrt and imgtool
  utilities, and an executable that runs pbrt's unit tests will be built.
* To make an XCode project file on OS X, run `cmake -G Xcode [path to pbrt-v3]`.
* Finally, on Windows, the cmake GUI will create MSVC solution files that
  you can load in MSVC.

If you plan to edit the lexer and parser for pbrt's input files
(`src/core/pbrtlex.ll` and `src/core/pbrtparase.y`), you'll also want to
have [bison](https://www.gnu.org/software/bison/) and
[flex](http://flex.sourceforge.net/) installed. On OS X, note that the
version of flex that ships with the developer tools is extremely old and is
unable to process `pbrtlex.ll`; you'll need to install a more recent
version of flex in this case.

### Debug and Release Builds ###

By default, the build files that are created that will compile an optimized
release build of pbrt. These builds give the highest performance when
rendering, but many runtime checks are disabled in these builds and
optimized builds are generally difficult to trace in a debugger.

To build a debug version of pbrt, set the `CMAKE_BUILD_TYPE` flag to
`Debug` when you run cmake to create build files to make a debug build. For
example, when running cmake from the command lne, provide it with the
argument `-DCMAKE_BUILD_TYPE=Debug`. Then build pbrt using the resulting
build files. (You may want to keep two build directories, one for release
builds and one for debug builds, so that you don't need to switch back and
forth.)

Debug versions of the system run much more slowly than release
builds. Therefore, in order to avoid surprisingly slow renders when
debugging support isn't desired, debug versions of pbrt print a banner
message indicating that they were built for debugging at startup time.

### Build Configurations ###

There are two configuration settings that must be set at compile time. The
first controls whether pbrt uses 32-bit or 64-bit values for floating-point
computation, and the second controls whether tristimulus RGB values or
sampled spectral values are used for rendering.  (Both of these aren't
amenable to being chosen at runtime, but must be determined at compile time
for efficiency).

To change them from their defaults (respectively, 32-bit
and RGB.), edit the file `src/core/pbrt.h`.

To select 64-bit floating point values, remove the comment symbol before
the line:
```
//#define PBRT_FLOAT_AS_DOUBLE
```
and recompile the system.

To select full-spectral rendering, comment out the first of these two
typedefs and remove the comment from the second one:
```
typedef RGBSpectrum Spectrum;
// typedef SampledSpectrum Spectrum;
```
Again, don't forget to recompile after making this change.
