# WaveOpticsFiberRendering

## Overview
This repository holds the rendering code of the SIGGRAPH Asia 2020 paper
**A Wave Optics Based Fiber Scattering Model**, by Mengqi (Mandy) Xia, Bruce Walter, Eric Michielssen, David Bindel and Steve Marschner. More information about this project can be found at the [project website](https://mandyxmq.github.io/research/wavefiber.html). This rendering code is based on [pbrt-v3](https://github.com/mmp/pbrt-v3.git).

Building WaveOpticsFiberRendering
-------------

To check the code together with all dependencies, be sure to use the
`--recursive` flag when cloning the repository, i.e.
```bash
$ git clone --recursive https://github.com/mandyxmq/WaveOpticsFiberRendering.git/
```
If you accidentally already cloned pbrt without this flag (or to update an
pbrt source tree after a new submodule has been added, run the following
command to also fetch the dependencies:
```bash
$ git submodule update --init --recursive
```

Create a new directory for the build, change to that directory, and run
`cmake [path to WaveOpticsFiberRendering]`. A Makefile will be created in the current
directory.  Next, run `make` to build pbrt, the obj2pbrt and imgtool
utilities, and an executable that runs pbrt's unit tests.  Depending on the
number of cores in your system, you will probably want to supply make with
the `-j` parameter to specify the number of compilation jobs to run in
parallel (e.g. `make -j8`).

Download data and scenes
--------------
Examples of fiber Scattering functions can be downloaded from this [shared folder](https://drive.google.com/drive/folders/12MIcSFucB0IbVl5Y41Lgt2_JoICuJozF?usp=sharing). Please place this data folder under repo's root directory. The example scene files can be downloaded from this [shared folder](https://drive.google.com/drive/folders/1I0wX_URJ2uFsw2NF3abgbP2EcVLHHKF6?usp=sharing). The render images are located in the directory render/. Usage: In WaveOpticsFiberRendering/build, run ./pbrt ../scenes/1layer/1e-6_forward.pbrt.
