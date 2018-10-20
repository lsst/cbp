# Code for the [collimated beam projector](https://arxiv.org/abs/1805.05867) (CBP)

Includes:

* `lsst.cbp.CoordinateConverter`: Compute the telescope and CBP pointing that will give you a desired
  beam arrangement, such as placing beam B at point P on the pupil and point D on a specified detector.
* `lsst.cbp.computeHolePositions`: compute hole positions for a CBP mask.

The package is designed for use with LSST DM's `scons` build system and `eups` package management system.
Assuming you have the basic LSST DM stack installed you can do the following, from within the package directory:

- `setup -r .` to setup the package and dependencies.
- `scons` to build the package and run unit tests.
- `scons install declare` to install the package and declare it to eups.
- `package-docs build` to build the documentation.
  This requires optional [dependencies](https://developer.lsst.io/stack/building-single-package-docs.html)
  beyond those required to build and use the package.
