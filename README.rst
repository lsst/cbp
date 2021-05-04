###
cbp
###

Code for the `collimated beam projector <https://arxiv.org/abs/1805.05867>`_ (CBP)

`Documentation <https://pipelines.lsst.io/modules/lsst.cbp/index.html>`_

The package is designed for use with the Vera Rubin LSST DM's `scons` build system and `eups` package management system.
Assuming you have the basic DM stack installed you can do the following, from within the package directory:

* ``setup -r .`` to setup the package and dependencies, at which point the unit tests can be run and the package can be used "in place".
* ``pytest`` or ``scons`` to run the unit tests.
* ``scons install declare`` to install the package and declare it to eups.
* ``package-docs build`` to `build the documentation <https://developer.lsst.io/stack/building-single-package-docs.html>`_
