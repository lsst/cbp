# Code for the collimated beam projector (CBP)

Includes:
`lsst.cbp.CoordinateConverter`: Compute the telescope and CBP pointing that will give you
    a desired beam arrangement, such as placing beam B at point P on the pupil
    and point D on a specified detector.

`lsst.cbp.computeHolePositions`: compute hole positions for a CBP mask.

Notes
-----

scons will not yet build the docs; to build them manually:

    cd doc
    sphinx-build . _build/html
