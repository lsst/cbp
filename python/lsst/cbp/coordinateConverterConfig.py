# This file is part of cbp.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""CoordinateConverterConfig class: configuration for
a CoordinateConverter
"""

__all__ = ["CoordinateConverterConfig"]

import numpy as np

from lsst.geom import degrees


class CoordinateConverterConfig:
    """Configuration for the CoordinateConverter.

    :ref:`Configuration <lsst.cbp.configuration>` for the
    lsst.cbp.CoordinateConverter.

    Parameters
    ----------
    telPupilOffset : `float`
        Offset of the telescope :ref:`pupil plane
        <lsst.cbp.pupil_position>` from the center of the telescope
        (mm, positive if closer to the sky).
    telPupilDiameter : `float`
        Diameter of telescope pupil (mm).
    telPupilObscurationDiameter : `float`
        Diameter of the telescope central obscuration at the
        `pupil plane <lsst.cbp.pupil_position>` (mm).
    telFocalPlaneDiameter : `float`
        Diameter of telescope flocal plane (mm).
    telFlipX : `bool`
        True if the x axis of the telescope focal plane
        is :ref:`flipped <lsst.cbp.flipped_x_axis>`
        with respect to the pupil frame.
    telAzimuthOffsetDeg : `float`
        Azimuth `offset`_ (degrees).
    telAzimuthScale : `float`
        Azimuth `scale`_; must be ±1
        (in order to handle wrap correctly).
    telAltitudeOffsetDeg : `float` (optional)
        Telescope altitude `offset`_ (degrees); defaults to 0.
    telAltitudeScale : `float` (optional)
        Telescope altitude `scale`_; defaults to 1.
    telAltitudeLimitsDeg : pair of `float`
        Telescope minimum, maximum allowed observed altitude (degrees).
    telRotOffsetDeg : `float`
        Telescope camera rotator offset (degrees).
    telRotScale : `float`
        Telescope camera rotator scale; must be ±1
        (in order to handle wrap correctly).
    defaultDetector : `str`
        Name of default detector.

    cbpPosition : triplet of `float`
        CBP x, y, z position of center of CBP relative to the center
        of the telescope, in the base frame (mm).
    cbpFocalLength : `float`
        Effective focal length of the CBP (mm);
        635 mm is an estimate for the LSST's CBP.
    cbpFlipX : `bool`
        True if the x axis of the CBP focal plane is flipped
        with respect to the pupil frame?
    cbpAzimuthOffsetDeg : `float`
        CBP azimuth `offset`_ (degrees).
    cbpAzimuthScale : `float`
        CBP azimuth `scale`_; must be ±1
        (in order to handle wrap correctly).
    cbpAltitudeOffsetDeg : `float` (optional)
        CBP altitude `offset`_ (degrees); defaults to 0.
    cbpAltitudeScale : `float` (optional)
        CBP altitude `scale`_; defaults to 1.
    cbpAltitudeLimitsDeg : pair of `float`
        CBP minimum, maximum allowed observed altitude (degrees).

    Raises
    ------
    ValueError
        Raised if ``telAzimuthScale``, ``cbpAzimuthScale``
        and/or ``telRotScale`` is not ±1.
    ValueError
        Raised if items with multiple values have the wrong length.

    Notes
    -----
    .. _offset:

    .. _scale:

    **Offset and Scale:**

    Azimuth, altitude and rotator offset and scale define the mapping
    between :ref:`internal angle <lsst.cbp.internal_angles>`
    and :ref:`observed angle <lsst.cbp.observed_angles>` as follows::

        observed angle = internal angle * scale + offset

    .. _fields:

    **Fields:**

    telPupilOffset : `float`
        Offset of the telescope :ref:`pupil plane
        <lsst.cbp.pupil_position>` from the center of the telescope
        (mm, + if closer to the sky).
    telPupilDiameter : `float`
        Diameter of telescope pupil (mm).
    telPupilObscurationDiameter : `float`
        Diameter of the telescope central obscuration at the
        `pupil plane <lsst.cbp.pupil_position>` (mm).
    telFocalPlaneDiameter : `float`
        Diameter of telescope flocal plane (mm).
    telFlipX : `bool`
        True if the x axis of the telescope focal plane
        is :ref:`flipped <lsst.cbp.flipped_x_axis>`
        with respect to the pupil frame.
    telAzAltOffset : pair of `lsst.geom.Angle`
        Telescope azimuth and altitude `offset`_ (degrees).
    telAzAltScale : pair of `float`
        Telescope azimuth and altitude `scale`_;
        azimuth scale is ±1.
    telAltitudeLimits : pair of `lsst.geom.Angle`
        Telescope minimum, maximum allowed observed altitude.
    telRotOffset : `lsst.geom.Angle`
        Telescope camera rotator offset.
    telRotScale : `float`
        Telescope camera rotator scale; must be ±1.
    defaultDetector : `str`
        Name of default detector.

    cbpFocalLength : `float`
        Effective focal length of the CBP (mm);
        635 mm is an estimate for the LSST's CBP.
    cbpFlipX : `bool`
        True if the x axis of the CBP focal plane is flipped
        with respect to the pupil frame?
    cbpAzAltOffset : pair of `lsst.geom.Angle`
        CBP azimuth and altitude `offset`_ (degrees).
    cbpAzAltScale : pair of `float`
        CBP azimuth and altitude `scale`_; azimuth scale is ±1.
    cbpAltitudeLimits : pair of `lsst.geom.Angle`
        CBP minimum, maximum allowed observed altitude.
    """

    def __init__(self, *, telPupilOffset,
                 telPupilDiameter, telPupilObscurationDiameter, telFocalPlaneDiameter, telFlipX,
                 telAzimuthOffsetDeg, telAzimuthScale,
                 telAltitudeOffsetDeg=0, telAltitudeScale=1, telAltitudeLimitsDeg,
                 telRotOffsetDeg, telRotScale,
                 defaultDetector,
                 cbpPosition, cbpFocalLength, cbpFlipX,
                 cbpAzimuthOffsetDeg, cbpAzimuthScale,
                 cbpAltitudeOffsetDeg=0, cbpAltitudeScale=1, cbpAltitudeLimitsDeg):
        self.telPupilOffset = telPupilOffset
        self.telPupilDiameter = telPupilDiameter
        self.telPupilObscurationDiameter = telPupilObscurationDiameter
        self.telFocalPlaneDiameter = telFocalPlaneDiameter
        self.telFlipX = telFlipX
        self.telAzAltOffset = (telAzimuthOffsetDeg*degrees, telAltitudeOffsetDeg*degrees)
        if abs(telAzimuthScale) != 1:
            raise ValueError("telAzimuthScale={} must be +/-1".format(telAzimuthScale))
        self.telAzAltScale = (telAzimuthScale, telAltitudeScale)
        if len(telAltitudeLimitsDeg) != 2:
            raise ValueError("telAltitudeLimitsDeg={!r} must have length 2".format(telAltitudeLimitsDeg))
        self.telAltitudeLimits = tuple(val*degrees for val in telAltitudeLimitsDeg)
        self.telRotOffset = telRotOffsetDeg*degrees
        if abs(telRotScale) != 1:
            raise ValueError("telRotScale={} must be +/-1".format(telRotScale))
        self.telRotScale = telRotScale

        self.defaultDetector = defaultDetector

        self.cbpPosition = cbpPosition
        self.cbpFocalLength = cbpFocalLength
        self.cbpFlipX = cbpFlipX
        self.cbpAzAltOffset = (cbpAzimuthOffsetDeg*degrees, cbpAltitudeOffsetDeg*degrees)
        if abs(cbpAzimuthScale) != 1:
            raise ValueError("cbpAzimuthScale={} must be +/-1".format(cbpAzimuthScale))
        self.cbpAzAltScale = (cbpAzimuthScale, cbpAltitudeScale)
        if len(cbpAltitudeLimitsDeg) != 2:
            raise ValueError("cbpAltitudeLimitsDeg={!r} must have length 2".format(cbpAltitudeLimitsDeg))
        self.cbpAltitudeLimits = tuple(val*degrees for val in cbpAltitudeLimitsDeg)

    @property
    def cbpPosition(self):
        """The position of the CBP in the base frame (mm, read/write)."""
        return self._cbpPosition

    @property
    def cbpDistance(self):
        """The distance from the telescope to the CBP (mm, read only)."""
        return self._cbpDistance

    @cbpPosition.setter
    def cbpPosition(self, cbpPosition):
        """Set the position of the CBP.

        Parameters
        ----------
        cbpPosition : triplet of `float`
            CBP x, y, z position of center of CBP relative to the center
            of the telescope, in the base frame (mm).
        """
        if len(cbpPosition) != 3:
            raise ValueError("cbpPosition={!r} must have length 2".format(cbpPosition))
        self._cbpPosition = np.array(cbpPosition)
        self._cbpDistance = np.linalg.norm(self.cbpPosition)
