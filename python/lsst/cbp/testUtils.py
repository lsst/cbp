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
"""SampleCoordinateConverter class: make a CoordinateConverter for tests"""

__all__ = ["SampleCoordinateConverter"]

import numpy as np

import lsst.geom
from lsst.afw.geom import makeRadialTransform
from lsst.afw.cameraGeom import (Camera, Amplifier, FIELD_ANGLE, ReadoutCorner,
                                 addDetectorBuilderFromConfig, DetectorConfig)
from .coordinateConverterConfig import CoordinateConverterConfig
from .coordinateConverter import CoordinateConverter
from .computeHolePositions import computeHolePositions
from .maskInfo import MaskInfo


class SampleCoordinateConverter:
    """An object containing a CoordinateConverter and the information
    used to create it.

    Parameters
    ----------
    detectorFracPosList : `iterable` of pair of `float` (optional)
        Position of the center of each detector, as a fraction of
        the width and height of the detector.
        The first element must have value (0, 0).
        See the field of the same name for more information.
        Defaults to::

            (
                (0, 0),
                (1.01, 0),  # 1.01: leave a 1% gap
                (-4, 7),    # a corner detector in the LSST camera
            )

    holeFracPosList : `iterable` of pair of `float` (optional)
        Positions of holes on a given detector,
        as a fraction of the distance from lower left corner
        to upper right corner. Thus (0.5, 0.5) is centered
        on the detector.
        Defaults to `((0, 0), (0.75, 0.75))`.

    Notes
    -----
    **Attributes**

    detectorWidthPix : `int`
        Width of each detector, in pixels.
    detectorHeightPix : `int`
        Height of each detector, in pixels.
    pixelSizeMm : `float`
        Width = height of each pixel, in mm.
    plateScale : `lsst.geom.Angle`
        Plate scale: in angle on the sky per mm on the focal plane.
    detectorFracPosList : `iterable` of pair of `float`
        Position of the center of each detector, as a fraction of
        the width and height of the detector. For instance
        (0, 0) is a detector centered on the focal plane
        and (1, 0) is adjacent to a centered detector,
        in the direction of increasing focal plane x.
    holeFracPosList : `iterable` of pair of `float`
        Positions of holes on a given detector,
        as a fraction of the distance from lower left corner
        to upper right corner. Thus (0.5, 0.5) is centered
        on the detector.
    cameraGeom : `lsst.afw.cameraGeom.Camera`
        Camera geometry. There will be one detector per entry in
        detectorFracPosList with names "D0", "D1", ...
        Detector "D0" is centered on the focal plane.
    config : `lsst.cbp.CoordinateConverterConfig`
        Basic configuration for ``coordinateConverter``.
    maskInfo : `lsst.cbp.MaskInfo`
        CBP mask information.
    coordinateConverter : `lsst.cbp.CoordinateConverter`
        The test coordinate converter.
    """

    def __init__(self, detectorFracPosList=None, holeFracPosList=None,
                 telFlipX=False, cbpFlipX=False):
        # these value are close to LSST and are rectangular
        # in order to catch axis transposition errors
        self.detectorWidthPix = 4000
        self.detectorHeightPix = 4095
        self.pixelSizeMm = 0.01
        self.plateScale = 20 * lsst.geom.arcseconds
        if holeFracPosList is None:
            holeFracPosList = ((0.5, 0.5), (0.75, 0.75))
        self.holeFracPosList = holeFracPosList
        if detectorFracPosList is None:
            detectorFracPosList = (
                (0, 0),
                (1.01, 0),  # 1.01: leave a 1% gap
                (-4.04, 7.07),    # approximately a corner detector in the LSST camera
            )
        self.detectorFracPosList = detectorFracPosList

        self.cameraGeom = self.makeCameraGeom()
        self.config = self.makeCoordinateConverterConfig(
            telFlipX=telFlipX,
            cbpFlipX=cbpFlipX,
        )
        self.maskInfo = self.makeMaskInfo()
        self.coordinateConverter = CoordinateConverter(
            config=self.config,
            maskInfo=self.maskInfo,
            cameraGeom=self.cameraGeom,
        )

    def makeCoordinateConverterConfig(self, telFlipX, cbpFlipX):
        """Make a coordinate converter config.

        Parameters
        ----------
        telFlipX : `bool`
            :ref:`Flip <lsst.cbp.flipped_x_axis>` the
            telescope focal plane.
        cbpFlipX : `bool`
            :ref:`Flip <lsst.cbp.flipped_x_axis>` the CBP focal plane.

        Returns
        -------
        config : `lsst.cbp.CoordinateConverterConfig`
            Coordinate converter config.
        """
        return CoordinateConverterConfig(
            telPupilOffset=101,
            telPupilDiameter=8500,
            telPupilObscurationDiameter=1000,
            telFocalPlaneDiameter=3000,
            telFlipX=telFlipX,
            telAzimuthOffsetDeg=-180,  # offset=-180, scale=-1 for az=0 N, 90 E
            telAzimuthScale=-1,
            telAltitudeOffsetDeg=0,
            telAltitudeScale=1,
            telAltitudeLimitsDeg=(0, 89),
            telRotOffsetDeg=0,
            telRotScale=-1,
            defaultDetector="D0",
            cbpPosition=(10000, 3000, 5000),
            cbpFocalLength=635,  # nominal value for LSST CBP
            cbpFlipX=cbpFlipX,
            cbpAzimuthOffsetDeg=-180,
            cbpAzimuthScale=-1,
            cbpAltitudeOffsetDeg=0,
            cbpAltitudeScale=1,
            cbpAltitudeLimitsDeg=(-70, 70),
        )

    def makeMaskInfo(self):
        """Make mask information.

        Returns
        -------
        maskInfo : `lsst.cbp.MaskInfo`
            Mask info.

        Notes
        -----
        The mask will have one hole per entry in self.holeFracPosList
        per detector.

        self.cameraGeom and self.config must be set before calling this
        method.
        """
        detectorBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                                       lsst.geom.Extent2I(self.detectorWidthPix, self.detectorHeightPix))
        bboxd = lsst.geom.Box2D(detectorBBox)
        bboxdMin = np.array(bboxd.getMin())
        bboxdDim = np.array(bboxd.getDimensions())
        detectorPositions = [bboxdMin + bboxdDim*np.array(fracPos) for fracPos in self.holeFracPosList]
        holePositions = computeHolePositions(
            detectorNames=None,
            detectorPositions=detectorPositions,
            cameraGeom=self.cameraGeom,
            cbpFlipX=self.config.cbpFlipX,
            cbpFocalLength=self.config.cbpFocalLength,
        )
        holeNames = ["beam{}".format(i) for i in range(len(holePositions))]
        return MaskInfo(
            name="test",
            holePositions=holePositions,
            holeNames=holeNames,
            defaultHole=0,
        )

    def makeCameraGeom(self):
        """Make a camera geometry.

        Returns
        -------
        cameraGeom : `lsst.afw.cameraGeom.Camera`
            Camera geometry.

        Notes
        -----
        There is one field per entry in self.detectorFracPosList
        with specifications set by self.detectorWidthPix,
        self.detectorHeightPix, and self.pixelSizeMm.

        The plate scale is set by self.plateScale
        and the amount of optical distortion is fixed.

        All detectors have the same shape (unlike LSST) and orientation
        (unlike HSC). Varying these is not necessary for testing the CBP
        and having all detectors the same simplifies the code.
        """
        radialCoeff = np.array([0.0, 1.0, 0.0, 0.925]) / self.plateScale.asRadians()
        fieldAngleToFocalPlane = makeRadialTransform(radialCoeff)
        focalPlaneToFieldAngle = fieldAngleToFocalPlane.inverted()

        cameraBuilder = Camera.Builder("testCamera")
        cameraBuilder.setTransformFromFocalPlaneTo(FIELD_ANGLE, focalPlaneToFieldAngle)
        ampBuilder = self._makeAmpBuilder()

        for i, fpPos in enumerate(self.detectorFracPosList):
            detectorConfig = self._makeDetectorConfig(id=i, fpPos=fpPos)
            addDetectorBuilderFromConfig(cameraBuilder, detectorConfig, [ampBuilder], focalPlaneToFieldAngle)

        return cameraBuilder.finish()

    def _makeAmpBuilder(self):
        """Construct a trivial amplifier builder.

        The CBP code does not care about the details of the amplifier, so this
        builder is as simple as possible: one amplifier that covers the whole
        CCD, with no overscan, and semi-plausible valus for everything else.

        Returns
        -------
        ampBuilder : `lsst.afw.cameraGeom.Amplifier.Builder`
            Amplifier builder.
        """
        ampExtent = lsst.geom.Extent2I(self.detectorWidthPix, self.detectorHeightPix)
        ampBBox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0), ampExtent)

        ampBuilder = Amplifier.Builder()
        ampBuilder.setName("TestAmp")
        ampBuilder.setBBox(ampBBox)
        ampBuilder.setGain(1.8)
        ampBuilder.setReadNoise(3.9)
        ampBuilder.setReadoutCorner(ReadoutCorner.LL)

        return ampBuilder

    def _makeDetectorConfig(self, id, fpPos):
        """Make a detector config for one detector.

        Parameters
        ----------
        id : `int`
            Detector ID.
        fpPos : trio of `float`
            Focal plane position of detector, in units of detector
            width/height. For example:

            - (0, 0) is a detector centered on the focal plane
            - (1, 0) is adjacent to a centered detector,
              in the direction of increasing focal plane x

        Returns
        -------
        detectorConfig : `lsst.afw.cameraGeom.DetectorConfig`
            Detector configuration for one detector.
        """
        bbox = lsst.geom.Box2I(lsst.geom.Point2I(0, 0),
                               lsst.geom.Extent2I(self.detectorWidthPix, self.detectorHeightPix))
        ctr = lsst.geom.Box2D(bbox).getCenter()
        pixelSizeMm = 0.01
        config = DetectorConfig()
        config.name = "D{}".format(id)
        config.id = id
        # detector serial number is not used by the CBP code,
        # but some value is required to construct a Detector
        config.serial = '0000011'
        config.detectorType = 0
        config.bbox_x0 = bbox.getMinX()
        config.bbox_x1 = bbox.getMaxX()
        config.bbox_y0 = bbox.getMinY()
        config.bbox_y1 = bbox.getMaxY()
        config.pixelSize_x = pixelSizeMm
        config.pixelSize_y = pixelSizeMm
        config.transformDict.nativeSys = 'Pixels'
        config.transformDict.transforms = None
        config.refpos_x = ctr[0]
        config.refpos_y = ctr[1]
        config.offset_x = fpPos[0] * pixelSizeMm * self.detectorWidthPix
        config.offset_y = fpPos[1] * pixelSizeMm * self.detectorHeightPix
        config.offset_z = 0.0  # No displacement along optic axis
        config.transposeDetector = False
        config.pitchDeg = 0.0
        config.yawDeg = 0.0
        config.rollDeg = 0.0
        return config
