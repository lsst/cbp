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
"""CoordinateConverter class: coordinate conversions for the CBP."""

__all__ = ["CoordinateConverter"]

import math

import numpy as np
import scipy.optimize

from . import coordUtils
from lsst.geom import Extent2D, Point2D, SpherePoint, radians
from lsst.afw.cameraGeom import PIXELS, FOCAL_PLANE, FIELD_ANGLE
from .beamInfo import BeamInfo

# Record errors from setFocalFieldAngle?
_RecordErrors = False
# list of errors from the root finder: error (abs value, radians), pupilPos, focalFieldAngle, beam
_ErrorList = None


def startRecordingErrors():
    """Start recording setFocalFieldAngle errors and reset accumulators."""
    global _RecordErrors, _ErrorList
    _RecordErrors = True
    _ErrorList = []


def stopRecordingErrors():
    """Stop recording errors."""
    global _RecordErrors
    _RecordErrors = False


def getRecordedErrors():
    """Return setFocalFieldAngle errors.

    Returns
    -------
    rootFinderErrors : `list`
        Each element is a tuple of the following:
        - |error| in radians
        - pupilPos
        - focalFieldAngle
        - beam name
    """
    return _ErrorList


class CoordinateConverter:
    """Coordinate conversions for the collimated beam projector (CBP).

    This object supports the following tasks:

    - Given a desired CBP "beam arrangement" (e.g. a particular beam
      should fall on a particular spot on the pupil
      and be focused to spot on a particular position of a detector),
      compute the necessary telescope and CBP pointing
      and the information about where each beam is directed.
    - Given the current telescope and CBP pointing and a desired
      offset to the resulting CBP beam arrangement, compute the new
      telescope and CBP pointing.

    See :ref:`how to use this object
    <lsst.cbp.coordinateConverter.howToUse>`
    for a summary of how to use this object.

    Parameters
    ----------
    config : `lsst.cbp.CoordinateConverterConfig`
        Telescope and CBP :ref:`configuration <lsst.cbp.configuration>`.
    maskInfo : `lsst.cbp.MaskInfo`
        Information about the CBP mask.
    camera : `lsst.afw.cameraGeom.Camera`
        Camera geometry.

    Notes
    -----
    .. _lsst.cbp.coordinateConverter.howToUse:

    **How to Use This Object**

    Call a :ref:`set method <lsst.cbp.coordinateConverter.setMethods>`,
    such as `setDetectorPos` to specify a desired CBP beam arrangement, or
    an :ref:`offset method <lsst.cbp.coordinateConverter.offsetMethods>`,
    such as `offsetDetectorPos`, to offset the current beam arrangement.

    This computes new values for telescope and CBP pointing, as attributes
    `telAzAltObserved`, `telRotObserved` and `cbpAzAltObserved`.
    Read these attributes and move the telescope and CBP accordingly.

    Get information about the beams. There are two ways to do this:

    - To get information for all beams, iterate on this object
      to get one `lsst.cbp.BeamInfo` per beam. Also the length
      of this object is the number of beams.
    - To get information for a specific beam, use ``[beamNameOrIndex]``;
      see ``__getitem__`` for details.
      Also attribute `beamNames` provides an iterable of beam names.

    That is basically it. However, it may also help to know the following:

    Whenever the telescope, camera rotator or CBP are moved,
    you must update the appropriate attribute(s) accordingly.
    Otherwise the beam information will be incorrect
    and offset commands will not work as expected.
    After each such update you can read the new beam information
    as usual.

    You are free to update configuration information at any time,
    by setting the appropriate attribute(s).
    The two items you are most likely to update are:

    - ``maskInfo``: set this when you change the mask
    - ``config.telRotOffset``: set this if you want to correct
      the orientation of the spot pattern on the focal plane.

    After updating configuration information, read the new beam information
    (and possibly new telescope and/or CBP position information,
    though few configuration parameters directly affect those)
    to see the effect of the change.

    .. _lsst.cbp.coordinateConverter.setMethods:

    **Set Methods**

    The methods `setFocalPlanePos`, `setDetectorPos` and
    `setFocalFieldAngle` all allow you to specify the desired
    arrangement for one beam:

    - The position of the beam on the pupil.
    - The position of the spot on the focal plane, expressed in different
        ways depending on the method. In most cases you will probably
        want to specify the position of the spot in pixels on a sepecified
        detector, in which case call `setFocalPlanePos`.

    These set methods simply update the pointing of the telescope and CBP
    (`telAzAltObserved`, `telRotObserved` and `cbpAzAltObserved`).
    If you then move the telescope and CBP as suggested, the beam
    should have the arrangement you specified, and the spot pattern of all
    the beams should be aligned with the detectors.

    .. _lsst.cbp.coordinateConverter.offsetMethods:

    **Offset Methods**

    The methods `offsetFocalPlanePos`, `offsetDetectorPos` and
    `offsetFocalFieldAngle` all allow you to offset the arrangement
    for one beam:

    - Offset the position of the beam on the pupil.
    - Offset the position of the spot on the focal plane, expressed in different
        ways depending on the method. In most cases you will probably
        want to specify the offset of the spot in pixels on a sepecified
        detector, in which case call `offsetFocalPlanePos`.

    These offset methods simply update the pointing of the telescope and
    CBP (`telAzAltObserved`, `telRotObserved` and `cbpAzAltObserved`).
    If you then move the telescope and CBP as suggested, the beam
    should have the arrangement you specified, and the spot pattern of all
    the beams should be aligned with the detectors.

    """

    def __init__(self, config, maskInfo, cameraGeom):
        self.config = config
        self.maskInfo = maskInfo
        self.cameraGeom = cameraGeom
        self._fieldAngleToFocalPlane = cameraGeom.getTransform(FIELD_ANGLE, FOCAL_PLANE)
        # amount to add to default hole position to compute telescope rotator angle (pixels);
        # I found that a wide range of values works, with (1, 0) comfortably in that range
        self._holeDelta = Extent2D(1, 0)
        self._telAzAlt = SpherePoint(np.nan, np.nan, radians)
        self._telRot = np.nan*radians
        self._cbpAzAlt = SpherePoint(np.nan, np.nan, radians)

    def setFocalFieldAngle(self, pupilPos, focalFieldAngle=None, beam=None):
        """Set the focal plane field angle of a beam.

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm).
        focalFieldAngle : pair of `float` (optional)
            :ref:`Focal plane field angle <lsst.cbp.focal_plane_field_angle>`
            of the specified beam (x, y rad); defaults to (0, 0).
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to
            ``self.maskInfo.defaultBeam``.
        """
        beam = self.maskInfo.asHoleName(beam)
        if focalFieldAngle is None:
            focalFieldAngle = Point2D(0, 0)
        else:
            focalFieldAngle = Point2D(*focalFieldAngle)

        # if the field angle is too small to matter, and too small
        # to determine its orientation, treat it as 0,0 (no iteration required)
        if math.hypot(*focalFieldAngle) < 1e-10:
            self.setPupilFieldAngle(pupilPos, focalFieldAngle, beam)
            return

        # minimize the field angle error as a function of the orientation
        # of the field angle; start with rotation = 0 so pupil field angle
        # equals focal plane field angle

        # record initial conditions, in case the miminizer fails to converge
        telAzAlt = self._telAzAlt
        telRot = self._telRot
        cbpAzAlt = self._cbpAzAlt
        try:
            class TrialFunctor:
                """Functor to compute error in focal plane orientation.

                Parameters
                ----------
                cco : `lsst.cbp.CoordinateConverter`
                    A coordinate converter object.
                pupilPos : pair of `float`
                    Position of the specified beam on the
                    :ref:`telescope pupil <lsst.cbp.pupil_position>`
                    (x, y mm).
                focalFieldAngle : pair of `float` (optional)
                    :ref:`Focal plane field angle
                    <lsst.cbp.focal_plane_field_angle>`
                    of the specified beam (x, y rad); defaults to (0, 0).
                beam : `int` or `str` (optional)
                    Name or index of beam; defaults to
                    ``self.maskInfo.defaultBeam``.
                """

                def __init__(self, cco, pupilPos, focalFieldAngle, beam):
                    self.cco = cco
                    self.pupilPos = pupilPos
                    # desired focal plane field angle
                    self.focalFieldAngle = focalFieldAngle
                    self.beam = beam
                    self.focalFieldAngleOrientation = math.atan2(focalFieldAngle[1],
                                                                 focalFieldAngle[0])*radians

                def __call__(self, rotAngleRadArr):
                    """Compute the error in focal plane orientation
                    (in radians) at the specified camera rotation angle.

                    Parameters
                    ----------
                    rotAngleRadArr : sequence of one `float`
                        The internal camera rotator angle, in radians.
                        It is passed in as a sequence of one element,
                        as required by scipy.optimize.
                    """
                    rotAngleRad = rotAngleRadArr[0]
                    pupilFieldAngle = coordUtils.rotate2d(pos=self.focalFieldAngle, angle=rotAngleRad*radians)
                    self.cco.setPupilFieldAngle(pupilPos=self.pupilPos,
                                                pupilFieldAngle=pupilFieldAngle,
                                                beam=self.beam)
                    measuredFocalFieldAngle = self.cco[beam].focalFieldAngle
                    measuredFocalFieldAngleOrientation = math.atan2(measuredFocalFieldAngle[1],
                                                                    measuredFocalFieldAngle[0])*radians
                    return self.focalFieldAngleOrientation.separation(
                        measuredFocalFieldAngleOrientation).asRadians()

            funcToFindRoot = TrialFunctor(cco=self, pupilPos=pupilPos,
                                          focalFieldAngle=focalFieldAngle, beam=beam)

            # Fit the rotator angle using a root finder
            with np.errstate(divide="ignore", invalid="ignore"):
                iterResult = scipy.optimize.root(fun=funcToFindRoot, x0=np.array([0.0], dtype=float),
                                                 options=dict(fatol=1e-11), method="broyden1")

            if not iterResult.success:
                raise RuntimeError("Iteration failed to converge")
            # call the function again to make sure the final value found is the one that is used
            err = funcToFindRoot(iterResult.x)

            global _RecordErrors, _ErrorList
            if _RecordErrors:
                _ErrorList.append((abs(err), pupilPos, focalFieldAngle, beam))

        except Exception:
            self._telAzAlt = telAzAlt
            self._telRot = telRot
            self._cbpAzAlt = cbpAzAlt
            raise

    def setFocalPlanePos(self, pupilPos, focalPlanePos=None, beam=None):
        """Set the position of a spot on the focal plane.

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the :ref:`telescope pupil
            <lsst.cbp.pupil_position>` (x, y mm)
        focalPlanePos : pair of `float` (optional).
            :ref:`Focal plane position <lsst.cbp.focal_plane>` of the spot
            formed by the specified beam (x, y mm); defaults to (0, 0).
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam.
        """
        beam = self.maskInfo.asHoleName(beam)
        if focalPlanePos is None:
            focalPlanePos = Point2D(0, 0)
        else:
            focalPlanePos = Point2D(*focalPlanePos)
        focalFieldAngle = self._fieldAngleToFocalPlane.applyInverse(focalPlanePos)
        self.setFocalFieldAngle(pupilPos=pupilPos, focalFieldAngle=focalFieldAngle, beam=beam)

    def setDetectorPos(self, pupilPos, detectorPos=None, detector=None, beam=None):
        """Set the position of a spot on a detector.

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the :ref:`telescope pupil
            <lsst.cbp.pupil_position>` (x, y mm).
        detectorPos : pair of `float` (optional)
            Position of the spot formed by the specified beam
            on the specified detector (x, y pixels);
            defaults to the center of the detector.
        detector : `str` (optional
            Name of detector; defaults to self.config.defaultDetector.
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam.
        """
        beam = self.maskInfo.asHoleName(beam)
        if detector is None:
            detector = self.config.defaultDetector
        detectorInfo = self.cameraGeom[detector]
        if detectorPos is None:
            detectorPos = detectorInfo.getCenter(PIXELS)
        else:
            detectorPos = Point2D(*detectorPos)
        pixelSys = detectorInfo.makeCameraSys(PIXELS)
        pixelsToFieldAngle = self.cameraGeom.getTransform(pixelSys, FIELD_ANGLE)
        focalFieldAngle = pixelsToFieldAngle.applyForward(detectorPos)
        self.setFocalFieldAngle(pupilPos=pupilPos, focalFieldAngle=focalFieldAngle, beam=beam)

    def offsetDetectorPos(self, pupilOffset=None,
                          detectorOffset=None, beam=None):
        """Offset the detector position and/or pupil position of a beam.

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilOffset : pair of `float` (optional)
            Offset of the position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm);
            defaults to (0, 0).
        detectorOffset : pair of `float` (optional)
            Offset of the position of the specified spot
            on the detector it is presently on (x, y pixels);
            defaults to (0, 0).
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam.
        """
        beamInfo = self[beam]
        if not beamInfo.isOnDetector:
            raise RuntimeError("This beam is not on a detector")
        if pupilOffset is None:
            pupilOffset = (0, 0)
        pupilOffset = Extent2D(*pupilOffset)
        if detectorOffset is None:
            detectorOffset = (0, 0)
        detectorOffset = Extent2D(*detectorOffset)
        newPupilPos = beamInfo.pupilPos + pupilOffset
        newDetectorPos = beamInfo.detectorPos + detectorOffset
        self.setDetectorPos(pupilPos=newPupilPos,
                            detectorPos=newDetectorPos,
                            detector=beamInfo.detectorName,
                            beam=beam)

    def offsetFocalPlanePos(self, pupilOffset=None,
                            focalPlaneOffset=None, beam=None):
        """Offset the focal plane position and/or pupil position of a beam.

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilOffset : pair of `float` (optional)
            Offset of the position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm);
            defaults to (0, 0).
        focalPlaneOffset : pair of `float` (optional)
            Offset of the position of the specified spot
            on the :ref:`focal plane <lsst.cbp.focal_plane>` (x, y mm);
            defaults to (0, 0).
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam.
        """
        beamInfo = self[beam]
        if pupilOffset is None:
            pupilOffset = (0, 0)
        pupilOffset = Extent2D(*pupilOffset)
        if focalPlaneOffset is None:
            focalPlaneOffset = (0, 0)
        focalPlaneOffset = Extent2D(*focalPlaneOffset)
        newPupilPos = beamInfo.pupilPos + pupilOffset
        newFocalPlanePos = beamInfo.focalPlanePos + focalPlaneOffset
        self.setFocalPlanePos(pupilPos=newPupilPos,
                              focalPlanePos=newFocalPlanePos,
                              beam=beam)

    def offsetFocalFieldAngle(self, pupilOffset=None,
                              focalFieldAngleOffset=None, beam=None):
        """Offset the focal plane field angle and/or pupil position
        of a beam.

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        Parameters
        ----------
        pupilOffset : pair of `float` (optional)
            Offset of the position of the specified beam on the
            :ref:`telescope pupil <lsst.cbp.pupil_position>` (x, y mm);
            defaults to (0, 0).
        focalFieldAngleOffset : pair of `float` (optional)
            Offset of the :ref:`focal plane field angle
            <lsst.cbp.focal_plane_field_angle>` of the specified beam
            (x, y mm); defaults to (0, 0).
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam.
        """
        beamInfo = self[beam]
        if pupilOffset is None:
            pupilOffset = (0, 0)
        pupilOffset = Extent2D(*pupilOffset)
        if focalFieldAngleOffset is None:
            focalFieldAngleOffset = (0, 0)
        focalFieldAngleOffset = Extent2D(*focalFieldAngleOffset)
        newPupilPos = beamInfo.pupilPos + pupilOffset
        newFocalFieldAngle = beamInfo.focalFieldAngle + focalFieldAngleOffset
        self.setFocalFieldAngle(pupilPos=newPupilPos,
                                focalFieldAngle=newFocalFieldAngle,
                                beam=beam)

    def __getitem__(self, beam):
        """Dict-like access to beam information.

        Get a BeamInfo for a beam specified by integer index or name.
        """
        beam = self.maskInfo.asHoleName(beam)
        return self.getBeamInfo(beam=beam)

    def __iter__(self):
        """Iterator over beam information.

        Do not modify this object during iteration, e.g. by modifying
        attributes or calling set or offset methods.
        """
        for beam in self.beamNames:
            yield self[beam]

    def __len__(self):
        """The number of beams."""
        return self.maskInfo.numHoles

    @property
    def cbpInBounds(self):
        """True if CBP observed altitude is in bounds (read only)."""
        alt = self.cbpAzAltObserved[1]
        return self.config.cbpAltitudeLimits[0] <= alt <= self.config.cbpAltitudeLimits[1]

    @property
    def telInBounds(self):
        """True if telescope observed altitude is in bounds (read only)."""
        alt = self.telAzAltObserved[1]
        return self.config.telAltitudeLimits[0] <= alt <= self.config.telAltitudeLimits[1]

    @property
    def beamNames(self):
        """Beam names, in index order (read only)."""
        return self.maskInfo.holeNames

    @property
    def cbpAzAltObserved(self):
        """Observed az/alt of the CBP (read/write),
        as an lsst.geom.SpherePoint`.
        """
        return SpherePoint(*[self._internalToObserved(
            internal=self._cbpAzAlt[i],
            offset=self.config.cbpAzAltOffset[i],
            scale=self.config.cbpAzAltScale[i]) for i in range(2)])

    @cbpAzAltObserved.setter
    def cbpAzAltObserved(self, cbpObs):
        """Set the observed az/alt of the CBP.

        Parameters
        ----------
        cbpObs : `lsst.geom.SpherePoint`
            Observed az/alt of the CBP.
        """
        self._cbpAzAlt = SpherePoint(*[self._observedToInternal(
            observed=cbpObs[i],
            offset=self.config.cbpAzAltOffset[i],
            scale=self.config.cbpAzAltScale[i]) for i in range(2)])

    @property
    def telAzAltObserved(self):
        """Observed az/alt of the telescope (read/write),
        as an `lsst.geom.SpherePoint`.
        """
        return SpherePoint(*[self._internalToObserved(
            internal=self._telAzAlt[i],
            offset=self.config.telAzAltOffset[i],
            scale=self.config.telAzAltScale[i]) for i in range(2)])

    @telAzAltObserved.setter
    def telAzAltObserved(self, telObs):
        """Set the observed az/alt of the telescope.

        Parameters
        ----------
        telObs : `lsst.geom.SpherePoint`
            Observed az/alt of the telescope.
        """
        self._telAzAlt = SpherePoint(*[self._observedToInternal(
            observed=telObs[i],
            offset=self.config.telAzAltOffset[i],
            scale=self.config.telAzAltScale[i]) for i in range(2)])

    @property
    def telRotObserved(self):
        """Observed angle of the telescope camera rotator (read/write),
        as an `lsst.geom.Angle`.
        """
        return self._internalToObserved(
            internal=self._telRot,
            offset=self.config.telRotOffset,
            scale=self.config.telRotScale)

    @telRotObserved.setter
    def telRotObserved(self, telRotObserved):
        """Set the observed angle of the telescope camera rotator.

        Parameters
        ----------
        telRotObserved : `lsst.geom.Angle`
            The observed angle of the telescope camera rotator.
        """
        self._telRot = self._observedToInternal(
            observed=telRotObserved,
            offset=self.config.telRotOffset,
            scale=self.config.telRotScale)

    @property
    def cbpAzAltInternal(self):
        """Internal az/alt of the CBP (read only),
        as an `lsst.geom.SpherePoint`.

        Primarily intended for testing.
        """
        return self._cbpAzAlt

    @property
    def telAzAltInternal(self):
        """Internal az/alt of the telescope (read only),
        as an `lsst.geom.SpherePoint`.

        Primarily intended for testing.
        """
        return self._telAzAlt

    @property
    def telRotInternal(self):
        """Internal angle of the telescope camera rotator (read only),
        as an `lsst.geom.SpherePoint`.

        Primarily intended for testing.
        """
        return self._telRot

    def setPupilFieldAngle(self, pupilPos, pupilFieldAngle=None, beam=None):
        """Set the pupil field angle and pupil position of a beam.

        Compute new telescope, camera rotator and CBP positions
        and thus update beam info.

        This method is primarily intended for internal use,
        to support the other set methods. It is public so it can be
        unit-tested.

        Parameters
        ----------
        pupilPos : pair of `float`
            Position of the specified beam on the :ref:`telescope pupil
            <lsst.cbp.pupil_position>` (x, y mm).
        pupilFieldAngle : pair of `float` (optional)
            Pupil field angle of specified beam (x, y rad);
            defaults to (0, 0).
        beam : `int` or `str` (optional)
            Name or index of beam; defaults to self.maskInfo.defaultBeam.
        """
        beam = self.maskInfo.asHoleName(beam)
        if pupilFieldAngle is None:
            pupilFieldAngle = Point2D(0, 0)
        else:
            pupilFieldAngle = Point2D(*pupilFieldAngle)
        beamPosAtCtr = coordUtils.computeShiftedPlanePos(pupilPos, pupilFieldAngle,
                                                         -self.config.telPupilOffset)
        beamVectorInCtrPupil = self._computeBeamVectorInCtrPupilFrame(
            beamPosAtCtr=beamPosAtCtr, pupilFieldAngle=pupilFieldAngle)
        cbpVectorInCtrPupil = self._computeCbpVectorInCtrPupilFrame(
            beamPosAtCtr=beamPosAtCtr,
            beamVectorInCtrPupil=beamVectorInCtrPupil)

        telAzAlt = coordUtils.computeAzAltFromBasePupil(
            vectorBase=self.config.cbpPosition,
            vectorPupil=cbpVectorInCtrPupil)

        beamVectorBase = coordUtils.convertVectorFromPupilToBase(
            vectorPupil=beamVectorInCtrPupil,
            pupilAzAlt=telAzAlt,
        )
        beamFieldAngleCbp = self._getBeamCbpFieldAngle(beam)
        beamUnitVectorCbpPupil = coordUtils.fieldAngleToVector(beamFieldAngleCbp, self.config.cbpFlipX)
        cbpAzAlt = coordUtils.computeAzAltFromBasePupil(
            vectorBase=-beamVectorBase,
            vectorPupil=beamUnitVectorCbpPupil,
        )

        self._telRot = self._computeCameraRotatorAngle(telAzAlt=telAzAlt, cbpAzAlt=cbpAzAlt)
        self._telAzAlt = telAzAlt
        self._cbpAzAlt = cbpAzAlt

    def _computeCameraRotatorAngle(self, telAzAlt, cbpAzAlt):
        """Compute the internal camera rotator angle needed for a given
        telescope and CBP pointing.

        Parameters
        ----------
        telAzAlt : `lsst.geom.SpherePoint`
            Telescope internal azimuth and altitude.
        cbpAzAlt : `lsst.geom.SpherePoint`
            CBP internal azimuth and altitude.

        Returns
        -------
        rotatorAangle : `lsst.geom.Angle`
            Internal camera rotator angle.
        """
        # compute focal plane position, ignoring self._telRot,
        # for two holes separated by x in the CBP equidistant from the center;
        # compute the angle that would make the spots line up with the x axis in the focal plane
        ctrHolePos = Point2D(0, 0)
        holeDelta = Extent2D(*coordUtils.getFlippedPos(self._holeDelta, flipX=self.config.cbpFlipX))
        holePos1 = ctrHolePos - holeDelta
        holePos2 = ctrHolePos + holeDelta
        pupilUnitVector1 = self._computeTelPupilUnitVectorFromHolePos(holePos1, telAzAlt=telAzAlt,
                                                                      cbpAzAlt=cbpAzAlt)
        pupilUnitVector2 = self._computeTelPupilUnitVectorFromHolePos(holePos2, telAzAlt=telAzAlt,
                                                                      cbpAzAlt=cbpAzAlt)
        # Rotation is done in a right-handed system, regardless of telFlipX
        pupilFieldAngle1 = coordUtils.vectorToFieldAngle(pupilUnitVector1, flipX=False)
        pupilFieldAngle2 = coordUtils.vectorToFieldAngle(pupilUnitVector2, flipX=False)
        focalPlane1 = self._fieldAngleToFocalPlane.applyForward(Point2D(*pupilFieldAngle1))
        focalPlane2 = self._fieldAngleToFocalPlane.applyForward(Point2D(*pupilFieldAngle2))
        deltaFocalPlane = np.subtract(focalPlane2, focalPlane1)
        return -math.atan2(deltaFocalPlane[1], deltaFocalPlane[0])*radians

    def _getBeamCbpFieldAngle(self, beam):
        """Return the field angle of the specified beam
        in the CBP pupil frame.

        Parameters
        ----------
        beam : `int`, `str` or None
            Name or index of beam; if None then
            ``self.maskInfo.defaultBeam``.

        Returns
        -------
        fieldAngle : a pair of floats, in radians
            Field angle of the specified beam in the CBP pupil frame
            (x, y radians).
        """
        holePos = self.maskInfo.getHolePos(beam)
        return tuple(math.atan(pos / self.config.cbpFocalLength) for pos in holePos)

    def _computeBeamVectorInCtrPupilFrame(self, beamPosAtCtr, pupilFieldAngle):
        """Compute the beam vector to the CBP in the centered pupil frame.

        Parameters
        ----------
        beamPosAtCtr : pair of `float`
            Position of beam on centered pupil (x, y mm).
        pupilFieldAngle : pair of `float`
            incident angle of beam on pupil (x, y rad).

        Returns
        -------
        beamVectorinCtrPupil : `numpy.array` of 3 `float`
            Beam vector in telescope centered pupil frame (mm).
        """
        # beamPosVec is a vector in the telescope pupil frame from the center of the centered pupil plane
        # to the point on that plane specified by beamPosAtCtr; this vector lies in the centered pupil plane
        beamPosVec = coordUtils.pupilPositionToVector(beamPosAtCtr, self.config.telFlipX)
        beamUnitVec = coordUtils.fieldAngleToVector(pupilFieldAngle, self.config.telFlipX)
        abyz = beamPosVec[1]*beamUnitVec[1] + beamPosVec[2]*beamUnitVec[2]
        beamPosMag = np.linalg.norm(beamPosAtCtr)
        cbpDistance = self.config.cbpDistance
        beamLength = -abyz + math.sqrt(math.fsum((cbpDistance**2, abyz**2, -beamPosMag**2)))
        return beamLength*np.array(beamUnitVec)

    def _computeCbpVectorInCtrPupilFrame(self, beamPosAtCtr, beamVectorInCtrPupil):
        """Compute a vector from telescope to CBP in the telescope's
        centered pupil frame.

        Parameters
        ----------
        beamPosAtCtr : pair of `float`
            Position of beam on centered pupil (x, y mm).
        beamVectorInCtrPupil : `numpy.array` of 3 `float`
            Vector in the telescope pupil frame from ``beamPosAtCtr`
            to the center of the CBP (mm).

        Returns
        -------
        cbpVectorInCtrPupilFrame : `numpy.array` of 3 `float`
            Vector in the telescope pupil frame from the center
            of the pupil frame to the center of the CBP (mm).
        """
        # beamPosVec is a vector in the telescope pupil frame from the center of the centered pupil plane
        # to the point on that plane specified by beamPosAtCtr; this vector lies in the centered pupil plane
        beamPosVec = coordUtils.pupilPositionToVector(beamPosAtCtr, self.config.telFlipX)
        return beamVectorInCtrPupil + beamPosVec

    def _rotateFocalPlaneToPupil(self, focalPlane):
        """Rotate a position or field angle in telescope focal plane
        orientation to telescope pupil orientation.

        Parameters
        ----------
        focalPlane : pair of `float`
            Focal plane position or field angle (x, y any units).

        Returns
        -------
        pupil : pair of `float`
            ``focalPlane`` rotated to the pupil plane (x, y, same units).
        """
        return coordUtils.rotate2d(pos=focalPlane, angle=-self._telRot)

    def _rotatePupilToFocalPlane(self, pupil):
        """Rotate a position or field angle in telescope pupil orientation
        to one in telescope focal plane orientation.

        Parameters
        ----------
        pupil : pair of `float`
            Pupil plane position or field angle (x, y any units)

        Returns
        -------
        focalPlane : pair of `float`
            ``pupil`` rotated to the focal plane (x, y same units)
        """
        return coordUtils.rotate2d(pos=pupil, angle=self._telRot)

    def getBeamInfo(self, beam, *, holePos=None):
        """Get beam information for a beam from a specified CBP
        beam or hole position.

        You may specify a hole position. This can be useful for unit
        tests and "what if" scenarios.

        Parameters
        ----------
        beam : `str`
            Beam name; ignored if holePos specified.
        holePos : pair of `float` (optional)
            Hole position on CBP mask (x, y mm);
            defaults to the actual hole position of the named beam.

        Returns
        -------
        beamInfo : `lsst.cbp.BeamInfo`
            Information about the specified beam.
        """
        if holePos is None:
            holePos = self.maskInfo.getHolePos(beam)

        # compute focal plane field angle of the beam
        beamPupilUnitVector = self._computeTelPupilUnitVectorFromHolePos(holePos, telAzAlt=self._telAzAlt,
                                                                         cbpAzAlt=self._cbpAzAlt)
        pupilFieldAngle = coordUtils.vectorToFieldAngle(beamPupilUnitVector, self.config.telFlipX)
        focalFieldAngle = self._rotatePupilToFocalPlane(pupilFieldAngle)

        # compute focal plane position of the beam
        focalPlanePos = self._fieldAngleToFocalPlane.applyForward(Point2D(*pupilFieldAngle))
        isOnFocalPlane = math.hypot(*focalPlanePos) < self.config.telFocalPlaneDiameter

        # vector from telescope centered pupil to actual pupil
        pupilOffset = np.array((self.config.telPupilOffset, 0, 0), dtype=float)

        # compute the pupil position of the beam on the actual telescope pupil
        # (first compute on the centered pupil, then subtract the pupil offset)
        cbpPositionPupil = coordUtils.convertVectorFromBaseToPupil(
            vectorBase=self.config.cbpPosition,
            pupilAzAlt=self._telAzAlt,
        ) - pupilOffset
        pupilNormalVector = np.array((1, 0, 0), dtype=float)
        pupilDistance = np.dot(cbpPositionPupil, pupilNormalVector) / \
            np.dot(beamPupilUnitVector, pupilNormalVector)
        pupilPosVector = cbpPositionPupil - beamPupilUnitVector*pupilDistance
        # the x component should be zero; y, z -> plane position ±x, y
        pupilPos = coordUtils.getFlippedPos((pupilPosVector[1], pupilPosVector[2]),
                                            flipX=self.config.telFlipX)
        pupilPosDiameter = math.hypot(*pupilPos)
        isOnPupil = self.config.telPupilObscurationDiameter <= pupilPosDiameter and \
            pupilPosDiameter <= self.config.telPupilDiameter
        return BeamInfo(
            cameraGeom=self.cameraGeom,
            name=beam,
            holePos=holePos,
            isOnPupil=isOnPupil,
            isOnFocalPlane=isOnFocalPlane,
            focalPlanePos=focalPlanePos,
            pupilPos=pupilPos,
            focalFieldAngle=focalFieldAngle,
            pupilFieldAngle=pupilFieldAngle,
        )

    def _computeCbpPupilUnitVectorFromHolePos(self, holePos):
        """Compute the CBP pupil unit vector of a beam.

        Parameters
        ----------
        holePos : pair of `float`
            Hole position on CBP mask (x, y mm).

        Returns
        -------
        cbpPupilUnitVector : `numpy.array` of 3 `float`
            The direction of the beam emitted by the CBP
            as a unit vector in the CBP pupil frame.
        """
        beamCbpFieldAngle = [math.atan(pos/self.config.cbpFocalLength) for pos in holePos]
        return coordUtils.fieldAngleToVector(beamCbpFieldAngle, self.config.cbpFlipX)

    def _computeTelPupilUnitVectorFromHolePos(self, holePos, telAzAlt, cbpAzAlt):
        """Compute the telescope pupil unit vector of a beam
        from its CBP hole position.

        Parameters
        ----------
        holePos : pair of `float`
            Hole position on CBP mask (x, y mm).
        telAzAlt : `lsst.geom.SpherePoint`
            Telescope internal azimuth and altitude.
        cbpAzAlt : `lsst.geom.SpherePoint`
            CBP internal azimuth and altitude.

        Returns
        -------
        telPupilUnitVector : `numpy.array` of 3 `float`
            The direction of the beam received by the telescope
            as a unit vector in the telescope pupil frame.
        """
        beamCbpPupilUnitVector = self._computeCbpPupilUnitVectorFromHolePos(holePos)
        beamCbpBaseUnitVector = coordUtils.convertVectorFromPupilToBase(
            vectorPupil=beamCbpPupilUnitVector,
            pupilAzAlt=cbpAzAlt,
        )
        beamTelBaseUnitVector = -beamCbpBaseUnitVector
        return coordUtils.convertVectorFromBaseToPupil(
            vectorBase=beamTelBaseUnitVector,
            pupilAzAlt=telAzAlt,
        )

    def _observedToInternal(self, observed, offset, scale):
        """Convert an observed angle into an internal angle.

        Computes ``internal angle = (observed angle / scale) - offset``.

        Parameters
        ----------
        observed : `lsst.geom.Angle`
            Observed angle.
        offset : `lsst.geom.Angle`
            Offset.
        scale : `float`
            Scale.

        Returns
        -------
        internal : `lsst.geom.Angle`
            Internal angle.
        """
        return (observed - offset)/scale

    def _internalToObserved(self, internal, offset, scale):
        """Convert an internal angle into an observed angle.

        Computes ``observed angle = (internal angle + offset) * scale``.

        Parameters
        ----------
        internal : `lsst.geom.Angle`
            Internal angle
        offset : `lsst.geom.Angle`
            Offset
        scale : `float`
            Scale

        Returns
        -------
        observed : `lsst.geom.Angle`
            Observed angle.
        """
        return internal*scale + offset
