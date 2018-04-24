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

import math
import unittest

import numpy as np

import lsst.utils.tests
from lsst.geom import degrees
from lsst.cbp import CoordinateConverterConfig


class CoordinateConverterTestCase(lsst.utils.tests.TestCase):

    def setUp(self):
        self.argDict = dict(
            telPupilDiameter=3000,
            telPupilOffset=-0.134,
            telPupilObscurationDiameter=333,
            telFocalPlaneDiameter=550,
            telFlipX=True,
            telAzimuthOffsetDeg=1.1,
            telAzimuthScale=-1,  # must be ±1
            telAltitudeOffsetDeg=2.2,
            telAltitudeScale=-2.3,
            telAltitudeLimitsDeg=(-3.1, 3.2),
            telRotOffsetDeg=4.1,
            telRotScale=-1,  # must be ±1
            defaultDetector=5,
            cbpPosition=(6.1, 6.2, 6.3),
            cbpFocalLength=7,
            cbpFlipX=False,
            cbpAzimuthOffsetDeg=8.1,
            cbpAzimuthScale=-1,  # must be ±1
            cbpAltitudeOffsetDeg=9.1,
            cbpAltitudeScale=9.2,
            cbpAltitudeLimitsDeg=(10.1, 10.2),
        )
        # names of arguments that show up as a field of the same name and type
        # (ignoring list vs. tuple)
        self.matchedNames = [
            "telPupilDiameter",
            "telPupilOffset",
            "telPupilObscurationDiameter",
            "telFocalPlaneDiameter",
            "telFlipX",
            "telRotScale",
            "defaultDetector",
            "cbpPosition",
            "cbpFocalLength",
            "cbpFlipX",
        ]

    def testBasics(self):
        argDict = self.argDict
        config = CoordinateConverterConfig(**argDict)
        self.checkValues(argDict, config)

    def testDefaults(self):
        argDict = self.argDict
        minimalArgDict = argDict.copy()
        defaultArgDict = dict(
            telAltitudeOffsetDeg=0,
            telAltitudeScale=1,
            cbpAltitudeOffsetDeg=0,
            cbpAltitudeScale=1,
        )
        for name in defaultArgDict:
            del minimalArgDict[name]
            argDict[name] = defaultArgDict[name]
        config = CoordinateConverterConfig(**minimalArgDict)
        self.checkValues(argDict, config)

    def testLengthErrors(self):
        for name, value in self.argDict.items():
            if not isinstance(value, tuple):
                continue
            argDict = self.argDict.copy()
            argDict[name] = value[0:-1]  # too few values
            with self.assertRaises(ValueError):
                CoordinateConverterConfig(**argDict)
            argDict[name] = value + (3.1,)  # too many values
            with self.assertRaises(ValueError):
                CoordinateConverterConfig(**argDict)

    def testScaleErrors(self):
        """Rotator and azimuth scales must be ±1 degree
        """
        for name in ("telAzimuthScale", "telRotScale", "cbpAzimuthScale"):
            argDict = self.argDict.copy()
            argDict[name] = 1.01
            with self.assertRaises(ValueError):
                CoordinateConverterConfig(**argDict)

    def checkValues(self, argDict, config):
        """Check the values in a CoordinateConverterConfig

        Parameters
        ----------
        argDict : `dict`
            Dictionary of arguments used to construct ``config``,
            plus values for defaulted arguments (if any)
        config : `lsst.cbp.CoordinateConverterConfig`
            The configuration to check
        """
        for name in self.matchedNames:
            argValue = argDict[name]
            fieldValue = getattr(config, name)
            if isinstance(argValue, tuple):
                self.assertEqual(len(argValue), len(fieldValue))
                for a, f in zip(argValue, fieldValue):
                    self.assertEqual(a, f)
            else:
                self.assertEqual(argDict[name], getattr(config, name))
        self.assertEqual(config.telRotOffset, argDict["telRotOffsetDeg"]*degrees)
        desiredAzAltScale = (argDict["telAzimuthScale"],
                             argDict["telAltitudeScale"])
        self.assertEqual(config.telAzAltScale, desiredAzAltScale)
        desiredAzAltOffset = (argDict["telAzimuthOffsetDeg"]*degrees,
                              argDict["telAltitudeOffsetDeg"]*degrees)
        self.assertEqual(config.telAzAltOffset, desiredAzAltOffset)
        desiredAzAltScale = (argDict["cbpAzimuthScale"],
                             argDict["cbpAltitudeScale"])
        self.assertEqual(config.cbpAzAltScale, desiredAzAltScale)
        desiredAzAltOffset = (argDict["cbpAzimuthOffsetDeg"]*degrees,
                              argDict["cbpAltitudeOffsetDeg"]*degrees)
        self.assertEqual(config.cbpAzAltOffset, desiredAzAltOffset)

        self.checkCbpDistance(config)

        # try setting a different cbpPosition; cbpPosition and cbpDistance should update accordingly
        cbpPosition2 = config.cbpPosition + np.array((1.5, -3.2, 9.4))
        config.cbpPosition = cbpPosition2
        np.testing.assert_equal(config.cbpPosition, cbpPosition2)
        self.checkCbpDistance(config)

    def checkCbpDistance(self, config):
        """Check that cbpPosition and cbpDistance are as desired

        Parameters
        ----------
        config : `lsst.cbp.CoordinateConverterConfig`
            The configuration to check
        """

        # this is less accurate than np.linalg.norm, as used by CoordinateConverterConfig,
        # but is fine for unit tests and it's nice to have a different implementation
        cbpPosition = config.cbpPosition
        desiredCbpDistance = math.sqrt(cbpPosition[0]**2 + cbpPosition[1]**2 + cbpPosition[2]**2)
        self.assertAlmostEqual(config.cbpDistance, desiredCbpDistance)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
