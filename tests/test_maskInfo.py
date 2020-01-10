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

import unittest

import lsst.utils.tests
import lsst.geom
from lsst.cbp import MaskInfo


class MaskInfoTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.holePositions = (
            (0, 0),
            (1, -2),
            (-3, 4),
            (5, 6),
        )
        # A useful set of hole names that are not the default.
        self.holeNames = ["beam{}".format(i) for i in range(len(self.holePositions))]

    def testBasics(self):
        name = "test mask"
        defaultHole = 2
        self.assertEqual(len(self.holeNames), len(self.holePositions))
        mi = MaskInfo(
            name=name,
            defaultHole=defaultHole,
            holePositions=self.holePositions,
            holeNames=self.holeNames,
        )
        self.assertEqual(mi.name, name)
        self.assertEqual(mi.numHoles, len(self.holePositions))
        self.assertEqual(list(mi.holeNames), self.holeNames)
        self.assertEqual(mi.asHoleName(None), self.holeNames[defaultHole])
        self.assertEqual(mi.getHolePos(None), lsst.geom.Point2D(*self.holePositions[defaultHole]))
        for i in range(-mi.numHoles, mi.numHoles):
            holeName = mi.asHoleName(i)
            desiredHoleName = self.holeNames[i]
            self.assertEqual(holeName, desiredHoleName)
            desiredHolePosition = lsst.geom.Point2D(*self.holePositions[i])
            self.assertEqual(mi.getHolePos(holeName), desiredHolePosition)
            self.assertEqual(mi.getHolePos(i), desiredHolePosition)

        with self.assertRaises(LookupError):
            mi.asHoleName("invalidName")
        with self.assertRaises(LookupError):
            mi.asHoleName(len(self.holePositions))
        with self.assertRaises(LookupError):
            mi.asHoleName(-len(self.holePositions) - 1)

    def testConstructorErrors(self):
        # defaultHole must not be None
        with self.assertRaises(ValueError):
            MaskInfo(name="test", defaultHole=None, holePositions=self.holePositions)
        with self.assertRaises(ValueError):
            MaskInfo(name="test", defaultHole=None, holePositions=self.holePositions,
                     holeNames=self.holeNames)

        # defaultHole must be a valid name or index,
        # but unlike None, this raises LookupError
        # so it has to be a separate test.
        for invalidHole in (len(self.holePositions), "bad"):
            with self.assertRaises(LookupError):
                MaskInfo(name="test", defaultHole=invalidHole, holePositions=self.holePositions)
            with self.assertRaises(LookupError):
                MaskInfo(name="test", defaultHole=invalidHole, holePositions=self.holePositions,
                         holeNames=self.holeNames)

        # The number of hole names must match the number of hole positions.
        holeNamesTooShort = self.holeNames[0: -1]
        with self.assertRaises(ValueError):
            MaskInfo(name="test", defaultHole=0, holePositions=self.holePositions,
                     holeNames=holeNamesTooShort)

        holeNamesTooLong = self.holeNames + ["extra"]
        with self.assertRaises(ValueError):
            MaskInfo(name="test", defaultHole=0, holePositions=self.holePositions,
                     holeNames=holeNamesTooLong)

    def testDefaultNames(self):
        mi = MaskInfo(
            name="test",
            defaultHole=1,
            holePositions=self.holePositions,
        )
        self.assertEqual(mi.numHoles, len(self.holePositions))
        self.assertEqual(list(mi.holeNames), [str(i) for i in range(len(self.holePositions))])


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
