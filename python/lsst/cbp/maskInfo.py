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
"""MaskInfo class: information about a CBP mask"""

__all__ = ["MaskInfo"]

from collections import OrderedDict

from lsst.geom import Point2D


class MaskInfo:
    """Information about a CBP mask.

    Parameters
    ----------
    name : `str`
        Name of mask.
    defaultHole : `str` or `int`
        Name or index of default hole.
    holePositions : `iterable` of (`float`, `float`)
        Position of a fiducial point for each hole in the mask (x, y mm).
    holeNames : `iterable` of `str` or None
        Name of each hole in the mask, in the same order as
        ``holePositions``.
        If None then set name = str(index) for each hole.

    Raises
    ------
    ValueError
        If ``defaultHole`` is `None`.
    LookupError
        If ``defaultHole`` is not a valid name or index.
    ValueError
        If ``holeNames`` is not `None` and has a different length
        than ``holePositions``.

    Notes
    -----
    **Attributes**

    name : `str`
        Name of mask.
    defaultBeam : `str`
        Name of default hole or beam.
    """

    def __init__(self, name, defaultHole, holePositions, holeNames=None):
        self.name = name
        if holeNames is not None:
            if len(holePositions) != len(holeNames):
                raise ValueError("Number of hole positions = {} != Number of hole names = {}".format(
                    len(holePositions), len(holeNames)
                ))
        else:
            holeNames = [str(i) for i in range(len(holePositions))]
        self._holePosDict = OrderedDict()
        for holeName, holePos in zip(holeNames, holePositions):
            self._holePosDict[holeName] = Point2D(*holePos)
        # parse "defaultBeam" after "holes", so we can check the default
        if defaultHole is None:
            raise ValueError("defaultHole cannot be None")
        self.defaultBeam = self.asHoleName(defaultHole)

    def asHoleName(self, hole):
        """Read a hole index, name or None as a name, and validate it.

        Parameters
        ----------
        hole : `int`, `str` or `None`
            If `None`  then return the default beam.
            If an integer index then return the corresponding beam name.
            If a string then return unchanged.

        Raises
        ------
        `LookupError` if beam is an integer and is out of range
        or if beam is a string and the name is unknown.
        """
        if hole is None:
            return self.defaultBeam
        if isinstance(hole, int):
            return [name for name in self._holePosDict][hole]
        if hole not in self._holePosDict:
            raise KeyError("Unknown name {!r}".format(hole))
        return hole

    def getHolePos(self, beam):
        """Return the position of a hole in focal plane x,y mm.

        Parameters
        ----------
        hole : `int`, `str` or `None`
            If `None`  then use the default beam.
            If an integer index then use the corresponding beam name.
        """
        return self._holePosDict[self.asHoleName(beam)]

    @property
    def numHoles(self):
        """The number of holes (read only)"""
        return len(self._holePosDict)

    @property
    def holeNames(self):
        """An iterable of hole names, in index order (read only)"""
        return self._holePosDict.keys()
