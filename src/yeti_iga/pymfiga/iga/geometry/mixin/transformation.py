from .base import BasePatchIGA
from typing import Optional
import numpy as np
import logging

logger = logging.getLogger("SRC.IGA.GEOMETRY")


class TransformationMixin(BasePatchIGA):
    def translate(self, vector: np.ndarray):
        """Translate control points by a vector in-place."""
        logger.info("Translating single patch")
        v = np.asarray(vector, dtype=float)
        for i in range(min(self.ndim, v.size)):
            self.ctrlpts[i] += v[i]

    def reflect(self, plane: str):
        """Reflect across planes: 'xy', 'xz', 'yz' (order-insensitive)."""
        logger.info("Reflecting single patch")
        plane = plane.lower()
        assert plane in {"xy", "yx", "xz", "zx", "yz", "zy"}
        setofaxis = {"x": 0, "y": 1, "z": 2}
        missing = (set("xyz") - set(plane)).pop()  # the normal axis
        axis = setofaxis[missing]
        self.ctrlpts[axis] *= -1.0

    def rotate(self, axis: str, angle: float, point: Optional[np.ndarray] = None):
        """
        Rotate the control points around a given axis ('x', 'y', or 'z') by a given angle (in radians),
        around a given point (default is the origin).
        """
        logger.info("Single patch has been rotated")
        assert self.ctrlpts.shape[0] >= 3, "Rotation only in 3D"
        assert axis in ["x", "y", "z"], "Axis must be 'x', 'y', or 'z'"
        idx = {"x": 0, "y": 1, "z": 2}[axis]
        if point is None:
            point = np.zeros(self._ndim_default)
        # Move to rotation center
        self.ctrlpts[:3, :] -= point[:3, np.newaxis]
        # Build rotation matrix
        c, s = np.cos(angle), np.sin(angle)
        R = np.eye(3)
        if idx == 0:  # x-axis
            R[1, 1] = c
            R[1, 2] = -s
            R[2, 1] = s
            R[2, 2] = c
        elif idx == 1:  # y-axis
            R[0, 0] = c
            R[0, 2] = s
            R[2, 0] = -s
            R[2, 2] = c
        elif idx == 2:  # z-axis
            R[0, 0] = c
            R[0, 1] = -s
            R[1, 0] = s
            R[1, 1] = c
        # Apply rotation
        self.ctrlpts[:3, :] = R @ self.ctrlpts[:3, :]
        # Move back
        self.ctrlpts[:3, :] += point[:3, np.newaxis]
