from typing import Sequence
import numpy as np


def kron_axis0(arrays: Sequence[np.ndarray]) -> np.ndarray:
    assert isinstance(arrays, (list, tuple))
    assert len(arrays) > 0, ValueError("Empty list")

    if len(arrays) == 1:
        return arrays[0]

    m = arrays[0].shape[1]
    for curr in arrays:
        assert curr.shape[1] == m

    result = arrays[0]
    for curr in arrays[1:]:
        result = np.reshape(curr[None, :, :] * result[:, None, :], (-1, m), order="F")

    return result
