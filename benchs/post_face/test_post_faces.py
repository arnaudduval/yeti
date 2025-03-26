# Copyright 2021 Arnaud Duval

# This file is part of Yeti.
#
# Yeti is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# Yeti is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Yeti. If not, see <https://www.gnu.org/licenses/>

# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
From a 3D geometry, a VTU file is generated with a field indicating face index
"""

import os

import numpy as np

# pylint: disable=c-extension-no-member


from yeti_iga.preprocessing.igaparametrization import IGAparametrization
import yeti_iga.postprocessing.postproc as pp


def test_post_faces(tmp_path):
    """
    Generate VTU file
    """
    filename = 'testU5'

    # Read data and create IGAparametrization object
    script_dir = os.path.dirname(os.path.realpath(__file__))
    iga_model = IGAparametrization(filename=f'{script_dir}/{filename}')

    nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
    nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

    nb_ref[:, 0] = np.array([3, 3, 3])
    nb_deg[:, 0] = np.array([1, 1, 1])

    nb_ref[:, 1] = np.array([3, 4, 4])
    nb_deg[:, 1] = np.array([1, 1, 1])

    nb_ref[:, 2] = np.array([3, 3, 3])
    nb_deg[:, 2] = np.array([1, 1, 1])

    iga_model.refine(nb_ref, nb_deg)

    # Output for faces
    pp.generate_faces_vtu(**iga_model.get_inputs4postproc_faces_vtu(
        filename + 'faces', nb_ref=np.array([1, 1, 1]), output_path=tmp_path))


if __name__ == '__main__':
    test_post_faces('results')
