import numpy as np

from yeti_iga import IgaModel
from yeti_iga.exchange import load_geomdl_json

def test_import_multipatch_volume(tmp_path):

    filename = 'mvol.json'

    patchs = load_geomdl_json(f'{tmp_path}/{filename}')

    model = IgaModel("3D solid")
    for patch in patchs:
        model.add_patch(patch)

    model.write_solution_vtu(np.zeros_like(model.idof_free), 'geometry.vtu')
    model.write_control_mesh_vtk('control_net.vtk')

if __name__ == '__main__':
    test_import_multipatch_volume('.')