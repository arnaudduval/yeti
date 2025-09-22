"""
Definition of the geometry of a stiffened cylinder
"""

import numpy as np

RADIUS = 318
HEIGHT = 600
TEMPLATE_INP_FILE = 'stiff_geom_2xhalf.inp.template'
NB_STIFFENERS_BY_HALF = 18
STIFF_THICKNESS = 1.6

# Set design variables values for X distributed stiffeners
def set_X_stiffeners(nb_stiffeners):
    if nb_stiffeners == 8:
        return np.array([0.0,     0.5,    0.0,   0.5,
                        0.0,   0.5,  0.5,   0.0,
                        0.5,    1.0,   0.0,   0.5,
                        0.5,   1.0,  0.5,   0.0,
                        0.0,     0.5,    1.0,   0.5,
                        0.0,   0.5,  0.5,   1.0,
                        0.5,    1.0,   1.0,   0.5,
                        0.5,   1.0,  0.5,   1.0])

    if nb_stiffeners == 24:
        return np.array([ 0.05, 0.35, 0.05, 0.05,
                            0.35, 0.65, 0.05, 0.05,
                            0.65, 0.95, 0.05, 0.05,
                            0.05, 0.35, 0.35, 0.35,
                            0.35, 0.65, 0.35, 0.35,
                            0.65, 0.95, 0.35, 0.35,
                            0.05, 0.35, 0.65, 0.65,
                            0.35, 0.65, 0.65, 0.65,
                            0.65, 0.95, 0.65, 0.65,
                            0.05, 0.35, 0.95, 0.95,
                            0.35, 0.65, 0.95, 0.95,
                            0.65, 0.95, 0.95, 0.95,
                            0.05, 0.05, 0.05, 0.35,
                            0.05, 0.05, 0.35, 0.65,
                            0.05, 0.05, 0.65, 0.95,
                            0.35, 0.35, 0.05, 0.35,
                            0.35, 0.35, 0.35, 0.65,
                            0.35, 0.35, 0.65, 0.95,
                            0.65, 0.65, 0.05, 0.35,
                            0.65, 0.65, 0.35, 0.65,
                            0.65, 0.65, 0.65, 0.95,
                            0.95, 0.95, 0.05, 0.35,
                            0.95, 0.95, 0.35, 0.65,
                            0.95, 0.95, 0.65, 0.95])

    if nb_stiffeners == 18:
        return np.array([0.0,   1./3.,    0.0,   1./3.,
                        0.0,   1./3.,    1./3.,   0.0,
                        0.0,   1./3.,    1./3.,   2./3.,
                        0.0,   1./3.,    2./3.,   1./3.,
                        0.0,   1./3.,    2./3.,   1.0,
                        0.0,   1./3.,    1.0,   2./3.,

                        1./3., 2./3.,    0.0,   1./3.,
                        1./3., 2./3.,    1./3.,   0.0,
                        1./3., 2./3.,    1./3.,   2./3.,
                        1./3., 2./3.,    2./3.,   1./3.,
                        1./3., 2./3.,    2./3.,   1.0,
                        1./3., 2./3.,    1.0,   2./3.,

                        2./3., 1.0,    0.0,   1./3.,
                        2./3., 1.0,    1./3.,   0.0,
                        2./3., 1.0,    1./3.,   2./3.,
                        2./3., 1.0,    2./3.,   1./3.,
                        2./3., 1.0,    2./3.,   1.0,
                        2./3., 1.0,    1.0,   2./3.,
                        ])

    raise Exception(f"Can not handle {nb_stiffeners} stiffeners")

# Set design variables values for H distributed stiffeners
def set_H_stiffeners(nb_stiffeners):
    # TODO Mapping is non linear
    print(nb_stiffeners)
    if nb_stiffeners % 2 != 0:
        raise Exception("Can only handle even number of stiffeners")

    x = np.array([], dtype=float)

    # Vertical stiffeners
    for istiff in range(int(nb_stiffeners/2)):
        interval = 1./(nb_stiffeners/2)
        u = interval/2. + istiff*interval
        x = np.append(x, np.array([u, u, 0., 1.]))

    # Horizontal stiffeners
    for istiff in range(int(nb_stiffeners/2)):
        interval = 1./(nb_stiffeners/2)
        v = interval/2. + istiff*interval
        x = np.append(x, np.array([0., 1., v, v]))

    return x

# Write inp file from a template
def write_inp_file(filename):
# Write inp file
    with open(TEMPLATE_INP_FILE, 'r') as template:
        with open(filename, 'w') as inpfile:
            for line in template.readlines():
                if line.startswith("#include nodes"):
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        # 28 nodes
                        # U30s
                        inpfile.write(f'{69+istiff*28}, 0.1, 0.4, 0.5\n')
                        inpfile.write(f'{70+istiff*28}, 0.9, 0.6, 0.5\n')
                        inpfile.write(f'{71+istiff*28}, 0.1, 0.4, 1.0\n')
                        inpfile.write(f'{72+istiff*28}, 0.9, 0.6, 1.0\n')

                        inpfile.write(f'{73+istiff*28}, 0.1, 0.4, 0.5\n')
                        inpfile.write(f'{74+istiff*28}, 0.9, 0.6, 0.5\n')
                        inpfile.write(f'{75+istiff*28}, 0.1, 0.4, 1.0\n')
                        inpfile.write(f'{76+istiff*28}, 0.9, 0.6, 1.0\n')
                        # U00
                        inpfile.write(f'{77+istiff*28}, 0.1, 0.4, 0.0\n')
                        inpfile.write(f'{78+istiff*28}, 0.9, 0.6, 0.0\n')
                        inpfile.write(f'{79+istiff*28}, 0.0, 1.0, 0.0\n')
                        inpfile.write(f'{80+istiff*28}, 1.0, 1.0, 0.0\n')
                        inpfile.write(f'{81+istiff*28}, 0.1, 0.4, 0.0\n')
                        inpfile.write(f'{82+istiff*28}, 0.9, 0.6, 0.0\n')
                        inpfile.write(f'{83+istiff*28}, 0.0, 1.0, 0.0\n')
                        inpfile.write(f'{84+istiff*28}, 1.0, 1.0, 0.0\n')

                        inpfile.write(f'{85+istiff*28}, 0.1, 0.4, 0.0\n')
                        inpfile.write(f'{86+istiff*28}, 0.9, 0.6, 0.0\n')
                        inpfile.write(f'{87+istiff*28}, 0.0, 1.0, 0.0\n')
                        inpfile.write(f'{88+istiff*28}, 1.0, 1.0, 0.0\n')
                        inpfile.write(f'{89+istiff*28}, 0.1, 0.4, 0.0\n')
                        inpfile.write(f'{90+istiff*28}, 0.9, 0.6, 0.0\n')
                        inpfile.write(f'{91+istiff*28}, 0.0, 1.0, 0.0\n')
                        inpfile.write(f'{92+istiff*28}, 1.0, 1.0, 0.0\n')

                        #U4
                        inpfile.write(f'{93+istiff*28}, 0.0, 0.0, 0.0\n')
                        inpfile.write(f'{94+istiff*28}, 0.0, 0.0, 0.0\n')
                        inpfile.write(f'{95+istiff*28}, 0.0, 0.0, 0.0\n')
                        inpfile.write(f'{96+istiff*28}, 0.0, 0.0, 0.0\n')
                elif line.startswith('#include curve elements'):
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        inpfile.write(f'{23+istiff*14}, {78+istiff*28}, {77+istiff*28}\n')
                        inpfile.write(f'{24+istiff*14}, {80+istiff*28}, {79+istiff*28}\n')
                        inpfile.write(f'{25+istiff*14}, {82+istiff*28}, {81+istiff*28}\n')
                        inpfile.write(f'{26+istiff*14}, {84+istiff*28}, {83+istiff*28}\n')

                        inpfile.write(f'{27+istiff*14}, {86+istiff*28}, {85+istiff*28}\n')
                        inpfile.write(f'{28+istiff*14}, {88+istiff*28}, {87+istiff*28}\n')
                        inpfile.write(f'{29+istiff*14}, {90+istiff*28}, {89+istiff*28}\n')
                        inpfile.write(f'{30+istiff*14}, {92+istiff*28}, {91+istiff*28}\n')

                elif line.startswith('#include lagrange elements'):
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        inpfile.write(f'{31+istiff*14}, {93+istiff*28}\n')
                        inpfile.write(f'{32+istiff*14}, {94+istiff*28}\n')
                        inpfile.write(f'{33+istiff*14}, {95+istiff*28}\n')
                        inpfile.write(f'{34+istiff*14}, {96+istiff*28}\n')
                elif line.startswith('#include stiffeners elements'):
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        # U30
                        inpfile.write(f'{21+istiff*14}, {72+istiff*28}, {71+istiff*28}, {70+istiff*28}, {69+istiff*28}\n')
                        inpfile.write(f'{22+istiff*14}, {76+istiff*28}, {75+istiff*28}, {74+istiff*28}, {73+istiff*28}\n')
                elif line.startswith('#include elsets'):
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        inpfile.write(f'*ELSET,ELSET=EltPatch{17+istiff*14}\n')
                        inpfile.write(f'{21+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{18+istiff*14}\n')
                        inpfile.write(f'{22+istiff*14}\n')

                        inpfile.write(f'*ELSET,ELSET=EltPatch{19+istiff*14}\n')
                        inpfile.write(f'{23+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{20+istiff*14}\n')
                        inpfile.write(f'{24+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{21+istiff*14}\n')
                        inpfile.write(f'{25+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{22+istiff*14}\n')
                        inpfile.write(f'{26+istiff*14}\n')

                        inpfile.write(f'*ELSET,ELSET=EltPatch{23+istiff*14}\n')
                        inpfile.write(f'{27+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{24+istiff*14}\n')
                        inpfile.write(f'{28+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{25+istiff*14}\n')
                        inpfile.write(f'{29+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{26+istiff*14}\n')
                        inpfile.write(f'{30+istiff*14}\n')

                        inpfile.write(f'*ELSET,ELSET=EltPatch{27+istiff*14}\n')
                        inpfile.write(f'{31+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{28+istiff*14}\n')
                        inpfile.write(f'{32+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{29+istiff*14}\n')
                        inpfile.write(f'{33+istiff*14}\n')
                        inpfile.write(f'*ELSET,ELSET=EltPatch{30+istiff*14}\n')
                        inpfile.write(f'{34+istiff*14}\n')
                elif line.startswith('#include nsets'):
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        inpfile.write(f'*NSET,NSET=CPcurveP{19+istiff*14}\n')
                        inpfile.write(f'{78+istiff*28}, {77+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPcurveP{20+istiff*14}\n')
                        inpfile.write(f'{80+istiff*28}, {79+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPcurveP{21+istiff*14}\n')
                        inpfile.write(f'{82+istiff*28}, {81+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPcurveP{22+istiff*14}\n')
                        inpfile.write(f'{84+istiff*28}, {83+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPcurveP{23+istiff*14}\n')
                        inpfile.write(f'{86+istiff*28}, {85+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPcurveP{24+istiff*14}\n')
                        inpfile.write(f'{88+istiff*28}, {87+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPcurveP{25+istiff*14}\n')
                        inpfile.write(f'{90+istiff*28}, {89+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPcurveP{26+istiff*14}\n')
                        inpfile.write(f'{92+istiff*28}, {91+istiff*28}\n')

                        inpfile.write(f'*NSET,NSET=CPlgrgeP{27+istiff*14}\n')
                        inpfile.write(f'{93+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPlgrgeP{28+istiff*14}\n')
                        inpfile.write(f'{94+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPlgrgeP{29+istiff*14}\n')
                        inpfile.write(f'{95+istiff*28}\n')
                        inpfile.write(f'*NSET,NSET=CPlgrgeP{30+istiff*14}\n')
                        inpfile.write(f'{96+istiff*28}\n')

                elif line.startswith('#include uel properties'):
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{17+istiff*14}, MATERIAL=Aluminium\n')
                        inpfile.write(f'{17+istiff*14}, 15, {STIFF_THICKNESS}\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{18+istiff*14}, MATERIAL=Aluminium\n')
                        inpfile.write(f'{18+istiff*14}, 16, {STIFF_THICKNESS}\n')

                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{19+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{19+istiff*14}, 1, {27+istiff*14}, 0\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{20+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{20+istiff*14}, {17+istiff*14}, {27+istiff*14}, 1\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{21+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{21+istiff*14}, 1, {28+istiff*14}, 0\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{22+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{22+istiff*14}, {17+istiff*14}, {28+istiff*14}, 1\n')

                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{23+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{23+istiff*14}, 2, {29+istiff*14}, 0\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{24+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{24+istiff*14}, {18+istiff*14}, {29+istiff*14}, 1\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{25+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{23+istiff*14}, 2, {30+istiff*14}, 0\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{26+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{24+istiff*14}, {18+istiff*14}, {30+istiff*14}, 1\n')

                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{27+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{27+istiff*14}, 0\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{28+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{28+istiff*14}, 1\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{29+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{29+istiff*14}, 0\n')
                        inpfile.write(f'*UEL PROPERTY, ELSET=EltPatch{30+istiff*14}, MATERIAL=Void\n')
                        inpfile.write(f'{30+istiff*14}, 1\n')

                elif line.startswith('#include boundary'):
                    # print(line)
                    for istiff in range(NB_STIFFENERS_BY_HALF):
                        inpfile.write(f'I1.CPcurveP{19+istiff*14}, 1, 3, 0.\n')
                        inpfile.write(f'I1.CPcurveP{20+istiff*14}, 1, 3, 0.\n')
                        inpfile.write(f'I1.CPcurveP{21+istiff*14}, 1, 3, 0.\n')
                        inpfile.write(f'I1.CPcurveP{22+istiff*14}, 1, 3, 0.\n')
                        inpfile.write(f'I1.CPcurveP{23+istiff*14}, 1, 3, 0.\n')
                        inpfile.write(f'I1.CPcurveP{24+istiff*14}, 1, 3, 0.\n')
                        inpfile.write(f'I1.CPcurveP{25+istiff*14}, 1, 3, 0.\n')
                        inpfile.write(f'I1.CPcurveP{26+istiff*14}, 1, 3, 0.\n')
                        # inpfile.write(f'I1.CPlgrgeP{27+istiff*14}, 2, 3, 0.\n')
                        inpfile.write(f'I1.CPlgrgeP{28+istiff*14}, 2, 3, 0.\n')
                        # inpfile.write(f'I1.CPlgrgeP{29+istiff*14}, 2, 3, 0.\n')
                        inpfile.write(f'I1.CPlgrgeP{30+istiff*14}, 2, 3, 0.\n')

                else:
                    inpfile.write(line)


# Write NB file
def write_nb_file(filename):
    # write NB file
    with open(filename, 'w') as nbfile:
        nbfile.write('*Dimension\n')
        nbfile.write('2,2,1,1,1,1,1,1,1,1,1,1,1,1,3,3')
        for istiff in range(NB_STIFFENERS_BY_HALF):
            nbfile.write(',2,2,1,1,1,1,1,1,1,1,1,1,1,1')
        nbfile.write('\n')

        nbfile.write('*Number of CP by element\n')
        nbfile.write('6,6,2,2,2,2,2,2,2,2,1,1,1,1,12,12')
        for istiff in range(NB_STIFFENERS_BY_HALF):
            nbfile.write(',4,4,2,2,2,2,2,2,2,2,1,1,1,1')
        nbfile.write('\n')

        nbfile.write('*Number of patch\n')
        nbfile.write(f'{16+NB_STIFFENERS_BY_HALF*14}\n')

        nbfile.write('*Total number of element\n')
        nbfile.write(f'{20+NB_STIFFENERS_BY_HALF*14}\n')

        nbfile.write('*Number of element by patch\n')
        nbfile.write('2,2,1,1,1,1,1,1,1,1,1,1,1,1,2,2')
        for istiff in range(NB_STIFFENERS_BY_HALF):
            nbfile.write(',1,1,1,1,1,1,1,1,1,1,1,1,1,1')
        nbfile.write('\n')

        nbfile.write('*Patch(1)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(2)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(3)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(4)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(5)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(6)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(7)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(8)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(9)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(10)\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(11)\n2\n0.0, 1.0\n')
        nbfile.write('*Patch(12)\n2\n0.0, 1.0\n')
        nbfile.write('*Patch(13)\n2\n0.0, 1.0\n')
        nbfile.write('*Patch(14)\n2\n0.0, 1.0\n')
        nbfile.write('*Patch(15)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
        nbfile.write('*Patch(16)\n7\n0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
        for istiff in range(NB_STIFFENERS_BY_HALF):
            nbfile.write(f'*Patch({17+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({18+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n4\n0.0, 0.0, 1.0, 1.0\n')

            nbfile.write(f'*Patch({19+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({20+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({21+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({22+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({23+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({24+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({25+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')
            nbfile.write(f'*Patch({26+istiff*14})\n4\n0.0, 0.0, 1.0, 1.0\n')

            nbfile.write(f'*Patch({27+istiff*14})\n2\n0.0, 1.0\n')
            nbfile.write(f'*Patch({28+istiff*14})\n2\n0.0, 1.0\n')
            nbfile.write(f'*Patch({29+istiff*14})\n2\n0.0, 1.0\n')
            nbfile.write(f'*Patch({30+istiff*14})\n2\n0.0, 1.0\n')

        nbfile.write('*Jpqr\n')
        nbfile.write('2,1\n2,1\n1\n1\n1\n1\n1\n1\n1\n1\n0\n0\n0\n0\n2,1,1\n2,1,1\n')
        for istiff in range(NB_STIFFENERS_BY_HALF):
            nbfile.write('1,1\n1,1\n1\n1\n1\n1\n1\n1\n1\n1\n0\n0\n0\n0\n')
        nbfile.write('*Nijk\n')
        nbfile.write('1, 3, 2\n2, 4, 2\n3, 3, 2\n4, 4, 2\n5, 2\n6, 2\n7, 2\n8, 2\n9, 2\n10, 2\n11, 2\n12, 2\n13, 1\n14, 1\n15, 1\n16, 1\n17, 3, 2, 2\n18, 4, 2, 2\n19, 3, 2, 2\n20, 4, 2, 2\n')
        for istiff in range(NB_STIFFENERS_BY_HALF):
            nbfile.write(f'{21+istiff*14}, 2, 2\n')
            nbfile.write(f'{22+istiff*14}, 2, 2\n')

            nbfile.write(f'{23+istiff*14}, 2\n')
            nbfile.write(f'{24+istiff*14}, 2\n')
            nbfile.write(f'{25+istiff*14}, 2\n')
            nbfile.write(f'{26+istiff*14}, 2\n')
            nbfile.write(f'{27+istiff*14}, 2\n')
            nbfile.write(f'{28+istiff*14}, 2\n')
            nbfile.write(f'{29+istiff*14}, 2\n')
            nbfile.write(f'{30+istiff*14}, 2\n')

            nbfile.write(f'{31+istiff*14}, 1\n')
            nbfile.write(f'{32+istiff*14}, 1\n')
            nbfile.write(f'{33+istiff*14}, 1\n')
            nbfile.write(f'{34+istiff*14}, 1\n')

        nbfile.write('*Weight\n')
        nbfile.write('1, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
        nbfile.write('2, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
        nbfile.write('3, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
        nbfile.write('4, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
        nbfile.write('5, 1.0, 1.0\n')
        nbfile.write('6, 1.0, 1.0\n')
        nbfile.write('7, 1.0, 1.0\n')
        nbfile.write('8, 1.0, 1.0\n')
        nbfile.write('9, 1.0, 1.0\n')
        nbfile.write('10, 1.0, 1.0\n')
        nbfile.write('11, 1.0, 1.0\n')
        nbfile.write('12, 1.0, 1.0\n')
        nbfile.write('13, 1.0\n')
        nbfile.write('14, 1.0\n')
        nbfile.write('15, 1.0\n')
        nbfile.write('16, 1.0\n')
        nbfile.write('17, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
        nbfile.write('18, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
        nbfile.write('19, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0\n')
        nbfile.write('20, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5, 1.0, 0.5, 0.5\n')
        for istiff in range(NB_STIFFENERS_BY_HALF):
            nbfile.write(f'{21+istiff*14}, 1.0, 1.0, 1.0, 1.0\n')
            nbfile.write(f'{22+istiff*14}, 1.0, 1.0, 1.0, 1.0\n')
            nbfile.write(f'{23+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{24+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{25+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{26+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{27+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{28+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{29+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{30+istiff*14}, 1.0, 1.0\n')
            nbfile.write(f'{31+istiff*14}, 1.0\n')
            nbfile.write(f'{32+istiff*14}, 1.0\n')
            nbfile.write(f'{33+istiff*14}, 1.0\n')
            nbfile.write(f'{34+istiff*14}, 1.0\n')


