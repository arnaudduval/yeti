"""
TODO
"""

import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.optimize import minimize

from preprocessing.igaparametrization import IGAparametrization
from stiffmtrx_elemstorage import sys_linmat_lindef_static as build_stiffmatrix
import reconstructionSOL as rsol
import postprocessing.postproc as pp
from preprocessing.igaparametrization import IGAmanip as manip
from coupling.cplgmatrix import cplg_matrix
from preprocessing.igaparametrization import OPTmodelling



NB_STIFFENERS_BY_HALF = 8
ALLOWED_VOLUME_VAR = 0.05
STIFF_THICKNESS = 1.7

with open('stiff_geom_2xhalf.inp.template', 'r') as template:
    with open('stiff_geom_2xhalf.inp', 'w') as inpfile:
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
                print(line)
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



with open('stiff_geom_2xhalf.NB', 'w') as nbfile:
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




iga_model = IGAparametrization(filename='stiff_geom_2xhalf')
iga_model._NBPINT[
    np.where(iga_model._ELT_TYPE == 'U00')] = 6

# Refine model
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

print('Raffinement initial + découpage patch 1')
nb_deg[:, 0] = np.array([0, 1, 0])
nb_deg[:, 1] = np.array([0, 1, 0])
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.43, 0.57]),
                    '2': np.array([0.1, 0.3]),
                    '3': np.array([])}
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

print('découpage patch 2')
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
additional_knots = {'patches': np.array([1]),
                    '1': np.array([0.43, 0.57]),
                    '2': np.array([0.1, 0.3]),
                    '3': np.array([])}
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

# Add loading on a single element
iga_model._indDLoad = np.array([ np.array([6,7]), np.array([6,7]), np.array([18,19]), np.array([18,19]) ])
iga_model._JDLType = np.array([62,63,62,63])
iga_model._ADLMAG = np.array([1.190, 7.143, -1.190, 7.143])
iga_model._load_target_nbelem = np.array([2, 2, 2, 2])
iga_model._nb_load = 4
iga_model._indDLoad_flat = np.array([], dtype=np.intp)
for load in iga_model._indDLoad:
    iga_model._indDLoad_flat = np.hstack((iga_model._indDLoad_flat, load))

additional_knots = {"patches": np.array([]),
                    "1": np.array([]), "2": np.array([]), "3": np.array([])}

# Refinement to get elements with approx the same size
print("Uniformisation maillage patch 1")
additional_knots = {'patches': np.array([0]),
                    '1': np.array([0.1075, 0.215, 0.3325, 0.6775, 0.785, 0.8925]),
                    '2': np.array([0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
                    '3': np.array([])}
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)
print("Uniformisation maillage patch 2")
additional_knots = {'patches': np.array([1]),
                    '1': np.array([0.1075, 0.215, 0.3325, 0.6775, 0.785, 0.8925]),
                    '2': np.array([0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]),
                    '3': np.array([])}
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
iga_model.refine(nb_ref, nb_deg, additional_knots=additional_knots)

def movestiffeners(coords0, igapara, var):
    igapara._COORDS[:, :] = coords0[:, :]

    assert var.shape[0] == NB_STIFFENERS_BY_HALF*4

    for istiff in range(NB_STIFFENERS_BY_HALF):
        igapara._COORDS[0, igapara._indCPbyPatch[16+istiff*14][:]-1] = var[[0 + istiff*4, 1 + istiff*4, 0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[17+istiff*14][:]-1] = var[[0 + istiff*4, 1 + istiff*4, 0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[18+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[20+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[22+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]
        igapara._COORDS[0, igapara._indCPbyPatch[24+istiff*14][[0,1]]-1] = var[[0 + istiff*4, 1 + istiff*4]]

        igapara._COORDS[1, igapara._indCPbyPatch[16+istiff*14][:]-1] = var[[2 + istiff*4, 3 + istiff*4, 2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[17+istiff*14][:]-1] = var[[2 + istiff*4, 3 + istiff*4, 2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[18+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[20+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[22+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]
        igapara._COORDS[1, igapara._indCPbyPatch[24+istiff*14][[0,1]]-1] = var[[2 + istiff*4, 3 + istiff*4]]

    return None

### Cas 24 raidisseurs
# x0 = np.array([ 0.05, 0.35, 0.05, 0.05,
#                 0.35, 0.65, 0.05, 0.05,
#                 0.65, 0.95, 0.05, 0.05,
#                 0.05, 0.35, 0.35, 0.35,
#                 0.35, 0.65, 0.35, 0.35,
#                 0.65, 0.95, 0.35, 0.35,
#                 0.05, 0.35, 0.65, 0.65,
#                 0.35, 0.65, 0.65, 0.65,
#                 0.65, 0.95, 0.65, 0.65,
#                 0.05, 0.35, 0.95, 0.95,
#                 0.35, 0.65, 0.95, 0.95,
#                 0.65, 0.95, 0.95, 0.95,
#                 0.05, 0.05, 0.05, 0.35,
#                 0.05, 0.05, 0.35, 0.65,
#                 0.05, 0.05, 0.65, 0.95,
#                 0.35, 0.35, 0.05, 0.35,
#                 0.35, 0.35, 0.35, 0.65,
#                 0.35, 0.35, 0.65, 0.95,
#                 0.65, 0.65, 0.05, 0.35,
#                 0.65, 0.65, 0.35, 0.65,
#                 0.65, 0.65, 0.65, 0.95,
#                 0.95, 0.95, 0.05, 0.35,
#                 0.95, 0.95, 0.35, 0.65,
#                 0.95, 0.95, 0.65, 0.95])

x0 = np.array([0.0,     0.5,    0.0,   0.5,
               0.0,   0.5,  0.5,   0.0,
               0.5,    1.0,   0.0,   0.5,
               0.5,   1.0,  0.5,   0.0,
               0.0,     0.5,    1.0,   0.5,
               0.0,   0.5,  0.5,   1.0,
               0.5,    1.0,   1.0,   0.5,
               0.5,   1.0,  0.5,   1.0])

movestiffeners(iga_model._COORDS.copy(), iga_model, x0)

# OTHER REFINEMENTS For optim problem
nb_deg = np.zeros((3, iga_model.nb_patch), dtype=np.intp)
nb_ref = np.zeros((3, iga_model.nb_patch), dtype=np.intp)

nb_ref[:, 0] = [1, 1, 0]
nb_ref[:, 1] = [1, 1, 0]

# Hull
nb_deg[:, 14] = [0, 1, 1]
nb_deg[:, 15] = [0, 1, 1]

# Curves 2 halfs
nb_ref[0, 2:10] = 3
# Lagrange patches 2 halfs
nb_deg[0, 10] = 1
nb_deg[0, 11] = 0
nb_deg[0, 12] = 1
nb_deg[0, 13] = 0

nb_ref[0, 10] = 3
nb_ref[0, 11] = 3
nb_ref[0, 12] = 3
nb_ref[0, 13] = 3

# Stiffeners list
for istiff in range(NB_STIFFENERS_BY_HALF):
    nb_deg[:, 16+istiff*14] = [1, 1, 0]
    nb_deg[:, 17+istiff*14] = [1, 1, 0]
    nb_ref[:, 16+istiff*14] = [4, 2, 0]
    nb_ref[:, 17+istiff*14] = [4, 2, 0]

    nb_ref[0, 18+istiff*14] = 4
    nb_ref[0, 19+istiff*14] = 4
    nb_ref[0, 20+istiff*14] = 4
    nb_ref[0, 21+istiff*14] = 4
    nb_ref[0, 22+istiff*14] = 4
    nb_ref[0, 23+istiff*14] = 4
    nb_ref[0, 24+istiff*14] = 4
    nb_ref[0, 25+istiff*14] = 4

    nb_deg[0, 26+istiff*14] = 1
    nb_deg[0, 27+istiff*14] = 0
    nb_deg[0, 28+istiff*14] = 1
    nb_deg[0, 29+istiff*14] = 0

    nb_ref[0, 26+istiff*14] = 3
    nb_ref[0, 27+istiff*14] = 3
    nb_ref[0, 28+istiff*14] = 3
    nb_ref[0, 29+istiff*14] = 3

# END OF OTHER REFINEMENTS
# NOTE : on peut faire un denier raffinement ici et enchainer sur un calcul linéaire

opt_pb = OPTmodelling(iga_model, NB_STIFFENERS_BY_HALF*4, movestiffeners,
                     nb_degreeElevationByDirection=nb_deg,
                     nb_refinementByDirection=nb_ref)

# List of patches for volume computation
listpatch = np.zeros(opt_pb._coarseParametrization._nb_patch, dtype=np.intp)
listpatch[0] = 1
listpatch[1] = 1
for istiff in range(NB_STIFFENERS_BY_HALF):
    listpatch[16+istiff*14] = 1
    listpatch[17+istiff*14] = 1

V0 = opt_pb.compute_volume(x0, listpatch)
C0 = opt_pb.compute_compliance_discrete(x0)



def var_vol(xk):
    return ALLOWED_VOLUME_VAR - abs(opt_pb.compute_volume(xk, listpatch) - V0) / V0

def grad_var_vol(xk):
    vk = opt_pb.compute_volume(xk, listpatch)
    return np.sign(vk - V0)*opt_pb.compute_gradVolume_AN(xk, listpatch)/V0

def compliance(xk):
    return opt_pb.compute_compliance_discrete(xk)/C0

def grad_compliance(xk):
    return (opt_pb.compute_gradCompliance_AN(xk)+ opt_pb.compute_gradCompliance_cplgOnly_AN(xk))/C0

iopt = 0

def save_xk(xk):
    global iopt
    print((f'\nIteration {iopt:03}'))
    pp.generatevtu(*opt_pb.coarse_parametrization.get_inputs4postprocVTU(
        f'opt-coarse{iopt:02}', np.zeros_like(opt_pb.coarse_parametrization.coords),
        nb_ref=4*np.array([1, 1, 0]),Flag=np.array([False]*3)))
    opt_pb.coarse_parametrization.generate_vtk4controlMeshVisu(f'opt-coarse{iopt:02}',0)

    sol, _ = rsol.reconstruction(
             **opt_pb.fine_parametrization.get_inputs4solution(opt_pb.save_sol_fine))
    pp.generatevtu(*opt_pb.fine_parametrization.get_inputs4postprocVTU(
                    f'opt-fine{iopt:02}',  sol.transpose(),
                    nb_ref=4*np.array([1, 1, 0]),
                    Flag=np.array([True, False, False])))
    iopt += 1

save_xk(x0)

# Bounds for design variables
bounds = ((0., 1.),) * NB_STIFFENERS_BY_HALF*4
constraints = ({'type': 'ineq', 'fun': var_vol, 'jac': grad_var_vol})
res = minimize(compliance, x0, method='SLSQP',
               jac=grad_compliance, bounds=bounds,
               constraints=constraints, callback=save_xk,
               tol=1.E-2)


print('Optimization succes : ', res['success'])
print(res['message'])
print('Result : ', res['x'])
print('Objective function value : ', res['fun'])
print('# of evaluations of objective function : ', res['nfev'])
print('# of evaluations of jacobian : ', res['njev'])
v_f = opt_pb.compute_volume(res['x'], listpatch)
c_f = opt_pb.compute_compliance_discrete(res['x'])
print('Volume: ', V0, '->', v_f, (100.*(v_f - V0)/V0), ' %')
print('Compliance: ', C0, '->',  c_f, (100.*(c_f - C0)/C0), ' %')


exit()

# Matrix assembly
ndof = iga_model._nb_dof_free
idof = iga_model._ind_dof_free[:ndof]-1



data, row, col, Fb = build_stiffmatrix(
                            *iga_model.get_inputs4system_elemStorage())
Kside = sp.coo_matrix((data, (row, col)),
                      shape=(iga_model.nb_dof_tot, iga_model.nb_dof_tot),
                      dtype='float64').tocsc()
Ktot = Kside + Kside.transpose()
K2solve = Ktot[idof, :][:, idof]
Fbsolve = Fb[idof]
del Kside, data, row, col

Cdata,Crow,Ccol = cplg_matrix( *iga_model.get_inputs4cplgmatrix() )
Cside = sp.coo_matrix((Cdata,(Crow,Ccol)), shape=(iga_model._nb_dof_tot, iga_model._nb_dof_tot),
                      dtype='float64').tocsc()
Ctot  = Cside + Cside.transpose()
del Cdata,Crow,Ccol,Cside


C2solve = Ctot[idof,:][:,idof] * K2solve.max()


# plt.spy(K2solve + C2solve)
# plt.show()

# Resolution
x = sp.linalg.spsolve(K2solve + C2solve, Fbsolve)

# Solution reconstruction
SOL, u = rsol.reconstruction(**iga_model.get_inputs4solution(x))
pp.generatevtu(*iga_model.get_inputs4postprocVTU(
    'Ariane',
    SOL.transpose(),
    nb_ref=np.array([4, 4, 0]),
    Flag=np.array([True, False, False])))
