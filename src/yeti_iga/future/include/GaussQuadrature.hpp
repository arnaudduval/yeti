#pragma once
#include <vector>
#include <cassert>
#include <stdexcept>

// Supports 1-20 points (degree up to 19 via IGABasis1D's degree+1 convention).
// Extend further the same way if needed (e.g. numpy.polynomial.legendre.leggauss).
inline void gauss_legendre_table(int ngauss, std::vector<double>& points, std::vector<double>& weights) {
    points.resize(ngauss);
    weights.resize(ngauss);

    if (ngauss == 1) {
        points[0] = 0.0;
        weights[0] = 2.0;
    } else if (ngauss == 2) {
        points[0] = -0.5773502691896257;
        points[1] = -points[0];
        weights[0] = weights[1] = 1.0;
    } else if (ngauss == 3) {
        points[0] = -0.7745966692414834;
        points[1] = 0.0;
        points[2] = -points[0];
        weights[0] = weights[2] = 0.5555555555555556;
        weights[1] = 0.8888888888888888;
    } else if (ngauss == 4) {
        points[0] = -0.8611363115940526;
        points[1] = -0.3399810435848563;
        points[2] = -points[1];
        points[3] = -points[0];
        weights[0] = weights[3] = 0.3478548451374538;
        weights[1] = weights[2] = 0.6521451548625461;
    } else if (ngauss == 5) {
        points[0] = -0.9061798459386640;
        points[1] = -0.5384693101056831;
        points[2] = 0.0;
        points[3] = -points[1];
        points[4] = -points[0];
        weights[0] = weights[4] = 0.2369268850561891;
        weights[1] = weights[3] = 0.4786286704993665;
        weights[2] = 0.5688888888888889;
    } else if (ngauss == 6) {
        points[0] = -0.9324695142031521;
        points[1] = -0.6612093864662645;
        points[2] = -0.2386191860831969;
        points[3] = -points[2];
        points[4] = -points[1];
        points[5] = -points[0];
        weights[0] = weights[5] = 0.1713244923791704;
        weights[1] = weights[4] = 0.3607615730481386;
        weights[2] = weights[3] = 0.4679139345726910;
    }
    else if (ngauss == 7) {
        points[0] = -0.9491079123427586;
        points[1] = -0.7415311855993945;
        points[2] = -0.4058451513773972;
        points[3] = 0.0;
        points[6] = -points[0];
        points[5] = -points[1];
        points[4] = -points[2];
        weights[0] = weights[6] = 0.1294849661688697;
        weights[1] = weights[5] = 0.2797053914892769;
        weights[2] = weights[4] = 0.3818300505051187;
        weights[3] = 0.4179591836734693;
    }
    else if (ngauss == 8) {
        points[0] = -0.9602898564975362;
        points[1] = -0.7966664774136267;
        points[2] = -0.525532409916329;
        points[3] = -0.1834346424956498;
        points[7] = -points[0];
        points[6] = -points[1];
        points[5] = -points[2];
        points[4] = -points[3];
        weights[0] = weights[7] = 0.1012285362903771;
        weights[1] = weights[6] = 0.2223810344533744;
        weights[2] = weights[5] = 0.3137066458778869;
        weights[3] = weights[4] = 0.3626837833783617;
    }
    else if (ngauss == 9) {
        points[0] = -0.9681602395076261;
        points[1] = -0.8360311073266358;
        points[2] = -0.6133714327005904;
        points[3] = -0.3242534234038089;
        points[4] = 0.0;
        points[8] = -points[0];
        points[7] = -points[1];
        points[6] = -points[2];
        points[5] = -points[3];
        weights[0] = weights[8] = 0.08127438836157416;
        weights[1] = weights[7] = 0.1806481606948574;
        weights[2] = weights[6] = 0.2606106964029356;
        weights[3] = weights[5] = 0.3123470770400029;
        weights[4] = 0.3302393550012598;
    }
    else if (ngauss == 10) {
        points[0] = -0.9739065285171717;
        points[1] = -0.8650633666889845;
        points[2] = -0.6794095682990244;
        points[3] = -0.4333953941292472;
        points[4] = -0.1488743389816312;
        points[9] = -points[0];
        points[8] = -points[1];
        points[7] = -points[2];
        points[6] = -points[3];
        points[5] = -points[4];
        weights[0] = weights[9] = 0.06667134430868814;
        weights[1] = weights[8] = 0.1494513491505804;
        weights[2] = weights[7] = 0.219086362515982;
        weights[3] = weights[6] = 0.2692667193099965;
        weights[4] = weights[5] = 0.2955242247147528;
    }
    else if (ngauss == 11) {
        points[0] = -0.978228658146057;
        points[1] = -0.8870625997680953;
        points[2] = -0.7301520055740494;
        points[3] = -0.5190961292068118;
        points[4] = -0.269543155952345;
        points[5] = 0.0;
        points[10] = -points[0];
        points[9] = -points[1];
        points[8] = -points[2];
        points[7] = -points[3];
        points[6] = -points[4];
        weights[0] = weights[10] = 0.05566856711617393;
        weights[1] = weights[9] = 0.1255803694649043;
        weights[2] = weights[8] = 0.1862902109277343;
        weights[3] = weights[7] = 0.2331937645919906;
        weights[4] = weights[6] = 0.2628045445102465;
        weights[5] = 0.2729250867779005;
    }
    else if (ngauss == 12) {
        points[0] = -0.9815606342467192;
        points[1] = -0.9041172563704748;
        points[2] = -0.7699026741943047;
        points[3] = -0.5873179542866175;
        points[4] = -0.3678314989981802;
        points[5] = -0.1252334085114689;
        points[11] = -points[0];
        points[10] = -points[1];
        points[9] = -points[2];
        points[8] = -points[3];
        points[7] = -points[4];
        points[6] = -points[5];
        weights[0] = weights[11] = 0.04717533638651141;
        weights[1] = weights[10] = 0.1069393259953191;
        weights[2] = weights[9] = 0.1600783285433464;
        weights[3] = weights[8] = 0.2031674267230657;
        weights[4] = weights[7] = 0.2334925365383546;
        weights[5] = weights[6] = 0.2491470458134027;
    }
    else if (ngauss == 13) {
        points[0] = -0.9841830547185881;
        points[1] = -0.9175983992229779;
        points[2] = -0.8015780907333099;
        points[3] = -0.6423493394403402;
        points[4] = -0.4484927510364469;
        points[5] = -0.2304583159551348;
        points[6] = 0.0;
        points[12] = -points[0];
        points[11] = -points[1];
        points[10] = -points[2];
        points[9] = -points[3];
        points[8] = -points[4];
        points[7] = -points[5];
        weights[0] = weights[12] = 0.04048400476531557;
        weights[1] = weights[11] = 0.0921214998377288;
        weights[2] = weights[10] = 0.1388735102197873;
        weights[3] = weights[9] = 0.1781459807619455;
        weights[4] = weights[8] = 0.2078160475368885;
        weights[5] = weights[7] = 0.2262831802628974;
        weights[6] = 0.2325515532308739;
    }
    else if (ngauss == 14) {
        points[0] = -0.9862838086968124;
        points[1] = -0.9284348836635735;
        points[2] = -0.827201315069765;
        points[3] = -0.6872929048116855;
        points[4] = -0.5152486363581541;
        points[5] = -0.3191123689278897;
        points[6] = -0.1080549487073437;
        points[13] = -points[0];
        points[12] = -points[1];
        points[11] = -points[2];
        points[10] = -points[3];
        points[9] = -points[4];
        points[8] = -points[5];
        points[7] = -points[6];
        weights[0] = weights[13] = 0.03511946033175175;
        weights[1] = weights[12] = 0.08015808715975999;
        weights[2] = weights[11] = 0.1215185706879033;
        weights[3] = weights[10] = 0.1572031671581936;
        weights[4] = weights[9] = 0.1855383974779379;
        weights[5] = weights[8] = 0.2051984637212957;
        weights[6] = weights[7] = 0.2152638534631579;
    }
    else if (ngauss == 15) {
        points[0] = -0.9879925180204854;
        points[1] = -0.9372733924007058;
        points[2] = -0.8482065834104272;
        points[3] = -0.7244177313601701;
        points[4] = -0.5709721726085388;
        points[5] = -0.3941513470775634;
        points[6] = -0.2011940939974345;
        points[7] = 0.0;
        points[14] = -points[0];
        points[13] = -points[1];
        points[12] = -points[2];
        points[11] = -points[3];
        points[10] = -points[4];
        points[9] = -points[5];
        points[8] = -points[6];
        weights[0] = weights[14] = 0.0307532419961172;
        weights[1] = weights[13] = 0.0703660474881084;
        weights[2] = weights[12] = 0.1071592204671714;
        weights[3] = weights[11] = 0.1395706779261544;
        weights[4] = weights[10] = 0.166269205816994;
        weights[5] = weights[9] = 0.1861610000155622;
        weights[6] = weights[8] = 0.1984314853271116;
        weights[7] = 0.2025782419255613;
    }
    else if (ngauss == 16) {
        points[0] = -0.9894009349916499;
        points[1] = -0.9445750230732326;
        points[2] = -0.8656312023878318;
        points[3] = -0.755404408355003;
        points[4] = -0.6178762444026438;
        points[5] = -0.4580167776572274;
        points[6] = -0.2816035507792589;
        points[7] = -0.09501250983763744;
        points[15] = -points[0];
        points[14] = -points[1];
        points[13] = -points[2];
        points[12] = -points[3];
        points[11] = -points[4];
        points[10] = -points[5];
        points[9] = -points[6];
        points[8] = -points[7];
        weights[0] = weights[15] = 0.02715245941175418;
        weights[1] = weights[14] = 0.06225352393864746;
        weights[2] = weights[13] = 0.09515851168249261;
        weights[3] = weights[12] = 0.1246289712555341;
        weights[4] = weights[11] = 0.1495959888165767;
        weights[5] = weights[10] = 0.1691565193950026;
        weights[6] = weights[9] = 0.1826034150449236;
        weights[7] = weights[8] = 0.1894506104550686;
    }
    else if (ngauss == 17) {
        points[0] = -0.9905754753144174;
        points[1] = -0.9506755217687677;
        points[2] = -0.8802391537269859;
        points[3] = -0.7815140038968014;
        points[4] = -0.6576711592166907;
        points[5] = -0.5126905370864769;
        points[6] = -0.3512317634538763;
        points[7] = -0.1784841814958479;
        points[8] = 0.0;
        points[16] = -points[0];
        points[15] = -points[1];
        points[14] = -points[2];
        points[13] = -points[3];
        points[12] = -points[4];
        points[11] = -points[5];
        points[10] = -points[6];
        points[9] = -points[7];
        weights[0] = weights[16] = 0.02414830286854758;
        weights[1] = weights[15] = 0.05545952937398796;
        weights[2] = weights[14] = 0.08503614831717912;
        weights[3] = weights[13] = 0.111883847193404;
        weights[4] = weights[12] = 0.1351363684685255;
        weights[5] = weights[11] = 0.1540457610768103;
        weights[6] = weights[10] = 0.1680041021564499;
        weights[7] = weights[9] = 0.1765627053669925;
        weights[8] = 0.1794464703562065;
    }
    else if (ngauss == 18) {
        points[0] = -0.991565168420931;
        points[1] = -0.9558239495713978;
        points[2] = -0.8926024664975557;
        points[3] = -0.8037049589725231;
        points[4] = -0.6916870430603532;
        points[5] = -0.5597708310739475;
        points[6] = -0.4117511614628426;
        points[7] = -0.2518862256915055;
        points[8] = -0.08477501304173529;
        points[17] = -points[0];
        points[16] = -points[1];
        points[15] = -points[2];
        points[14] = -points[3];
        points[13] = -points[4];
        points[12] = -points[5];
        points[11] = -points[6];
        points[10] = -points[7];
        points[9] = -points[8];
        weights[0] = weights[17] = 0.0216160135264815;
        weights[1] = weights[16] = 0.04971454889496939;
        weights[2] = weights[15] = 0.07642573025488957;
        weights[3] = weights[14] = 0.1009420441062873;
        weights[4] = weights[13] = 0.1225552067114787;
        weights[5] = weights[12] = 0.1406429146706509;
        weights[6] = weights[11] = 0.1546846751262655;
        weights[7] = weights[10] = 0.164276483745833;
        weights[8] = weights[9] = 0.1691423829631439;
    }
    else if (ngauss == 19) {
        points[0] = -0.9924068438435844;
        points[1] = -0.96020815213483;
        points[2] = -0.9031559036148179;
        points[3] = -0.8227146565371428;
        points[4] = -0.7209661773352294;
        points[5] = -0.600545304661681;
        points[6] = -0.4645707413759609;
        points[7] = -0.3165640999636298;
        points[8] = -0.1603586456402254;
        points[9] = 0.0;
        points[18] = -points[0];
        points[17] = -points[1];
        points[16] = -points[2];
        points[15] = -points[3];
        points[14] = -points[4];
        points[13] = -points[5];
        points[12] = -points[6];
        points[11] = -points[7];
        points[10] = -points[8];
        weights[0] = weights[18] = 0.01946178822972652;
        weights[1] = weights[17] = 0.04481422676569959;
        weights[2] = weights[16] = 0.06904454273764117;
        weights[3] = weights[15] = 0.09149002162244999;
        weights[4] = weights[14] = 0.111566645547334;
        weights[5] = weights[13] = 0.1287539625393362;
        weights[6] = weights[12] = 0.1426067021736065;
        weights[7] = weights[11] = 0.1527660420658596;
        weights[8] = weights[10] = 0.1589688433939543;
        weights[9] = 0.1610544498487836;
    }
    else if (ngauss == 20) {
        points[0] = -0.993128599185095;
        points[1] = -0.9639719272779138;
        points[2] = -0.9122344282513259;
        points[3] = -0.8391169718222188;
        points[4] = -0.7463319064601508;
        points[5] = -0.636053680726515;
        points[6] = -0.5108670019508271;
        points[7] = -0.3737060887154195;
        points[8] = -0.2277858511416451;
        points[9] = -0.07652652113349734;
        points[19] = -points[0];
        points[18] = -points[1];
        points[17] = -points[2];
        points[16] = -points[3];
        points[15] = -points[4];
        points[14] = -points[5];
        points[13] = -points[6];
        points[12] = -points[7];
        points[11] = -points[8];
        points[10] = -points[9];
        weights[0] = weights[19] = 0.01761400713915089;
        weights[1] = weights[18] = 0.04060142980038645;
        weights[2] = weights[17] = 0.06267204833410879;
        weights[3] = weights[16] = 0.08327674157670471;
        weights[4] = weights[15] = 0.1019301198172407;
        weights[5] = weights[14] = 0.1181945319615186;
        weights[6] = weights[13] = 0.1316886384491769;
        weights[7] = weights[12] = 0.1420961093183824;
        weights[8] = weights[11] = 0.1491729864726042;
        weights[9] = weights[10] = 0.1527533871307263;
    }
    else {
        throw std::runtime_error("Number of gauss points not supported");
    }

}
