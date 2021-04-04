/*
 * Copyright (c) 2021
 *      Nakata, Maho
 *      All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */

#include <mpblas.h>
#include <mplapack.h>

void Rlaruv(common &cmn, arr_ref<INTEGER> iseed, INTEGER const &n, REAL *x) {
    FEM_CMN_SVE(Rlaruv);
    const INTEGER lv = 128;
    INTEGER j = 0;
    if (is_called_first_time) {
        {
            static const INTEGER values[] = {494, 322, 2508, 2549};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(1, j);
            }
        }
        {
            static const INTEGER values[] = {2637, 789, 3754, 1145};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(2, j);
            }
        }
        {
            static const INTEGER values[] = {255, 1440, 1766, 2253};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(3, j);
            }
        }
        {
            static const INTEGER values[] = {2008, 752, 3572, 305};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(4, j);
            }
        }
        {
            static const INTEGER values[] = {1253, 2859, 2893, 3301};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(5, j);
            }
        }
        {
            static const INTEGER values[] = {3344, 123, 307, 1065};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(6, j);
            }
        }
        {
            static const INTEGER values[] = {4084, 1848, 1297, 3133};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(7, j);
            }
        }
        {
            static const INTEGER values[] = {1739, 643, 3966, 2913};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(8, j);
            }
        }
        {
            static const INTEGER values[] = {3143, 2405, 758, 3285};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(9, j);
            }
        }
        {
            static const INTEGER values[] = {3468, 2638, 2598, 1241};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(10, j);
            }
        }
        {
            static const INTEGER values[] = {688, 2344, 3406, 1197};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(11, j);
            }
        }
        {
            static const INTEGER values[] = {1657, 46, 2922, 3729};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(12, j);
            }
        }
        {
            static const INTEGER values[] = {1238, 3814, 1038, 2501};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(13, j);
            }
        }
        {
            static const INTEGER values[] = {3166, 913, 2934, 1673};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(14, j);
            }
        }
        {
            static const INTEGER values[] = {1292, 3649, 2091, 541};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(15, j);
            }
        }
        {
            static const INTEGER values[] = {3422, 339, 2451, 2753};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(16, j);
            }
        }
        {
            static const INTEGER values[] = {1270, 3808, 1580, 949};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(17, j);
            }
        }
        {
            static const INTEGER values[] = {2016, 822, 1958, 2361};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(18, j);
            }
        }
        {
            static const INTEGER values[] = {154, 2832, 2055, 1165};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(19, j);
            }
        }
        {
            static const INTEGER values[] = {2862, 3078, 1507, 4081};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(20, j);
            }
        }
        {
            static const INTEGER values[] = {697, 3633, 1078, 2725};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(21, j);
            }
        }
        {
            static const INTEGER values[] = {1706, 2970, 3273, 3305};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(22, j);
            }
        }
        {
            static const INTEGER values[] = {491, 637, 17, 3069};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(23, j);
            }
        }
        {
            static const INTEGER values[] = {931, 2249, 854, 3617};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(24, j);
            }
        }
        {
            static const INTEGER values[] = {1444, 2081, 2916, 3733};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(25, j);
            }
        }
        {
            static const INTEGER values[] = {444, 4019, 3971, 409};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(26, j);
            }
        }
        {
            static const INTEGER values[] = {3577, 1478, 2889, 2157};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(27, j);
            }
        }
        {
            static const INTEGER values[] = {3944, 242, 3831, 1361};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(28, j);
            }
        }
        {
            static const INTEGER values[] = {2184, 481, 2621, 3973};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(29, j);
            }
        }
        {
            static const INTEGER values[] = {1661, 2075, 1541, 1865};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(30, j);
            }
        }
        {
            static const INTEGER values[] = {3482, 4058, 893, 2525};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(31, j);
            }
        }
        {
            static const INTEGER values[] = {657, 622, 736, 1409};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(32, j);
            }
        }
        {
            static const INTEGER values[] = {3023, 3376, 3992, 3445};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(33, j);
            }
        }
        {
            static const INTEGER values[] = {3618, 812, 787, 3577};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(34, j);
            }
        }
        {
            static const INTEGER values[] = {1267, 234, 2125, 77};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(35, j);
            }
        }
        {
            static const INTEGER values[] = {1828, 641, 2364, 3761};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(36, j);
            }
        }
        {
            static const INTEGER values[] = {164, 4005, 2460, 2149};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(37, j);
            }
        }
        {
            static const INTEGER values[] = {3798, 1122, 257, 1449};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(38, j);
            }
        }
        {
            static const INTEGER values[] = {3087, 3135, 1574, 3005};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(39, j);
            }
        }
        {
            static const INTEGER values[] = {2400, 2640, 3912, 225};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(40, j);
            }
        }
        {
            static const INTEGER values[] = {2870, 2302, 1216, 85};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(41, j);
            }
        }
        {
            static const INTEGER values[] = {3876, 40, 3248, 3673};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(42, j);
            }
        }
        {
            static const INTEGER values[] = {1905, 1832, 3401, 3117};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(43, j);
            }
        }
        {
            static const INTEGER values[] = {1593, 2247, 2124, 3089};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(44, j);
            }
        }
        {
            static const INTEGER values[] = {1797, 2034, 2762, 1349};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(45, j);
            }
        }
        {
            static const INTEGER values[] = {1234, 2637, 149, 2057};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(46, j);
            }
        }
        {
            static const INTEGER values[] = {3460, 1287, 2245, 413};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(47, j);
            }
        }
        {
            static const INTEGER values[] = {328, 1691, 166, 65};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(48, j);
            }
        }
        {
            static const INTEGER values[] = {2861, 496, 466, 1845};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(49, j);
            }
        }
        {
            static const INTEGER values[] = {1950, 1597, 4018, 697};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(50, j);
            }
        }
        {
            static const INTEGER values[] = {617, 2394, 1399, 3085};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(51, j);
            }
        }
        {
            static const INTEGER values[] = {2070, 2584, 190, 3441};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(52, j);
            }
        }
        {
            static const INTEGER values[] = {3331, 1843, 2879, 1573};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(53, j);
            }
        }
        {
            static const INTEGER values[] = {769, 336, 153, 3689};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(54, j);
            }
        }
        {
            static const INTEGER values[] = {1558, 1472, 2320, 2941};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(55, j);
            }
        }
        {
            static const INTEGER values[] = {2412, 2407, 18, 929};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(56, j);
            }
        }
        {
            static const INTEGER values[] = {2800, 433, 712, 533};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(57, j);
            }
        }
        {
            static const INTEGER values[] = {189, 2096, 2159, 2841};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(58, j);
            }
        }
        {
            static const INTEGER values[] = {287, 1761, 2318, 4077};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(59, j);
            }
        }
        {
            static const INTEGER values[] = {2045, 2810, 2091, 721};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(60, j);
            }
        }
        {
            static const INTEGER values[] = {1227, 566, 3443, 2821};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(61, j);
            }
        }
        {
            static const INTEGER values[] = {2838, 442, 1510, 2249};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(62, j);
            }
        }
        {
            static const INTEGER values[] = {209, 41, 449, 2397};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(63, j);
            }
        }
        {
            static const INTEGER values[] = {2770, 1238, 1956, 2817};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(64, j);
            }
        }
        {
            static const INTEGER values[] = {3654, 1086, 2201, 245};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(65, j);
            }
        }
        {
            static const INTEGER values[] = {3993, 603, 3137, 1913};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(66, j);
            }
        }
        {
            static const INTEGER values[] = {192, 840, 3399, 1997};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(67, j);
            }
        }
        {
            static const INTEGER values[] = {2253, 3168, 1321, 3121};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(68, j);
            }
        }
        {
            static const INTEGER values[] = {3491, 1499, 2271, 997};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(69, j);
            }
        }
        {
            static const INTEGER values[] = {2889, 1084, 3667, 1833};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(70, j);
            }
        }
        {
            static const INTEGER values[] = {2857, 3438, 2703, 2877};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(71, j);
            }
        }
        {
            static const INTEGER values[] = {2094, 2408, 629, 1633};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(72, j);
            }
        }
        {
            static const INTEGER values[] = {1818, 1589, 2365, 981};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(73, j);
            }
        }
        {
            static const INTEGER values[] = {688, 2391, 2431, 2009};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(74, j);
            }
        }
        {
            static const INTEGER values[] = {1407, 288, 1113, 941};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(75, j);
            }
        }
        {
            static const INTEGER values[] = {634, 26, 3922, 2449};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(76, j);
            }
        }
        {
            static const INTEGER values[] = {3231, 512, 2554, 197};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(77, j);
            }
        }
        {
            static const INTEGER values[] = {815, 1456, 184, 2441};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(78, j);
            }
        }
        {
            static const INTEGER values[] = {3524, 171, 2099, 285};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(79, j);
            }
        }
        {
            static const INTEGER values[] = {1914, 1677, 3228, 1473};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(80, j);
            }
        }
        {
            static const INTEGER values[] = {516, 2657, 4012, 2741};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(81, j);
            }
        }
        {
            static const INTEGER values[] = {164, 2270, 1921, 3129};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(82, j);
            }
        }
        {
            static const INTEGER values[] = {303, 2587, 3452, 909};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(83, j);
            }
        }
        {
            static const INTEGER values[] = {2144, 2961, 3901, 2801};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(84, j);
            }
        }
        {
            static const INTEGER values[] = {3480, 1970, 572, 421};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(85, j);
            }
        }
        {
            static const INTEGER values[] = {119, 1817, 3309, 4073};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(86, j);
            }
        }
        {
            static const INTEGER values[] = {3357, 676, 3171, 2813};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(87, j);
            }
        }
        {
            static const INTEGER values[] = {837, 1410, 817, 2337};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(88, j);
            }
        }
        {
            static const INTEGER values[] = {2826, 3723, 3039, 1429};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(89, j);
            }
        }
        {
            static const INTEGER values[] = {2332, 2803, 1696, 1177};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(90, j);
            }
        }
        {
            static const INTEGER values[] = {2089, 3185, 1256, 1901};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(91, j);
            }
        }
        {
            static const INTEGER values[] = {3780, 184, 3715, 81};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(92, j);
            }
        }
        {
            static const INTEGER values[] = {1700, 663, 2077, 1669};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(93, j);
            }
        }
        {
            static const INTEGER values[] = {3712, 499, 3019, 2633};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(94, j);
            }
        }
        {
            static const INTEGER values[] = {150, 3784, 1497, 2269};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(95, j);
            }
        }
        {
            static const INTEGER values[] = {2000, 1631, 1101, 129};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(96, j);
            }
        }
        {
            static const INTEGER values[] = {3375, 1925, 717, 1141};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(97, j);
            }
        }
        {
            static const INTEGER values[] = {1621, 3912, 51, 249};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(98, j);
            }
        }
        {
            static const INTEGER values[] = {3090, 1398, 981, 3917};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(99, j);
            }
        }
        {
            static const INTEGER values[] = {3765, 1349, 1978, 2481};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(100, j);
            }
        }
        {
            static const INTEGER values[] = {1149, 1441, 1813, 3941};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(101, j);
            }
        }
        {
            static const INTEGER values[] = {3146, 2224, 3881, 2217};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(102, j);
            }
        }
        {
            static const INTEGER values[] = {33, 2411, 76, 2749};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(103, j);
            }
        }
        {
            static const INTEGER values[] = {3082, 1907, 3846, 3041};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(104, j);
            }
        }
        {
            static const INTEGER values[] = {2741, 3192, 3694, 1877};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(105, j);
            }
        }
        {
            static const INTEGER values[] = {359, 2786, 1682, 345};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(106, j);
            }
        }
        {
            static const INTEGER values[] = {3316, 382, 124, 2861};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(107, j);
            }
        }
        {
            static const INTEGER values[] = {1749, 37, 1660, 1809};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(108, j);
            }
        }
        {
            static const INTEGER values[] = {185, 759, 3997, 3141};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(109, j);
            }
        }
        {
            static const INTEGER values[] = {2784, 2948, 479, 2825};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(110, j);
            }
        }
        {
            static const INTEGER values[] = {2202, 1862, 1141, 157};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(111, j);
            }
        }
        {
            static const INTEGER values[] = {2199, 3802, 886, 2881};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(112, j);
            }
        }
        {
            static const INTEGER values[] = {1364, 2423, 3514, 3637};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(113, j);
            }
        }
        {
            static const INTEGER values[] = {1244, 2051, 1301, 1465};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(114, j);
            }
        }
        {
            static const INTEGER values[] = {2020, 2295, 3604, 2829};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(115, j);
            }
        }
        {
            static const INTEGER values[] = {3160, 1332, 1888, 2161};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(116, j);
            }
        }
        {
            static const INTEGER values[] = {2785, 1832, 1836, 3365};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(117, j);
            }
        }
        {
            static const INTEGER values[] = {2772, 2405, 1990, 361};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(118, j);
            }
        }
        {
            static const INTEGER values[] = {1217, 3638, 2058, 2685};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(119, j);
            }
        }
        {
            static const INTEGER values[] = {1822, 3661, 692, 3745};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(120, j);
            }
        }
        {
            static const INTEGER values[] = {1245, 327, 1194, 2325};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(121, j);
            }
        }
        {
            static const INTEGER values[] = {2252, 3660, 20, 3609};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(122, j);
            }
        }
        {
            static const INTEGER values[] = {3904, 716, 3285, 3821};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(123, j);
            }
        }
        {
            static const INTEGER values[] = {2774, 1842, 2046, 3537};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(124, j);
            }
        }
        {
            static const INTEGER values[] = {997, 3987, 2107, 517};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(125, j);
            }
        }
        {
            static const INTEGER values[] = {2573, 1368, 3508, 3017};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(126, j);
            }
        }
        {
            static const INTEGER values[] = {1148, 1848, 3525, 2141};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(127, j);
            }
        }
        {
            static const INTEGER values[] = {545, 2366, 3801, 1537};
            data_of_type<INTEGER> data(FEM_VALUES_AND_SIZE);
            for (j = 1; j <= 4; j = j + 1) {
                data, mm(128, j);
            }
        }
    }
    INTEGER i1 = 0;
    INTEGER i2 = 0;
    INTEGER i3 = 0;
    INTEGER i4 = 0;
    INTEGER i = 0;
    INTEGER it4 = 0;
    const INTEGER ipw2 = 4096;
    INTEGER it3 = 0;
    INTEGER it2 = 0;
    INTEGER it1 = 0;
    const REAL one = 1.0;
    const REAL r = one / ipw2;
    //
    //  -- LAPACK auxiliary routine --
    //  -- LAPACK is a software package provided by Univ. of Tennessee,    --
    //  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
    //
    //     .. Scalar Arguments ..
    //     ..
    //     .. Array Arguments ..
    //     ..
    //
    //  =====================================================================
    //
    //     .. Parameters ..
    //     ..
    //     .. Local Scalars ..
    //     ..
    //     .. Local Arrays ..
    //     ..
    //     .. Intrinsic Functions ..
    //     ..
    //     .. Data statements ..
    //     ..
    //     .. Executable Statements ..
    //
    i1 = iseed[1 - 1];
    i2 = iseed[2 - 1];
    i3 = iseed[3 - 1];
    i4 = iseed[4 - 1];
    //
    for (i = 1; i <= min(n, lv); i = i + 1) {
    //
    statement_20:
        //
        //        Multiply the seed by i-th power of the multiplier modulo 2**48
        //
        it4 = i4 * mm[(i - 1) + (4 - 1) * ldmm];
        it3 = it4 / ipw2;
        it4 = it4 - ipw2 * it3;
        it3 += i3 * mm[(i - 1) + (4 - 1) * ldmm] + i4 * mm[(i - 1) + (3 - 1) * ldmm];
        it2 = it3 / ipw2;
        it3 = it3 - ipw2 * it2;
        it2 += i2 * mm[(i - 1) + (4 - 1) * ldmm] + i3 * mm[(i - 1) + (3 - 1) * ldmm] + i4 * mm[(i - 1) + (2 - 1) * ldmm];
        it1 = it2 / ipw2;
        it2 = it2 - ipw2 * it1;
        it1 += i1 * mm[(i - 1) + (4 - 1) * ldmm] + i2 * mm[(i - 1) + (3 - 1) * ldmm] + i3 * mm[(i - 1) + (2 - 1) * ldmm] + i4 * mm[(i - 1)];
        it1 = mod(it1, ipw2);
        //
        //        Convert 48-bit INTEGEReger to a real number in the INTEGERerval (0,1)
        //
        x[i - 1] = r * (it1.real() + r * (it2.real() + r * (it3.real() + r * it4.real())));
        //
        if (x[i - 1] == 1.0) {
            //           If a real number has n bits of precision, and the first
            //           n bits of the 48-bit INTEGEReger above happen to be all 1 (which
            //           will occur about once every 2**n calls), then X( I ) will
            //           be rounded to exactly 1.0.
            //           Since X( I ) is not supposed to return exactly 0.0 or 1.0,
            //           the statistically correct thing to do in this situation is
            //           simply to iterate again.
            //           N.B. the case X( I ) = 0.0 should not be possible.
            i1 += 2;
            i2 += 2;
            i3 += 2;
            i4 += 2;
            goto statement_20;
        }
        //
    }
    //
    //     Return final value of seed
    //
    iseed[1 - 1] = it1;
    iseed[2 - 1] = it2;
    iseed[3 - 1] = it3;
    iseed[4 - 1] = it4;
    //
    //     End of Rlaruv
    //
}
