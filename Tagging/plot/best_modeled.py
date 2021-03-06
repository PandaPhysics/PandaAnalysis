#!/usr/bin/env python

formulae = [ 
            ('fj1ECFN_1_2_05/pow(fj1ECFN_1_2_40,0.12)',0,.6),
            ('fj1ECFN_1_2_10/pow(fj1ECFN_1_2_20,0.50)',0,1),
            ('fj1ECFN_1_2_10/pow(fj1ECFN_1_2_40,0.25)',0,.6),
            ('fj1ECFN_1_2_20/pow(fj1ECFN_1_2_10,2.00)',2,10),
            ('fj1ECFN_1_2_20/pow(fj1ECFN_1_2_40,0.50)',0,.8),
            ('fj1ECFN_1_2_40/pow(fj1ECFN_1_2_20,2.00)',3,12),
            ('fj1ECFN_1_3_05/fj1ECFN_1_2_05',0,.5),
            ('fj1ECFN_1_3_20/pow(fj1ECFN_1_2_05,4.00)',0,4),
            ('fj1ECFN_1_3_20/fj1ECFN_2_3_10',0,1),
            ('fj1ECFN_1_3_20/pow(fj1ECFN_3_3_05,1.33)',0,2),
            ('fj1ECFN_1_3_40/pow(fj1ECFN_1_2_20,2.00)',0,.25),
            ('fj1ECFN_1_3_40/pow(fj1ECFN_2_3_10,2.00)',2,17),
            ('fj1ECFN_1_3_40/fj1ECFN_2_3_20',0,1),
            ('fj1ECFN_1_3_40/pow(fj1ECFN_2_3_40,0.50)',0,.1),
            ('fj1ECFN_1_3_40/pow(fj1ECFN_3_3_10,1.33)',0,2),
            ('fj1ECFN_1_3_40/pow(fj1ECFN_3_3_20,0.67)',0,.5),
            ('fj1ECFN_2_3_10/pow(fj1ECFN_1_2_05,4.00)',1,7),
            ('fj1ECFN_2_3_10/fj1ECFN_1_3_20',1,4),
            ('fj1ECFN_2_3_10/pow(fj1ECFN_3_3_05,1.33)',1,2.5),
            ('fj1ECFN_2_3_10/pow(fj1ECFN_3_3_40,0.17)',0,.1),
            ('fj1ECFN_2_3_20/pow(fj1ECFN_1_2_10,4.00)',0,10),
            ('fj1ECFN_2_3_20/pow(fj1ECFN_1_2_20,2.00)',0,.5),
            ('fj1ECFN_2_3_20/fj1ECFN_1_3_40',1,5),
            ('fj1ECFN_2_3_40/pow(fj1ECFN_3_3_20,1.33)',0,3.5),
            ('fj1ECFN_3_3_05/pow(fj1ECFN_1_2_05,3.00)',.5,3.5),
            ('fj1ECFN_3_3_05/pow(fj1ECFN_1_3_20,0.75)',.5,3),
            ('fj1ECFN_3_3_10/pow(fj1ECFN_1_2_10,3.00)',0,5),
            ('fj1ECFN_3_3_10/pow(fj1ECFN_1_3_40,0.75)',0.5,4),
            ('fj1ECFN_3_3_10/pow(fj1ECFN_2_3_20,0.75)',0,2),
            ('fj1ECFN_3_3_20/pow(fj1ECFN_2_3_40,0.75)',0,2),
            ('fj1ECFN_3_3_20/pow(fj1ECFN_3_3_40,0.50)',0,.5),
            ('fj1ECFN_3_3_40/pow(fj1ECFN_1_2_40,3.00)',0,10),
            ('fj1ECFN_1_4_20/pow(fj1ECFN_1_3_10,2.00)',0,2),
            ('fj1ECFN_1_4_20/pow(fj1ECFN_2_3_05,2.00)',0,1),
            ('fj1ECFN_1_4_40/pow(fj1ECFN_1_3_20,2.00)',0,2.5),
            ('fj1ECFN_2_4_05/pow(fj1ECFN_1_3_05,2.00)',1,3),
            ('fj1ECFN_2_4_10/pow(fj1ECFN_1_3_10,2.00)',1,4),
            ('fj1ECFN_2_4_10/pow(fj1ECFN_2_3_05,2.00)',0,1.5),
            ('fj1ECFN_2_4_20/pow(fj1ECFN_1_2_05,8.00)',0,10),
            ('fj1ECFN_2_4_20/pow(fj1ECFN_1_3_20,2.00)',0,5),
            ('fj1ECFN_2_4_20/pow(fj1ECFN_3_3_05,2.67)',0,4),
            ('fj1ECFN_2_4_20/pow(fj1ECFN_2_4_40,0.50)',0,.1),

    ]

