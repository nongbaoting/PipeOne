#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2020/7/2 8:58
# @author  : Baoting Nong'
# @email   : 523135753@qq.com'
import fire


def myrun(model_trained, testX, testy, testfeature_name, Imp_f ):

    w_idx_test = [i for i, ft in enumerate(testfeature_name) if ft in list(Imp_f)]
    testX_sub, testy_sub, testfeature_name_sub = testX[:, w_idx_test], testy, [ testfeature_name[i] for i in w_idx_test ]

    y_test_pred = model_trained.predict(testX_sub)
    test_sensitivity, test_specificity = sensitivity_specificity(testy, y_test_pred)
    print(test_sensitivity, test_specificity)

    return  test_sensitivity, test_specificity

