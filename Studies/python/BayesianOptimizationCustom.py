# Definition of the model optimization using Bayesian optimization as a method
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

from bayes_opt import UtilityFunction
from bayes_opt import BayesianOptimization
from bayes_opt.observer import JSONLogger
from bayes_opt.event import Events

import json
import copy
import numpy as np

class BayesianOptimizationCustom:
    def __init__(self, param_ranges_json, params_init_probe_json, fn, probed_points_target_output,
                 probed_points_opt_output, num_iter, num_init_points, **kwargs):
        with open(param_ranges_json) as json_file:
            self.params_typed_range = json.load(json_file)

        self.params_ranges = {}
        self.probed_points = []
        self.probed_points_target_output =  probed_points_target_output
        self.probed_points_opt_output = probed_points_opt_output

        for key,values in self.params_typed_range.items():
            if values['type'] in ["int", "float"]:
                self.params_ranges[key] = tuple(values['range'])
            elif  values['type'] ==  "list" :
                self.params_ranges[key] = (0, len(values['range']) -1 )

        with open(params_init_probe_json) as json_file:
            self.init_points_to_probe = json.load(json_file)

        self.fn = fn
        self.optimizer = BayesianOptimization(
            f=None,
            pbounds=self.params_ranges,
            random_state=1,
            verbose=1
        )

    def TransformParams(self, params):
        new_params = {}
        for key, value in params.items():
            if self.params_typed_range[key]['type'] == 'int':
                new_params[key] = int(round(value))
            elif self.params_typed_range[key]['type'] == 'float':
                new_params[key] = float(value)
            elif self.params_typed_range[key]['type']== 'list':
                new_params[key] = self.params_typed_range[key]['range'][int(round(value))]
            else:
                raise Exception('type {} not allowed'.format(self.params_typed_range[key]['type']))

            if self.params_typed_range[key]['depends_on'] != 'none':
                if params[self.params_typed_range[key]['depends_on'] ] == 0:
                    new_params[key] = self.params_typed_range[key]['range'][0]
        return new_params

    def InverseTransformParams(self, params):
        new_params = params.copy()
        for key, values in params.items():
            if self.params_typed_range[key]['type'] == 'list':
                new_params[key] = self.params_typed_range[key]['range'].index(values)
        return new_params

    def ComparePoints(self, p1, p2):
        are_equal = True
        for key, values in p1.items():
            if p1[key] != p2[key]:
                return False
        return True

    def FindResult(self, p):
        for p2 in self.probed_points:
            if self.ComparePoints(p, p2["params"]):
                return p2["target"]
        return None

    def MaximizeStep(self, p, is_target_point=False):
        if is_target_point:
            p_target = p
            p_opt = self.InverseTransformParams(p)
        else:
            p_target = self.TransformParams(p)
            p_opt = p

        result = self.FindResult(p_target)
        #If the point haven't been tested, save it
        if result is None:
            result = self.fn(**p_target)
            has_been_tested = 1
        else:
            has_been_tested = 0

        entry_target = {"target": float(result), "params": p_target}
        entry_opt = {
            "target": float(result),
            "params": { key: float(value) for key,value in p_opt.items() }
        }

        with open(self.probed_points_opt_output, 'a') as f:
            f.write(json.dumps(entry_opt) + '\n')

        with open(self.probed_points_target_output, 'a') as f:
            f.write(json.dumps(entry_target) + '\n')

        self.probed_points.append(entry_target)

        equal = True
        for p in self.probed_points:
            for key, values in p['params'].items():
                if p['params'][key] != entry_target['params'][key]:
                    equal = False
        if equal == False:
            self.optimizer.register(params=p_opt,target=result)

        return has_been_tested

    def maximize(self,
                 init_points=num_init_points,
                 n_iter=num_iter,
                 acq='ucb',
                 kappa=2.576,
                 xi=0.0,
                 **gp_params):

        #Probe with all known good points
        for p in self.init_points_to_probe:
            self.MaximizeStep(p, is_target_point=True)

        # Suggest-Evaluate-Register paradigm
        utility = UtilityFunction(kind=acq, kappa=kappa, xi=xi)

        n = 0
        while n < n_iter:
            #create point to probe
            p_opt = self.optimizer.suggest(utility)
            n += self.MaximizeStep(p_opt, is_target_point=False)

      # for each param var up and down and try with edges parameter
        for key, values in params_typed_range.items():
            for n in self.params_range[key]:
                new_p_opt = copy.deepcopy(p_opt)
                p_opt[key] = n
                self.MaximizeStep(p_opt, is_target_point=True)

        p_opt_max = p_opt.max
        p_target_max = self.TransformParams(p_opt_max)

        return p_target_max
