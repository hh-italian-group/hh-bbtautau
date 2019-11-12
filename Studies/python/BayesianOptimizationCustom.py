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
                 probed_points_opt_output, num_iter, random_state):
        with open(param_ranges_json) as json_file:
            self.params_typed_range = json.load(json_file)
        self.num_iter = num_iter
        self.random_state = random_state
        self.params_ranges = {}
        self.probed_points = []
        self.probed_points_opt = []
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
            random_state=self.random_state,
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

        for key, value in params.items():
            #Has some dependency wiht other variables
            if 'depends_on' in self.params_typed_range[key]:
                # the variable on which it depends has a value of zero
                if new_params[self.params_typed_range[key]['depends_on'] ] == 0:
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

    def FindResult(self, p, is_target_point):
        points = self.probed_points if is_target_point else self.probed_points_opt
        for p2 in points:
            if self.ComparePoints(p, p2["params"]):
                return p2["target"]
        return None

    def LoadPoints(self, probed_points_target_output_prev, probed_points_opt_output_prev):
        target_points = []
        for line in open(probed_points_target_output_prev, 'r'):
            target_points.append(json.loads(line))

        opt_points = []
        for line in open(probed_points_opt_output_prev, 'r'):
            opt_points.append(json.loads(line))

        for p in target_points:
            result = self.FindResult(p['params'], True)
            # if has not been tested before
            if result is None:
                self.probed_points.append(p)
                with open(self.probed_points_target_output, 'a') as f:
                    f.write(json.dumps(p) + '\n')

        for p in opt_points:
            result = self.FindResult(p['params'], False)
            # if has not been tested before
            if result is None:
                self.probed_points_opt.append(p)

                with open(self.probed_points_opt_output, 'a') as f:
                    f.write(json.dumps(p) + '\n')

                self.optimizer.register(params=p['params'],target=p['target'])

    def MaximizeStep(self, p, is_target_point=False):
        if is_target_point:
            p_target = p
            p_opt = self.InverseTransformParams(p)
        else:
            p_target = self.TransformParams(p)
            p_opt = p

        if self.FindResult(p_opt, False) is not None:
            return 0

        result = self.FindResult(p_target, True)

        #If the point haven't been tested, save it

        if result is None:
            print('starting evaluating for the params:\n{}'.format(json.dumps(p_target, indent=4)))
            result = self.fn(**p_target)
            print('result: ',result)
            has_been_tested = 1
        else:
            has_been_tested = 0

        entry_target = {"target": float(result), "params": p_target}
        entry_opt = {
            "target": float(result),
            "params": { key: float(value) for key,value in p_opt.items() }
        }

        if has_been_tested == 1:
            with open(self.probed_points_target_output, 'a') as f:
                f.write(json.dumps(entry_target) + '\n')

            self.probed_points.append(entry_target)

        self.probed_points_opt.append(entry_opt)

        self.optimizer.register(params=p_opt,target=result)
        with open(self.probed_points_opt_output, 'a') as f:
            f.write(json.dumps(entry_opt) + '\n')

        return has_been_tested

    def maximize(self, n_iter, acq='ucb', kappa=2.576, xi=0.0):

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

        p_opt_max = self.optimizer.max
        p_target_max = self.TransformParams(p_opt_max['params'])

        # for each param var up and down and try with edges parameter
        for key, values in self.params_typed_range.items():
            new_target_point = copy.deepcopy(p_target_max)
            for v in values['range']:
                new_target_point[key] = v
                self.MaximizeStep(new_target_point, is_target_point=True)

        return self.TransformParams(self.optimizer.max['params']), self.optimizer.max['target']
