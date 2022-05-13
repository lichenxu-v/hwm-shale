#!/usr/bin/env python

from sympy import *
from random import random
from copy import copy

class Supply:
    def __init__(self, file):
        self.pair = {}
        self.satisfy_demand = {}
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                i, s = line.split('\t')
                self.pair[i] = int(s)

    def get_supply(self, i):
        return self.pair[i]

    def get_satisfy_demand(self, i):
        return self.satisfy_demand[i]

    def get_all_i(self):
        return list(self.pair.keys())


class Demand:
    def __init__(self, file, supply):
        self.demand = {}
        self.penalty = {}
        self.target_supply = {}
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                j, d, p, ii = line.split('\t')
                self.demand[j] = int(d)
                self.penalty[j] = float(p)
                self.target_supply[j] = ii.split(',')
        self._set_supply_satisfy_demand(supply)

    def _set_supply_satisfy_demand(self, supply):
        for (j, ii) in self.target_supply.items():
            for i in ii:
                if i not in supply.satisfy_demand:
                    supply.satisfy_demand[i] = []
                supply.satisfy_demand[i].append(j)

    def get_demand(self, j):
        return self.demand[j]

    def get_penalty(self, j):
        return self.penalty[j]

    def get_target_supply(self, j):
        return self.target_supply[j]
    
    def get_all_j(self):
        return list(self.demand.keys())

    def get_v(self, j):
        return 1.0


class HWM:
    def __init__(self, supply, demand):
        self.supply = supply
        self.demand = demand
        self.solve = Solve()

    def initialize(self):
        self.r_i = {}
        self.s_i = {}
        self.d_j = {}
        self.alpha_j = {}
        self.z_j = {}
        for j in self.demand.get_all_j():
            sum = 0.0
            for i in self.demand.get_target_supply(j):
                sum += self.supply.get_supply(i)
            self.d_j[j] = self.demand.get_demand(j) 
            self.z_j[j] = self.demand.get_demand(j) / sum
            print("demand order --->", j, self.demand.get_demand(j), sum, self.z_j[j])
        

    def update_alpha(self):
        for i in self.supply.get_all_i():
            self.r_i[i] = self.supply.get_supply(i)
            self.s_i[i] = self.r_i[i]
        z_j_l = sorted(self.z_j.items(), key=lambda d:d[0], reverse=True)

        print(z_j_l)
        for j, z_j in z_j_l:
            d_j = self.demand.get_demand(j)
            print("z_j:",j, z_j, d_j)
            alpha = Symbol('alpha')
            flag = True
            for i in self.demand.get_target_supply(j):
                print("\tz_i:",i, self.r_i[i], self.s_i[i])
                a = self.r_i[i] / self.s_i[i]
                if flag:
                    sum = Piecewise((self.r_i[i], alpha >= a), (self.s_i[i] * alpha, alpha < a)) 
                    flag = False
                else:
                    sum += Piecewise((self.r_i[i], alpha >= a), (self.s_i[i] * alpha, alpha < a)) 
            print(sum-d_j+1, alpha)
            result = solve(sum - d_j + 1, alpha)

            if len(result) == 0 or result[0] < 0.0:
                self.alpha_j[j] = 1.0
            else:
                self.alpha_j[j] = result[0]

            print("alpha_j:",self.alpha_j[j])

            for i in self.demand.get_target_supply(j):
                self.r_i[i] -= min(self.r_i[i], self.s_i[i] * self.alpha_j[j])


    def output(self):
        print('output alpha_j --->', sorted(self.alpha_j.items(), key=lambda d:d[0]))
        print(self.alpha_j)



class Online:
    def __init__(self, supply, demand, alpha_j):
        self.supply = supply
        self.demand = demand
        self.z_j = {}
        self.alpha_j = alpha_j
        self.allocation_j = {}
        self.remaind_i = {}
        for i in self.supply.get_all_i():
            self.remaind_i[i] = supply.get_supply(i)
        for j in self.demand.get_all_j():
            sum = 0.0
            for i in self.demand.get_target_supply(j):
                sum += self.supply.get_supply(i)
            self.z_j[j] = self.demand.get_demand(j) / sum
            self.allocation_j[j] = 0
        self.solve = Solve()


    def allocation(self, i):
        x_ij = {}
        for j in self.supply.get_satisfy_demand(i):
            x_ij[j] = self.alpha_j[j]
        
        #sum = 0.0
        #for (j, p) in x_ij.items():
        #    sum += p
        #if sum < 1.0:
        #    print('there is %f chance that no conract is selected' % (1.0 - sum))

        r = random()
        sum = 0.0
        for (j, p) in x_ij.items():
            sum += p
            if r < sum:
                self.allocation_j[j] += 1
                self.remaind_i[i] -= 1
                break


class Solve:
    # http://math.stackexchange.com/questions/145458/solve-equations-using-the-max-function
    def max(self, coef, y):
        solutions = []
        sum_cons = 0.0
        sum_coef = 0.0
        for i in range(len(coef)):
            sum_cons += coef[i][1]
            sum_coef += coef[i][2]
        coef = sorted(coef, key=lambda t:t[0])
        res = (sum_cons - y) / sum_coef
        if res <= coef[0][0]:
            solutions.append(res)
        for i in range(1, len(coef)):
            sum_cons -= coef[i-1][1]
            sum_coef -= coef[i-1][2]
            res = (sum_cons - y) / sum_coef
            if res <= coef[i][0]:
                solutions.append(res)
                break
        return solutions

    def max2(self):
        pass

    #http://www.bitsofpancake.com/math/minimum-and-maximum-of-two-functions/
    #http://math.stackexchange.com/questions/602553/how-to-invert-max-and-min-operators-in-equations
    def minmax(self):
        pass
                


class Debug:
    def __init__(self, hwm, online):
        self.hwm = hwm
        self.online = online

    def print_supply(self):
        ii = self.hwm.supply.get_all_i()
        ii.sort()
        print("\nsupply:")
        print("supply_node\tinventory\tsatisfy_demand")
        for i in ii:
            print('%s\t\t%d\t\t%s' % (i, self.hwm.supply.get_supply(i), ','.join(self.hwm.supply.get_satisfy_demand(i))))

    def print_demand(self):
        jj = self.hwm.demand.get_all_j()
        jj.sort()
        print("\ndemand:")
        print("demand_node\tdemand\tpenalty\ttarget_supply")
        for j in jj:
            print('%s\t\t%d\t%f\t%s' % (j, self.hwm.demand.get_demand(j), self.hwm.demand.get_penalty(j), \
                                        ','.join(self.hwm.demand.get_target_supply(j))))

    def print_online_allocation(self):
        jj = self.online.demand.get_all_j()
        jj.sort()
        print("\nallocation:")
        print("demand_node\tdemand\t\tallocation")
        for j in jj:
            print('%s\t\t%d\t\t%d' % (j, self.online.demand.get_demand(j), self.online.allocation_j[j]))

    def print_online_remaind(self):
        ii = self.online.supply.get_all_i()
        ii.sort()
        print("\nremaind:")
        print("supply_node\tinventory\tremaind")
        for i in ii:
            print('%s\t\t%d\t\t%s' % (i, self.online.supply.get_supply(i), self.online.remaind_i[i]))


def main():
    supply = Supply('./supply.txt')
    demand = Demand('./demand.txt', supply)

    hwm = HWM(supply, demand)
    hwm.initialize()
    hwm.update_alpha()
    hwm.output()

    online = Online(supply, demand, hwm.alpha_j)
    for i in supply.get_all_i():
        inventory = supply.get_supply(i)
        while inventory > 0:
            online.allocation(i)
            inventory -= 1

    debug = Debug(hwm, online)
    debug.print_supply()
    debug.print_demand()
    debug.print_online_allocation()
    debug.print_online_remaind()



if __name__ == '__main__' :
    main()
