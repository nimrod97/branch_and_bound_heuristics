from scipy.optimize import linprog
import numpy as np
import math


class BNB():
    def __init__(self, c, A_ub, b_ub, A_eq, b_eq, bounds):
        self.c = c
        self.A_ub = A_ub
        self.b_ub = b_ub
        self.A_eq = A_eq
        self.b_eq = b_eq
        self.bounds = bounds
        self.res = linprog(self.c, self.A_ub, self.b_ub, self.A_eq, self.b_eq, self.bounds, method='revised simplex')
        self.count = 0
        self.optimal = 0
        self.vars = []
        self.stop = False

    """ 
    ordering the vars according to their appearances in the constraints equations before solving
    the problem
    """

    def solve_by_most_appearances_in_constraints(self):
        # calculating for each variable in the problem, its number of appearances in the constraints
        # equations
        num_of_vars = len(self.A_ub[0]) if self.A_ub is not None else len(self.A_eq[0])
        appearances = {}
        for i in range(num_of_vars):
            appearances[i] = 0
        try:
            for eq in self.A_ub:
                for i in range(len(eq)):
                    if eq[i] != 0:
                        appearances[i] += 1
        except:
            {}
        try:
            if self.A_eq.all() != None:
                for eq in self.A_eq:
                    for i in range(len(eq)):
                        if eq[i] != 0:
                            appearances[i] += 1
        except:
            {}

        # running the algorithm according to the order of appearences of each var
        order = dict(sorted(appearances.items(), key=lambda x: x[1], reverse=True))
        order = list(order)
        while not self.are_all_integers(self.res.x):
            if self.stop:
                break
            for i in range(len(order)):
                x = self.res.x[order[i]]
                if not x.is_integer():
                    self.helper(order[i])

        # the solution is not feasible
        if self.stop or self.res.success != True:
            print(self.res)
            return

        if self.optimal == 0:
            self.optimal = self.res.fun
            self.vars = self.res.x

        print(f"optimal value: {self.optimal}")
        print(f"variables: {self.vars}")
        print(f"num of branches: {self.count}")

    """
    solving the problem according to the regular order of the variables: x1,x2,...,x_n
    """

    def solve_by_var_order(self):
        while not self.are_all_integers(self.res.x):
            if self.stop == True:
                break
            for i in range(len(self.res.x)):
                x = self.res.x[i]
                if not x.is_integer():
                    self.helper(i)

        # the solution is not feasible
        if self.stop or self.res.success != True:
            print(self.res)
            return

        if self.optimal == 0:
            self.optimal = self.res.fun
            self.vars = self.res.x
        print(f"optimal value: {self.optimal}")
        print(f"variables: {self.vars}")
        print(f"num of branches: {self.count}")

    """
    every iteration, the algorithm will choose the vars order according to their distance from
    their nearest integer
    """

    def solve_by_farest_num_from_integer_first(self):
        # compute distance of each number in the result from the closest integer number
        # the max distance is 0.5
        distances = {}
        for i in range(len(self.res.x)):
            x = self.res.x[i]
            distances[i] = (min(abs(x - math.ceil(x)), abs(x - math.floor(x))))
        distances = dict(sorted(distances.items(), key=lambda x: x[1], reverse=True))
        if distances[list(distances)[0]] == 0:
            print(self.res.fun)
            print(self.res.x)
            print(self.count)
            return
        while distances[list(distances)[0]] != 0:
            if self.stop:
                break
            key = list(distances)[0]
            x = self.res.x[key]
            if not x.is_integer():
                self.helper(key)
                # update the distances dict for finding new float values
                for i in range(len(self.res.x)):
                    x = self.res.x[i]
                    distances[i] = (min(abs(x - math.ceil(x)), abs(x - math.floor(x))))
                distances = dict(sorted(distances.items(), key=lambda x: x[1], reverse=True))

        # the solution is not feasible
        if self.stop or self.res.success != True:
            print(self.res)
            return

        if self.optimal == 0:
            self.optimal = self.res.fun
            self.vars = self.res.x

        print(f"optimal value: {self.optimal}")
        print(f"variables: {self.vars}")
        print(f"num of branches: {self.count}")

    """
    checking if all the variables in the problem are integers
    """

    def are_all_integers(self, vars):
        for x in vars:
            if not x.is_integer():
                return False
        return True

    """
    recursive function that uses the branch and bound if there is a need to.
    every method is using this function throughout it's solution
    """

    def helper(self, i):
        x = self.res.x[i]
        if x.is_integer():
            return

        # new branching to 2 subproblems
        self.count += 1
        last_bound = self.bounds[i]
        bound1 = (self.bounds[i][0], math.floor(x))
        bound2 = (math.ceil(x), self.bounds[i][1])
        self.bounds[i] = bound1
        try:
            res1 = linprog(self.c, self.A_ub, self.b_ub, self.A_eq, self.b_eq, self.bounds,
                           method='revised simplex')
        except:
            self.res.x[i] = self.bounds[i][0]
            bound1 = (self.bounds[i][0], self.bounds[i][0])
            self.bounds[i] = bound1
            res1 = self.res

        if self.are_all_integers(res1.x) and res1.success == True:
            if res1.fun < self.optimal:
                self.optimal = res1.fun
                self.vars = res1.x
        self.bounds[i] = bound2
        try:
            res2 = linprog(self.c, self.A_ub, self.b_ub, self.A_eq, self.b_eq, self.bounds,
                           method='revised simplex')
        except:
            self.res.x[i] = self.bounds[i][1]
            bound2 = (self.bounds[i][1], self.bounds[i][1])
            self.bounds[i] = bound2
            res2 = self.res

        if self.are_all_integers(res2.x) and res2.success == True:
            if res2.fun < self.optimal:
                self.optimal = res2.fun
                self.vars = res2.x

        if res1.fun < res2.fun:
            if res1.success == True:
                self.res = res1
                self.bounds[i] = bound1
                self.stop = False
            elif res2.success == True:
                self.res = res2
                self.bounds[i] = bound2
                self.stop = False
            else:
                self.bounds[i] = last_bound
                self.res = res1
                self.stop = True
                return
        elif res1.fun > res2.fun:
            if res2.success == True:
                self.res = res2
                self.bounds[i] = bound2
                self.stop = False
            elif res1.success == True:
                self.res = res1
                self.bounds[i] = bound1
                self.stop = False
            else:
                self.bounds[i] = last_bound
                self.res = res2
                self.stop = True
                return

        else:  # dealing with edge cases
            if bound1[0] != bound1[1]:
                if bound2[0] != bound2[1]:
                    if res1.success == True:
                        self.res = res1
                        self.bounds[i] = bound1
                        self.stop = False
                    elif res2.success == True:
                        self.res = res2
                        self.bounds[i] = bound2
                        self.stop = False
                    else:
                        self.bounds[i] = last_bound
                        self.res = res2
                        self.stop = True
                        return
                else:
                    self.res = res2
                    self.bounds[i] = bound2
                    self.stop = False
            else:
                self.res = res1
                self.bounds[i] = bound1
                self.stop = False

        self.helper(i)


def test1():
    c = np.array([1, 1, 1, 1, 1, 1, 1])
    A_ub = np.array(
        [[-1, 0, 0, -1, -1, -1, -1], [-1, -1, 0, 0, -1, -1, -1], [-1, -1, -1, 0, 0, -1, -1], [-1, -1, -1, -1, 0, 0, -1],
         [-1, -1, -1, -1, -1, 0, 0],
         [0, -1, -1, -1, -1, -1, 0], [0, 0, -1, -1, -1, -1, -1]])
    b_ub = np.array([-8, -6, -5, -4, -6, -7, -9])
    A_eq = None
    b_eq = None
    bounds = []
    for i in range(7):
        bounds.append((0, None))
    solver = BNB(c, A_ub, b_ub, A_eq, b_eq, bounds)

    solver.solve_by_var_order()
    # solver.solve_by_farest_num_from_integer_first()
    # solver.solve_by_most_appearances_in_constraints()


def test2():
    c = np.array([-600, -700])
    A = np.array([[2, 3], [6, 5]])
    b = np.array([12, 30])
    Aeq = None
    beq = None
    bounds = []
    for i in range(2):
        bounds.append((0, None))
    solver = BNB(c, A, b, Aeq, beq, bounds)

    solver.solve_by_var_order()
    # solver.solve_by_farest_num_from_integer_first()
    # solver.solve_by_most_appearances_in_constraints()


def test3():
    c = np.array([-1, 3, 4, -2, 5, -2, 1])
    A = np.array([[1, 2, -1, 1.5, -1, 1, 3.75], [0, -2, 0, 1, 0, 0, 0], [0, 0, 1, 3, 0, 0, 0], [4.5, 0, 0, -5, 0, 2, 0],
                  [0, 0, 3, -7.75, 0, 0, 0],
                  [0, 0, 0, 0, 0, -1, -1], [0, 0, 0, -3, -2, 0, 0]])
    b = np.array([-10, -7, 15, -5, 0, -8, -10])
    A_eq = None
    b_eq = None
    bounds = []
    bounds.append((1, None))
    for i in range(6):
        bounds.append((0, None))
    solver = BNB(c, A, b, A_eq, b_eq, bounds)

    solver.solve_by_var_order()
    # solver.solve_by_farest_num_from_integer_first()
    # solver.solve_by_most_appearances_in_constraints()


def test4():
    c = np.array([-5, -5, -6, -7])
    A = np.array([[2, 4, 0, 5], [0, -4, 0, 3], [0, 2, 5, 3]])
    b = np.array([10, 6, 11])
    A_eq = None
    b_eq = None
    bounds = []
    for i in range(4):
        bounds.append((-5, None))
    solver = BNB(c, A, b, A_eq, b_eq, bounds)

    solver.solve_by_var_order()
    # solver.solve_by_farest_num_from_integer_first()
    # solver.solve_by_most_appearances_in_constraints()


if __name__ == '__main__':
    print("test 1 results")
    test1()
    print("test 2 results")
    test2()
    print("test 3 results")
    test3()
    print("test 4 results")
    test4()
