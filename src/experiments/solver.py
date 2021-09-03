import abc
import tempfile
import itertools
from typing import List, Tuple

import dwave_networkx.algorithms.tsp as dnx
import networkx as nx
import numpy as np
from dwave.system import LeapHybridSampler
from ortools.constraint_solver import pywrapcp, routing_enums_pb2
from input import *
from vrp_problem import VRPProblem
from vrp_solvers import FullQuboSolver, SolutionPartitioningSolver, DBScanSolver


class Solver(metaclass=abc.ABCMeta):
    def __init__(self, distance_matrix: np.array) -> None:
        self.distance_matrix = distance_matrix

    @abc.abstractmethod
    def solveTSP() -> Tuple[List[int], int]:
        pass

    def name(self):
        return type(self).__name__

    def calculate_cost(self, path):
        path_cost = 0
        for i in range(len(path)):
            a = i % len(path)
            b = (i+1) % len(path)
            path_cost += self.distance_matrix[path[a]][path[b]]
        return int(path_cost)


class GoogleORToolsSolver(Solver):
    """ Code taken from https://developers.google.com/optimization/routing/tsp,
    with minor alterations. """

    def __init__(self, distance_matrix: np.array) -> None:
        super().__init__(distance_matrix)

    def solveTSP(self) -> Tuple[List[int], int]:
        manager = pywrapcp.RoutingIndexManager(len(self.distance_matrix), 1, 0)
        routing = pywrapcp.RoutingModel(manager)

        # Define cost of each arc.
        routing.SetArcCostEvaluatorOfAllVehicles(routing.RegisterTransitCallback(
            lambda from_idx, to_idx: self.distance_matrix[manager.IndexToNode(from_idx)][manager.IndexToNode(to_idx)]))

        # Setting first solution heuristic.
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = (
            routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

        # Solve the problem.
        solution = routing.SolveWithParameters(search_parameters)

        path = self.__get_routes(solution, routing, manager)[0]

        return path, super().calculate_cost(path)

    def __get_routes(self, solution, routing, manager):
        """Get vehicle routes from a solution and store them in an array."""
        routes = []
        for route_nbr in range(routing.vehicles()):
            index = routing.Start(route_nbr)
            route = [manager.IndexToNode(index)]
            while not routing.IsEnd(index):
                index = solution.Value(routing.NextVar(index))
                route.append(manager.IndexToNode(index))
            routes.append(route)
        return routes


class DWaveLeapHybridSolver(Solver):
    def __init__(self, distance_matrix: np.array, max_time: int, max_attempts: int) -> None:
        super().__init__(distance_matrix)
        self.max_time = max_time
        self.max_attempts = max_attempts

    def __distance_matrix2complete_graph(self) -> nx.Graph:
        N = len(self.distance_matrix)
        G = nx.complete_graph(N)
        edges = [(i, j, self.distance_matrix[i][j]) for i, j in filter(
            lambda x: x[0] < x[1], list(itertools.product(range(N), range(N))))]
        G.add_weighted_edges_from(edges)
        return G

    def solveTSP(self) -> Tuple[List[int], int]:
        routes = []
        G = self.__distance_matrix2complete_graph()
        for attempt in range(self.max_attempts):
            routes.append(dnx.traveling_salesperson(
                G, LeapHybridSampler(), start=0, time_limit=self.max_time))
        path = routes[0]
        return path, super().calculate_cost(path)


class DWaveVRPSolver(Solver):
    def __init__(self, distance_matrix: np.array, solver_type: str) -> None:
        super().__init__(distance_matrix)
        self.solver_type = solver_type

    def __createVRPProblem(self) -> VRPProblem:
        magazines_num = 1
        nodes_num = len(self.distance_matrix)
        vehicles = 1
        capacities = np.ones((vehicles), dtype=int)
        sources = [i for i in range(magazines_num)]
        dests = [i for i in range(magazines_num, nodes_num)]
        weights = np.zeros((nodes_num), dtype=int)
        return VRPProblem(sources, self.distance_matrix,
                          capacities, dests, weights)

    def solveTSP(self) -> Tuple[List[int], int]:
        only_one_const = 10000000.
        order_const = 1.
        problem = self.__createVRPProblem()
        solver = SolutionPartitioningSolver(
            problem, DBScanSolver(problem, anti_noiser=True))
        solution = solver.solve(
            only_one_const, order_const, solver_type=self.solver_type)

        if solution == None or solution.check() == False:
            print("Solver hasn't find solution.\n")
            return [], -1

        path = solution.solution[0]
        return path, super().calculate_cost(path)

