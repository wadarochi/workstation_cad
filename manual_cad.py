# -*- coding: utf-8 -*-
import numpy as np

from scipy import linalg
from scipy.spatial.transform import Rotation as R


class WorkstationConfig(object):
    # Length
    TotalWidth = 60
    TotalHeight = 40
    BackDepth = 30
    FrontDepth = 20
    BottomHeight = 20
    TopHeight = 20

    # Degree
    SideTheta = 45
    TopTheta = 30

    # hole
    HoleX = 10
    HoleZ = 10
    HoleRadius = 3


def normalize_array(original_array: np.array) -> np.array:
    return original_array / np.sqrt((original_array ** 2).sum())


def rotate_around_vector(point: np.array, axes: np.array, theta: float) -> np.array:
    unit_axes = normalize_array(axes)


class Plane(object):

    def __init__(self, a: float, b: float, c: float, d: float):
        """
        A plane is defined by the equation:
        ax+by+cz=d
        """
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def rotate_around_vector(self, axes: np.array, theta: float, point: np.array) -> '__class__':
        normal = np.array([self.a, self.b, self.c])
        normalized_axes = normalize_array(axes)
        r = R.from_rotvec(theta / 180.0 * np.pi * normalized_axes)
        new_normal = r.apply(normal)
        new_d = np.dot(new_normal, point)
        return Plane(new_normal[0], new_normal[1], new_normal[2], new_d)

    def get_point_with_xz(self, x: float, z: float) -> np.array:
        return np.array([x, (self.d - self.a * x - self.c * z) / self.b, z])


class Polygon(object):

    def __init__(self, tag: str):
        self.points = []
        self.tag = tag

    def add_new_point(self, point: np.array) -> None:
        self.points.append(point)

    def add_existed_point(self, point: np.array) -> None:
        self.points.append(np.array(point))

    def manim(self) -> None:
        manim_code = "(["
        for p in self.points:
            manim_code += "({}, {}, {}), ".format(*list(p / 10.0))
        return manim_code + "], PINK),"

    def print(self) -> None:
        indent = " " * 12
        print(indent + "# " + self.tag)
        print(indent + self.manim())


def find_conjuction_point(plane_a: Plane, plane_b: Plane, plane_c: Plane) -> np.array:
    a = np.array(
        [
            [plane_a.a, plane_a.b, plane_a.c],
            [plane_b.a, plane_b.b, plane_b.c],
            [plane_c.a, plane_c.b, plane_c.c]
        ]
    )
    b = np.array([plane_a.d, plane_b.d, plane_c.d])
    return linalg.solve(a, b)


class Workstation(object):

    def __init__(
        self, bottom_board: Polygon, top_board: Polygon, left_board: Polygon,
        right_board: Polygon, front_left_board: Polygon, front_board: Polygon,
        front_right_board: Polygon, front_top_board: Polygon, left_hole_center: np.array,
        right_hole_center: np.array
    ):
        self.bottom_board = bottom_board
        self.top_board = top_board
        self.left_board = left_board
        self.right_board = right_board
        self.front_left_board = front_left_board
        self.front_board = front_board
        self.front_right_board = front_right_board
        self.front_top_board = front_top_board
        self.left_hole_center = left_hole_center
        self.right_hole_center = right_hole_center

    def print(self):
        self.bottom_board.print()
        self.top_board.print()

        self.left_board.print()
        self.right_board.print()

        self.front_left_board.print()
        self.front_board.print()
        self.front_right_board.print()

        self.front_top_board.print()

        indent = " " * 12
        l = self.left_hole_center / 10.0
        r = self.right_hole_center / 10.0
        print(indent + "# left_hole")
        print(indent + f"([{l[0]}, {l[1]}, {l[2]}], {WorkstationConfig.HoleRadius / 10.0}, PINK),")
        print(indent + "# right")
        print(indent + f"([{r[0]}, {r[1]}, {r[2]}], {WorkstationConfig.HoleRadius / 10.0}, PINK),")


def solve_all_polygon():
    back = Plane(1, 0, 0, -WorkstationConfig.BackDepth)
    bottom = Plane(0, 0, 1, 0)
    left = Plane(0, 1, 0, -WorkstationConfig.TotalWidth / 2)
    right = Plane(0, 1, 0, WorkstationConfig.TotalWidth / 2)
    front = Plane(1, 0, 0, WorkstationConfig.FrontDepth)
    top = Plane(0, 0, 1, WorkstationConfig.TopHeight)

    back_board = Polygon("back")
    back_board.add_new_point(find_conjuction_point(back, left, bottom))
    back_board.add_new_point(find_conjuction_point(back, right, bottom))
    back_board.add_new_point(find_conjuction_point(back, right, top))
    back_board.add_new_point(find_conjuction_point(back, left, top))

    bottom_board = Polygon("bottom")
    bottom_board.add_existed_point(back_board.points[1])
    bottom_board.add_existed_point(back_board.points[0])
    # To be continue

    left_board = Polygon("left")
    left_board.add_existed_point(back_board.points[0])
    left_board.add_existed_point(back_board.points[3])
    left_board.add_existed_point(left_board.points[1])
    left_board.points[2][0] += WorkstationConfig.BackDepth
    left_board.add_existed_point(np.array(left_board.points[2]))
    left_board.points[3][2] = left_board.points[0][2]

    right_board = Polygon("right")
    right_board.add_existed_point(back_board.points[1])
    right_board.add_existed_point(back_board.points[2])
    right_board.add_existed_point(right_board.points[1])
    right_board.points[2][0] += WorkstationConfig.BackDepth
    right_board.add_existed_point(right_board.points[2])
    right_board.points[3][2] = right_board.points[0][2]

    top_board = Polygon("top")
    top_board.add_existed_point(back_board.points[2])
    top_board.add_existed_point(back_board.points[3])
    top_board.add_existed_point(left_board.points[2])
    top_board.add_existed_point(right_board.points[2])

    front_left = left.rotate_around_vector(left_board.points[2] - left_board.points[3], WorkstationConfig.SideTheta, left_board.points[2])
    front_right = right.rotate_around_vector(right_board.points[3] - right_board.points[2], WorkstationConfig.SideTheta, right_board.points[3])

    bottom_board.add_existed_point(left_board.points[3])
    left_bottom_of_front_board = find_conjuction_point(bottom, front_left, front)
    bottom_board.add_new_point(left_bottom_of_front_board)
    right_bottom_of_front_board = find_conjuction_point(bottom, front, front_right)
    bottom_board.add_new_point(right_bottom_of_front_board)
    bottom_board.add_existed_point(right_board.points[3])
    # bottom_board Done

    front_top = top.rotate_around_vector(right_board.points[2] - left_board.points[2], WorkstationConfig.TopTheta, right_board.points[2])

    front_left_board = Polygon("front_left")
    front_left_board.add_existed_point(left_board.points[3])
    front_left_board.add_existed_point(left_board.points[2])
    front_left_board.add_new_point(find_conjuction_point(front_top, front_left, front))
    front_left_board.add_existed_point(left_bottom_of_front_board)

    front_right_board = Polygon("front_right")
    front_right_board.add_existed_point(right_board.points[3])
    front_right_board.add_existed_point(right_board.points[2])
    front_right_board.add_new_point(find_conjuction_point(front_top, front_right, front))
    front_right_board.add_existed_point(right_bottom_of_front_board)

    front_top_board = Polygon("front_top")
    front_top_board.add_existed_point(top_board.points[2])
    front_top_board.add_existed_point(top_board.points[3])
    front_top_board.add_existed_point(front_right_board.points[2])
    front_top_board.add_existed_point(front_left_board.points[2])

    front_board = Polygon("front")
    front_board.add_existed_point(left_bottom_of_front_board)
    front_board.add_existed_point(right_bottom_of_front_board)
    front_board.add_existed_point(front_right_board.points[2])
    front_board.add_existed_point(front_left_board.points[2])

    # Two holes
    left_hole_center = front_left.get_point_with_xz(WorkstationConfig.HoleX, WorkstationConfig.HoleZ)
    right_hole_center = front_right.get_point_with_xz(WorkstationConfig.HoleX, WorkstationConfig.HoleZ)
    # Done all

    return Workstation(
        bottom_board, top_board, left_board, right_board, front_left_board,
        front_board, front_right_board, front_top_board, left_hole_center,
        right_hole_center
    )


def main():
    return solve_all_polygon()


if __name__ == "__main__":
    w = main()
    w.print()
