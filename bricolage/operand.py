from enum import IntEnum

# Taken from the table here
# http://mathworld.wolfram.com/BooleanFunction.html

# A  B | 0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
# 0  0 | 0  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1
# 0  1 | 0  0  0  0  1  1  1  1  0  0  0  0  1  1  1  1
# 1  0 | 0  0  1  1  0  0  1  1  0  0  1  1  0  0  1  1
# 1  1 | 0  1  0  1  0  1  0  1  0  1  0  1  0  1  0  1
#
# 4 bits can capture all Binary Boolean functions
class Operand(IntEnum):
    FALSE = 0
    AND = 1
    A_AND_NOT_B = 2
    A = 3
    NOT_A_AND_B = 4
    B = 5
    XOR = 6
    OR = 7
    NOR = 8
    XNOR = 9
    NOT_B = 10
    A_OR_NOT_B = 11
    NOT_A = 12
    NOT_A_OR_B = 13
    NAND = 14
    TRUE = 15


def calculate(op, a, b):
    return op & (8 >> ((a << 1) | b)) != 0


def calc_all(op):
    return [calculate(op, a, b) for a, b in ((0, 0), (0, 1), (1, 0), (1, 1))]


# TODO: Sort out short descriptions
# _short_desc =
# {
#     FALSE = 0
#     AND = 1
#     A_AND_NOT_B = 2
#     A = 3
#     NOT_A_AND_B = 4
#     B = 5
#     XOR = 6
#     OR = 7
#     NOR = 8
#     XNOR = 9
#     NOT_B = 10
#     A_OR_NOT_B = 11
#     NOT_A = 12
#     NOT_A_OR_B = 13
#     NAND = 14
#     TRUE = 15
# }
