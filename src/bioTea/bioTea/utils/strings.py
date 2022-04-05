from itertools import zip_longest
from colorama import Fore
import re

# 7-bit C1 ANSI sequences
ansi_escape = re.compile(
    r"""
    \x1B  # ESC
    (?:   # 7-bit C1 Fe (except CSI)
        [@-Z\\-_]
    |     # or [ for CSI, followed by a control sequence
        \[
        [0-?]*  # Parameter bytes
        [ -/]*  # Intermediate bytes
        [@-~]   # Final byte
    )
""",
    re.VERBOSE,
)


def strip_colors(string):
    return ansi_escape.sub("", string)


def combine_linewise(a, b, padding="", strip=False, align_bottom=False):
    lines_a = a.split("\n")
    lines_b = b.split("\n")

    if align_bottom:
        lines_a = list(reversed(lines_a))
        lines_b = list(reversed(lines_b))

    if len(lines_a) > len(lines_b):
        fill = " " * len(lines_a)
    elif len(lines_b) > len(lines_a):
        fill = " " * len(lines_b)
    else:
        fill = ""

    result = []
    for line_a, line_b in zip_longest(lines_a, lines_b, fillvalue=fill):
        if strip:
            line_a, line_b = line_a.strip(), line_b.strip()
        line_a = line_a + padding
        result.append(line_a + line_b)

    if align_bottom:
        result = reversed(result)

    return "\n".join(result)


def make_square(logo, side="left"):
    assert side in ["left", "right"]
    lines = logo.split("\n")
    print(lines)
    longest = max([len(strip_colors(x)) for x in lines])
    print(longest)

    if side == "left":
        res = [line.ljust(longest) for line in lines]
    elif side == "right":
        res = [line.rjust(longest) for line in lines]

    return "\n".join(res)


TEA = """                                                    __/\__
             ;,'                               . _  \\\\''//
     _o_    ;:;' __    _     _______________   -( )-/_||_\\
 ,-.'---`.__ ;  / /_  (_)___/_  __/ ____/   |   .'. \_()_/
((j`=====',-'  / __ \/ / __ \/ / / __/ / /| |    |   | . \\
 `-\     /    / /_/ / / /_/ / / / /___/ ___ |    ϕ---| .  \\
    `-=-'    /_.___/_/\____/_/ /_____/_/  |_|   .'. ,\_____'.
                              W I Z A R D
"""

TEAPOT = (Fore.RESET + "\n").join(
    [
        Fore.LIGHTBLACK_EX + "             ;,'",
        Fore.BLUE + "     _o_    " + Fore.LIGHTBLACK_EX + ";:;'",
        Fore.BLUE + " ,-.'---`.__ " + Fore.LIGHTBLACK_EX + ";",
        Fore.BLUE + "((j`=====',-'",
        Fore.BLUE + " `-\     /",
        Fore.BLUE + "    `-=-'",
    ]
)

BIOTEA = (Fore.RESET + "\n").join(
    [
        "    __    _     " + Fore.LIGHTGREEN_EX + "_______________ ",
        "   / /_  (_)___" + Fore.LIGHTGREEN_EX + "/_  __/ ____/   |",
        "  / __ \/ / __ \\" + Fore.LIGHTGREEN_EX + "/ / / __/ / /| |",
        " / /_/ / / /_/ " + Fore.LIGHTGREEN_EX + "/ / / /___/ ___ |",
        "/_.___/_/\____" + Fore.LIGHTGREEN_EX + "/_/ /_____/_/  |_|",
    ]
)

TEA_LOGO = combine_linewise(TEAPOT)
