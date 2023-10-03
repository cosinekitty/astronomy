#!/usr/bin/env python3
#
#   goplan.py  -  Don Cross  -  2023-10-03
#
#   Using the dependencies of function calls in the C# version of Astronomy Engine,
#   make a plan for the order we should complete functions in the Go version.
#
import sys
from typing import Set, Dict, List, Optional


class Symbol:
    def __init__(self, name: str) -> None:
        self.name = name
        self.children: Set[str] = set()
        self.visited = False

    def insert(self, child: str) -> None:
        self.children.add(child)

    def __str__(self) -> str:
        text = self.name
        for child in sorted(self.children):
            text += ',' + child
        return text


def LoadDependencies(inFileName: str) -> Dict[str, Symbol]:
    with open(inFileName, 'rt') as infile:
        all: Dict[str, Symbol] = {}
        parent: Optional[Symbol] = None
        lnum = 0
        for line in infile:
            lnum += 1
            line = line.rstrip()
            name = line.lstrip()
            if name == '':
                continue    # ignore blank lines
            if name == line:
                # There is no whitespace to the left of the name.
                # This name is a parent (caller).
                parent = Symbol(name)
                if name in all:
                    print('FATAL({} line {}): parent {} was already defined.'.format(inFileName, lnum, name))
                    sys.exit(1)
                all[name] = parent
            else:
                # Whitespace exists to the left of the name, so this is a child (called function).
                if parent is None:
                    print('FATAL({} line {}): child {} has no parent.'.format(inFileName, lnum, name))
                    sys.exit(1)
                parent.insert(name)
        # Reduce verbosity in the input file by allowing implicit creation
        # of a symbol that has no children.
        undef: Set[str] = set()
        for parent in all.values():
            for child in parent.children:
                if child not in all:
                    undef.add(child)
        for name in undef:
            all[name] = Symbol(name)
        return all


def MakePlan(inFileName: str) -> int:
    all = LoadDependencies(inFileName)

    # Sweep through all symbols repeatedly.
    level = 0
    while True:
        # On each pass, print all unvisited parents that contain no unvisited children.
        # Mark each parent visited each time we print a parent.
        # Keep going until we fail to print something.

        namelist: List[str] = []

        for parent in all.values():
            if not parent.visited:
                blockers = 0
                for child in parent.children:
                    if not all[child].visited:
                        blockers += 1
                if blockers == 0:
                    namelist.append(parent.name)
        if len(namelist) == 0:
            break
        for name in sorted(namelist):
            parent = all[name]
            parent.visited = True
            print('{},{}'.format(level, parent))
        level += 1

    # Report an error if any unvisited symbols remain.
    undefCount = 0
    for symbol in all.values():
        if not symbol.visited:
            print('UNVISITED: {}'.format(symbol.name))
            undefCount += 1
    if undefCount > 0:
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(MakePlan('dependencies.txt'))
