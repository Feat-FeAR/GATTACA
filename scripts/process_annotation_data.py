"""Process a file with raw annotations from R to a processed file

The input file MUST have the following columns in this order:
  ID, SYMBOL, GENENAME, ENSEMBL, package_name, version
And adhere to the .csv standard.

The script performs the following steps:
    - Read in the data;
    - Find eventual collisions in the data;
    - Remove duplicate entries from the same packages (if any);
    - Fuse together entries for the same probe from different packages;
    - Write out the result.

This script uses the following external packages:
    - tqdm
    - jellyfish
    -
"""

from copy import copy
from itertools import tee
from pathlib import Path
from functools import cached_property
import csv

from tqdm import tqdm
from jellyfish import jaro_distance


def window(iterable, size):
    iters = tee(iterable, size)
    for i in range(1, size):
        for each in iters[i:]:
            next(each, None)
    return zip(*iters)


class Collision:
    def __init__(self, probe, gene1, gene2) -> None:
        self.probe = probe
        self.gene1 = gene1
        self.gene2 = gene2

    @cached_property
    def distance(self):
        return jaro_distance(self.gene1, self.gene2)

    def __str__(self) -> str:
        return f"Collision: {self.probe} -- {self.gene1}, {self.gene2} -- dist: {self.distance}"

    def __eq__(self, __o: object) -> bool:
        if type(__o) is not Collision:
            raise NotImplementedError()

        return self.probe == __o.probe and (
            all(elem in [self.gene1, self.gene2] for elem in [__o.gene1, __o.gene2])
        )

    def has_loc(self) -> bool:
        return self.gene1.startswith("LOC") or self.gene2.startswith("LOC")


# ID,SYMBOL,GENENAME,ENSEMBL,package_name,version
class Annotation:
    def __init__(self, pieces: str) -> None:
        pieces = [x if x != "NA" else None for x in pieces]

        self.id = pieces[0]
        self.symbol = pieces[1]
        self.description = pieces[2]
        self.ensembl_ids = pieces[3]
        if self.ensembl_ids:
            self.ensembl_ids = self.ensembl_ids.split(" /// ")
        else:
            self.ensembl_ids = [None]
        self.packages = {pieces[4]: pieces[5]}

    def __str__(self) -> str:

        if all([x is not None for x in self.ensembl_ids]):
            ensembl = " /// ".join(self.ensembl_ids)
        else:
            ensembl = "NA"
        return '{},{},{},{},{},{}'.format(
            f'"{self.id}"',
            f'"{self.symbol or "NA"}"',
            f'"{self.description or "NA"}"',
            f'"{ensembl}"',
            f'"{" /// ".join(self.packages.keys())}"',
            f'"{" /// ".join(self.packages.values())}"',
        )

    def add_to_packages(self, other) -> None:
        if type(other) != type(self):
            raise NotImplementedError(f"Cannot compare {type(self)} with {type(other)}")
        length_before = len(self.packages)
        old = self.packages
        self.packages.update(other.packages)
        if length_before == len(self.packages):
            print(old)
            print(other.packages)
            raise ValueError("The update updated some identical keys...???")

    def __eq__(self, __o: object) -> bool:
        if type(__o) != type(self):
            raise NotImplementedError(f"Cannot compare {type(self)} with {type(__o)}")

        return (self.id == __o.id and self.symbol == __o.symbol)

    @property
    def is_leaky(self):
        return self.symbol is None and self.description is None and \
            all([x is None for x in self.ensembl_ids])

    def identical_to(self, __o: object) -> bool:
        if type(__o) != type(self):
            raise NotImplementedError(f"Cannot compare {type(self)} with {type(__o)}")

        return self.__dict__ == __o.__dict__

    def merge_with(self, __o: object) -> bool:
        if type(__o) != type(self):
            raise NotImplementedError(f"Cannot compare {type(self)} with {type(__o)}")

        self.add_to_packages(__o)
        self.ensembl_ids.extend(__o.ensembl_ids)


def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)


def sort_file(file_in, file_out) -> None:
    with Path(file_in).open("r") as file_obj_in, Path(file_out).open("w+") as file_obj_out:
        lines = csv.reader(file_obj_in)

        print("Sorting lines...")
        lines = list(lines)
        lines.sort(key = lambda x: [1])

        writer = csv.writer(file_obj_out, quotechar='"')
        writer.writerows(lines)

    return None


def find_collisions(file_in, file_out) -> None:

    collisions = []

    with Path(file_in).open("r") as file_obj_in:
        next(file_obj_in) # Skip the csv header

        for row1, row2 in tqdm(window(csv.reader(file_obj_in), 2)):
            annot1 = Annotation(row1)
            annot2 = Annotation(row2)

            if annot1 == annot2:

                if annot1.is_leaky ^ annot2.is_leaky:
                    print(f"Found a leaky vs nonleaky annotation: {annot1} vs {annot2}")

                collisions.append(Collision(annot1.id, annot1.symbol, annot2.symbol))

    if len(collisions) == 0:
        print("No collisions found.")
        return None

    print(f"Found {len(collisions)} collisions.")

    print("Cleaning up collisions...")
    collisions.sort(key = lambda x: x.probe)
    clean = copy(collisions)
    cleanup = True
    done_cleanup = False
    cycle = 0

    while cleanup:
        print("Cleanup -- cycle {cycle}...")
        for col1, col2 in tqdm(pairwise(collisions)):
            if col1 == col2:
                clean.append(col1)
                done_cleanup = True
            else:
                clean.extend([col1, col2])

        if not done_cleanup:
            cleanup = False
        else:
            done_cleanup = False

        collisions = copy(clean)
        collisions.sort(key = lambda x: x.probe)
        clean = []
        cycle += 1

    print("Removing LOCs...")
    len_prior = len(collisions)
    collisions = [x for x in collisions if not x.has_loc()]
    len_after = len(collisions)
    print(f"Removed {len_prior - len_after} collisions (bef: {len_prior}, aft: {len_after}).")

    print("Sorting based on dist...")
    collisions.sort(key = lambda x: x.distance)

    with Path(file_out).open("w+") as file:
        file.writelines([f"{str(x)}\n" for x in collisions])

    print(f"Found {len(collisions)} collisions.")

    return True


def remove_duplicates(file_in, file_out):
    with Path(file_in).open("r") as file, Path(file_out).open("w+") as fileout:
        skipped = 0
        merged = 0
        fileout.write(f"{next(file)}\n") # Skip the header

        annotations = [Annotation(x) for x in csv.reader(file)]
        annotations.sort(key=lambda x: x.id)

        annot_old = annotations[0]
        for annot in tqdm(annotations[1:]):
            written = False

            if annot_old.identical_to(annot):
                # This line is the same as the previous line. Skip it.
                skipped += 1
                continue


            if annot.id == annot_old.id:
                # This line is the same as the previous line. Merge them.
                merged += 1
                annot_old.merge_with(annot)
                continue

            fileout.write(f"{annot_old}\n")
            written = True
            annot_old = annot

        if written:
            fileout.write(f"{annot}\n")

    skip_perc = round(skipped / len(annotations) * 100, 4)
    merged_perc = round(merged / len(annotations) * 100, 4)
    print(f"Read {len(annotations)} lines, merged {merged} lines ({merged_perc}%) and skipped {skipped} lines ({skip_perc}%)")


def main(input_path, output_path = None):
    if not output_path:
        output_path = f"{input_path}_processed.csv"

    print("Sorting file...")
    sort_file(input_path, Path(input_path).parent / "sortedinput.csv")
    print("Searching for collisions...")
    collisions = find_collisions(
        Path(input_path).parent / "sortedinput.csv",
        Path(input_path).parent / "collisions.txt"
    )
    if collisions:
        print("Found collisions. Stopping.")

    print("Removing duplicates...")
    remove_duplicates(
        Path(input_path).parent / "sortedinput.csv",
        output_path
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument("file_in", help="Input file with annotations")
    parser.add_argument(
        "file_out",
        nargs="?",
        help="Optional output file (if unspecified, creates a file alongside the input file)."
    )

    args = parser.parse_args()

    main(
        input_path=args.file_in,
        output_path=args.file_out
    )
