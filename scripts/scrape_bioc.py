"""Scrape Bioconductor for annotation package names and versions

Give the output CSV path to dump the scraped data to
"""
from pathlib import Path
import bs4
from bs4 import BeautifulSoup
import urllib.request
import re


def parse_row(row: str) -> str:
    cell_re = re.compile('<td>(.*?)<\/td>')
    strip_a_href = re.compile('<a href=".*?">(.*?)<\/a>')

    matches = cell_re.findall(row)

    if not matches:
        return None

    if strip_a_href.findall(matches[0]):
        matches[0] = strip_a_href.findall(matches[0])[0]
    return tuple(matches)

def main(output_path):
    with urllib.request.urlopen("https://bioconductor.org/packages/release/data/annotation/") as file:
        data = file.read().decode("utf8")

    soup = BeautifulSoup(data, features="html5lib")

    table_rows = soup.find("table").find_all("tr")

    table_rows = set([parse_row(str(row)) for row in table_rows])
    table_rows = set([x for x in table_rows if x is not None])

    with Path(output_path).open("w+") as file:
        file.writelines("package_name,maintainer,info\n")
        for x in table_rows:
            file.writelines('"{}","{}","{}"\n'.format(*x))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument("file_out")

    args = parser.parse_args()

    main(
        output_path=args.file_out
    )
