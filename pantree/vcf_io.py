import gzip
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

from Bio import bgzf


def is_compressed_vcf_path(path: str | Path) -> bool:
    return str(path).endswith(".gz")


@contextmanager
def open_vcf_text_for_read(path: str | Path) -> Iterator[TextIO]:
    if is_compressed_vcf_path(path):
        with gzip.open(path, "rt") as handle:
            yield handle
    else:
        with open(path, "r") as handle:
            yield handle


@contextmanager
def open_vcf_text_for_write(path: str | Path) -> Iterator[TextIO]:
    if is_compressed_vcf_path(path):
        with bgzf.open(path, "wt") as handle:
            yield handle
    else:
        with open(path, "w") as handle:
            yield handle
