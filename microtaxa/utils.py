import os
from typing import Optional, Tuple


class FastaParser:

    def __init__(self, file: str):
        self.__fasta = open(file, 'r')

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return

    def __iter__(self):
        self.__fasta.seek(0)
        return self

    def __next__(self):
        r = self.next()
        if r:
            return r
        else:  # r is None
            raise StopIteration

    def next(self) -> Optional[Tuple[str, str]]:
        """
        Returns the next read of the fasta file
        If it reaches the end of the file, return None
        """
        header = self.__fasta.readline().rstrip()[1:]
        if header == '':
            return None

        seq_lines = []
        while True:
            pos = self.__fasta.tell()
            line = self.__fasta.readline().rstrip()

            if line.startswith('>'):
                self.__fasta.seek(pos)
                return header, ''.join(seq_lines)

            end_of_file = line == ''
            if end_of_file:
                return header, ''.join(seq_lines)

            seq_lines.append(line)

    def close(self):
        self.__fasta.close()


def get_temp_path(
        prefix: str = 'temp',
        suffix: str = '') -> str:

    i = 1
    while True:
        fpath = f'{prefix}{i:03}{suffix}'
        if not os.path.exists(fpath):
            return fpath
        i += 1
