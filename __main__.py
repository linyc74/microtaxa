import argparse
import microtaxa


__VERSION__ = '1.0.1-beta'


PROG = 'python microtaxa'
DESCRIPTION = f'Microbial taxonomy analysis by reference sequence alignment (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['-r', '--ref-fa'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the reference fasta file, e.g. SILVA database',
        }
    },
    {
        'keys': ['-s', '--sample-sheet'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the sample sheet (CSV, TSV, or XLSX format), "Sample" and "Group" cloumns are required',
        }
    },
    {
        'keys': ['-f', '--fq-dir'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the directory containing all input fastq files',
        }
    },
    {
        'keys': ['-1', '--fq1-suffix'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'suffix of read 1 fastq files',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-2', '--fq2-suffix'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'suffix of read 2 fastq files, "None" for single end reads (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'qiime2_pipeline_outdir',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['-i', '--min-percent-identity'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.97,
            'help': 'minimum percent identity (range 0, 1) for sequence alignment (default: %(default)s)',
        }
    },
    {
        'keys': ['--clip-r1-5-prime'],
        'properties': {
            'type': int,
            'required': False,
            'default': 0,
            'help': 'hard clip <int> bp from 5\' end of read 1 (default: %(default)s)',
        }
    },
    {
        'keys': ['--clip-r2-5-prime'],
        'properties': {
            'type': int,
            'required': False,
            'default': 0,
            'help': 'hard clip <int> bp from 5\' end of read 2 (default: %(default)s)',
        }
    },
    {
        'keys': ['--colormap'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'Set1',
            'help': 'matplotlib colormap for plotting, or comma-separated color names, e.g. "darkred,lightgreen,skyblue" (default: %(default)s)',
        }
    },
    {
        'keys': ['--invert-colors'],
        'properties': {
            'action': 'store_true',
            'help': 'invert the order of colors',
        }
    },
    {
        'keys': ['--publication-figure'],
        'properties': {
            'action': 'store_true',
            'help': 'plot figures in the form and quality for paper publication',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __VERSION__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running MicroTaxa {__VERSION__}\n', flush=True)
        microtaxa.entrypoint(
            ref_fa=args.ref_fa,
            sample_sheet=args.sample_sheet,
            fq_dir=args.fq_dir,
            fq1_suffix=args.fq1_suffix,
            fq2_suffix=args.fq2_suffix,
            clip_r1_5_prime=args.clip_r1_5_prime,
            clip_r2_5_prime=args.clip_r2_5_prime,
            min_percent_identity=args.min_percent_identity,
            colormap=args.colormap,
            invert_colors=args.invert_colors,
            publication_figure=args.publication_figure,
            outdir=args.outdir,
            threads=args.threads,
            debug=args.debug)


if __name__ == '__main__':
    EntryPoint().main()
