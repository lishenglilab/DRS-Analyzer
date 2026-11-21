#!/usr/bin/env python3
"""Check input structure and exact pairing for the F1-1 Real-mini example."""
import argparse
import gzip
from pathlib import Path
from uuid import UUID

import h5py


def read_fastq(path, count):
    records = {}
    with gzip.open(path, 'rb') as stream:
        while True:
            header = stream.readline().rstrip(b'\r\n')
            if not header:
                if stream.peek(1):
                    raise ValueError('empty FASTQ header')
                break
            sequence, plus, quality = (stream.readline().rstrip(b'\r\n') for _ in range(3))
            if not header.startswith(b'@') or plus != b'+' or not sequence or len(sequence) != len(quality):
                raise ValueError('invalid FASTQ record')
            uid = header[1:].decode('ascii')
            if str(UUID(uid)) != uid or uid in records:
                raise ValueError('noncanonical or duplicate FASTQ read ID')
            records[uid] = (sequence, quality)
    if len(records) != count:
        raise ValueError('unexpected FASTQ record count')
    return records


def verify_paired_fast5(path, fastq):
    basecall = 'Analyses/Basecall_1D_000/BaseCalled_template/Fastq'
    raw_count = 0
    with h5py.File(path, 'r') as doc:
        if set(doc) != {'read_' + uid for uid in fastq}:
            raise ValueError('FASTQ and FAST5 read ID sets differ')
        for uid, (sequence, quality) in fastq.items():
            group = doc['read_' + uid]
            raw = group['Raw']
            signal = raw['Signal']
            if raw.attrs['read_id'] != uid.encode('ascii'):
                raise ValueError('FAST5 read index differs from Raw read ID')
            if (signal.ndim != 1 or signal.dtype.kind != 'i' or signal.dtype.itemsize != 2
                    or not signal.size or len(signal[:]) != int(raw.attrs['duration'])):
                raise ValueError('invalid FAST5 raw signal or duration')
            sample_rate = group['channel_id'].attrs['sampling_rate']
            if not 0 < sample_rate < float('inf'):
                raise ValueError('invalid sampling rate')
            if group[basecall][()].splitlines() != [b'@' + uid.encode('ascii'), sequence, b'+', quality]:
                raise ValueError('embedded FASTQ differs from paired FASTQ')
            raw_count += signal.size
    return raw_count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('example_dir', type=Path)
    parser.add_argument('--include-qc', action='store_true', help='also check the optional independent QC FASTQ')
    args = parser.parse_args()
    example = args.example_dir
    try:
        fastq = read_fastq(example / 'F1-1-paired-mini.fastq.gz', 12)
        raw_count = verify_paired_fast5(example / 'fast5/F1-1-paired-mini.fast5', fastq)
        if args.include_qc:
            qc = read_fastq(example / 'F1-1-QC-mini.fastq.gz', 1500)
            if qc.keys() & fastq.keys():
                raise ValueError('QC and paired FASTQ read ID sets intersect')
    except (OSError, KeyError, TypeError, ValueError, IndexError) as error:
        message = str(error)
        safe_messages = {
            'empty FASTQ header', 'invalid FASTQ record', 'noncanonical or duplicate FASTQ read ID',
            'unexpected FASTQ record count', 'FASTQ and FAST5 read ID sets differ',
            'FAST5 read index differs from Raw read ID', 'invalid FAST5 raw signal or duration',
            'invalid sampling rate', 'embedded FASTQ differs from paired FASTQ',
            'QC and paired FASTQ read ID sets intersect',
        }
        parser.exit(1, 'Structure check failed: ' + (message if message in safe_messages else 'unreadable or malformed input') + '\n')
    print(f'PASS: {len(fastq)} exact FASTQ/FAST5 pairs; {raw_count} int16 signal samples readable.')
    if args.include_qc:
        print(f'PASS: {len(qc)} independent QC records; {len(qc.keys() & fastq.keys())} shared IDs.')


if __name__ == '__main__':
    main()
