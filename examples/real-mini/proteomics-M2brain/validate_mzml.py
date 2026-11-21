#!/usr/bin/env python3
"""Validate the indexed M2brain DDA example using the Python standard library."""
import argparse
import base64
from collections import Counter
import hashlib
import math
from pathlib import Path
import re
import struct
from xml.etree import ElementTree as ET
from xml.parsers import expat
import zlib

NS = '{http://psi.hupo.org/ms/mzml}'
NATIVE = re.compile(r'controllerType=0 controllerNumber=1 scan=([1-9][0-9]*)')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def cv(node, accession):
    return [x for x in node.iter(NS+'cvParam') if x.get('accession') == accession]


def one_cv(node, accession):
    fields = cv(node, accession)
    require(len(fields) == 1, 'Missing or duplicate required PSI-MS field: ' + accession)
    return fields[0]


def peak_arrays(spectrum):
    length = int(spectrum.get('defaultArrayLength'))
    require(length > 0, 'Empty spectrum')
    arrays = {}
    for array in spectrum.iter(NS+'binaryDataArray'):
        kinds = [a for a in ('MS:1000514', 'MS:1000515') if cv(array, a)]
        require(len(kinds) == 1 and kinds[0] not in arrays, 'Unrecognized or duplicate peak array')
        one_cv(array, 'MS:1000523')
        one_cv(array, 'MS:1000574')
        encoded = array.find(NS+'binary').text
        require(len(encoded) == int(array.get('encodedLength')), 'Encoded binary array length mismatch')
        data = zlib.decompress(base64.b64decode(encoded, validate=True))
        require(len(data) == 8*length, 'Decoded binary array length mismatch')
        require(all(math.isfinite(p[0]) for p in struct.iter_unpack('<d', data)), 'Nonfinite peak value')
        arrays[kinds[0]] = data
    require(set(arrays) == {'MS:1000514','MS:1000515'}, 'Required m/z and intensity arrays not both present')
    return arrays


def check_index(raw, root, spectra):
    index = root.find(NS+'indexList')
    require(index.get('count') == '1' and len(index) == 1 and index[0].get('name') == 'spectrum',
            'Unexpected indexed mzML index list')
    offsets = [(entry.get('idRef'), int(entry.text)) for entry in index[0]]
    require([ident for ident,_ in offsets] == [x.get('id') for x in spectra], 'Spectrum index IDs differ')
    require(len({pos for _,pos in offsets}) == len(spectra), 'Repeated index positions')
    offsets_map = {}
    parser = expat.ParserCreate(namespace_separator='}')
    def record(name, attrs):
        tag = name.rsplit('}',1)[-1]
        if tag == 'spectrum':
            offsets_map[attrs['id']] = parser.CurrentByteIndex
        elif tag == 'indexList':
            offsets_map['_index'] = parser.CurrentByteIndex
    parser.StartElementHandler = record
    parser.Parse(raw, True)
    require(all(offsets_map.get(ident) == pos for ident,pos in offsets), 'Spectrum index byte offsets differ')
    require(offsets_map['_index'] == int(root.findtext(NS+'indexListOffset')), 'Index list byte offset differs')
    checksum = re.search(rb'<fileChecksum>([a-f0-9]{40})</fileChecksum>\s*</indexedmzML>\s*$', raw)
    require(checksum is not None and hashlib.sha1(raw[:checksum.start(1)]).hexdigest().encode() == checksum[1],
            'Internal indexed mzML SHA-1 does not match')
    return checksum[1].decode()


def validate(path):
    raw = path.read_bytes()
    require(b'<!DOCTYPE' not in raw and b'<!ENTITY' not in raw, 'Document type and external entities are unsupported')
    root = ET.fromstring(raw)
    require(root.tag == NS+'indexedmzML', 'Expected an indexed mzML file')
    mzml = root.find(NS+'mzML')
    require(mzml.get('version') == '1.1.0', 'Expected mzML 1.1.0')
    run = mzml.find(NS+'run')
    spectra = run.find(NS+'spectrumList')
    require(int(spectra.get('count')) == len(spectra) == 1196, 'Unexpected count of spectra')
    require(not list(run.iter(NS+'chromatogram')), 'Unexpected chromatogram entries')
    counts = Counter()
    seen = {}
    reference_count = 0
    last_time = -1.0
    for i,spectrum in enumerate(spectra):
        native = spectrum.get('id')
        match = NATIVE.fullmatch(native)
        require(match is not None and int(spectrum.get('index')) == i and native not in seen,
                'Malformed or repeated native spectrum ID')
        level = int(one_cv(spectrum, 'MS:1000511').get('value'))
        require(level in (1,2), 'Unexpected spectrum MS level')
        rt = one_cv(spectrum, 'MS:1000016')
        unit = rt.get('unitAccession')
        require(unit in ('UO:0000031','UO:0000010'), 'Unexpected retention time unit')
        seconds = float(rt.get('value')) * (60 if unit=='UO:0000031' else 1)
        require(math.isfinite(seconds) and seconds >= last_time, 'Retention time is invalid or unordered')
        last_time = seconds
        precursors = list(spectrum.iter(NS+'precursor'))
        charge = None
        if level == 2:
            require(len(precursors) == 1 and precursors[0].get('spectrumRef') in seen
                    and seen[precursors[0].get('spectrumRef')]['level'] == 1, 'MS2 precursor is not a preceding MS1')
            charge = int(one_cv(precursors[0], 'MS:1000041').get('value'))
            mass = float(one_cv(precursors[0], 'MS:1000744').get('value'))
            require(charge > 0 and math.isfinite(mass) and mass > 0, 'Invalid precursor charge or m/z')
            reference_count += 1
        else:
            require(not precursors, 'Unexpected MS1 precursor')
        peaks = peak_arrays(spectrum)
        seen[native] = {'level':level,'scan':int(match[1]),'rt_s':seconds,'charge':charge,
                        'mzml_node':spectrum,'peaks':peaks}
        counts[level] += 1
    require(counts == {1:52,2:1144} and reference_count == 1144, 'Unexpected MS levels or precursor closure')
    checksum = check_index(raw,root,spectra)
    return {'raw':raw,'root':root,'mzml':mzml,'run':run,'spectra':seen,'counts':dict(counts),
            'precursor_references':reference_count,'index_sha1':checksum,'sha256':hashlib.sha256(raw).hexdigest()}


def validate_target_only(path):
    accession = r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})'
    header = re.compile(r'sp\|('+accession+r')\|[A-Z0-9_]+_MOUSE OX=10090 SV=[0-9]+')
    identifiers = set()
    sequence_count = 0
    has_sequence = False
    for row in path.read_text(encoding='ascii').splitlines():
        if row.startswith('>'):
            require(not identifiers or has_sequence, 'Empty FASTA record')
            match = header.fullmatch(row[1:])
            require(match is not None and match[1] not in identifiers, 'Invalid or repeated Swiss-Prot target header')
            identifiers.add(match[1])
            sequence_count += 1
            has_sequence = False
        else:
            require(bool(identifiers) and re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYXUBZ]+', row),
                    'Invalid FASTA sequence record')
            has_sequence = True
    require(has_sequence and sequence_count == 57, 'Expected 57 nonempty target-only FASTA records')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mzml', type=Path)
    parser.add_argument('--fasta', type=Path)
    parser.add_argument('--sha256', help='expected mzML SHA-256')
    opts = parser.parse_args()
    try:
        proof = validate(opts.mzml)
        if opts.sha256:
            require(proof['sha256'] == opts.sha256.lower(), 'mzML SHA-256 does not match')
        if opts.fasta:
            validate_target_only(opts.fasta)
    except (ValueError, TypeError, KeyError, AttributeError, OSError, ET.ParseError, expat.ExpatError, zlib.error) as error:
        # Avoid exposing paths or source metadata from failed inputs.
        if type(error) is ValueError and re.fullmatch(r'[A-Za-z0-9 :/.-]+', str(error)):
            parser.exit(1, 'Invalid example: '+str(error)+'\n')
        parser.exit(1, 'Invalid or unreadable example file.\n')
    print('PASS: 52 MS1, 1144 MS2, 1196 indexed spectra, 1144 precursor references, strict binary arrays and internal SHA-1.')
    if opts.fasta:
        print('PASS: 57 canonical target-only FASTA records.')


if __name__ == '__main__':
    main()
