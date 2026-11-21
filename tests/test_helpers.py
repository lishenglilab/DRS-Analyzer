"""Run the real QC processes and helpers against small synthetic inputs."""

import csv
import gzip
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile


REPO = Path(__file__).resolve().parent.parent
TESTS = REPO / "tests"
DATA = TESTS / "data"
QC_COLUMNS = [
    "sample_id", "Size_GB", "Total_reads", "pass_percent",
    "Median_pass_read_length", "Mean_pass_read_length", "Mapped_read", "Mapped_ratio",
]
M5C_HEADER = ["transcript", "position", "base", "fraction"]


def equal(actual, expected):
    if actual != expected:
        raise AssertionError(f"expected {expected!r}, observed {actual!r}")


def invoke(directory, command, log="console.log"):
    args = [str(value) for value in command]
    result = subprocess.run(
        args, cwd=directory, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        encoding="utf-8", errors="replace", check=False,
        env={**os.environ, "LC_ALL": "C.UTF-8", "NXF_DISABLE_CHECK_LATEST": "true"},
    )
    with (directory / log).open("a", encoding="utf-8") as stream:
        stream.write("$ " + " ".join(shlex.quote(arg) for arg in args) + "\n")
        stream.write(result.stdout)
        stream.write(f"\n[exit {result.returncode}]\n")
    return result


def success(result):
    if result.returncode != 0:
        raise AssertionError(f"expected exit 0, observed exit {result.returncode}:\n{result.stdout[-1800:]}")


def rejected(result, pattern):
    if result.returncode == 0 or not re.search(pattern, result.stdout, re.I):
        raise AssertionError(f"expected nonzero exit with {pattern!r}, observed exit {result.returncode}:\n"
                             f"{result.stdout[-1800:]}")
    return f"exit {result.returncode}, matched {pattern}"


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.reader(stream, delimiter="\t"))


def qc_table(path, expected):
    rows = read_tsv(path)
    equal(len(rows), 2)
    equal(rows[0], QC_COLUMNS)
    equal(len(rows[1]), len(QC_COLUMNS))
    equal(rows[1][0], "toy")
    actual = [float(value) for value in rows[1][1:]]
    for name, observed, wanted in zip(QC_COLUMNS[1:], actual, expected):
        if not math.isclose(observed, wanted, rel_tol=0, abs_tol=1e-8):
            raise AssertionError(f"{name}: expected {wanted}, observed {observed}")
    return f"{path.name}: expected toy,{expected}; observed toy,{actual}"


def build_fixtures(directory):
    fq = (DATA / "toy.fastq").read_bytes()
    lines = fq.splitlines()
    equal(len(lines), 12)
    for header, sequence, separator, quality in zip(lines[::4], lines[1::4], lines[2::4], lines[3::4]):
        equal((header[:1], separator[:1], len(sequence)), (b"@", b"+", len(quality)))
    fastq = directory / "toy.clean.fastq.gz"
    fastq.write_bytes(gzip.compress(fq, mtime=0))
    bam = directory / "toy.bam"
    success(invoke(directory, ["samtools", "view", "-b", "-o", bam, DATA / "toy.sam"]))
    return fastq, bam


def nextflow(case, script, *args, config=TESTS / "ci.config"):
    return invoke(case, [
        "nextflow", "-C", config, "-log", case / "nextflow.log",
        "run", script, "-ansi-log", "false", "-w", case / "work",
        "--outdir", case / "out", "--qc_dir", case / "qc", *args,
    ])


def emitted_file(case, name, expected_basename):
    rows = read_tsv(case / "out" / name)
    equal(len(rows), 1)
    equal(len(rows[0]), 2)
    sample, output = rows[0]
    path = Path(output)
    equal((sample, path.name), ("toy", expected_basename))
    return path


def seqkit(case, fastq):
    success(nextflow(case, TESTS / "qc_processes.nf", "-entry", "LENGTHS", "--fastq", fastq))
    path = emitted_file(case, "lengths.tsv", "toy_Length_pass.txt")
    expected = [["#id", "length"], ["read_a", "4"], ["read_b", "6"], ["read_c", "11"]]
    actual = read_tsv(path)
    equal(actual, expected)
    return f"lengths emit (toy, {path.name}); expected {expected}; observed {actual}"


def bad_fastq(case):
    malformed = case / "malformed.fastq.gz"
    malformed.write_bytes(gzip.compress(b"@bad_read\nACGT\n+\nII\n", mtime=0))
    result = nextflow(case, TESTS / "qc_processes.nf", "-entry", "LENGTHS", "--fastq", malformed)
    message = rejected(result, r"bad fastq format|unequal sequence and quality|invalid FASTQ")
    equal((case / "out" / "lengths.tsv").exists(), False)
    return f"expected rejection of mismatched sequence/quality and no emit; observed {message}"


def bam_flagstat(case, bam):
    success(nextflow(case, TESTS / "qc_processes.nf", "-entry", "FLAGSTAT", "--bam", bam))
    stat = emitted_file(case, "flagstat.tsv", "toy.map.stat").read_text(encoding="utf-8")
    total = re.findall(r"^(\d+) \+ \d+ in total\b", stat, re.M)
    mapped = re.findall(r"^(\d+) \+ \d+ mapped \(", stat, re.M)
    stats = read_tsv(case / "out" / "stats.tsv")
    equal((total, mapped, stats), (["3"], ["2"], [["toy", "2"]]))
    return f"expected total=3, mapped=2, stats=[toy,2]; observed total={total}, mapped={mapped}, stats={stats}"


def r_qc(case, fail_contents, expected):
    (case / "toy_Length_pass.txt").write_text("#id\tlength\nread_a\t4\nread_b\t6\nread_c\t11\n", encoding="utf-8")
    (case / "toy_mapping.txt").write_text("toy 2\n", encoding="utf-8")
    # This Size_GB is injected independently of the integration test's FASTQ.
    (case / "toy_size.txt").write_text("toy 1.25\n", encoding="utf-8")
    (case / "read_number_fail.txt").write_text(fail_contents, encoding="utf-8")
    success(invoke(case, ["Rscript", "--vanilla", REPO / "bin" / "combine_qc_stats.R"]))
    return qc_table(case / "QC.xls", expected)


def m5c(case, reference, plus, minus, expected):
    result = invoke(case, [
        "python3", REPO / "bin" / "mod" / "m5c_res_process.py", "--out_dir", case,
        "--ref_tr_fa", reference, "--plus_file", DATA / plus, "--minus_file", DATA / minus,
    ])
    report = case / "m5c.prediction.txt"
    if expected is None:
        message = rejected(result, r"KeyError: 'unknown'")
        equal(report.exists(), False)
        return f"expected unknown ID rejection and no TSV; observed {message}"
    success(result)
    actual = read_tsv(report)
    equal(actual, [M5C_HEADER, *expected])
    return f"expected {expected}; observed {actual[1:]}"


def quality_report(case, fastq, bam, layout, fail_mode, missing=None):
    bam_dir = case / "input_bam"
    bam_dir.mkdir()
    fastq_dir = case / "input_fastq"
    fastq_dir.mkdir()
    if missing != "bam":
        bam_target = bam_dir / "sample" / "align" if layout == "nested" else bam_dir
        bam_target.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(bam, bam_target / "toy.bam")
    if missing != "fastq":
        filename = "other.clean.fastq.gz" if missing == "mismatch" else fastq.name
        shutil.copyfile(fastq, fastq_dir / filename)
    if fail_mode != "absent":
        (case / "qc").mkdir()
        (case / "qc" / "read_number_fail.txt").write_text(
            "toy 4\n" if fail_mode == "nonempty" else "", encoding="utf-8",
        )

    result = nextflow(
        case, REPO / "main.nf", "--module", "quality_report",
        "--bam_dir", bam_dir, "--clean_fastq_dir", fastq_dir,
    )
    report = case / "out" / "table" / "QC.xls"
    if missing:
        pattern = {
            "bam": "No BAM files found",
            "fastq": "No FASTQ files found",
            "mismatch": r"BAM and FASTQ sample IDs do not match.*BAM only: \[toy\].*FASTQ only: \[other\]",
        }[missing]
        message = rejected(result, pattern)
        equal(report.exists(), False)
        return f"expected rejection and no QC table; observed {message}"

    success(result)
    total = 4 if fail_mode == "nonempty" else 3
    pass_percent = 75 if fail_mode == "nonempty" else 100
    return qc_table(report, [0, total, pass_percent, 6, 7, 2, 66.67])


def write_versions(directory):
    commands = [
        ["nextflow", "-version"],
        ["java", "-version"],
        ["seqkit", "version"],
        ["samtools", "--version"],
        ["Rscript", "--vanilla", "-e",
         'cat(R.version.string, "\n"); for (p in c("tidyverse", "data.table")) '
         'cat(p, as.character(packageVersion(p)), "\n")'],
        ["python3", "-c", 'import sys, Bio; print(sys.version); print("Biopython", Bio.__version__)'],
    ]
    for command in commands:
        success(invoke(directory, command, log="versions.log"))


def main():
    output = Path(tempfile.mkdtemp(prefix="drs-tests-")).resolve()
    print(f"Test artifacts: {output}", flush=True)
    write_versions(output)
    fixtures = output / "fixtures"
    fixtures.mkdir()
    fastq, bam = build_fixtures(fixtures)
    ref = DATA / "toy.fa"
    with_collision = fixtures / "with-collision.fa"
    with_collision.write_text(ref.read_text(encoding="utf-8") + ">edge\nTTTTTT\n", encoding="utf-8")
    normal = [["toy_F", "2", "C", "0.25"], ["toy_F", "6", "C", "0.75"], ["toy_R", "3", "G", "0.5"]]
    boundaries = [["Fedge_F_F", "2", "C", "0.25"], ["Redge_R_R", "3", "G", "0.5"]]

    cases = [
        ("T1_SEQKIT_01", seqkit, fastq),
        ("T1_SEQKIT_02", bad_fastq),
        ("T1_BAM_01", bam_flagstat, bam),
        ("T1_R_01", r_qc, "toy 4\n", [1.25, 4, 75, 6, 7, 2, 66.67]),
        ("T1_R_02", r_qc, "", [1.25, 3, 100, 6, 7, 2, 66.67]),
        ("T1_M5C_01", m5c, ref, "plus.wig", "minus.wig", normal),
        ("T1_M5C_02", m5c, ref, "plus-edge.wig", "minus-edge.wig", boundaries),
        ("T1_M5C_03", m5c, ref, "plus-unknown.wig", "minus.wig", None),
        ("T1_M5C_04", m5c, with_collision, "plus-edge.wig", "minus-edge.wig", boundaries),
        ("T1_CLI_01", quality_report, fastq, bam, "flat", "absent"),
        ("T1_CLI_02", quality_report, fastq, bam, "nested", "empty"),
        ("T1_CLI_03", quality_report, fastq, bam, "nested", "nonempty"),
        ("T1_CLI_04", quality_report, fastq, bam, "flat", "absent", "bam"),
        ("T1_CLI_05", quality_report, fastq, bam, "flat", "absent", "fastq"),
        ("T1_CLI_06", quality_report, fastq, bam, "nested", "absent", "mismatch"),
    ]

    failures = 0
    with (output / "results.log").open("w", encoding="utf-8") as summary:
        for name, test, *args in cases:
            directory = output / name
            directory.mkdir()
            try:
                details = test(directory, *args)
                record = f"PASS {name} | {details}"
            except (AssertionError, OSError, ValueError) as error:
                failures += 1
                record = f"FAIL {name} | {error}"
            line = f"{record}\n  {directory}"
            print(line, flush=True)
            summary.write(line + "\n")
            summary.flush()

        totals = f"{len(cases) - failures}/{len(cases)} passed, {failures} failed"
        print(totals, flush=True)
        summary.write(totals + "\n")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
