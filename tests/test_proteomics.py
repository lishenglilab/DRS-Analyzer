"""Exercise the Nextflow FragPipe boundary with local tools and synthetic inputs."""
import argparse
from pathlib import Path
import tempfile

from test_helpers import DATA, REPO, equal, nextflow, read_tsv, rejected, success


PSM = "Protein\tIntensity\nMOCK_TARGET\t42\n"

PHILOSOPHER = """#!/usr/bin/env python3
from pathlib import Path
import sys

args = sys.argv[1:]
if args[0] == "database":
    mode = "@MODE@"
    if mode != "missing":
        target = Path(args[-1]).read_text()
        decoy = ">rev_mock\\nPEPTIDE\\n" if mode != "bad_prefix" else ">decoy_mock\\nPEPTIDE\\n"
        Path("date-decoys-mock.fas").write_text(target + decoy)
        if mode == "multiple":
            Path("second-decoys-mock.fas").write_text(target + decoy)
elif args[0] != "workspace":
    sys.exit(13)
"""

FRAGPIPE = """#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--headless", action="store_true")
for flag in ["workflow", "manifest", "workdir", "config-msfragger",
             "config-ionquant", "config-philosopher", "threads"]:
    parser.add_argument("--" + flag, required=True)
args = parser.parse_args()

def require(condition, message, code=19):
    if not condition:
        print(message, file=sys.stderr)
        sys.exit(code)

manifest = Path(args.manifest)
require(manifest.is_file(), "mock manifest argv does not exist", 17)
require(args.headless and args.threads == "1", "headless or thread argument mismatch")
require(Path(args.workdir).is_dir(), "workdir unavailable")
require(Path(args.config_msfragger).read_text() == "mock MSFragger\\n", "wrong MSFragger path")
require(Path(args.config_ionquant).read_text() == "mock IonQuant\\n", "wrong IonQuant path")
require(os.access(args.config_philosopher, os.X_OK), "Philosopher unavailable")
props = dict(line.split("=", 1) for line in Path(args.workflow).read_text().splitlines() if "=" in line)
require(props["database.decoy-tag"] == "rev_", "wrong decoy configuration")
db = Path(props["database.db-path"])
require(db.is_absolute() and db.is_file(), "decoy DB unavailable")
require(">rev_mock\\n" in db.read_text(), "missing decoy")

paths = set()
for line in manifest.read_text().splitlines():
    filename, experiment, replicate, mode = line.split("\\t")
    spectrum = Path(filename)
    require(spectrum.is_absolute() and not spectrum.is_symlink() and
            (Path.cwd().parent / "spectra") in spectrum.parents,
            "spectrum was not physically staged to the work directory")
    require(spectrum.read_text() == "<mzML>" + experiment + "</mzML>\\n" and mode == "DDA",
            "wrong spectrum contents or mode")
    require(spectrum not in paths, "spectrum staging collision")
    paths.add(spectrum)

print("mock FragPipe received staged spectra and database", flush=True)
mode = "@MODE@"
require(mode != "failure", "mock FragPipe intentional failure")
if mode != "missing_output":
    Path("combined_protein.tsv").write_text("Protein\\tIntensity\\nMOCK_TARGET\\t42\\n")
"""

CASES = (
    ("FP_01", "canonical", None),
    ("FP_02", "alias", None),
    ("FP_03", "spaces", None),
    ("FP_04", "multi", None),
    ("FP_05", "no_fasta", r"--protein_ref .*required"),
    ("FP_06", "stale_fasta", r"--protein_ref .*readable regular file"),
    ("FP_07", "no_workflow", r"--proteom_workflow .*required"),
    ("FP_08", "no_manifest", r"--proteomics_manifest .*required"),
    ("FP_09", "stale_manifest", r"--proteomics_manifest .*readable regular file"),
    ("FP_10", "bad_columns", r"manifest line 1.*four tab-separated columns"),
    ("FP_11", "no_mzml", r"manifest line 1.*mzML path is required"),
    ("FP_12", "stale_mzml", r"manifest line 1.*readable mzML file"),
    ("FP_13", "dia", r"manifest line 1.*DDA"),
    ("FP_14", "no_exec", r"--fragpipe_path .*executable"),
    ("FP_15", "no_jar", r"--ionquant_jar .*readable regular file"),
    ("FP_16", "tool_failure", r"mock FragPipe intentional failure"),
    ("FP_17", "missing_output", r"no nonempty TSV or CSV output"),
    ("FP_18", "decoy_many", r"must produce exactly one target-decoy FASTA"),
    ("FP_19", "decoy_prefix", r"target and rev_ decoy headers"),
    ("FP_20", "decoy_missing", r"must produce exactly one target-decoy FASTA"),
    ("FP_21", "bad_workflow_prefix", r"exactly one database.db-path and database.decoy-tag=rev_"),
)


def make_executable(path, contents):
    path.write_text(contents, encoding="utf-8")
    path.chmod(0o755)


def run_case(case, script, kind, expected_error):
    spaced = kind == "spaces"
    inputs = case / ("inputs with ' and $ spaces" if spaced else "inputs")
    inputs.mkdir()
    fasta = inputs / "target protein.fasta" if spaced else inputs / "target.fasta"
    fasta.write_bytes((DATA / "toy.fa").read_bytes())
    workflow = inputs / ("LFQ template.workflow" if spaced else "LFQ.workflow")
    prefix = "other_" if kind == "bad_workflow_prefix" else "rev_"
    workflow.write_text("database.db-path=/archive/stale.fasta\ndatabase.decoy-tag=" + prefix + "\n")
    philo = inputs / ("mock philosopher" if spaced else "philosopher")
    fragpipe = inputs / ("mock fragpipe" if spaced else "fragpipe")
    msfragger, ionquant = inputs / "MS Fragger.jar", inputs / "Ion Quant.jar"
    if not spaced:
        msfragger, ionquant = inputs / "MSFragger.jar", inputs / "IonQuant.jar"
    msfragger.write_text("mock MSFragger\n")
    ionquant.write_text("mock IonQuant\n")

    runs = [("first", "")]
    if kind == "multi":
        runs.append(("second", "2"))
    rows = []
    for i, (experiment, replicate) in enumerate(runs):
        directory = inputs / (f"batch {i}" if spaced else f"batch{i}")
        directory.mkdir()
        mzml = directory / ("spectra ' $.mzML" if spaced else "same.mzML")
        mzml.write_text(f"<mzML>{experiment}</mzML>\n")
        spectrum = str(mzml) if not spaced else str(mzml.relative_to(inputs))
        rows.append([spectrum, experiment, replicate, "DDA"])

    if kind == "bad_columns":
        rows[0].pop(2)
    if kind == "no_mzml":
        rows[0][0] = ""
    if kind == "stale_mzml":
        rows[0][0] = "old/path.mzML"
    if kind == "dia":
        rows[0][-1] = "DIA"
    name = "manifest.fp-manifest"
    if kind == "alias":
        name = "renamed.fp-manifest"
    elif spaced:
        name = "manifest alias with spaces.fp-manifest"
    manifest = inputs / name
    manifest.write_text("".join("\t".join(row) + "\n" for row in rows))
    if kind == "no_jar":
        ionquant.unlink()
    if kind == "stale_fasta":
        fasta.unlink()
    if kind == "stale_manifest":
        manifest.unlink()

    philo_modes = {"decoy_many": "multiple", "decoy_prefix": "bad_prefix", "decoy_missing": "missing"}
    fragpipe_modes = {"tool_failure": "failure", "missing_output": "missing_output"}
    make_executable(philo, PHILOSOPHER.replace("@MODE@", philo_modes.get(kind, "single")))
    make_executable(fragpipe, FRAGPIPE.replace("@MODE@", fragpipe_modes.get(kind, "success")))
    if kind == "no_exec":
        fragpipe.chmod(0o644)

    args = ["--module", "proteomics", "--philosopher_path", philo, "--fragpipe_path", fragpipe,
            "--msfragger_jar", msfragger, "--ionquant_jar", ionquant]
    for missing, option, value in (("no_fasta", "--protein_ref", fasta),
                                   ("no_workflow", "--proteom_workflow", workflow),
                                   ("no_manifest", "--proteomics_manifest", manifest)):
        if kind != missing:
            args.extend((option, value))

    result = nextflow(case, script, *args)
    output = case / "out/protein/fragpipe_output"
    if expected_error:
        outcome = rejected(result, expected_error)
        equal(output.exists(), False)
        if kind in {"no_fasta", "stale_fasta", "no_workflow", "no_manifest", "stale_manifest", "bad_columns",
                    "no_mzml", "stale_mzml", "dia", "no_exec", "no_jar"}:
            equal((case / "out/protein_ref").exists(), False)
        return f"{outcome}; no FragPipe result published"

    success(result)
    equal((output / "combined_protein.tsv").read_text(), PSM)
    equal((output / "fragpipe.log").read_text(), "mock FragPipe received staged spectra and database\n")
    rows = read_tsv(output / "manifest.fp-manifest")
    equal([row[1:] for row in rows], [[exp, rep, "DDA"] for exp, rep in runs])
    for row in rows:
        equal(Path(row[0]).read_text(), f"<mzML>{row[1]}</mzML>\n")
    return f"exit 0; {len(runs)} staged mzML inputs, manifest fields, log and exact synthetic result verified"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=[name for name, _, _ in CASES])
    parser.add_argument("--pipeline", type=Path, default=REPO / "main.nf")
    args = parser.parse_args()
    output = Path(tempfile.mkdtemp(prefix="drs-fragpipe-tests-")).resolve()
    print(f"FragPipe test artifacts: {output}", flush=True)
    count = failures = 0
    with (output / "results.log").open("w", encoding="utf-8") as log:
        for name, kind, pattern in CASES:
            if args.case and args.case != name:
                continue
            directory = output / name
            directory.mkdir()
            count += 1
            try:
                record = f"PASS {name} | " + run_case(directory, args.pipeline.resolve(), kind, pattern)
            except (AssertionError, OSError, ValueError) as error:
                failures += 1
                record = f"FAIL {name} | {error}"
            print(record, flush=True)
            log.write(record + "\n")
        totals = f"{count - failures}/{count} FragPipe interface tests passed, {failures} failed"
        print(totals, flush=True)
        log.write(totals + "\n")
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
