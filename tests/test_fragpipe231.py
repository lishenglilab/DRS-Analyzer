"""Run deterministic FragPipe 23.1 CLI cases without the external software."""
import argparse
import json
from pathlib import Path
import tempfile

from test_helpers import DATA, REPO, TESTS, equal, nextflow, read_tsv, rejected, success
from test_proteomics import PHILOSOPHER, PSM, make_executable


ENGINES = (
    "MSFragger-4.3.jar", "IonQuant-1.11.11.jar", "diaTracer-1.3.3.jar",
    "ext/bruker/libtimsdata-2-21-0-4.so",
)
DRY_MARKER = "It's a dry-run, not running the commands."
MOCK = r"""#!/usr/bin/env python3
import argparse
import json
import os
from pathlib import Path
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--headless", action="store_true")
parser.add_argument("--dry-run", action="store_true")
parser.add_argument("--ram")
for flag in ("workflow", "manifest", "workdir", "config-tools-folder", "config-python", "threads"):
    parser.add_argument("--" + flag, required=True)
args = parser.parse_args()
assert args.headless and args.threads == "1" and Path(args.workdir) == Path(".")
assert sys.argv.count("--ram") == 1 and args.ram == "@RAM@"
assert args.manifest == "manifest.fp-manifest"
assert os.access(args.config_python, os.X_OK) and Path(args.config_python).is_file()
assert Path(args.config_python).read_text().endswith("# Python path token\n")
assert (Path(args.config_tools_folder) / "MSFragger-4.3.jar").read_text() == "MSFragger-4.3.jar"
assert (Path(args.config_tools_folder) / "ext/bruker/libtimsdata-2-21-0-4.so").is_file()
manifest = Path(args.manifest)
row, = manifest.read_text().splitlines()
spectrum, label, replicate, acquisition = row.split("\t")
spectrum = Path(spectrum)
assert (Path.cwd().parent / "spectra") in spectrum.parents and not spectrum.is_symlink()
assert spectrum.read_text() == "<mzML>mock single run</mzML>\n"
assert (label, replicate, acquisition) == ("mock-single", "", "DDA")
workflow = Path(args.workflow)
settings = dict(l.split("=", 1) for l in workflow.read_text().splitlines())
assert settings["database.decoy-tag"] == "rev_"
assert ">rev_mock\n" in Path(settings["database.db-path"]).read_text()
Path("mock_argv.json").write_text(json.dumps(sys.argv[1:]))
behavior = "@BEHAVIOR@"
if behavior == "external_failure":
    print("modern mock external failure", file=sys.stderr)
    sys.exit(31)
if args.dry_run:
    assert settings["ionquant.mbr"] == "0"
    print("FragPipe version " + ("24.0" if behavior == "wrong_version_log" else "23.1"))
    if behavior != "missing_marker":
        print("It's a dry-run, not running the commands.")
    if behavior == "log_error":
        print("ERROR - Cannot execute Python: Permission denied")
    Path("fragpipe.workflow").write_text(workflow.read_text())
    Path("fragpipe-files.fp-manifest").write_text(manifest.read_text())
    if behavior != "missing_plan":
        Path("fragpipe.job").write_text("synthetic MSFragger and IonQuant plan\n")
else:
    print("modern mock search interface")
    Path("combined_protein.tsv").write_text("Protein\tIntensity\nMOCK_TARGET\t42\n")
"""

CASES = (
    ("F231_01", "normal", None),
    ("F231_02", "dry", None),
    ("F231_03", "spaces", None),
    ("F231_04", "no_tools_option", r"--fragpipe_tools_dir is required"),
    ("F231_05", "tools_file", r"--fragpipe_tools_dir must point to a readable directory"),
    ("F231_06", "no_ionquant", r"--fragpipe_tools_dir member IonQuant-1\.11\.11\.jar.*readable regular file"),
    ("F231_07", "no_native", r"--fragpipe_tools_dir member ext/bruker/.*readable regular file"),
    ("F231_08", "no_python_option", r"--fragpipe_python is required"),
    ("F231_09", "python_dir", r"--fragpipe_python must point to a readable regular file"),
    ("F231_10", "python_not_executable", r"--fragpipe_python must point to an executable file"),
    ("F231_11", "unknown_cli", r"--fragpipe_cli must be legacy or 23.1"),
    ("F231_12", "wrong_package", r"--fragpipe_cli 23.1 does not match the FragPipe launcher package"),
    ("F231_13", "default_legacy_on_231", r"--fragpipe_cli legacy does not match the FragPipe launcher package"),
    ("F231_14", "legacy_jars", r"FragPipe 23.1 uses --fragpipe_tools_dir"),
    ("F231_15", "legacy_modern_args", r"require --fragpipe_cli 23.1"),
    ("F231_16", "external_failure", r"modern mock external failure"),
    ("F231_17", "log_error", r"FragPipe dry-run reported an error"),
    ("F231_18", "missing_plan", r"FragPipe dry-run did not create a nonempty fragpipe.job"),
    ("F231_19", "wrong_version_log", r"FragPipe did not confirm a 23.1 dry-run"),
    ("F231_20", "missing_marker", r"FragPipe did not confirm a 23.1 dry-run"),
    ("F231_21", "ram_real", None),
    ("F231_22", "ram_dry", None),
    ("F231_23", "ram_missing", r"--fragpipe_ram_gb must be a positive integer"),
    ("F231_24", "ram_text", r"--fragpipe_ram_gb must be a positive integer"),
    ("F231_25", "ram_exceeds", r"--fragpipe_ram_gb.*exceeds the FRAGPIPE task memory"),
    ("F231_26", "ram_duplicate", r"--fragpipe_ram_gb may be specified only once"),
    ("F231_27", "ram_duplicate_equals", r"--fragpipe_ram_gb may be specified only once"),
    ("F231_28", "ram_in_path", None),
    ("F231_29", "ram_legacy", r"--fragpipe_ram_gb requires --fragpipe_cli 23.1"),
    ("F231_30", "ram_zero", r"--fragpipe_ram_gb must be a positive integer"),
    ("F231_31", "ram_fraction", r"--fragpipe_ram_gb must be a positive integer"),
    ("F231_32", "ram_task_default", None),
)


def run_case(case, kind, expected):
    space = kind == "spaces"
    fixtures = case / ("fixtures ' and $ spaces" if space else "fixtures")
    if kind == "ram_in_path":
        fixtures = case / "fixtures ' --fragpipe_ram_gb 3 and $ spaces"
    install = fixtures / "install"
    (install / "bin").mkdir(parents=True)
    (install / "lib").mkdir()
    packaged = "20.0" if kind in {"wrong_package", "legacy_modern_args", "ram_legacy"} else "23.1"
    (install / "lib" / f"fragpipe-{packaged}.jar").touch()
    launcher = install / "bin" / "fragpipe"
    expected_ram = "2" if kind == "ram_task_default" else "1"
    make_executable(launcher, MOCK.replace("@BEHAVIOR@", kind).replace("@RAM@", expected_ram))
    philosopher = fixtures / "philosopher"
    make_executable(philosopher, PHILOSOPHER.replace("@MODE@", "single"))
    source = fixtures / "mock input.mzML"
    source.write_text("<mzML>mock single run</mzML>\n")
    manifest = fixtures / "aliased manifest.fp-manifest"
    manifest.write_text(f"{source.name}\tmock-single\t\tDDA\n")
    workflow = fixtures / "single.workflow"
    workflow.write_text("database.db-path=old.fasta\ndatabase.decoy-tag=rev_\nionquant.mbr=0\n")
    tools = fixtures / "engine bundle"
    for name in ENGINES:
        member = tools / name
        member.parent.mkdir(parents=True, exist_ok=True)
        if (kind, name) not in {("no_ionquant", ENGINES[1]), ("no_native", ENGINES[3])}:
            member.write_text(name)
    python = fixtures / "pinned python"
    make_executable(python, "#!/usr/bin/env bash\n# Python path token\n")
    if kind == "python_not_executable":
        python.chmod(0o644)

    args = ["--module", "proteomics", "--protein_ref", DATA / "toy.fa",
            "--proteomics_manifest", manifest, "--proteom_workflow", workflow,
            "--philosopher_path", philosopher, "--fragpipe_path", launcher]
    modes = {"unknown_cli": "22.0", "legacy_modern_args": "legacy", "ram_legacy": "legacy"}
    if kind != "default_legacy_on_231":
        args += ["--fragpipe_cli", modes.get(kind, "23.1")]
    if kind not in {"no_tools_option", "ram_legacy"}:
        args += ["--fragpipe_tools_dir", workflow if kind == "tools_file" else tools]
    if kind not in {"no_python_option", "ram_legacy"}:
        args += ["--fragpipe_python", install / "bin" if kind == "python_dir" else python]
    if kind == "legacy_jars":
        args += ["--msfragger_jar", tools / ENGINES[0]]
    ram_flags = {
        "ram_real": ["--fragpipe_ram_gb", "1"],
        "ram_dry": ["--fragpipe_ram_gb=1"],
        "ram_missing": ["--fragpipe_ram_gb"],
        "ram_text": ["--fragpipe_ram_gb", "large"],
        "ram_exceeds": ["--fragpipe_ram_gb", "2"],
        "ram_duplicate": ["--fragpipe_ram_gb", "1", "--fragpipe_ram_gb", "1"],
        "ram_duplicate_equals": ["--fragpipe_ram_gb=1", "--fragpipe_ram_gb=1"],
        "ram_in_path": ["--fragpipe_ram_gb", "1"],
        "ram_legacy": ["--fragpipe_ram_gb", "1"],
        "ram_zero": ["--fragpipe_ram_gb", "0"],
        "ram_fraction": ["--fragpipe_ram_gb", "1.5"],
    }
    args += ram_flags.get(kind, [])
    dry = kind in {"dry", "spaces", "log_error", "missing_plan", "wrong_version_log", "missing_marker", "ram_dry"}
    if dry:
        args += ["--fragpipe_dry_run", "true"]
    output = case / "out/protein/fragpipe_output"
    config = TESTS / "ci.config"
    if kind in {"ram_real", "ram_dry", "ram_task_default"}:
        config = case / "ram-test.config"
        config.write_text((TESTS / "ci.config").read_text() +
                          "\nprocess { withName: FRAGPIPE { memory = '2 GB' } }\n")
    result = nextflow(case, REPO / "main.nf", *args, config=config)
    if expected:
        observed = rejected(result, expected)
        equal(result.returncode, 1)
        equal(output.exists(), False)
        if kind in {"external_failure", "log_error", "missing_plan", "wrong_version_log", "missing_marker"}:
            if kind == "external_failure":
                rejected(result, r"exit status \(31\)")
        elif kind != "ram_exceeds":
            equal((case / "out/protein_ref").exists(), False)
        if kind.startswith("ram_"):
            equal(list((case / "work").rglob("mock_argv.json")), [])
        return observed + "; no FragPipe output published"

    success(result)
    actual = json.loads((output / "mock_argv.json").read_text())
    expected_argv = ["--headless", "--workflow", "workflow.workflow", "--manifest", "manifest.fp-manifest",
                     "--workdir", ".", "--config-tools-folder", str(tools), "--config-python", str(python)]
    if dry:
        expected_argv += ["--dry-run"]
    equal(actual.count("--ram"), 1)
    equal(actual, expected_argv + ["--ram", expected_ram, "--threads", "1"])
    if dry:
        equal((output / "fragpipe.job").read_text(), "synthetic MSFragger and IonQuant plan\n")
        equal(read_tsv(output / "fragpipe-files.fp-manifest")[0][1:], ["mock-single", "", "DDA"])
        equal((output / "combined_protein.tsv").exists(), False)
        return "exit 0; exact 23.1 argv including --dry-run; one staged DDA row and plan published"
    equal((output / "combined_protein.tsv").read_text(), PSM)
    return "exit 0; exact 23.1 argv; staged manifest and mock result verified"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=[name for name, _, _ in CASES])
    params = parser.parse_args()
    output = Path(tempfile.mkdtemp(prefix="drs-fragpipe231-tests-")).resolve()
    print(f"FragPipe 23.1 test artifacts: {output}", flush=True)
    count = failed = 0
    with (output / "results.log").open("w", encoding="utf-8") as log:
        for name, kind, expected in CASES:
            if params.case and name != params.case:
                continue
            case = output / name
            case.mkdir()
            count += 1
            try:
                line = f"PASS {name} | " + run_case(case, kind, expected)
            except (AssertionError, OSError, ValueError) as error:
                failed += 1
                line = f"FAIL {name} | {error}"
            print(line, flush=True)
            log.write(line + "\n")
        totals = f"{count - failed}/{count} FragPipe 23.1 tests passed, {failed} failed"
        print(totals, flush=True)
        log.write(totals + "\n")
    return int(bool(failed))


if __name__ == "__main__":
    raise SystemExit(main())
