import subprocess
import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--run-dir", required=True,
        help="Nextflow launch directory containing a completed mergeExperiments run"
    )


@pytest.fixture(scope="session")
def run_dir(request):
    return request.config.getoption("--run-dir")


@pytest.fixture(scope="session")
def work_dirs(run_dir):
    """Return dict mapping short process name -> work dir path for the most recent run."""
    r = subprocess.run(
        "nextflow log | tail -n1 | cut -f3",
        shell=True, cwd=run_dir, capture_output=True, text=True, check=True
    )
    run_name = r.stdout.strip()
    assert run_name, f"No nextflow runs found in {run_dir}"

    r = subprocess.run(
        ["nextflow", "log", run_name, "-f", "process,workdir"],
        cwd=run_dir, capture_output=True, text=True, check=True
    )
    dirs = {}
    for line in r.stdout.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) != 2:
            continue
        process_name = parts[0].strip().split(':')[-1]
        dirs[process_name] = parts[1].strip()
    return dirs
