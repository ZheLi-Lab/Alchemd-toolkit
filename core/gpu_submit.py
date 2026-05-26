#!/usr/bin/env python
"""
gpu_submit.py - GPU-aware queue submission tool.

Reads a queue.lst file, monitors GPU availability via pynvml,
and submits jobs to free GPUs one job per GPU. Tracks progress
in a .state.lst file for resumable operation.

State file format (space-separated):
  STATUS work_dir submit_sh pid gpu_id
  STATUS is one of: PENDING RUNNING DONE FAILED
  pid and gpu_id are integers or 'None'

Usage:
  python gpu_submit.py run_all_queue.lst
  python gpu_submit.py run_all_queue.lst --gpus 0 1 2
  python gpu_submit.py run_all_queue.lst --mem_threshold 5000 --interval 120
"""

import argparse
import os
import sys
import time
import logging
import subprocess
import signal

# ---------------------------------------------------------------------------
# Dependency check: must happen before any other local imports
# ---------------------------------------------------------------------------
_missing = []
try:
    import pynvml
except ImportError:
    _missing.append("pynvml")

try:
    import psutil
except ImportError:
    _missing.append("psutil")

if _missing:
    print(f"ERROR: Missing required package(s): {', '.join(_missing)}")
    print(f"       Please install with: pip install {' '.join(_missing)}")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
FOREIGN_PROCESS_NAMES = (
    "python", "bash", "sh",
    "pmemd", "pmemd.cuda",          # Amber
    "gmx", "gmx_mpi", "mdrun",      # GROMACS
    "g09", "g16", "gaussian",       # Gaussian
)

STATUS_PENDING = "PENDING"
STATUS_RUNNING = "RUNNING"
STATUS_DONE    = "DONE"
STATUS_FAILED  = "FAILED"


# ---------------------------------------------------------------------------
# State file helpers
# ---------------------------------------------------------------------------

def get_state_file_path(queue_file: str) -> str:
    """run_all_queue.lst  ->  run_all_queue.state.lst"""
    base, ext = os.path.splitext(queue_file)
    return base + ".state" + ext


def load_state(state_file: str) -> list:
    jobs = []
    with open(state_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            status    = parts[0]
            work_dir  = parts[1]
            submit_sh = parts[2]
            pid    = int(parts[3]) if len(parts) > 3 and parts[3] != "None" else None
            gpu_id = int(parts[4]) if len(parts) > 4 and parts[4] != "None" else None
            jobs.append({
                "status":    status,
                "work_dir":  work_dir,
                "submit_sh": submit_sh,
                "pid":       pid,
                "gpu_id":    gpu_id,
            })
    return jobs


def save_state(state_file: str, jobs: list) -> None:
    """Atomically write state file using a temp file + os.replace."""
    tmp = state_file + ".tmp"
    with open(tmp, "w") as f:
        for job in jobs:
            pid_s    = str(job["pid"])    if job["pid"]    is not None else "None"
            gpu_id_s = str(job["gpu_id"]) if job["gpu_id"] is not None else "None"
            f.write(f"{job['status']} {job['work_dir']} {job['submit_sh']} {pid_s} {gpu_id_s}\n")
    os.replace(tmp, state_file)


def init_state_from_queue(queue_file: str, state_file: str) -> list:
    jobs = []
    with open(queue_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            jobs.append({
                "status":    STATUS_PENDING,
                "work_dir":  parts[0],
                "submit_sh": parts[1],
                "pid":       None,
                "gpu_id":    None,
            })
    save_state(state_file, jobs)
    logger.info(f"Initialized state file with {len(jobs)} jobs: {state_file}")
    return jobs


# ---------------------------------------------------------------------------
# Process helpers
# ---------------------------------------------------------------------------

def check_process_alive(pid: int) -> bool:
    """Return True if pid refers to a live, non-zombie process we launched."""
    if not psutil.pid_exists(pid):
        return False
    try:
        proc = psutil.Process(pid)
        if proc.status() == psutil.STATUS_ZOMBIE:
            # Reap the zombie if we are still its parent; ignore if not.
            try:
                os.waitpid(pid, os.WNOHANG)
            except ChildProcessError:
                pass
            return False
        name = proc.name().lower()
        # Use exact token matching to avoid 'sh' matching 'ssh', 'fish', etc.
        return name in ("bash", "python", "python3", "sh")
    except psutil.NoSuchProcess:
        return False


def update_running_jobs(jobs: list) -> list:
    """Transition RUNNING jobs to DONE when their process is no longer alive."""
    for job in jobs:
        if job["status"] != STATUS_RUNNING:
            continue
        pid = job["pid"]
        if pid is None or not check_process_alive(pid):
            logger.info(f"Job finished: {job['submit_sh']} (PID={pid})")
            job["status"]  = STATUS_DONE
            job["pid"]     = None
            job["gpu_id"]  = None
    return jobs


def get_managed_pids(jobs: list) -> set:
    return {j["pid"] for j in jobs if j["status"] == STATUS_RUNNING and j["pid"] is not None}


# ---------------------------------------------------------------------------
# GPU helpers
# ---------------------------------------------------------------------------

def check_gpu_availability(allowed_gpus, managed_pids: set, mem_threshold_mb: int) -> list:
    """
    Return list of GPU IDs that are free.

    A GPU is considered free when:
      - free memory >= mem_threshold_mb
      - no foreign (non-managed) python / bash / sh compute processes are running on it
    """
    available = []
    try:
        pynvml.nvmlInit()
        device_count = pynvml.nvmlDeviceGetCount()

        for gpu_id in range(device_count):
            if allowed_gpus is not None and gpu_id not in allowed_gpus:
                continue

            handle = pynvml.nvmlDeviceGetHandleByIndex(gpu_id)

            # Memory check
            mem_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
            free_mb  = mem_info.free / (1024 * 1024)
            if free_mb < mem_threshold_mb:
                logger.debug(f"GPU {gpu_id}: low free memory ({free_mb:.0f} MB), skipping")
                continue

            # Compute process check
            try:
                compute_procs = pynvml.nvmlDeviceGetComputeRunningProcesses(handle)
            except pynvml.NVMLError:
                compute_procs = []

            occupied_by_others = False
            for proc in compute_procs:
                pid = proc.pid
                if pid in managed_pids:
                    continue
                try:
                    name = psutil.Process(pid).name().lower()
                    if any(tok in name for tok in FOREIGN_PROCESS_NAMES):
                        logger.debug(f"GPU {gpu_id}: foreign process PID={pid} ({name})")
                        occupied_by_others = True
                        break
                except psutil.NoSuchProcess:
                    pass

            if not occupied_by_others:
                available.append(gpu_id)

    finally:
        try:
            pynvml.nvmlShutdown()
        except Exception:
            pass

    return available


# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------

def submit_job(job: dict, gpu_id: int) -> int:
    """
    Launch submit.sh on the given GPU in a detached process session.
    The child process survives parent termination.
    Returns the PID of the launched bash process.
    """
    proc = subprocess.Popen(
        ["/bin/bash", job["submit_sh"], str(gpu_id)],
        cwd=job["work_dir"],
        start_new_session=True,   # setsid() - detach from parent session
        stdin=subprocess.DEVNULL,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return proc.pid


# ---------------------------------------------------------------------------
# Progress reporting
# ---------------------------------------------------------------------------

def log_progress(jobs: list) -> None:
    counts = {STATUS_PENDING: 0, STATUS_RUNNING: 0, STATUS_DONE: 0, STATUS_FAILED: 0}
    for job in jobs:
        counts[job["status"]] = counts.get(job["status"], 0) + 1
    logger.info(
        f"Progress | Pending={counts[STATUS_PENDING]}  "
        f"Running={counts[STATUS_RUNNING]}  "
        f"Done={counts[STATUS_DONE]}  "
        f"Failed={counts[STATUS_FAILED]}  "
        f"Total={len(jobs)}"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="GPU-aware queue submission script.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python gpu_submit.py run_all_queue.lst
  python gpu_submit.py run_all_queue.lst --gpus 0 1
  python gpu_submit.py run_all_queue.lst --gpus 2 --mem_threshold 5000 --sleep 60
        """,
    )
    parser.add_argument(
        "queue",
        type=str,
        help="Path to the queue.lst file (two columns: work_dir submit_sh).",
    )
    parser.add_argument(
        "--gpus",
        type=int,
        nargs="+",
        default=None,
        help="GPU IDs to submit to (e.g. --gpus 0 1 2). Default: all detected GPUs.",
    )
    parser.add_argument(
        "--mem_threshold",
        type=int,
        default=3000,
        help="Minimum free GPU memory in MB to consider a GPU available. Default: 3000.",
    )
    parser.add_argument(
        "--sleep",
        type=int,
        default=300,
        help="Sleep duration in seconds between polling cycles. Default: 300.",
    )

    args = parser.parse_args()

    queue_file = os.path.abspath(args.queue)
    if not os.path.exists(queue_file):
        logger.error(f"Queue file not found: {queue_file}")
        sys.exit(1)

    # Verify NVML is functional before doing anything
    try:
        pynvml.nvmlInit()
        device_count = pynvml.nvmlDeviceGetCount()
        pynvml.nvmlShutdown()
        logger.info(f"NVML OK. {device_count} GPU(s) detected.")
    except pynvml.NVMLError as exc:
        logger.error(f"NVML initialization failed: {exc}")
        logger.error("Ensure NVIDIA drivers are installed and GPUs are accessible.")
        sys.exit(1)

    allowed_gpus = set(args.gpus) if args.gpus is not None else None
    if allowed_gpus is not None:
        logger.info(f"Restricting submission to GPU(s): {sorted(allowed_gpus)}")

    state_file = get_state_file_path(queue_file)

    # ------------------------------------------------------------------
    # Load or initialize state
    # ------------------------------------------------------------------
    if os.path.exists(state_file):
        jobs = load_state(state_file)
        # On resume: reconcile RUNNING jobs whose process may have finished
        jobs = update_running_jobs(jobs)
        save_state(state_file, jobs)
        logger.info(f"Resumed from state file: {state_file}")
        log_progress(jobs)

        if all(j["status"] in (STATUS_DONE, STATUS_FAILED) for j in jobs):
            logger.info("All jobs already completed. Nothing to do.")
            return
    else:
        jobs = init_state_from_queue(queue_file, state_file)

    # ------------------------------------------------------------------
    # Main scheduling loop
    # ------------------------------------------------------------------
    while True:
        # Update finished jobs
        jobs = update_running_jobs(jobs)
        save_state(state_file, jobs)
        log_progress(jobs)

        pending = sum(1 for j in jobs if j["status"] == STATUS_PENDING)
        running = sum(1 for j in jobs if j["status"] == STATUS_RUNNING)

        if pending == 0 and running == 0:
            logger.info("All jobs completed. Exiting.")
            break

        # GPUs locked by our own managed jobs.
        # This intentionally includes jobs in analysis phases (CPU-only, no GPU
        # compute visible via NVML) since their bash process is still alive.
        # A GPU is only released when the entire bash script exits.
        occupied_gpus = {
            j["gpu_id"]
            for j in jobs
            if j["status"] == STATUS_RUNNING and j["gpu_id"] is not None
        }

        # Query GPU availability (ignores our managed PIDs as foreign)
        managed_pids = get_managed_pids(jobs)
        try:
            available_gpus = check_gpu_availability(allowed_gpus, managed_pids, args.mem_threshold)
        except Exception as exc:
            logger.warning(f"GPU availability check failed: {exc}")
            available_gpus = []

        # GPUs that are available AND not already running one of our jobs
        free_gpus = [g for g in available_gpus if g not in occupied_gpus]
        logger.info(f"Available GPUs: {available_gpus}  |  Free for new jobs: {free_gpus}")

        if free_gpus and pending > 0:
            pending_jobs = [j for j in jobs if j["status"] == STATUS_PENDING]
            for gpu_id in free_gpus:
                if not pending_jobs:
                    break
                job = pending_jobs.pop(0)
                try:
                    pid = submit_job(job, gpu_id)
                    job["status"]  = STATUS_RUNNING
                    job["pid"]     = pid
                    job["gpu_id"]  = gpu_id
                    occupied_gpus.add(gpu_id)
                    done_count = sum(1 for j in jobs if j["status"] in (STATUS_DONE, STATUS_FAILED))
                    logger.info(
                        f"Submitted -> GPU {gpu_id}: {job['submit_sh']} (PID={pid}) "
                        f"[{done_count}/{len(jobs)} completed]"
                    )
                except Exception as exc:
                    logger.error(f"Submission failed for {job['submit_sh']}: {exc}")
                    job["status"] = STATUS_FAILED
            save_state(state_file, jobs)
        elif pending > 0:
            logger.info(f"No free GPUs. Waiting {args.sleep}s ...")

        # Check again after submission before sleeping
        pending = sum(1 for j in jobs if j["status"] == STATUS_PENDING)
        running = sum(1 for j in jobs if j["status"] == STATUS_RUNNING)
        if pending == 0 and running == 0:
            logger.info("All jobs completed. Exiting.")
            break

        time.sleep(args.sleep)


if __name__ == "__main__":
    main()
