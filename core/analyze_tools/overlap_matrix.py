import os
import logging
import glob
import re
from typing import List, Tuple, Dict, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Module-level BAR function initialization for pymbar 3.x/4.x compatibility
_BAR_FUNC = None
_PYMBAR_V4 = None

def _init_bar():
    """Initialize BAR function once at module load."""
    global _BAR_FUNC, _PYMBAR_V4
    if _BAR_FUNC is not None:
        return
    try:
        from pymbar.other_estimators import bar
        _BAR_FUNC = bar
        _PYMBAR_V4 = True
        logger.debug("Using pymbar 4.x BAR function")
    except ImportError:
        from pymbar import BAR
        _BAR_FUNC = BAR
        _PYMBAR_V4 = False
        logger.debug("Using pymbar 3.x BAR function")

_init_bar()


def fermi(x: np.ndarray) -> np.ndarray:
    """Fermi function: f(x) = 1 / (1 + exp(x))"""
    return 1.0 / (1.0 + np.exp(np.clip(x, -700, 700)))


def compute_pairwise_overlap(u_ij: np.ndarray, u_ji: np.ndarray) -> float:
    """
    Compute overlap between two states using BAR + Fermi function acceptance rates.

    Uses pymbar.bar to solve for optimal delta_F, then computes overlap
    from Fermi function acceptance probabilities.

    Note: Cannot use bar_overlap() because it internally uses MBAR which
    requires a common reference point for all states. Our CSV data has
    each state referenced to itself, so we must use BAR (which handles
    two independent samples) and compute overlap from acceptance rates.

    Reference: Klimovich et al. (2015) J. Comput. Aided Mol. Des.

    Args:
        u_ij: [N_i] array - delta_u = u_j - u_i for samples from state i (forward work)
        u_ji: [N_j] array - delta_u = u_i - u_j for samples from state j (reverse work)

    Returns:
        overlap: Phase space overlap between 0 and 1
    """
    # Remove inf/nan values
    w_F = u_ij[np.isfinite(u_ij)]  # forward work: u_j - u_i from state i
    w_R = u_ji[np.isfinite(u_ji)]  # reverse work: u_i - u_j from state j

    if len(w_F) == 0 or len(w_R) == 0:
        return 0.0

    T_F = len(w_F)
    T_R = len(w_R)

    try:
        # Use BAR to solve for delta_F
        if _PYMBAR_V4:
            # pymbar 4.x: returns dict {'Delta_f': ..., 'dDelta_f': ...}
            result = _BAR_FUNC(w_F, w_R, compute_uncertainty=False)
            delta_F = result['Delta_f']
        else:
            # pymbar 3.x: returns tuple (df, ddf)
            delta_F = _BAR_FUNC(w_F, w_R, compute_uncertainty=False)

        # Compute shift C = ln(T_F/T_R) - delta_F (from bar implementation)
        M = np.log(T_F / T_R)
        C = M - delta_F

        # Compute Fermi function values (from pymbar uncertainty calculation)
        # fF = 1 / (1 + exp(w_F + C))
        # fR = 1 / (1 + exp(w_R - C))
        f_F = fermi(w_F + C)
        f_R = fermi(w_R - C)

        # Average acceptance rates
        avg_f_F = np.mean(f_F)
        avg_f_R = np.mean(f_R)

        # Overlap = 2 * min(avg_f_F, avg_f_R)
        # overlap = 2.0 * min(avg_f_F, avg_f_R)
        overlap = (np.sum(f_F) + np.sum(f_R)) / (T_F + T_R)

        return min(1.0, max(0.0, overlap))

    except Exception as e:
        logger.warning(f"BAR overlap calculation failed: {e}")
        return 0.0


def load_overlap_matrix_from_csvs(csv_files: List[str], window_radius: int = 5) -> Tuple[np.ndarray, List[str]]:
    """
    Load overlap matrix from multiple state CSV files.

    Each CSV contains relative energies (k_B*T units) for samples from that state
    evaluated at nearby states (current ± window_radius).

    Overlap is computed pairwise using BAR-style estimator for states within
    each other's windows.

    Args:
        csv_files: List of state_g*_s*.csv file paths (sorted by state index)
        window_radius: Radius of local energy evaluation window (default: 5)

    Returns:
        overlap_matrix: [K, K] overlap matrix (0 for non-overlapping windows)
        lambda_values: List of lambda state identifiers

    Raises:
        ValueError: If CSV files are empty or malformed
        FileNotFoundError: If required CSV files don't exist
    """
    if not csv_files:
        raise ValueError("No CSV files provided")

    K = len(csv_files)
    logger.info(f"Computing overlap matrix from {K} state files (window radius: ±{window_radius})")

    # Initialize overlap matrix
    overlap_matrix = np.zeros((K, K), dtype=np.float64)

    # Load all CSVs and extract energy data
    state_data = []
    for k, csv_file in enumerate(csv_files):
        if not os.path.exists(csv_file):
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        try:
            df = pd.read_csv(csv_file, delimiter='|')

            # Extract energy columns (exclude 'times(ps)')
            energy_cols = [col for col in df.columns
                          if '(' in col and ')' in col and 'time' not in col.lower()]

            if not energy_cols:
                raise ValueError(f"No lambda state columns found in {csv_file}")

            energies = df[energy_cols].values  # [N_k, num_evaluated_states]

            # Infer state indices from relative positions
            self_col_idx = min(k, window_radius)
            num_cols = len(energy_cols)
            first_state_idx = k - self_col_idx
            last_state_idx = first_state_idx + num_cols - 1

            state_data.append({
                'state_idx': k,
                'energies': energies,
                'first_state': first_state_idx,
                'last_state': last_state_idx,
                'self_col': self_col_idx,
                'num_samples': len(df)
            })

            logger.debug(f"State {k}: {len(df)} samples, window [{first_state_idx}, {last_state_idx}]")

        except Exception as e:
            logger.error(f"Error reading {csv_file}: {e}")
            raise

    # Compute pairwise overlaps
    logger.info("Computing pairwise overlaps...")

    for i in range(K):
        data_i = state_data[i]

        # Diagonal: self-overlap (always 1.0)
        overlap_matrix[i, i] = 0.5

        # Off-diagonal: compute overlap with states in window
        for j in range(i + 1, K):
            data_j = state_data[j]

            # Check if j is in i's window
            j_in_i_window = (data_i['first_state'] <= j <= data_i['last_state'])
            # Check if i is in j's window
            i_in_j_window = (data_j['first_state'] <= i <= data_j['last_state'])

            if j_in_i_window and i_in_j_window:
                # Both states evaluated each other

                # Get relative energies
                # u_ij: samples from i evaluated at j, relative to i
                j_col_in_i = j - data_i['first_state']
                i_col_in_i = data_i['self_col']  # i's self column
                u_ij = data_i['energies'][:, j_col_in_i] - data_i['energies'][:, i_col_in_i]

                # u_ji: samples from j evaluated at i, relative to j (reverse work)
                i_col_in_j = i - data_j['first_state']
                j_col_in_j = data_j['self_col']  # j's self column
                u_ji = data_j['energies'][:, i_col_in_j] - data_j['energies'][:, j_col_in_j]

                # Compute overlap
                overlap = compute_pairwise_overlap(u_ij, u_ji)
                overlap_matrix[i, j] = overlap
                overlap_matrix[j, i] = overlap  # Symmetric

                logger.debug(f"Overlap[{i},{j}] = {overlap:.4f}")

            elif j_in_i_window or i_in_j_window:
                # Only one-way evaluation available - estimate with less accuracy
                logger.debug(f"States {i} and {j}: asymmetric window coverage")

                if j_in_i_window:
                    # Only have u_ij
                    j_col_in_i = j - data_i['first_state']
                    i_col_in_i = data_i['self_col']
                    u_ij = data_i['energies'][:, j_col_in_i] - data_i['energies'][:, i_col_in_i]
                    u_ij_finite = u_ij[np.isfinite(u_ij)]
                    if len(u_ij_finite) > 0:
                        # Rough estimate: exp(-mean(u_ij)/2)
                        overlap = np.exp(-np.mean(np.abs(u_ij_finite)) / 2.0)
                        overlap_matrix[i, j] = overlap
                        overlap_matrix[j, i] = overlap

                elif i_in_j_window:
                    # Only have u_ji
                    i_col_in_j = i - data_j['first_state']
                    j_col_in_j = data_j['self_col']
                    u_ji = data_j['energies'][:, i_col_in_j] - data_j['energies'][:, j_col_in_j]
                    u_ji_finite = u_ji[np.isfinite(u_ji)]
                    if len(u_ji_finite) > 0:
                        overlap = np.exp(-np.mean(np.abs(u_ji_finite)) / 2.0)
                        overlap_matrix[i, j] = overlap
                        overlap_matrix[j, i] = overlap

    # Generate lambda state labels
    lambda_values = [f"State_{i}" for i in range(K)]

    logger.info(f"Overlap matrix computed: shape={overlap_matrix.shape}")
    logger.info(f"  Non-zero entries: {np.count_nonzero(overlap_matrix)}/{K*K}")
    logger.info(f"  Min overlap: {np.min(overlap_matrix):.4f}")
    logger.info(f"  Max overlap: {np.max(overlap_matrix):.4f}")
    logger.info(f"  Mean overlap: {np.mean(overlap_matrix):.4f}")

    return overlap_matrix, lambda_values


def generate_overlap_report(overlap_matrix: np.ndarray,
                           lambda_values: List[str],
                           threshold: float = 0.03) -> Dict:
    """
    Generate diagnostic report for overlap quality.

    Args:
        overlap_matrix: [K, K] overlap matrix
        lambda_values: List of lambda state identifiers
        threshold: Minimum acceptable overlap between adjacent states

    Returns:
        report_dict: Dictionary with warnings and statistics
    """
    K = overlap_matrix.shape[0]
    report = {
        'num_states': K,
        'min_overlap': float(np.min(overlap_matrix)),
        'max_overlap': float(np.max(overlap_matrix)),
        'mean_overlap': float(np.mean(overlap_matrix)),
        'diagonal_mean': float(np.mean(np.diag(overlap_matrix))),
        'poor_overlaps': [],
        'warnings': []
    }

    # Check adjacent state overlaps
    adjacent_overlaps = []
    for i in range(K - 1):
        overlap = overlap_matrix[i, i+1]
        adjacent_overlaps.append(overlap)

        if overlap < threshold:
            report['poor_overlaps'].append({
                'state_i': i,
                'state_j': i + 1,
                'overlap': float(overlap),
                'lambda_i': lambda_values[i] if i < len(lambda_values) else f"State {i}",
                'lambda_j': lambda_values[i+1] if i+1 < len(lambda_values) else f"State {i+1}"
            })

    report['min_adjacent_overlap'] = float(min(adjacent_overlaps)) if adjacent_overlaps else 0.0
    report['mean_adjacent_overlap'] = float(np.mean(adjacent_overlaps)) if adjacent_overlaps else 0.0

    if report['poor_overlaps']:
        report['warnings'].append(
            f"Found {len(report['poor_overlaps'])} adjacent state pairs with overlap < {threshold}"
        )

    if report['min_overlap'] < 0:
        report['warnings'].append(f"Negative overlap values detected (min: {report['min_overlap']:.4f})")

    logger.info(f"Overlap report: {len(report['poor_overlaps'])} poor overlaps, "
                f"min adjacent overlap: {report['min_adjacent_overlap']:.4f}")

    return report


def analyze_overlap_from_edge_folder(edge_folder: str,
                                     group_num: int = 1,
                                     window_radius: int = 5) -> Optional[Tuple[np.ndarray, List[str], Dict]]:
    """
    Convenience function to analyze overlap from an edge calculation folder.

    Args:
        edge_folder: Path to edge calculation folder containing state CSV files
        group_num: Group number for CSV pattern matching
        window_radius: Radius of local energy evaluation window (default: 5)

    Returns:
        Tuple of (overlap_matrix, lambda_values, report_dict) or None if failed
    """
    try:
        pattern = os.path.join(edge_folder, f'state_g{group_num}_s*.csv')
        csv_files = sorted(
            glob.glob(pattern),
            key=lambda x: int(re.search(r's(\d+)\.csv', x).group(1))
        )

        if not csv_files:
            logger.warning(f"No state CSV files found in {edge_folder}")
            return None

        logger.info(f"Found {len(csv_files)} state CSV files in {edge_folder}")

        overlap_matrix, lambda_values = load_overlap_matrix_from_csvs(csv_files, window_radius=window_radius)
        report = generate_overlap_report(overlap_matrix, lambda_values)

        return overlap_matrix, lambda_values, report

    except Exception as e:
        logger.error(f"Failed to analyze overlap from {edge_folder}: {e}")
        return None

if __name__ == "__main__":
    # Test case 1: Perfect overlap (same distribution)
    u_ij = np.random.randn(1000)  # random noise around 0
    u_ji = np.random.randn(1000)
    overlap = compute_pairwise_overlap(u_ij, u_ji)
    print(f"Perfect overlap: {overlap}")  # Should be close to 0.5

    # Test case 2: No overlap
    u_ij = np.ones(1000) * 10  # all very positive
    u_ji = np.ones(1000) * (-10)  # all very negative
    overlap = compute_pairwise_overlap(u_ij, u_ji)
    print(f"No overlap: {overlap}")  # Should be close to 0
