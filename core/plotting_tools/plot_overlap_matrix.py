import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import logging

logger = logging.getLogger(__name__)


def plot_overlap_heatmap(overlap_matrix, lambda_values, output_png,
                         cmap='viridis', show_values=False, figsize=(12, 10),
                         threshold=0.03):
    """
    Generate overlap matrix heatmap visualization.

    Args:
        overlap_matrix: [K, K] overlap matrix
        lambda_values: Lambda state labels (list of strings)
        output_png: Output PNG file path
        cmap: Colormap name (default: 'viridis' - purple to yellow)
        show_values: Whether to annotate cells with values (default: False for large matrices)
        figsize: Figure size tuple
        threshold: Overlap threshold for marking poor regions

    Returns:
        None (saves figure to file)
    """
    K = overlap_matrix.shape[0]
    logger.info(f"Generating overlap matrix heatmap: {K}x{K}, output={output_png}")

    fig, ax = plt.subplots(figsize=figsize)

    mask = np.zeros_like(overlap_matrix, dtype=bool)

    sns.heatmap(
        overlap_matrix,
        mask=mask,
        cmap=cmap,
        vmin=0,
        vmax=0.5,
        square=True,
        linewidths=0.5,
        cbar_kws={'label': 'Overlap', 'shrink': 0.8},
        ax=ax,
        annot=show_values if K <= 20 else False,
        fmt='.2f' if show_values else '',
        annot_kws={'fontsize': 6} if show_values else None
    )

    if K <= 51:
        tick_indices = range(0, K, max(1, K // 10))
    else:
        tick_indices = range(0, K, max(1, K // 20))

    tick_labels = [f"{i}" for i in tick_indices]

    ax.set_xticks([i + 0.5 for i in tick_indices])
    ax.set_yticks([i + 0.5 for i in tick_indices])
    ax.set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(tick_labels, rotation=0, fontsize=8)

    ax.set_xlabel('Lambda State Index', fontsize=12, fontweight='bold')
    ax.set_ylabel('Lambda State Index', fontsize=12, fontweight='bold')
    ax.set_title('Phase Space Overlap Matrix', fontsize=14, fontweight='bold', pad=20)

    poor_overlaps = []
    for i in range(K - 1):
        if overlap_matrix[i, i+1] < threshold:
            poor_overlaps.append((i, i+1, overlap_matrix[i, i+1]))

    if poor_overlaps:
        for i, j, val in poor_overlaps[:10]:
            rect = plt.Rectangle((j, i), 1, 1, fill=False, edgecolor='red', linewidth=2)
            ax.add_patch(rect)

    nonzero_overlaps = overlap_matrix[overlap_matrix > 0]
    min_overlap = np.min(nonzero_overlaps) if len(nonzero_overlaps) > 0 else 0.0
    max_overlap = np.max(overlap_matrix)
    mean_overlap = np.mean(overlap_matrix)

    adjacent_overlaps = [overlap_matrix[i, i+1] for i in range(K-1)]
    min_adjacent = min(adjacent_overlaps) if adjacent_overlaps else 0

    stats_text = f"Min: {min_overlap:.3f} | Min Adjacent: {min_adjacent:.3f}"

    ax.text(0.5, 1.05, stats_text,
            transform=ax.transAxes,
            ha='center', va='bottom',
            fontsize=10,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight', transparent=False)
    plt.close()

    logger.info(f"Overlap heatmap saved to {output_png}")


def plot_overlap_diagonal(overlap_matrix, lambda_values, output_png, figsize=(10, 6)):
    """
    Plot diagonal and off-diagonal elements of overlap matrix.

    Useful for quick visualization of adjacent state overlaps.

    Args:
        overlap_matrix: [K, K] overlap matrix
        lambda_values: Lambda state labels
        output_png: Output PNG file path
        figsize: Figure size tuple
    """
    K = overlap_matrix.shape[0]

    diagonal = np.diag(overlap_matrix)
    off_diagonal = np.array([overlap_matrix[i, i+1] for i in range(K-1)])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)

    ax1.plot(range(K), diagonal, 'o-', linewidth=2, markersize=4)
    ax1.set_xlabel('Lambda State Index', fontweight='bold')
    ax1.set_ylabel('Self Overlap', fontweight='bold')
    ax1.set_title('Diagonal Elements', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([0, 1.05])

    ax2.plot(range(K-1), off_diagonal, 'o-', linewidth=2, markersize=4, color='orange')
    ax2.axhline(y=0.03, color='r', linestyle='--', label='Threshold (0.03)')
    ax2.set_xlabel('Lambda State Index', fontweight='bold')
    ax2.set_ylabel('Adjacent State Overlap', fontweight='bold')
    ax2.set_title('Off-Diagonal Elements (i, i+1)', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_ylim([0, max(1.05, max(off_diagonal) * 1.1)])

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.close()

    logger.info(f"Overlap diagonal plot saved to {output_png}")
