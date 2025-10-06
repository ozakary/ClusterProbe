#!/usr/bin/env python3
"""
Xenon Environment Analysis Script for Trajectory Files

This script analyzes the local environment around Xenon atoms in trajectory structures,
flags anomalous snapshots based on coordination number thresholds.

Usage:
    python xe_trajectory_analyzer.py trajectory.xyz [options]

Options:
    --rcut FLOAT          Cutoff radius for neighbor search (default: 6.0 Å)
    --num_atoms INT       Minimum number of neighbors for good snapshots (default: 10)

Files analyzed:
    - Multi-frame XYZ file containing trajectory snapshots
"""

import sys
import argparse
import numpy as np
from tqdm import tqdm
from ase.io import read


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze local environment around Xenon atoms in trajectory snapshots",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        "xyz_file",
        type=str,
        help="XYZ trajectory file containing multiple snapshots"
    )
    
    parser.add_argument(
        "--rcut", 
        type=float, 
        default=6.0,
        help="Cutoff radius for neighbor search in Angstroms (default: 6.0)"
    )
    
    parser.add_argument(
        "--num_atoms", 
        type=int, 
        default=10,
        help="Minimum number of neighbors for good snapshots (default: 10)"
    )
    
    return parser.parse_args()


def analyze_snapshot(atoms, snapshot_idx, rcut, min_neighbors):
    """
    Analyze a single snapshot.
    
    Args:
        atoms: ASE Atoms object
        snapshot_idx: Index of the snapshot
        rcut: Cutoff radius for neighbors
        min_neighbors: Minimum number of neighbors required
        
    Returns:
        dict: Analysis results containing snapshot info and neighbor counts
    """
    try:
        # Find Xenon atoms
        xe_indices = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == 'Xe']
        
        if len(xe_indices) == 0:
            return {
                'snapshot': snapshot_idx,
                'success': False,
                'error': 'No Xenon atoms found',
                'n_atoms': len(atoms),
                'n_xe': 0,
                'xe_neighbors': [],
                'is_anomalous': True
            }
        
        # Calculate neighbors for each Xe atom
        xe_neighbor_counts = []
        xe_detailed_info = []
        
        for xe_idx in xe_indices:
            xe_pos = atoms.positions[xe_idx]
            
            # Calculate distances to all other atoms
            distances = []
            neighbor_indices = []
            
            for i, pos in enumerate(atoms.positions):
                if i != xe_idx:  # Exclude the Xe atom itself
                    dist = np.linalg.norm(pos - xe_pos)
                    if dist <= rcut:
                        distances.append(dist)
                        neighbor_indices.append(i)
            
            neighbor_count = len(distances)
            xe_neighbor_counts.append(neighbor_count)
            
            # Store detailed information
            xe_info = {
                'xe_index': xe_idx,
                'position': xe_pos,
                'neighbor_count': neighbor_count,
                'neighbor_indices': neighbor_indices,
                'distances': distances,
                'min_distance': min(distances) if distances else None,
                'max_distance': max(distances) if distances else None,
                'avg_distance': np.mean(distances) if distances else None
            }
            xe_detailed_info.append(xe_info)
        
        # Check if snapshot is anomalous (any Xe with too few neighbors)
        is_anomalous = any(count < min_neighbors for count in xe_neighbor_counts)
        
        return {
            'snapshot': snapshot_idx,
            'success': True,
            'error': None,
            'n_atoms': len(atoms),
            'n_xe': len(xe_indices),
            'xe_neighbors': xe_neighbor_counts,
            'xe_detailed_info': xe_detailed_info,
            'is_anomalous': is_anomalous,
            'atoms': atoms
        }
        
    except Exception as e:
        return {
            'snapshot': snapshot_idx,
            'success': False,
            'error': str(e),
            'n_atoms': 0,
            'n_xe': 0,
            'xe_neighbors': [],
            'is_anomalous': True
        }


def print_summary_statistics(results, rcut, min_neighbors):
    """Print comprehensive summary statistics."""
    print("\n" + "="*80)
    print("SNAPSHOT ANALYSIS SUMMARY")
    print("="*80)
    
    print(f"Analysis Parameters:")
    print(f"  - Cutoff radius: {rcut:.2f} Å")
    print(f"  - Minimum neighbors threshold: {min_neighbors}")
    print()
    
    # Basic statistics
    total_snapshots = len(results)
    successful_analyses = sum(1 for r in results if r['success'])
    failed_analyses = total_snapshots - successful_analyses
    
    print(f"Snapshot Statistics:")
    print(f"  - Total snapshots found: {total_snapshots}")
    print(f"  - Successfully analyzed: {successful_analyses}")
    print(f"  - Failed analyses: {failed_analyses}")
    print()
    
    if failed_analyses > 0:
        print("Failed analyses:")
        for r in results:
            if not r['success']:
                print(f"  - Snapshot {r['snapshot']}: {r['error']}")
        print()
    
    # Filter successful results for further analysis
    good_results = [r for r in results if r['success']]
    
    if not good_results:
        print("No successful analyses to report.")
        return
    
    # Atom count statistics
    atom_counts = [r['n_atoms'] for r in good_results]
    xe_counts = [r['n_xe'] for r in good_results]
    
    print(f"Snapshot Composition:")
    print(f"  - Total atoms per snapshot: {min(atom_counts)} - {max(atom_counts)} "
          f"(avg: {np.mean(atom_counts):.1f})")
    print(f"  - Xe atoms per snapshot: {min(xe_counts)} - {max(xe_counts)} "
          f"(avg: {np.mean(xe_counts):.1f})")
    print()
    
    # Coordination analysis
    all_neighbor_counts = []
    for r in good_results:
        all_neighbor_counts.extend(r['xe_neighbors'])
    
    print(f"Xenon Coordination Analysis:")
    print(f"  - Total Xe atoms analyzed: {len(all_neighbor_counts)}")
    print(f"  - Coordination numbers: {min(all_neighbor_counts)} - {max(all_neighbor_counts)} "
          f"(avg: {np.mean(all_neighbor_counts):.1f} ± {np.std(all_neighbor_counts):.1f})")
    print()
    
    # Anomaly statistics
    anomalous_results = [r for r in good_results if r['is_anomalous']]
    good_snapshots = len(good_results) - len(anomalous_results)
    
    print(f"Quality Assessment:")
    print(f"  - Good snapshots: {good_snapshots} ({100*good_snapshots/len(good_results):.1f}%)")
    print(f"  - Anomalous snapshots: {len(anomalous_results)} "
          f"({100*len(anomalous_results)/len(good_results):.1f}%)")
    print()
    
    if anomalous_results:
        print("Anomalous Snapshots Details:")
        for r in anomalous_results:
            snapshot_name = f"Snapshot {r['snapshot']}"
            neighbor_counts = r['xe_neighbors']
            min_neighbors_in_snapshot = min(neighbor_counts) if neighbor_counts else 0
            
            print(f"  - {snapshot_name}: Xe coordination = {neighbor_counts} "
                  f"(min: {min_neighbors_in_snapshot})")
    
    print("\n" + "="*80)


def main():
    """Main function."""
    args = parse_arguments()
    
    print("Xenon Environment Analysis for Trajectory")
    print("="*50)
    print(f"Parameters:")
    print(f"  - Input file: {args.xyz_file}")
    print(f"  - Cutoff radius: {args.rcut} Å")
    print(f"  - Minimum neighbors: {args.num_atoms}")
    print()
    
    # Read trajectory file
    print("Reading trajectory file...")
    try:
        # Read all snapshots from the XYZ file
        trajectory = read(args.xyz_file, index=':')
        
        # Handle single snapshot case
        if not isinstance(trajectory, list):
            trajectory = [trajectory]
            
        print(f"Found {len(trajectory)} snapshots.")
        print()
        
    except Exception as e:
        print(f"Error reading trajectory file: {str(e)}")
        sys.exit(1)
    
    # Analyze snapshots
    print("Analyzing snapshots...")
    results = []
    
    for idx, atoms in enumerate(tqdm(trajectory, desc="Processing snapshots")):
        result = analyze_snapshot(atoms, idx, args.rcut, args.num_atoms)
        results.append(result)
    
    # Print detailed results
    print("\nDetailed Results:")
    print("-" * 80)
    print(f"{'Snapshot':<12} | {'Atoms':<5} | {'Xe':<2} | {'Neighbors':<20} | {'Status'}")
    print("-" * 80)
    for result in results:
        snapshot_name = f"Snapshot {result['snapshot']}"
        
        if result['success']:
            status = "ANOMALOUS" if result['is_anomalous'] else "GOOD"
            neighbors_str = ", ".join(map(str, result['xe_neighbors']))
            
            print(f"{snapshot_name:<12} | {result['n_atoms']:4d} | {result['n_xe']} | "
                  f"[{neighbors_str:<18}] | {status}")
        else:
            print(f"{snapshot_name:<12} | ERROR: {result['error']}")
    print("-" * 80)
    
    # Print summary statistics
    print_summary_statistics(results, args.rcut, args.num_atoms)
    
    print("\nAnalysis completed!")


if __name__ == "__main__":
    main()
