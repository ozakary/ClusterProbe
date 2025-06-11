#!/usr/bin/env python3
"""
Cluster Sanity Check Script

This script analyzes the local environment around Xenon atoms in cluster structures
and flags anomalous clusters based on coordination number thresholds.

Usage:
    python cluster_sanity_check.py [options]

Options:
    --rcut FLOAT          Cutoff radius for neighbor search (default: 6.0 Å)
    --num_atoms INT       Minimum number of neighbors for good clusters (default: 10)
    --sort_out_anomalies  Whether to move anomalous clusters to bad_seeds folder
                         (choices: yes/no, default: no)
    --base_dir PATH       Base directory containing cluster folders (default: .)
"""

import os
import sys
import argparse
import shutil
import glob
from pathlib import Path
import numpy as np
from tqdm import tqdm
from ase import Atoms
from ase.io import read
from ase.neighborlist import neighbor_list


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze local environment around Xenon atoms in clusters",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
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
        help="Minimum number of neighbors for good clusters (default: 10)"
    )
    
    parser.add_argument(
        "--sort_out_anomalies", 
        choices=["yes", "no"], 
        default="no",
        help="Move anomalous clusters to bad_seeds folder (default: no)"
    )
    
    parser.add_argument(
        "--base_dir", 
        type=str, 
        default=".",
        help="Base directory containing cluster folders (default: current directory)"
    )
    
    return parser.parse_args()


def find_cluster_files(base_dir):
    """Find all cluster XYZ files in the directory structure."""
    pattern = os.path.join(base_dir, "cluster_*/coord_*ClusterAroundXeNew.xyz")
    cluster_files = glob.glob(pattern)
    
    # Sort by cluster number for consistent ordering
    def extract_cluster_number(filepath):
        try:
            # Extract cluster number from path
            cluster_dir = os.path.basename(os.path.dirname(filepath))
            return int(cluster_dir.split('_')[1])
        except (IndexError, ValueError):
            return 0
    
    cluster_files.sort(key=extract_cluster_number)
    return cluster_files


def analyze_cluster(xyz_file, rcut, min_neighbors):
    """
    Analyze a single cluster file.
    
    Returns:
        dict: Analysis results containing cluster info and neighbor counts
    """
    try:
        # Read the cluster structure
        atoms = read(xyz_file)
        
        # Find Xenon atoms
        xe_indices = [i for i, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == 'Xe']
        
        if len(xe_indices) == 0:
            return {
                'file': xyz_file,
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
        
        # Check if cluster is anomalous (any Xe with too few neighbors)
        is_anomalous = any(count < min_neighbors for count in xe_neighbor_counts)
        
        return {
            'file': xyz_file,
            'success': True,
            'error': None,
            'n_atoms': len(atoms),
            'n_xe': len(xe_indices),
            'xe_neighbors': xe_neighbor_counts,
            'xe_detailed_info': xe_detailed_info,
            'is_anomalous': is_anomalous,
            'atoms': atoms  # Keep for potential further analysis
        }
        
    except Exception as e:
        return {
            'file': xyz_file,
            'success': False,
            'error': str(e),
            'n_atoms': 0,
            'n_xe': 0,
            'xe_neighbors': [],
            'is_anomalous': True
        }


def move_anomalous_cluster(cluster_file, bad_seeds_dir):
    """Move an anomalous cluster folder to the bad_seeds directory."""
    cluster_dir = os.path.dirname(cluster_file)
    cluster_name = os.path.basename(cluster_dir)
    
    source_path = cluster_dir
    target_path = os.path.join(bad_seeds_dir, cluster_name)
    
    try:
        # Create bad_seeds directory if it doesn't exist
        os.makedirs(bad_seeds_dir, exist_ok=True)
        
        # Move the entire cluster folder
        if os.path.exists(target_path):
            shutil.rmtree(target_path)  # Remove existing folder if present
        
        shutil.move(source_path, target_path)
        return True, f"Moved {cluster_name} to bad_seeds"
        
    except Exception as e:
        return False, f"Failed to move {cluster_name}: {str(e)}"


def print_summary_statistics(results, rcut, min_neighbors):
    """Print comprehensive summary statistics."""
    print("\n" + "="*80)
    print("CLUSTER ANALYSIS SUMMARY")
    print("="*80)
    
    print(f"Analysis Parameters:")
    print(f"  - Cutoff radius: {rcut:.2f} Å")
    print(f"  - Minimum neighbors threshold: {min_neighbors}")
    print()
    
    # Basic statistics
    total_clusters = len(results)
    successful_analyses = sum(1 for r in results if r['success'])
    failed_analyses = total_clusters - successful_analyses
    
    print(f"File Statistics:")
    print(f"  - Total cluster files found: {total_clusters}")
    print(f"  - Successfully analyzed: {successful_analyses}")
    print(f"  - Failed analyses: {failed_analyses}")
    print()
    
    if failed_analyses > 0:
        print("Failed analyses:")
        for r in results:
            if not r['success']:
                print(f"  - {os.path.basename(r['file'])}: {r['error']}")
        print()
    
    # Filter successful results for further analysis
    good_results = [r for r in results if r['success']]
    
    if not good_results:
        print("No successful analyses to report.")
        return
    
    # Atom count statistics
    atom_counts = [r['n_atoms'] for r in good_results]
    xe_counts = [r['n_xe'] for r in good_results]
    
    print(f"Cluster Composition:")
    print(f"  - Total atoms per cluster: {min(atom_counts)} - {max(atom_counts)} "
          f"(avg: {np.mean(atom_counts):.1f})")
    print(f"  - Xe atoms per cluster: {min(xe_counts)} - {max(xe_counts)} "
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
    good_clusters = len(good_results) - len(anomalous_results)
    
    print(f"Quality Assessment:")
    print(f"  - Good clusters: {good_clusters} ({100*good_clusters/len(good_results):.1f}%)")
    print(f"  - Anomalous clusters: {len(anomalous_results)} "
          f"({100*len(anomalous_results)/len(good_results):.1f}%)")
    print()
    
    if anomalous_results:
        print("Anomalous Clusters Details:")
        for r in anomalous_results:
            cluster_name = os.path.basename(os.path.dirname(r['file']))
            neighbor_counts = r['xe_neighbors']
            min_neighbors_in_cluster = min(neighbor_counts) if neighbor_counts else 0
            print(f"  - {cluster_name}: Xe coordination = {neighbor_counts} "
                  f"(min: {min_neighbors_in_cluster})")
    
    print("\n" + "="*80)


def main():
    """Main function."""
    args = parse_arguments()
    
    print("Cluster Sanity Check Analysis")
    print("="*50)
    print(f"Parameters:")
    print(f"  - Base directory: {args.base_dir}")
    print(f"  - Cutoff radius: {args.rcut} Å")
    print(f"  - Minimum neighbors: {args.num_atoms}")
    print(f"  - Sort out anomalies: {args.sort_out_anomalies}")
    print()
    
    # Find cluster files
    print("Searching for cluster files...")
    cluster_files = find_cluster_files(args.base_dir)
    
    if not cluster_files:
        print(f"No cluster files found in {args.base_dir}")
        print("Expected pattern: cluster_*/coord_*ClusterAroundXeNew.xyz")
        sys.exit(1)
    
    print(f"Found {len(cluster_files)} cluster files.")
    print()
    
    # Analyze clusters
    print("Analyzing clusters...")
    results = []
    
    for cluster_file in tqdm(cluster_files, desc="Processing clusters"):
        result = analyze_cluster(cluster_file, args.rcut, args.num_atoms)
        results.append(result)
    
    # Print detailed results
    print("\nDetailed Results:")
    print("-" * 80)
    for result in results:
        cluster_name = os.path.basename(os.path.dirname(result['file']))
        
        if result['success']:
            status = "ANOMALOUS" if result['is_anomalous'] else "GOOD"
            neighbors_str = ", ".join(map(str, result['xe_neighbors']))
            print(f"{cluster_name:15} | {result['n_atoms']:4d} atoms | "
                  f"{result['n_xe']} Xe | Neighbors: [{neighbors_str:15}] | {status}")
        else:
            print(f"{cluster_name:15} | ERROR: {result['error']}")
    
    # Print summary statistics
    print_summary_statistics(results, args.rcut, args.num_atoms)
    
    # Handle anomalous cluster sorting
    if args.sort_out_anomalies == "yes":
        anomalous_results = [r for r in results if r['success'] and r['is_anomalous']]
        
        if anomalous_results:
            print(f"\nMoving {len(anomalous_results)} anomalous clusters to bad_seeds folder...")
            
            bad_seeds_dir = os.path.join(args.base_dir, "bad_seeds")
            moved_count = 0
            
            for result in tqdm(anomalous_results, desc="Moving anomalous clusters"):
                success, message = move_anomalous_cluster(result['file'], bad_seeds_dir)
                if success:
                    moved_count += 1
                else:
                    print(f"Warning: {message}")
            
            print(f"Successfully moved {moved_count}/{len(anomalous_results)} anomalous clusters.")
        else:
            print("\nNo anomalous clusters found to move.")
    
    print("\nAnalysis completed!")


if __name__ == "__main__":
    main()
