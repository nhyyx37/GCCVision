#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GCCVision: Genome Contribution Calculator and Visualizer v1.0

An Integrated Toolkit for Calculating and Visualizing Parental Genome Contribution in Breeding Populations

Usage:
GCCVision.py -v input.vcf -a Parent1 -b Parent2 -s Sample1,Sample2 -o results -t 4 -f 2 -w 5

Or using existing informative sites file:
GCCVision.py -i results_informative_sites.tsv -v input.vcf -s Sample1,Sample2 -o results -t 4 -f 2 -w 5
"""

import sys
import gzip
import argparse
import os
import csv
import time
import pandas as pd
from collections import defaultdict, Counter
import concurrent.futures

def print_header(title):
    """Print main program header"""
    print(f"\n{title}")

def print_section(title):
    """Print section"""
    pass  # 简化：不显示

def print_step(step_num, description):
    """Print simplified step"""
    print(f"\n[{step_num}] {description}")

def print_info(key, value):
    """Print key-value information"""
    # 只显示特定的信息
    allowed_keys = ["Sample columns found", "Summary file"]
    if key in allowed_keys:
        print(f"{key}: {value}")

def print_result(title, data):
    """Print results"""
    if not data:
        return
    
    # CONFIGURATION 显示
    if title == "CONFIGURATION":
        for key, value in data.items():
            print(f"{key}: {value}")
        return
    
    # INFORMATIVE SITES SUMMARY 显示
    if title == "INFORMATIVE SITES SUMMARY":
        print("Informative Sites Summary")
        for key, value in data.items():
            print(f"{key}: {value}")
        return
    
    # 样本结果显示
    if "Parent A contribution" in data and "Parent B contribution" in data:
        sample_name = title.replace("SAMPLE ", "").replace(" RESULTS", "")
        print(f"{sample_name}: A={data['Parent A contribution']}, B={data['Parent B contribution']}")

def print_status(message, status_type="info"):
    """Print status message"""
    if status_type == "error":
        print(f"Error: {message}")
    elif status_type == "success":
        print(f"✓ {message}")

def is_gzipped(filename):
    """Check if file is gzip compressed"""
    try:
        with open(filename, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except IOError:
        print(f"Error: Cannot open file {filename}")
        sys.exit(1)

def open_file(filename):
    """Open file based on file type"""
    try:
        if is_gzipped(filename):
            return gzip.open(filename, 'rt')
        else:
            return open(filename, 'r')
    except Exception as e:
        print(f"Error: Exception occurred while opening file {filename}: {e}")
        sys.exit(1)



def parse_gt_field(gt_str):
    """Parse genotype field, return genotype index list and whether it's heterozygous"""
    # Unified processing of genotype separators (/ or |)
    genotype_parts = []
    separator = None
    
    if '/' in gt_str:
        genotype_parts = gt_str.split('/')
        separator = '/'
    elif '|' in gt_str:
        genotype_parts = gt_str.split('|')
        separator = '|'
    else:
        # Handle haploid calls
        if gt_str.count('.') == 0 and gt_str.isdigit():
            genotype_parts = [gt_str]
        else:
            return None, None, "bad_gt_format"
    
    # Check if there's missing data in genotype
    if '.' in genotype_parts:
        return None, None, "missing_gt_allele"
    
    try:
        # Ensure all parts are numeric
        if not all(part.isdigit() for part in genotype_parts):
            raise ValueError(f"Genotype contains non-numeric parts: {genotype_parts}")
        genotype_indices = [int(i) for i in genotype_parts]
    except ValueError:
        return None, None, "non_integer_gt"
    
    # Determine if heterozygous
    is_heterozygous = False
    if len(genotype_indices) > 1:
        is_heterozygous = len(set(genotype_indices)) > 1  # Has different allele indices
    
    return genotype_indices, separator, "success"

def is_homozygous(gt_indices):
    """Check if genotype is homozygous"""
    if not gt_indices:
        return False
    return len(set(gt_indices)) == 1

def get_window_representative(window, prev_repr=None):
    """
    Determine window representative value
    
    Args:
    window -- Window containing SNP sites
    prev_repr -- Previous window representative value
    
    Returns:
    Window representative value
    """
    # Calculate frequency of each marker value in window
    counter = Counter(window)
    total = len(window)
    
    # Find marker value with frequency >50%
    for mark, count in counter.items():
        if count / total > 0.5:
            return mark
    
    # If no marker value has frequency >50%, inherit previous window representative value
    return prev_repr

def filter_snps(input_file, output_file, filter_window_size=5):
    """Use sliding window method for SNP site noise reduction"""
    # Read input file
    df = pd.read_csv(input_file, sep='\t')
    
    # Process by chromosome grouping
    chromosomes = df['CHROM'].unique()
    filtered_dfs = []
    
    for chrom in chromosomes:
        chrom_df = df[df['CHROM'] == chrom].copy()
        chrom_df = chrom_df.sort_values('POS')
        
        # Extract site type list
        site_types = chrom_df['Site_Type'].tolist()
        positions = chrom_df['POS'].tolist()
        
        # Store window representative values for each site
        position_to_window_reprs = {pos: [] for pos in positions}
        
        # Sliding window processing
        prev_repr = None
        
        for i in range(len(site_types) - filter_window_size + 1):
            window = site_types[i:i+filter_window_size]
            window_repr = get_window_representative(window, prev_repr)
            
            # If window representative value is not None, update previous window representative value
            if window_repr is not None:
                prev_repr = window_repr
            
            # Assign window representative value to each site in window
            for j in range(filter_window_size):
                if i+j < len(positions):
                    position_to_window_reprs[positions[i+j]].append(window_repr)
        
        # Determine final SNP type for each site
        final_site_types = []
        
        for pos in positions:
            window_reprs = position_to_window_reprs[pos]
            if not window_reprs:
                final_site_types.append(site_types[positions.index(pos)])  # Keep original value
            else:
                # Calculate frequency of window representative values
                counter = Counter(window_reprs)
                total = len(window_reprs)
                
                # Find marker value with frequency >50%
                final_type = None
                for mark, count in counter.items():
                    if mark is not None and count / total > 0.5:
                        final_type = mark
                        break
                
                # If no marker value with frequency >50% found, keep original value
                if final_type is None:
                    final_type = site_types[positions.index(pos)]
                
                final_site_types.append(final_type)
        
        # Update chromosome DataFrame Site_Type column
        chrom_df['Site_Type'] = final_site_types
        filtered_dfs.append(chrom_df)
    
    # Merge results from all chromosomes
    result_df = pd.concat(filtered_dfs)
    
    # Write to output file
    result_df.to_csv(output_file, sep='\t', index=False)
    return result_df

def parse_vcf(vcf_file, parent_a, parent_b, progeny_samples=None, output_prefix=None):
    """Parse VCF file and find sites where both parents are homozygous but have different bases"""
    try:
        # Find indices of specified columns
        sample_indices = {}
        header_line = None
        
        # Read file header to get sample indices
        with open_file(vcf_file) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    header_line = line.strip()
                    headers = header_line.split('\t')
                    try:
                        sample_indices[parent_a] = headers.index(parent_a)
                        sample_indices[parent_b] = headers.index(parent_b)
                        
                        if progeny_samples:
                            for sample in progeny_samples:
                                if sample in headers:
                                    sample_indices[sample] = headers.index(sample)
                                else:
                                    print(f"Warning: Sample '{sample}' not found in VCF file")
                    except ValueError as e:
                        print(f"Error: Sample '{parent_a}' or '{parent_b}' not found in VCF file. Please ensure file format is correct.")
                        print(f"Sample list in VCF file: {headers[9:]}")
                        sys.exit(1)
                    break
        
        if not header_line:
            print("Error: Sample header line not found. Please ensure VCF file contains header line.")
            sys.exit(1)
            
        # print_info("Sample columns found", ', '.join(sample_indices.keys()))
            
        # Output informative sites file
        informative_file = f"{output_prefix}_informative_sites.tsv"
        
        # Find informative sites and output
        with open_file(vcf_file) as f, open(informative_file, 'w') as outfile:
            # Write header line
            outfile.write("CHROM\tPOS\tREF\tALT\tGT_A\tBASE_A\tGT_B\tBASE_B\tExpected_F1\n")
            
            informative_sites = []
            chromosome_counts = defaultdict(int)
            
            for line_num, line in enumerate(f, 1):
                # Skip comment lines
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) <= max(sample_indices.values()):
                    continue
                    
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                
                # Only process simple biallelic SNPs, skip complex variants
                if ',' in alt or '*' in alt:
                    continue
                
                # Get parent A and parent B data
                parent_a_data = fields[sample_indices[parent_a]].split(':')
                parent_b_data = fields[sample_indices[parent_b]].split(':')
                
                # Extract genotype
                format_field = fields[8]
                try:
                    gt_index = format_field.split(':').index('GT')
                    parent_a_gt = parent_a_data[gt_index]
                    parent_b_gt = parent_b_data[gt_index]
                except (ValueError, IndexError):
                    continue
                
                # Parse genotype
                parent_a_indices, parent_a_sep, status_a = parse_gt_field(parent_a_gt)
                parent_b_indices, parent_b_sep, status_b = parse_gt_field(parent_b_gt)
                
                if status_a != "success" or status_b != "success":
                    continue
                
                # Check if both parents are homozygous
                parent_a_is_homozygous = is_homozygous(parent_a_indices)
                parent_b_is_homozygous = is_homozygous(parent_b_indices)
                
                # Check if two parents have different genotypes
                if parent_a_is_homozygous and parent_b_is_homozygous:
                    parent_a_allele_idx = parent_a_indices[0]
                    parent_b_allele_idx = parent_b_indices[0]
                    
                    if parent_a_allele_idx != parent_b_allele_idx:
                        # Determine parent A base
                        if parent_a_allele_idx == 0:
                            base_a = ref
                        else:
                            base_a = alt
                        
                        # Determine parent B base
                        if parent_b_allele_idx == 0:
                            base_b = ref
                        else:
                            base_b = alt
                        
                        # Predict F1 genotype
                        expected_f1 = f"{base_a}/{base_b}"
                        
                        # Save informative site
                        informative_sites.append({
                            'chrom': chrom,
                            'pos': int(pos),
                            'ref': ref,
                            'alt': alt,
                            'gt_a': parent_a_gt,
                            'base_a': base_a,
                            'gt_b': parent_b_gt,
                            'base_b': base_b,
                            'expected_f1': expected_f1
                        })
                        
                        # Output to file
                        outfile.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{parent_a_gt}\t{base_a}\t{parent_b_gt}\t{base_b}\t{expected_f1}\n")
                        
                        # Count chromosome sites
                        chromosome_counts[chrom] += 1
        
        total_sites = len(informative_sites)
        
        # Display results
        print_result("INFORMATIVE SITES SUMMARY", {
            "Total sites found": total_sites,
            "Chromosomes processed": len(chromosome_counts),
            "Output file": informative_file
        })
        
        return informative_sites, chromosome_counts, sample_indices
    
    except Exception as e:
        print(f"Error: Exception occurred while processing VCF file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def load_informative_sites(informative_sites_file):
    """Load informative sites data from file"""
    try:

        
        # Store informative sites
        informative_sites = []
        chromosome_counts = defaultdict(int)
        
        with open(informative_sites_file, 'r') as f:
            # Skip header line
            header = f.readline()
            
            # Read data lines
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 9:  # Ensure at least 9 columns
                    continue
                
                chrom = fields[0]
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                
                ref = fields[2]
                alt = fields[3]
                gt_a = fields[4]
                base_a = fields[5]
                gt_b = fields[6]
                base_b = fields[7]
                expected_f1 = fields[8]
                
                # Save informative site
                informative_sites.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gt_a': gt_a,
                    'base_a': base_a,
                    'gt_b': gt_b,
                    'base_b': base_b,
                    'expected_f1': expected_f1
                })
                
                # Count chromosome sites
                chromosome_counts[chrom] += 1
        
        total_sites = len(informative_sites)
        
        print_result("INFORMATIVE SITES SUMMARY", {
            "Total sites found": total_sites,
            "Chromosomes processed": len(chromosome_counts)
        })
        
        return informative_sites, chromosome_counts
    
    except Exception as e:
        print(f"Error: Exception occurred while loading informative sites file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

def process_progeny_sample(informative_sites, vcf_file, sample_name, sample_index, output_prefix, filter_times=0, filter_window_size=5):
    """Process single sample, calculate genome contribution rate"""
    # 不显示processing信息
    
    # Statistical variables
    found_sites_count = 0
    
    # Convert informative sites to dictionary for fast lookup
    informative_dict = {}
    for site in informative_sites:
        informative_dict[(site['chrom'], site['pos'])] = site
    
    # Output detailed information file
    details_file = f"{output_prefix}_{sample_name}_details.tsv"
    
    # Save site information
    site_data = []
    
    with open_file(vcf_file) as infile:
        # Skip comment lines, find sample columns
        headers = None
        for line in infile:
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                break
            elif line.startswith('#'):
                continue
        
        if not headers:
            print(f"Error: Header line not found in file {vcf_file}")
            return None, None
        
        # Process data lines
        for line in infile:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) <= sample_index:
                continue
                
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue
                
            # Check if it's an informative site
            site_key = (chrom, pos)
            if site_key not in informative_dict:
                continue
                
            found_sites_count += 1
            
            # Get informative site data
            info = informative_dict[site_key]
            parent_allele_a = info['base_a']
            parent_allele_b = info['base_b']
            ref = fields[3]
            alt = fields[4]
            
            # Get sample genotype
            format_field = fields[8]
            sample_data = fields[sample_index].split(':')
            
            # Find GT index
            try:
                gt_index = format_field.split(':').index('GT')
                sample_gt = sample_data[gt_index]
            except (ValueError, IndexError):
                continue
            
            # Parse sample genotype
            sample_indices, separator, status = parse_gt_field(sample_gt)
            if status != "success":
                continue
            
            # Get sample actual alleles
            sample_alleles = []
            for idx in sample_indices:
                if idx == 0:
                    sample_alleles.append(ref)
                else:
                    sample_alleles.append(alt)
            
            if not sample_alleles:
                continue
            
            # Calculate matching with parents
            a_count = 0
            b_count = 0
            other_count = 0
            
            for allele in sample_alleles:
                if allele == parent_allele_a:
                    a_count += 1
                elif allele == parent_allele_b:
                    b_count += 1
                else:
                    other_count += 1
            
            # Determine site type
            site_type = "Unknown"
            
            if a_count == 1 and b_count == 1 and other_count == 0:
                site_type = "Het_AB"
            elif a_count == 2 and b_count == 0 and other_count == 0:
                site_type = "Hom_A"
            elif a_count == 0 and b_count == 2 and other_count == 0:
                site_type = "Hom_B"
            elif other_count > 0:
                # Ignore sites containing Other_Allele, don't add to site_data
                continue
            else:
                # Ignore other abnormal sites
                continue
            
            # Save site information
            site_data.append({
                'CHROM': chrom,
                'POS': pos,
                'Site_Type': site_type
            })
    
    # Create simplified detailed information file with only CHROM, POS, Site_Type columns
    df_simple = pd.DataFrame(site_data)
    
    # If no sites found, ensure creating an empty DataFrame with correct columns
    if len(df_simple) == 0:
        df_simple = pd.DataFrame(columns=['CHROM', 'POS', 'Site_Type'])
    
    df_simple.to_csv(details_file, sep='\t', index=False)
    
    # Apply filtering based on filter times
    filtered_df = df_simple
    if filter_times > 0:
        filtered_details_files = []
        
        for i in range(filter_times):
            filtered_details_file = f"{output_prefix}_{sample_name}_details_filtered_{i+1}.tsv"
            filtered_details_files.append(filtered_details_file)
            
            # Apply filtering
            filtered_df = filter_snps(
                details_file if i == 0 else filtered_details_files[i-1], 
                filtered_details_file,
                filter_window_size
            )
    
    # Use filtered data to calculate statistics
    filtered_sites_dict = {}
    for idx, row in filtered_df.iterrows():
        filtered_sites_dict[(row['CHROM'], row['POS'])] = row['Site_Type']
    
    # Calculate statistical data
    chr_stats = defaultdict(lambda: {
        'A': 0, 'B': 0, 
        'found': 0, 
        'effective': 0,
        'het': 0,
        'hom_a': 0,
        'hom_b': 0
    })
    
    sites_het = 0
    sites_hom_a = 0
    sites_hom_b = 0
    total_a_allele_count = 0
    total_b_allele_count = 0
    effective_sites_count = 0
    
    for (chrom, pos), site_type in filtered_sites_dict.items():
        chr_stats[chrom]['found'] += 1
        
        if site_type == "Het_AB":
            sites_het += 1
            chr_stats[chrom]['effective'] += 1
            chr_stats[chrom]['het'] += 1
            chr_stats[chrom]['A'] += 1
            chr_stats[chrom]['B'] += 1
            total_a_allele_count += 1
            total_b_allele_count += 1
            effective_sites_count += 1
            
        elif site_type == "Hom_A":
            sites_hom_a += 1
            chr_stats[chrom]['effective'] += 1
            chr_stats[chrom]['hom_a'] += 1
            chr_stats[chrom]['A'] += 2
            total_a_allele_count += 2
            effective_sites_count += 1
            
        elif site_type == "Hom_B":
            sites_hom_b += 1
            chr_stats[chrom]['effective'] += 1
            chr_stats[chrom]['hom_b'] += 1
            chr_stats[chrom]['B'] += 2
            total_b_allele_count += 2
            effective_sites_count += 1
    
    # Display processing results
    contribution_stats = {
        "Original sites": found_sites_count,
        "Effective sites": effective_sites_count,
        "Het_AB sites": sites_het,
        "Hom_A sites": sites_hom_a,
        "Hom_B sites": sites_hom_b
    }
    
    # Prepare results
    results = []
    
    # Add genome-wide statistics
    total_alleles = total_a_allele_count + total_b_allele_count
    if total_alleles > 0:
        contribution_a = (total_a_allele_count / total_alleles) * 100
        contribution_b = (total_b_allele_count / total_alleles) * 100
        
        contribution_stats["Parent A contribution"] = f"{contribution_a:.2f}%"
        contribution_stats["Parent B contribution"] = f"{contribution_b:.2f}%"
        
        print_result(f"SAMPLE {sample_name} RESULTS", contribution_stats)
        
        results.append({
            'sample': sample_name,
            'chrom': 'Whole_genome',
            'found': found_sites_count,
            'effective': effective_sites_count,
            'het': sites_het,
            'hom_a': sites_hom_a,
            'hom_b': sites_hom_b,
            'a_count': total_a_allele_count,
            'b_count': total_b_allele_count,
            'total': total_alleles,
            'contribution_a': contribution_a,
            'contribution_b': contribution_b
        })
    
    # Add chromosome-level statistics
    for chrom in sorted(chr_stats.keys(), key=lambda x: (len(x), x)):
        stats = chr_stats[chrom]
        chr_a_count = stats['A']
        chr_b_count = stats['B']
        chr_total = chr_a_count + chr_b_count
        
        if chr_total > 0:
            chr_contribution_a = (chr_a_count / chr_total) * 100
            chr_contribution_b = (chr_b_count / chr_total) * 100
            
            results.append({
                'sample': sample_name,
                'chrom': chrom,
                'found': stats['found'],
                'effective': stats['effective'],
                'het': stats['het'],
                'hom_a': stats['hom_a'],
                'hom_b': stats['hom_b'],
                'a_count': chr_a_count,
                'b_count': chr_b_count,
                'total': chr_total,
                'contribution_a': chr_contribution_a,
                'contribution_b': chr_contribution_b
            })
    
    return results, filtered_df

def calculate_contribution_rates(informative_sites, vcf_file, sample_indices, parent_a, parent_b, progeny_samples, output_prefix, filter_times=0, filter_window_size=5, chromosome_counts=None, max_workers=1):
    """Calculate genome contribution rates of progeny samples"""
    print_section("Calculating Contribution Rates")
    start_time = time.time()
    
    # Only process samples specified by -s parameter
    actual_progeny_samples = {}
    for sample_name in progeny_samples:
        if sample_name in sample_indices:
            actual_progeny_samples[sample_name] = sample_indices[sample_name]
    
    if not actual_progeny_samples:
        print_status("No valid samples to test found, skipping contribution rate calculation", "warning")
        return
    
    # Summary file
    summary_file = f"{output_prefix}_summary.tsv"
    with open(summary_file, 'w') as f_summary:
        f_summary.write("Sample\tChromosome\tFound_Sites\tEffective_Sites\tHet_Sites\tHom_A_Sites\tHom_B_Sites\tParent_A_Contribution(%)\tParent_B_Contribution(%)\n")
    print_info("Summary file", summary_file)
    
    # Process samples
    all_results = []
    
    # Determine number of processes/threads
    if max_workers is None or max_workers <= 0:
        max_workers = 1
    
    if max_workers == 1:
        # Single-threaded processing
        for sample_name, sample_index in actual_progeny_samples.items():
            results, df = process_progeny_sample(
                informative_sites, 
                vcf_file, 
                sample_name, 
                sample_index, 
                output_prefix,
                filter_times,
                filter_window_size
            )
            
            if results:
                all_results.extend(results)
                
                # Write results to summary file
                with open(summary_file, 'a') as f_summary:
                    for result in results:
                        if 'contribution_a' in result and 'contribution_b' in result:
                            f_summary.write(f"{result['sample']}\t{result['chrom']}\t"
                                          f"{result['found']}\t{result['effective']}\t{result['het']}\t"
                                          f"{result['hom_a']}\t{result['hom_b']}\t"
                                          f"{result['contribution_a']:.4f}\t{result['contribution_b']:.4f}\n")
    else:
        # Multi-threaded/multi-process processing
        executor_class = concurrent.futures.ProcessPoolExecutor if sys.platform != "win32" else concurrent.futures.ThreadPoolExecutor
        
        with executor_class(max_workers=max_workers) as executor:
            # Submit tasks
            future_to_sample = {
                executor.submit(
                    process_progeny_sample,
                    informative_sites,
                    vcf_file,
                    sample_name,
                    sample_index,
                    output_prefix,
                    filter_times,
                    filter_window_size
                ): sample_name for sample_name, sample_index in actual_progeny_samples.items()
            }
            
            # Process results
            for future in concurrent.futures.as_completed(future_to_sample):
                sample_name = future_to_sample[future]
                
                try:
                    results, df = future.result()
                    
                    if results:
                        all_results.extend(results)
                        
                        # Write results to summary file
                        with open(summary_file, 'a') as f_summary:
                            for result in results:
                                if 'contribution_a' in result and 'contribution_b' in result:
                                    f_summary.write(f"{result['sample']}\t{result['chrom']}\t"
                                                  f"{result['found']}\t{result['effective']}\t{result['het']}\t"
                                                  f"{result['hom_a']}\t{result['hom_b']}\t"
                                                  f"{result['contribution_a']:.4f}\t{result['contribution_b']:.4f}\n")
                except Exception as exc:
                    print(f"Error processing sample {sample_name}: {exc}")
                    import traceback
                    traceback.print_exc()
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    # Final summary
    print(f"✓ Processed {len(actual_progeny_samples)} samples in {processing_time:.1f}s")

def parse_samples(samples_arg):
    """Parse samples parameter, support comma-separated sample list or txt file (one sample per line)"""
    # Check if it's a txt file
    if samples_arg.endswith('.txt'):
        # Check if file exists
        if not os.path.exists(samples_arg):
            print(f"Error: Sample list file {samples_arg} does not exist")
            sys.exit(1)
            
        # Read samples from file
        samples = []
        with open(samples_arg, 'r') as f:
            for line in f:
                sample = line.strip()
                if sample:  # Skip empty lines
                    samples.append(sample)
        
        if not samples:
            print(f"Error: No valid samples found in sample list file {samples_arg}")
            sys.exit(1)
            

        return samples
    else:
        # Process comma-separated sample list
        return samples_arg.split(',')

def main():
    """Run genome contribution rate calculation tool"""
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(
        description='GCCVision: Genome Contribution Calculator and Visualizer v1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Description:
  This program specifically processes multi-sample merged VCF files processed by GATK, implementing the following workflow:
  1. Find informative sites from VCF file where both parents are homozygous but have different bases
  2. Record site types (Het_AB/Hom_A/Hom_B) for each sample (keeping only CHROM, POS, Site_Type)
  3. Use sliding window method for specified number of filtering rounds to remove possible error sites (optional)
  4. Calculate genome contribution rates for each sample based on filtered sites
  5. Visualize the results in a browser-based interactive HTML report
  
Examples:
  GCCVision.py -v input.vcf -a Parent1 -b Parent2 -s Sample1,Sample2 -o results -t 4 -f 2 -w 5
  GCCVision.py -v input.vcf -a Parent1 -b Parent2 -s samples.txt -o results -t 4 -f 2 -w 5
  
Using existing informative sites file:
  GCCVision.py -i results_informative_sites.tsv -v input.vcf -s Sample1,Sample2 -o results -t 4 -f 2 -w 5
  GCCVision.py -i results_informative_sites.tsv -v input.vcf -s samples.txt -o results -t 4 -f 2 -w 5
        """
    )
    parser.add_argument('-v', '--vcf', dest='vcf_file', required=True, help='Input multi-sample merged VCF file (can be .gz compressed file)')
    parser.add_argument('-i', '--informative-sites', dest='informative_sites_file', help='Existing informative sites file, if specified will skip finding informative sites step')
    parser.add_argument('-a', '--parent-a', help='Parent A sample name')
    parser.add_argument('-b', '--parent-b', help='Parent B sample name') 
    parser.add_argument('-s', '--samples', required=True, help='Samples to test. Can be: 1) comma-separated sample name list; 2) txt file containing one sample per line')
    parser.add_argument('-o', '--output-prefix', required=True, help='Output file prefix')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of parallel processing threads, default is 1')
    parser.add_argument('-f', '--filter-times', type=int, default=0, help='Number of times to filter SNP sites using sliding window method, default is 0 (no filtering)')
    parser.add_argument('-w', '--filter-window-size', type=int, default=5, help='Sliding window size, default is 5')
    
    args = parser.parse_args()
    
    # Program header
    print_header("GCCVision: Genome Contribution Calculator and Visualizer v1.0")
    
    # Check if VCF file exists
    if not os.path.exists(args.vcf_file):
        print_status(f"Input file {args.vcf_file} does not exist", "error")
        sys.exit(1)
    
    # Parse samples parameter and display configuration
    print_step(1, "Configuration")
    progeny_samples = parse_samples(args.samples)
    
    # Display configuration
    config_data = {
        "Input VCF": args.vcf_file,
        "Output prefix": args.output_prefix,
        "Threads": args.threads,
        "Filtering": f"{args.filter_times} times (window={args.filter_window_size})" if args.filter_times > 0 else "Disabled"
    }
    
    # Add parent info if available
    if args.parent_a and args.parent_b:
        config_data["Parent A"] = args.parent_a
        config_data["Parent B"] = args.parent_b
    
    print_result("CONFIGURATION", config_data)
    
    if args.informative_sites_file:
        # Use existing informative sites file mode
        if not os.path.exists(args.informative_sites_file):
            print_status(f"Informative sites file {args.informative_sites_file} does not exist", "error")
            sys.exit(1)
        
        print_step(2, "Loading Informative Sites")
        
        # Load informative sites
        informative_sites, chromosome_counts = load_informative_sites(args.informative_sites_file)
        
        print_step(3, "Calculating Contribution Rates")
        # Read VCF file header line to get sample indices
        sample_indices = {}
        with open_file(args.vcf_file) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    headers = line.strip().split('\t')
                    for i, header in enumerate(headers):
                        if header in progeny_samples:
                            sample_indices[header] = i
                    break
        
        if not sample_indices:
            print_status("No specified samples found in VCF file", "warning")
            sys.exit(1)
        
        # Calculate genome contribution rates
        calculate_contribution_rates(
            informative_sites,
            args.vcf_file,
            sample_indices,
            "ParentA",  # Use default names since we only care about genome contribution rates
            "ParentB",
            progeny_samples,
            args.output_prefix,
            args.filter_times,
            args.filter_window_size,
            chromosome_counts,
            args.threads
        )
    else:
        # Standard mode: analyze informative sites from VCF file
        if not args.parent_a or not args.parent_b:
            print_status("Parent A and B parameters required in standard mode", "error")
            sys.exit(1)
        
        print_step(2, "Identificating Informative Sites")
        
        # Find informative sites
        informative_sites, chromosome_counts, sample_indices = parse_vcf(
            args.vcf_file, 
            parent_a=args.parent_a,
            parent_b=args.parent_b,
            progeny_samples=progeny_samples,
            output_prefix=args.output_prefix
        )
        
        print_step(3, "Calculating Contribution Rates")
        # Calculate genome contribution rates
        calculate_contribution_rates(
            informative_sites,
            args.vcf_file,
            sample_indices,
            args.parent_a,
            args.parent_b,
            progeny_samples,
            args.output_prefix,
            args.filter_times,
            args.filter_window_size,
            chromosome_counts,
            args.threads
        )
    
    print_status("All processing completed successfully", "success")

if __name__ == "__main__":
    main() 