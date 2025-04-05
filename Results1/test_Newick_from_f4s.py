"""
==============================================================================
Title:          test_Newick_from_f4s.py
Description:    Find violations of clades in a Newick tree using a table of f4 statistics.

Author:         Javier Maravall-López
Affiliation:    Harvard University
Email:          fmaravalllopez@fas.harvard.edu
Date Created:   April 17, 2023
Last Modified:  October 24, 2024

Usage: 
    $ python test_Newick_from_f4s.py -t <treefile> -f <f4table> [-z <Zthreshold>] [-o <output>]

    Arguments:
        -t, --treefile       : Path to a tree in Newick format (required)
        -f, --f4table        : Path to a text f4 table with rows of the form 
                               [Outgroup {Population in the tree} {Population in the tree} 
                               {Population in the tree} Z-score] (required)
        -z, --Zthreshold     : Z-score absolute value threshold for significance (default: 3)
        -o, --output         : Path to an output file (optional)

Inputs:
    - treefile: 
        Path to a tree file in Newick format. This file defines the tree structure with nodes
        representing populations and branches representing cladality relationships between them.

    - f4table: 
        Path to a text file representing the f4 table. Each row should have the form: 
        [Outgroup {Population in the tree} {Population in the tree} {Population in the tree} Z-score].

    - Zthreshold: 
        Optional argument specifying the Z-score absolute value threshold for significance. 
        Default value is set to 3.

    - output: 
        Optional path to the output file. If not specified, results will be printed to the console.

Outputs:
    - The output file (or console output) contains cladality violations found for the different
    clades in the tree as this is traversed, and the statistics that support such violations.

Dependencies:
    - ete3          : For working with phylogenetic trees.
    - argparse      : For parsing command-line arguments.
    - itertools     : For efficient looping and iterator functions.
    - sys           : For system-specific parameters and functions.

Version information:
    - Python:        3.9.5
    - ete3:          3.1.3
    - argparse:      1.1
    - sys:           3.9.5 (default, Jun  4 2021, 12:28:51) [GCC 7.5.0]
        
Notes:
    - It is important that f4-statistics be of the form f4(Outgroup, P, P, P), where P is the set 
    of ALL populations present in the tree. 
    Symmetric statistics (i.e., f(O, A, B, C)=f(O, A, C, B)) will be automatically filled in 
    by the script and need not be present in the table. 

==============================================================================
"""




import ete3
import argparse
import itertools
import sys

def get_arguments():
    parser = argparse.ArgumentParser(description='Arguments')
    
    # Create an argument group for required arguments
    required = parser.add_argument_group('required arguments')
    
    required.add_argument('-t', '--treefile',
                           type=str,
                           required=True,
                           help='Path to a tree in Newick format')
    
    required.add_argument('-f', '--f4table',
                           type=str,
                           required=True,
                           help='Path to a text f4 table with rows of the form [Outgroup {Population in the tree} {Population in the tree} {Population in the tree} Z-score]')
    
    parser.add_argument('-z', '--Zthreshold',
                        type=float,
                        default=3,
                        help='Z-score absolute value threshold for significance (default: 3)')
                        
    parser.add_argument('-o', '--output',
                        type=str,
                        help='Path to an output file')
    
    args = parser.parse_args()
    return args

def test_tree(t, f4_fn, threshold):
    
    # Read f4_table from a file
    with open(f4_fn) as f:
        f4_table = [tuple(line.strip().split()) for line in f]
    
    # Create a set of tuples from the f4_table
    f4_set = set(f4_table)
    
    # Iterate over a copy of the f4_set
    for row in f4_set.copy():
        # Create a new row with A, B, D, C and opposite sign of the original row
        new_row = (row[0], row[1], row[3], row[2], str(-1*float(row[4])))
        # If the new row is not in the f4_set, add it to the f4_table
        if new_row not in f4_set:
            f4_set.add(new_row)
    
    # Convert the f4_set back to a list
    f4_table = list(f4_set)
    
    Outgroup = f4_table[0][0]
    
    print("\n We will test the following tree:")
    
    print(t)
    
    # Create an empty set to store the outermost leaves 
    outer_leaves = set()
    
    # Traverse the tree
    for node in t.traverse():
        # Check if the node is a leaf
        if node.is_leaf():
            
            # Check if the parent node has no internal nodes as children
            if all(child.is_leaf() for child in node.up.get_children()):
                outer_leaves.add(node)
    
    # Keep only one sister leaf per pair of sisters
    
    processed = set()
    target_leaves = outer_leaves.copy()
                
    for node in outer_leaves: 
        
        if not node in processed: 
            
            sister = node.get_sisters()[0]
            
            target_leaves.remove(sister)
            
            processed.add(node)
            processed.add(sister)
            
    print("\nA set of target (outermost, non-sister) leaves is: " + ", ".join(leaf.name for leaf in target_leaves))
            
    total = t.get_leaf_names()
    
    i=1
    
    for l in target_leaves:
        
        print("\n\n\n\n" + str(i) + ". Testing from leaf", l.name, "\n")
    
        p_1 = l.up
        
        j = 1
        
        # Plot clade
        
        print("\n" + str(i) + "." + str(j) + ". Testing the following clade:\n")
        
        print(p_1)
        
        # Compute sets of populations in the 2nd, 4th and 3rd positions respectively 
        
        v = set(l.get_leaf_names())
        
        u = set(p_1.get_leaf_names())-v
        
        t = set(total) - u - v
        
        print("\nThat is, looking for evidence to reject the hypothesis that the set of populations", u, end=" ")
        
        print("\nforms a clade with the set of populations", v, end=" ")
    
        print("\nwith respect to the rest of the populations in the tree...")
    
        ## Look for evidence 
    
        Subtable = []
        for p1, p2, p3 in itertools.product(t, u, v):
            Subtable.append((Outgroup, p1, p2, p3))
    
        results = []
        for row in f4_table:
            if (row[0], row[1], row[2], row[3]) in Subtable:
                results.append(row)
    
        evidence_found = False
    
        for row in results:
            if abs(float(row[4])) > threshold:
                if not evidence_found:
                    print("\nEvidence found in the form of the following results:")
                    evidence_found = True
        
                output = "\tThe statistic f({0}, {1}, {2}, {3}) has a Z-score of {4}.".format(row[0], row[1], row[2], row[3], row[4])
                print(output)
        
        if evidence_found:
            print("\nClade rejected.")
        
        if not evidence_found:
            print("\nNo evidence found.")
            print("Clade not rejected.")
    
        while not p_1.up.is_root():
        
            print("\nMoving one level up the tree...\n")
            
            j += 1
        
            p_2 = p_1.up
            
            print("\n" + str(i) + "." + str(j) + ". Testing the following clade:\n")
        
            print(p_2)
        
            v = set(p_1.get_leaf_names())
        
            aux = set(p_2.get_leaf_names())
    
            u = aux - v
        
            t = set(total) - aux
            
            print("\nThat is, looking for evidence to reject the hypothesis that the set of populations", u, end=" ")
            
            print("\nforms a clade with the set of populations", v, end=" ")
        
            print("\nwith respect to the rest of the populations in the tree...")
        
            ## Look for evidence 
        
            Subtable = []
            for p1, p2, p3 in itertools.product(t, u, v):
                Subtable.append((Outgroup, p1, p2, p3))
        
            results = []
            for row in f4_table:
                if (row[0], row[1], row[2], row[3]) in Subtable:
                    results.append(row)
        
            evidence_found = False
        
            for row in results:
                if abs(float(row[4])) > threshold:
                    if not evidence_found:
                        print("\nEvidence found in the form of the following results:")
                        evidence_found = True
            
                    output = "\tThe statistic f({0}, {1}, {2}, {3}) has a Z-score of {4}.".format(row[0], row[1], row[2], row[3], row[4])
                    print(output)
            
            if evidence_found:
                print("\nClade rejected.")
            
            if not evidence_found:
                print("\nNo evidence found.")
                print("Clade not rejected.")
        
            p_1 = p_2
        
        i += 1
        
    print("\n\n\n\nThe set of target leaves has been exhausted.")
    print("Tree testing completed.\n\n\n")
    
args = get_arguments()

# Import arguments 
tree_fn = args.treefile
f4_fn = args.f4table
threshold = args.Zthreshold
output_fn = args.output
    
# Import the tree
t = ete3.Tree(tree_fn)

if output_fn:
    with open(output_fn, 'w') as f:

        sys.stdout = f
        
        test_tree(t, f4_fn, threshold)
        
    sys.stdout = sys.__stdout__
        
else:
    
    test_tree(t, f4_fn, threshold)
